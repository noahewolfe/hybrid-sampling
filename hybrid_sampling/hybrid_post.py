#!/usr/bin/env python
""" Script to perform hybrid analysis based on previous result """
import os
import sys
from copy import deepcopy

import bilby
import numpy as np
from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_component_masses
from bilby_pipe.data_analysis import DataAnalysisInput
from bilby_pipe.utils import convert_string_to_dict, log_version_information, logger

# fmt: off
import matplotlib  # isort:skip
matplotlib.use("agg")
# fmt: on


def tempered_weights(nested_samples, beta):
    """

    Tempers the posterior weights from a nested sampling analysis. See Section B and
    Eq.~10 of https://arxiv.org/abs/2208.12872.

    Notes
    =====
    We use the notation of https://arxiv.org/abs/2208.12872. First we read in or compute
    the log prior volume weights, :math:`w_i`, up to an overall normalization. If those
    are not directly provided, we compute it as

    .. math::
        \ln w_i \propto \ln p_i - \ln \cal L_i.

    To get tempered posterior weights at an inverse temperature :math:`\beta_T`, we have

    .. math::
        \ln p_{i, \beta_T} \propto \ln w_i + \beta_T \ln \cal L_i.

    Finally, we normalize these tempered posterior weights by the tempered evidence,
    which is computed as

    .. math::
        \ln \cal Z_{\beta_T} = \ln \sum_i \exp \left( \ln w_i + \beta_T \ln \cal L_i
        \right).

    Note that we return :math:`p_{i, \beta_T}`, not :math:`\ln p_{i, \beta_T}`.

    Parameters
    ==========
    nested_samples: pandas.DataFrame
        Resultant dataframe from a `dynesty` run, or generally, the `nested_samples`
        member of a `bilby.core.result` object. (This code has not yet been tested
        with other nested samplers.)
    beta: float or array_like
        Inverse temperature(s) at which to compute tempered posterior weights.

    Returns
    =======
    array_like: tempered posterior weights, shaped like (number of nested samples, )
    or (number of nested samples, number of inverse temperatures).

    """

    from scipy.special import logsumexp

    if "weights" in nested_samples:
        ln_weights = (
            np.log(nested_samples["weights"]) - nested_samples["log_likelihood"]
        ).values
    else:
        ln_weights = nested_samples["log_prior_volume"].values

    if isinstance(beta, (list, np.ndarray)):
        ln_weights = np.tile(ln_weights, (len(beta), 1))
        ln_weights += (
            np.tile(nested_samples["log_likelihood"].values, (len(beta), 1)).T * beta
        ).T
    else:
        ln_weights += nested_samples["log_likelihood"].values * beta

    return np.exp(ln_weights - logsumexp(ln_weights))


def set_tempered_nested_samples(
    samples, nested_samples, parameter_keys, nwalkers, temperatures, set_idxs=None
):
    """

    Sets, in place, an array (or dictionary) of tempered nested samples.

    In particular, set an array (or dictionary) representing a subsampling of
    nested samples, where the probability of including a particular sample
    is determined by the associated nested sampling weights which have been
    raised to the power of a "temperature" beta.

    Parameters
    ==========
    samples: dict
        Dictionary of samples; the value at each key is an np.ndarray array
        shaped like (ntemps, nwalkers)
    nested_samples: pandas.DataFrame
        pandas DataFrame of nested samples, indexed by parameter names.
    parameter_keys: list
        List of strings, of the names of the parameters that we want to sample.
        This may represent a subset of the keys of the samples dictionary.
    nwalkers: int
        Number of walkers in our eventual ensemble of parallel-tempered MCMC chains.
    temperatures: np.ndarray or list
        List or array of temperatures beta (of length we call "ntemps").
    set_idxs: tuple
        Numpy-comptaible tuple of arrays for indexing the values at each key of samples;
        if we only want to set the value of samples at certain (temperature, walker)
        positions in the arrays at each key of samples.

    Returns
    =======
    None (this modifies the dictionary samples).

    """

    ndim = len(parameter_keys)
    ntemps = len(temperatures)

    if set_idxs is None:
        # indices for each position in an array (ntemps, nwalkers)
        set_idxs = tuple(np.mgrid[0:ntemps, 0:nwalkers])

    tempered_sample_idxs = np.array(
        [
            np.random.choice(
                len(nested_samples),
                (nwalkers, ndim),
                p=tempered_weights(nested_samples.copy(), temp),
                replace=True,
            )
            for temp in temperatures
        ]
    )
    # --> shape is (ntemps, nwalkers, ndim) b/c we generate an array
    # shapes like (nwalkers, ndim) at each temp

    # then, move the last axis to the front, so it's (ndim, ntemps, nwalkers)
    tempered_sample_idxs = np.moveaxis(tempered_sample_idxs, -1, 0)

    for k, key in enumerate(parameter_keys):
        samples[key][set_idxs] = nested_samples[key].values[
            tempered_sample_idxs[k, :, :]
        ]


class HybridInput(DataAnalysisInput):
    def __init__(self, args, unknown_args, test=False):
        super(HybridInput, self).__init__(
            args=args, unknown_args=unknown_args, test=test
        )
        self.extra_prior = args.extra_prior
        self.extra_initialization = args.extra_initialization
        self.extra_waveform_arguments = args.extra_waveform_arguments
        self.extra_label = args.extra_label

    def get_likelihood_and_priors(self):
        likelihood, priors = super(HybridInput, self).get_likelihood_and_priors()
        if self.extra_prior is not None:
            priors.update(bilby.core.prior.PriorDict(self.extra_prior))
        return likelihood, priors

    def _set_pos0_within_prior(
        self, pos0, nested_samples, likelihood, parameter_keys, temperatures, nwalkers
    ):
        """
        Modify pos0, in place, such that all initial points are within the
        prior support.

        Parameters
        ==========
        pos0: dict
            Dictionary of initial points; the value at each key is an np.ndarray array
            shaped like (ntemps, nwalkers)
        nested_samples: pandas.DataFrame
            pandas DataFrame of nested samples, indexed by parameter names.
        likelihood: bilby.gw.likelihood.GravitationalWaveTransient
            bilby gravitational wave likelihood object that can evaluate
            the log likelihood ratio at pos0 parameter values.
        parameter_keys: list
            List of strings, of the names of the parameters that we want to sample.
            This may represent a subset of the keys of the samples dictionary.
        temperatures: np.ndarray or list
            List or array of temperatures beta (of length we call "ntemps").
        nwalkers: int
            Number of walkers in our eventual ensemble of parallel-tempered MCMC chains.
        
        Returns
        =======
        None (this modifies the dictionary samples).

        """
        
        logger.info(
            "Adjusting pos0 seed samples to fit within prior/overlap constraints"
        )

        old_parameters = deepcopy(likelihood.parameters)

        ntemps = len(temperatures)
        extra_init_prior_dict = bilby.core.prior.PriorDict(self.extra_initialization)

        test_log_like_ratio = np.full(
            shape=(ntemps, nwalkers),
            fill_value=-np.inf
        )

        inf_indices = np.nonzero(test_log_like_ratio == -np.inf)
        
        tries = 0
        while inf_indices[0].shape[0] > 0:
            badtemps = np.unique(inf_indices[0])
            badwalkers = np.unique(inf_indices[1])

            n_badwalkers = len(badwalkers)

            set_tempered_nested_samples(
                pos0,
                nested_samples,
                parameter_keys,
                n_badwalkers,
                temperatures[badtemps],
                set_idxs=np.ix_(badtemps, badwalkers),
            )

            new_extra_init_prior_samples = extra_init_prior_dict.sample(
                size=(inf_indices[0].shape[0])
            )
            for k in extra_init_prior_dict.keys():
                pos0[k][inf_indices] = new_extra_init_prior_samples[k]

            for i,j in np.transpose(inf_indices):
                parameters = { k : pos0[k][i,j] for k in pos0.keys() }
                likelihood.parameters.update(parameters)
                test_log_like_ratio[i,j] = likelihood.log_likelihood_ratio()

            inf_indices = np.nonzero(test_log_like_ratio == -np.inf)
            tries += 1

        if old_parameters is not None:
            likelihood.parameters = old_parameters
        else:
            likelihood.parameters = None
        
        logger.info(f"Done, after {tries} tries.")

    def setup_initial_points(self):
        from ptemcee import default_beta_ladder

        initial_result = bilby.core.result.read_in_result(
            f"{self.outdir}/result/{self.label}_result.{self.result_format}"
        )
        ntemps = self.sampler_kwargs.get("ntemps", 10)
        nwalkers = self.sampler_kwargs.get("nwalkers", 5)
        temperatures = default_beta_ladder(
            ndim=len(initial_result.search_parameter_keys) + 1,
            ntemps=ntemps,
            Tmax=self.sampler_kwargs.get("Tmax", None),
        )

        gr_parameters = initial_result.search_parameter_keys

        likelihood, complete_priors = self.get_likelihood_and_priors()
        pos0 = complete_priors.sample((ntemps, nwalkers))
        set_tempered_nested_samples(
            pos0, initial_result.nested_samples, gr_parameters, nwalkers, temperatures
        )

        if self.extra_prior is not None:
            self._set_pos0_within_prior(
                pos0,
                initial_result.nested_samples,
                likelihood,
                gr_parameters,
                temperatures,
                nwalkers,
            )

        self.sampler_kwargs["pos0"] = pos0

    def get_default_waveform_arguments(self):
        wfa = super(HybridInput, self).get_default_waveform_arguments()
        if self.extra_waveform_arguments is not None:
            wfa.update(convert_string_to_dict(self.extra_waveform_arguments))
        inverse_asd = (
            np.sum(
                [1 / ifo.power_spectral_density_array for ifo in self.interferometers],
                axis=0,
            )
            ** 0.5
        )
        wfa["inverse_asd"] = np.nan_to_num(inverse_asd, nan=0, posinf=0)
        return wfa

    def run_sampler(self):
        self.setup_initial_points()
        self.label += f"_{self.extra_label}"
        self.sampler = "ptemcee"

        super(HybridInput, self).run_sampler()


def create_hybrid_parser():
    """Hybrid data analysis parser creation"""
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="bilby_pipe hybrid sampling stage")
    parser.add_argument("result", type=str)
    parser.add_argument("--extra-prior-file", type=str, default=None)
    parser.add_argument("--extra-initialization-prior", type=str, default=None)
    parser.add_argument("--waveform-arguments", type=str, default=None)
    parser.add_argument("--extra-label", type=str, default="hybrid")
    return parser


def parse_args(parser):
    from argparse import Namespace

    args = parser.parse_args()
    old_result = bilby.core.result.read_in_result(args.result)
    hybrid_args = old_result.meta_data["command_line_args"].copy()
    hybrid_args["data_dump_file"] = os.path.abspath(old_result.meta_data["data_dump"])
    hybrid_args["label"] = old_result.label
    hybrid_args["outdir"] = old_result.outdir[:-7]
    hybrid_args["extra_prior"] = args.extra_prior_file
    hybrid_args["extra_waveform_arguments"] = args.waveform_arguments
    hybrid_args["extra_label"] = args.extra_label
    if args.extra_initialization_prior is not None:
        hybrid_args["extra_initialization"] = args.extra_initialization_prior
    else:
        hybrid_args["extra_initialization"] = args.extra_prior_file
    return (Namespace(**hybrid_args), old_result.meta_data["unknown_command_line_args"])


def main():
    """Data analysis main logic"""
    parser = create_hybrid_parser()
    args, unknown_args = parse_args(parser)
    log_version_information()
    analysis = HybridInput(args, unknown_args)
    analysis.run_sampler()
    sys.exit(0)
