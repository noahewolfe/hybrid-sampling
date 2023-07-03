# hybrid-sampling
Code to implement the "hybrid sampling" method for initializing parallel-tempered ensemble Markov Chain Monte Carlo posterior sampling in bilby.
Currently, this is written to perform phenomenological tests of general relativity (GR), but we'd love to generalize the code to handle any data analysis problems involving model misspecification and/or upgrading from more simple to more complex models.
For conceptual details on the hybrid sampling method generally, and its particular application to tests of GR, see ["Accelerating Tests of General Relativity with Hybrid Sampling"](https://arxiv.org/abs/2208.12872).

# Installation

0. (For tests of GR) Install a [lightly-modified fork of the `bilby_tgr` package](https://git.ligo.org/noah.wolfe/bilby-tgr)
1. Clone this repository
2. Install via pip, e.g. inside the cloned directory:
`pip install .`

Note, the requirements list that `bilby_pipe` version 1.0.6 or later is required. It
likely also works with some earlier versions of `bilby_pipe`, but this has not been
explicitly tested.

# Quickstart

A good place to start is the hybrid analysis of an injected (simulated) beyond-GR binary black hole merger signal, in `examples/beyond_gr_injection`.
This assumes that you are working on a cluster with the `condor` submit system (e.g. the LIGO CIT cluster).

To run this example, once you're working in that directory, all you need to do is `bilby_pipe_build_hybrid config.ini --submit`.
This will first sample the posterior on the intrinsic + extrinsic binary black hole parameters assuming GR describes the waveform, and then launch 28 runs, each sampling in those parameters plus one of the phenomelogical deviations from GR (with either an overlap cut on the prior for the deviation, or not).

Conceptually, here are the things needed to setup a bilby_pipe-style hybrid sampling analysis configuration file for use with this code:
- Write a `bilby_pipe` configuration file to analyze a simulated signal / real data as you would normally (see [their docs](https://lscsoft.docs.ligo.org/bilby_pipe/0.3.12/index.html) for more information).
    - In the example, this is `config.ini`.
- Use a waveform model that assumes the most general form of the problem you want to study (e.g. to the `frequency-domain-source-model` argument).
    - For testing GR we use the frequency-domain source model `bilby_tgr.source.generic_non_gr_binary_black_hole` which allows for generic, phenomeological deviations from GR.
    - 
- Supply a `.prior` file to the `prior-file` argument under the "simple" waveform model you want to assume in the first step of sampling.
    - In the example, this is `modified.prior`, which is a standard set of priors for a binary black hole merger plus all of the phenomeological GR deviations set to zero.
    - To emphasize, **it is under these prior assumptions that `dynesty` samples in to generate the initial points for the second step of sampling with `ptemcee`.**
- Supply "seed" distributions from which to draw initial points for parameters that were not sampled in during the first step of sampling
    - Stored in a directory, specified by `hybrid-seed-priors-dir`.
    - In the example, this is `init_injection`.
- Supply prior distributions for parameters that were not sampled in during the first step of sampling.
    - Stored in a directory, specified by `hybrid-priors-dir`.
    - In the example, this is `new_injection`.
- List each of the extra parameters to sample in and the overlap cut to apply to their prior.
    - Listed in a file specified by `hybrid-runs`.
    - In the example, this is `queue.txt`.
    - **For each combination of parameter + overlap cut, you will get a new `ptemcee` analysis** (Currently, the code does not allow for specifying multiple new parameters to analyze simultaneously.)
- Finally, make sure to include the following arguments in the 

See `bilby_pipe_build_hybrid --help` for more detail.

Note that `bilby_pipe_build_hybrid` **accepts all of the same command-line arguments
and configuration files as `bilby_pipe`**; the main difference is that `bilby_pipe_build_hybrid` also looks for three additional command-line arguments,
- `hybrid-seed-priors-dir`: path to a directory containing distributions from which to draw initial points for any parameters not sampled in during the `dynesty` analysis step. In the beyond-GR injection example, this is `init_injection`.
- `hybrid-priors-dir`: path to a directory containing prior distributions for parameters not sampled in during the `dynesty` analysis step. In the example, this is `new_injection`.
- `hybrid-runs`: list of additional parameters to sample in during the `ptemcee` step, and the prior overlap cut to apply on that parameter. In the example, this is `queue.txt`.
- `hybrid-label`: custom label to differentiate the `ptemcee` analyses from the `dynesty` analysis.

# Citation

If you use this method in your work, please include the following citation:

```
@article{Wolfe:2022nkv,
    author = "Wolfe, Noah E. and Talbot, Colm and Golomb, Jacob",
    title = "{Accelerating tests of general relativity with gravitational-wave signals using hybrid sampling}",
    eprint = "2208.12872",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.107.104056",
    journal = "Phys. Rev. D",
    volume = "107",
    number = "10",
    pages = "104056",
    year = "2023"
}
```