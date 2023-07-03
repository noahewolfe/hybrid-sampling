# hybrid-sampling
Code to implement the "hybrid sampling" method for initializing parallel-tempered ensemble Markov Chain Monte Carlo posterior sampling in bilby.
Currently, this is written to perform phenomenological tests of general relativity (GR), but we'd love to generalize the code to handle any data analysis problems involving model misspecification and/or upgrading from more simple to more complex models.
For conceptual details on the hybrid sampling method generally, and its particular application to tests of GR, see ["Accelerating Tests of General Relativity with Hybrid Sampling"](https://arxiv.org/abs/2208.12872).

# Installation

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

Note that `bilby_pipe_build_hybrid` **accepts all of the same command-line arguments
and configuration files as `bilby_pipe`**; the main difference is that `bilby_pipe_build_hybrid` also looks for three additional command-line arguments,
- `hybrid-seed-priors-dir`: path to a directory containing distributions from which to draw initial points for any parameters not sampled in during the `dynesty` analysis step. In the beyond-GR injection example, this is `init_injection`.
- `hybrid-priors-dir`: path to a directory containing prior distributions for parameters not sampled in during the `dynesty` analysis step. In the example, this is `new_injection`.
- `hybrid-runs`: list of additional parameters to sample in during the `ptemcee` step, and the prior overlap cut to apply on that parameter. In the example, this is `queue.txt`.
- `hybrid-label`: custom label to differentiate the `ptemcee` analyses from the `dynesty` analysis.
See `bilby_pipe_build_hybrid --help` for more detail.

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