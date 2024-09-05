# A Hybrid solver for lasso and basis-pursuit denoise

This code accompanies the paper: E. van den Berg, "A hybrid quasi-Newton projected-gradient method with application to Lasso and basis-pursuit denoising," Math. Prog. Comp. 12, 1–38 (2020); see the [arXiv preprint](https://arxiv.org/abs/1611.05483) or the [published journal paper](https://doi.org/10.1007/s12532-019-00163-5).

#### Getting started

The Matlab scripts depend on several MEX files for speed. The source code for the extensions are located in the `private` directory and can be compiled using the `compile_projector_mex.m` and `compile_product_B_mex.m` scripts (this may require installation/configuation of supported compilers for Matlab).
Some of the scripts may require installation of [SPGL1](https://friedlander.io/spgl1/) and the [SPARCO](https://github.com/MPF-Optimization-Laboratory/Sparco) test problem set for sparse reconstruction.

#### The solver

The hybrid solver builds on [SPGL1](https://friedlander.io/spgl1/). Incremental versions of the solver were developed and used at different stages of the project. The latest version of the solver can be found in `solver_v05.m`. An earlier version `solver_v04.m` is also provided to generate results some of the figures and tables.


#### Reproducing figures and tables from the paper

There are several scripts that reproduce tables and figures from the paper. When available, these scripts use precomputed results that are stored in the `cache` directory. In order to regenerate the results, delete the appropriate files in the cache directory. The following table summarizes the relevant scripts and the associated tables and figures

| Script                    | generates          |
| ------------------------- | ------------------ |
| `plot_Pareto.m`           | Figure 1           |
| `plot_mutual_coherence.m` | Figure 5a          |
| `plot_Lipschitz.m`        | Figures 5a-c       |
| `plot_problem_duality.m`  | Figures 7a,b       |
| `plot_duality_gap.m`      | Figure 7c,d        |
| `table_sparse_x0.m`       | Table 1            |
| `table_correlated.m`      | Tables 2 and 3     |
| `table_sparco.m`          | Tables 4, 5, and 6 |
