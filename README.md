# Geometric Structures in Quantum Algorithms for Linear Differential Equations

Reproducible **classical Python experiments** supporting the numerical one-dimensional Dirichlet-Poisson case study in the manuscript *Spectral Operator Theory and Geometric Structures in Quantum Algorithms for Linear Differential Equations* (manuscript reference **PHYSSCR-156718.R1**).

The repository contains classical discretization, filtering, fidelity, condition-ratio, input-sensitivity, and runtime diagnostics. **No quantum circuit or quantum hardware experiment is implemented, and no quantum speedup is inferred from these simulations.**

## Software and dependencies

The revision calculations were checked with:

- Python 3.13.5
- NumPy 2.3.5
- Matplotlib 3.10.8

Install the declared dependencies with:

```bash
python -m pip install -r requirements.txt
```

## Reproduction commands

From the repository directory:

```bash
python poisson_reproduce.py
python mesh_cutoff_reproduce.py
python targeted_filter_study.py
```

No downloaded dataset is required. The exact solution and forcing profiles are manufactured analytically.

## Script-to-result map

### `poisson_reproduce.py`
Reproduces the original one-dimensional Dirichlet-Poisson discretization and hard-cutoff calculations.

Outputs include:

- `poisson_results.csv`
- mesh-convergence data used for the revised convergence figure
- hard-filter fidelity and projective-angle data for the original mixed input

The revised manuscript separately reports a least-squares convergence slope and an `O(h^2)` reference line using these discretization errors.

### `mesh_cutoff_reproduce.py`
Reproduces the fixed full-retained-subspace condition-budget example with `Gamma = 3000`.

Outputs include:

- `mesh_cutoff_results.csv`
- the critical-mesh calculation
- retained-mode and hard-cutoff fidelity diagnostics

For the stated manufactured forcing, the first mesh with no nonzero admissible hard-cutoff output is `N = 946`. This is a result for the stated unpreconditioned hard-cutoff setup, not a general impossibility result for quantum linear-system algorithms.

### `targeted_filter_study.py`
Reproduces the additional **classical** comparisons requested in the R1 review.

It generates:

- `filter_comparison_parameter_study.csv`
- `classical_runtime_memory.csv`
- `hard_vs_tikhonov_tradeoff_R1.pdf`
- `input_profile_sensitivity_R1.pdf`
- `mesh_filter_sensitivity_R1.pdf`

The script compares:

- hard low-eigenvalue inverse cutoffs;
- the smooth Tikhonov-type multiplier `lambda/(lambda^2 + eta^2)`;
- low-frequency input modes `(1,2,3)`;
- the mixed input modes `(1,5,11)`;
- higher-frequency input modes `(7,11,17)`;
- meshes `N = 63, 127, 255, 511`.

It records normalized-state fidelity, Fubini-Study angle, a unit-operator-norm observable-error bound, solution-norm fraction, the relevant spectral-multiplier ratio, and reproducibility-only classical dense-eigendecomposition timing/array footprint.

## Data files committed to the repository

- `poisson_results.csv`: numerical values for the original Poisson discretization/filter study.
- `mesh_cutoff_results.csv`: selected mesh-dependent hard-cutoff values reported in the manuscript.

The additional parameter-study CSV files are generated deterministically by `targeted_filter_study.py` when the script is run.

## Interpretation and scope

The Poisson experiment uses the standard second-difference Dirichlet operator. The original mixed manufactured solution contains sine modes 1, 5, and 11. The full-retained-subspace condition ratio and the input-active condition ratio answer different questions and are reported separately in the manuscript.

The Tikhonov comparison is a classical regularization comparison requested by the reviewer. It is not presented as a new quantum algorithm or as evidence of quantum computational advantage.

Any postselection quantity in the manuscript is explicitly algebraic/conditional on an assumed exact block encoding. This repository does not provide a block-encoding construction, QSVT polynomial, quantum circuit, query-complexity estimate, or hardware result.

## Repository scope

This repository contains numerical reproducibility materials only. It does not contain the unpublished manuscript or peer-review correspondence.
