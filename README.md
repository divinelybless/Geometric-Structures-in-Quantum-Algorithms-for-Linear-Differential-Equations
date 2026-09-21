# Geometric Structures in Quantum Algorithms for Linear Differential Equations

Reproducible **classical Python experiments** supporting the numerical one-dimensional Dirichlet–Poisson case study in the manuscript *Spectral Operator Theory and Geometric Structures in Quantum Algorithms for Linear Differential Equations* (manuscript reference **PHYSSCR-156718**).

The scripts examine finite-difference discretization, hard spectral filtering, the retained-subspace condition number, normalized-state fidelity, and a mesh-dependent hard-cutoff example. **No quantum circuit or quantum hardware experiment is implemented.** Any reported postselection probability is an idealized quantity conditional on an exact block encoding, not a measured hardware result.

## Files

- `poisson_reproduce.py`: reproduces the 1D Poisson discretization/filtering experiment, `poisson_results.csv`, and the `poisson_mesh_convergence` and `poisson_filter_fidelity` figures (PDF and PNG).
- `mesh_cutoff_reproduce.py`: reproduces the fixed retained-subspace condition-number budget example (Γ = 3000), `mesh_cutoff_results.csv`, and the `mesh_cutoff_fidelity` figure (PDF and PNG).
- `poisson_results.csv` and `mesh_cutoff_results.csv`: the numerical data accompanying the manuscript. The scripts overwrite these files when run.
- `requirements.txt`: Python dependencies.

## Reproduce the results

Use Python 3.10 or newer. From the repository directory:

```bash
python -m pip install -r requirements.txt
python poisson_reproduce.py
python mesh_cutoff_reproduce.py
```

The scripts save the CSV files and figures in their own directory. No downloaded data are needed: the forcing and exact solution are manufactured analytically. A verified run used Python 3.13.5, NumPy 2.3.5, and Matplotlib 3.10.8. The reproducibility scripts regenerated all 12 Poisson rows and 9 mesh-cutoff rows identically in that environment.

## Scope and interpretation

The Poisson experiment uses the standard second-difference Dirichlet operator and an exact solution with sine modes 1, 5, and 11. The full-retained-subspace condition-number constraint differs from a condition number restricted to the input-active modes. For the fixed Γ = 3000 example, the first grid size with no nonzero output after the admissible hard cutoff is N = 946; this is **not** a general impossibility result for quantum linear-system algorithms.

This repository contains only numerical reproducibility materials, not the unpublished manuscript, peer-review correspondence, or a quantum algorithm implementation. Please refer to the manuscript for the analytical assumptions and interpretation. No license is specified here; public visibility alone does not grant reuse rights.
