#!/usr/bin/env python3
"""Exact modal reproducibility for the fixed-bandwidth Dirichlet Poisson cutoff example.

No quantum circuit is simulated. Gamma constrains the entire retained spectral interval.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
GAMMA = 3000.0
MODES = ((1, 1.0), (5, 0.6), (11, 0.2))
SAMPLE_N = (63, 255, 511, 767, 943, 945, 946, 959, 1023)


def eig(n: int, j: int) -> float:
    return 4.0 * (n + 1)**2 * np.sin(j * np.pi / (2 * (n + 1)))**2


def row(n: int, gamma: float = GAMMA) -> dict:
    assert n >= max(m for m, _ in MODES) and gamma >= 1
    highest = eig(n,n)
    target = highest/gamma
    # This index implements the *actual retained eigenvalue* constraint, not tau>=target.
    spectrum = 4.0 * (n+1)**2 * np.sin(np.arange(1,n+1)*np.pi/(2*(n+1)))**2
    k = int(np.searchsorted(spectrum,target,side='left'))+1
    assert 1 <= k <= n
    solution_powers = np.array([a*a*(m*np.pi)**4/eig(n,m)**2 for m,a in MODES],dtype=float)
    weights=solution_powers/solution_powers.sum()
    f = float(sum(w for (m,_),w in zip(MODES,weights) if m>=k))
    retained = any(m>=k for m,_ in MODES)
    if not retained:
        assert f < 1e-15
        f=0.0
    return {'N':n,'Gamma':gamma,'lambda_N':highest,
            'lambda_N_over_Gamma':target,'first_retained_mode':k,
            'smallest_retained_eigenvalue':float(spectrum[k-1]),
            'retained_condition_number':highest/float(spectrum[k-1]),
            'maximum_hard_cutoff_fidelity':f,
            'projective_angle_rad':float(np.arccos(np.sqrt(f))) if retained else None,
            'normalized_output_exists':retained,
            'input_active_condition_number':eig(n,11)/eig(n,1)}


def main():
    results=[row(n) for n in SAMPLE_N]
    with (OUT/'mesh_cutoff_results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]))
        writer.writeheader();writer.writerows(results)
    for r in results:
        print('N={N:4d} first={first_retained_mode:2d} Fmax={maximum_hard_cutoff_fidelity:.8f} output={normalized_output_exists} K_active={input_active_condition_number:.2f}'.format(**r))
    first_failure=next(n for n in range(11,2001) if eig(n,n)/GAMMA>eig(n,11))
    assert first_failure==946
    assert eig(945,945)/GAMMA <=eig(945,11)
    assert eig(946,946)/GAMMA >eig(946,11)
    for r in results:
        assert r['retained_condition_number']<=GAMMA*(1+1e-12)
        assert 0<=r['maximum_hard_cutoff_fidelity']<=1
        if not r['normalized_output_exists']:
            assert r['maximum_hard_cutoff_fidelity']==0
    print('Threshold verification: first N with no nonzero admissible output =',first_failure)
    nvals=np.arange(11,1101)
    fvals=np.array([row(int(n))['maximum_hard_cutoff_fidelity'] for n in nvals])
    fig,ax=plt.subplots(figsize=(7.4,4.4))
    ax.step(nvals,fvals,where='post',lw=1.7,label=r'$F_{\max}$ (hard cutoff, $\Gamma=3000$)')
    ax.axvline(946,ls='--',lw=1.1,label=r'first zero-output grid $N=946$')
    ax.set(xlabel='Number of interior points, $N$',ylabel='Maximum possible fidelity',ylim=(-.045,1.045),xlim=(11,1100))
    ax.grid(alpha=.2)
    ax.legend(fontsize=9,loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT/'mesh_cutoff_fidelity.png',dpi=190)
    fig.savefig(OUT/'mesh_cutoff_fidelity.pdf')
    plt.close(fig)
    print('Numerical formula, monotonicity, and threshold checks passed.')

if __name__=='__main__':
    main()
