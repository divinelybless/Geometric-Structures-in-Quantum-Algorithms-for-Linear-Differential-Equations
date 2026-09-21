"""Reproduce the 1D Dirichlet-Poisson mesh/filter experiment.

Run: python poisson_reproduce.py
Uses NumPy and Matplotlib; plots and CSV saved beside this script.
All 'success probabilities' are CONDITIONAL on an exact block encoding of
f_tau(A_h)/alpha_f, not simulated circuits or measured quantum hardware data.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
GRID=(31,63,127,255)
TAUS=(0.0,50.0,300.0)

def problem(n):
    h=1.0/(n+1)
    grid=h*np.arange(1,n+1)
    u=np.sin(np.pi*grid)+0.6*np.sin(5*np.pi*grid)+0.2*np.sin(11*np.pi*grid)
    b=np.pi**2*np.sin(np.pi*grid)+0.6*(5*np.pi)**2*np.sin(5*np.pi*grid)+0.2*(11*np.pi)**2*np.sin(11*np.pi*grid)
    idx=np.arange(1,n+1)
    lam=4/h**2*np.sin(np.pi*idx/(2*(n+1)))**2
    V=np.sqrt(2/(n+1))*np.sin(np.pi*np.outer(np.arange(1,n+1),idx)/(n+1))
    beta=V.T@b
    x=V@(beta/lam)
    # This is the standard, symmetric second-difference Dirichlet matrix.
    A=(np.diag(2*np.ones(n))-np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/h**2
    assert np.allclose(A@x,b,atol=1e-8,rtol=1e-10)
    assert np.allclose(x, np.linalg.solve(A,b),atol=1e-8,rtol=1e-10)
    return h,grid,u,b,lam,V,beta,x

def run():
    rows=[]
    for n in GRID:
        h,grid,u,b,lam,V,beta,x=problem(n)
        err=np.sqrt(h)*np.linalg.norm(x-u)
        rel=err/(np.sqrt(h)*np.linalg.norm(u))
        full_fidelity=abs(np.vdot(x,u))**2/(np.vdot(x,x).real*np.vdot(u,u).real)
        for tau in TAUS:
            retained=lam>=tau
            assert retained.any()
            coeff=np.where(retained,beta/lam,0)
            xt=V@coeff
            F=abs(np.vdot(x,xt))**2/(np.vdot(x,x).real*np.vdot(xt,xt).real)
            F_spectral=np.sum(abs(coeff)**2)/np.sum(abs(beta/lam)**2)
            assert np.isclose(F,F_spectral,rtol=1e-11,atol=1e-11)
            angle=np.arccos(np.sqrt(np.clip(F,0,1)))
            alpha_f=1/np.min(lam[retained]); bnorm=np.linalg.norm(b)
            p=np.linalg.norm(xt)**2/(alpha_f**2*bnorm**2)
            assert 0<=p<=1+1e-12
            row=dict(n=n,h=h,tau=tau,lambda_min=float(lam[0]),lambda_max=float(lam[-1]),
                kappa_full=float(lam[-1]/lam[0]),lambda_min_retained=float(np.min(lam[retained])),
                kappa_retained=float(lam[-1]/np.min(lam[retained])),retained_dimension=int(retained.sum()),
                discretization_error_L2h=float(err),discretization_relative_L2h=float(rel),
                exact_vs_discrete_fidelity=float(full_fidelity),filtered_error_to_continuum_L2h=float(np.sqrt(h)*np.linalg.norm(xt-u)),
                filtered_vs_continuum_fidelity=float(abs(np.vdot(xt,u))**2/(np.vdot(xt,xt).real*np.vdot(u,u).real)),filtered_fidelity=float(F),
                fs_angle_radians=float(angle),postselection_probability_ideal_exact_block=float(p))
            rows.append(row)
        print(f'n={n:3d}, h={h:.7f}, error_L2h={err:.8f}, relative_error={rel:.7f}, kappa={lam[-1]/lam[0]:.3f}, fidelity_to_u={full_fidelity:.9f}')
    with (OUT/'poisson_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    unfiltered=[r for r in rows if r['tau']==0]
    orders=[np.log(unfiltered[i]['discretization_error_L2h']/unfiltered[i+1]['discretization_error_L2h'])/np.log(unfiltered[i]['h']/unfiltered[i+1]['h']) for i in range(len(unfiltered)-1)]
    print('observed convergence orders:', ', '.join(f'{p:.4f}' for p in orders))
    fig,ax=plt.subplots(figsize=(6.2,4.0))
    ax.loglog([r['h'] for r in unfiltered],[r['discretization_error_L2h'] for r in unfiltered],'o-',label=r'$\|x_h-u|_{grid}\|_{2,h}$')
    ax.loglog([r['h'] for r in unfiltered],[unfiltered[-1]['discretization_error_L2h']*(r['h']/unfiltered[-1]['h'])**2 for r in unfiltered],'--',label=r'$O(h^2)$ reference')
    ax.set(xlabel='Mesh spacing h',ylabel=r'Weighted discrete $L^2$ error',title='Dirichlet Poisson: mesh convergence')
    ax.grid(True,which='both',alpha=.3);ax.legend();fig.tight_layout();fig.savefig(OUT/'poisson_mesh_convergence.pdf');fig.savefig(OUT/'poisson_mesh_convergence.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.2,4.0))
    for tau in TAUS:
        r=[r for r in rows if r['tau']==tau]
        ax.plot([e['n'] for e in r],[e['filtered_fidelity'] for e in r],'o-',label=rf'$\tau={tau:g}$')
    ax.set(xlabel='Interior grid points N',ylabel='Fidelity with unfiltered discrete inverse',ylim=(0,1.04),title='Input-dependent cost of hard spectral filtering')
    ax.grid(True,alpha=.3);ax.legend();fig.tight_layout();fig.savefig(OUT/'poisson_filter_fidelity.pdf');fig.savefig(OUT/'poisson_filter_fidelity.png',dpi=170);plt.close(fig)
    print('N | tau | effective kappa | fidelity | FS angle | ideal p')
    for r in rows:
        print(f"{r['n']:3d} | {r['tau']:3.0f} | {r['kappa_retained']:11.2f} | {r['filtered_fidelity']:.6f} | {r['fs_angle_radians']:.6f} | {r['postselection_probability_ideal_exact_block']:.6f}")
    print('PASS: discrete linear solve, spectral fidelity, ideal probability, and convergence calculations')
    return rows,orders

if __name__=='__main__': run()
