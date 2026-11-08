#!/usr/bin/env python3
"""Targeted classical filter/input-sensitivity study requested in PHYSSCR-156718.R1.

No quantum circuit or quantum speedup is simulated or inferred.
The script compares a hard low-eigenvalue inverse cutoff with the smooth
Tikhonov-type inverse multiplier lambda/(lambda^2+eta^2) for manufactured
1D Dirichlet-Poisson inputs, and records classical dense-eigensolver costs.
"""
from pathlib import Path
import csv, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
MIXED=((1,1.0),(5,0.6),(11,0.2))
PROFILES={
    "low": ((1,1.0),(2,0.35),(3,0.15)),
    "mixed": MIXED,
    "high": ((7,1.0),(11,0.6),(17,0.2)),
}
MESHES=(63,127,255,511)
HARD_PARAMS=(0.0,25.0,50.0,100.0,200.0,300.0,500.0,1000.0)
TIK_PARAMS=(1.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0)

def spectrum(n):
    h=1/(n+1); j=np.arange(1,n+1)
    return h,4/h**2*np.sin(np.pi*j/(2*(n+1)))**2

def basis(n):
    i=np.arange(1,n+1); j=np.arange(1,n+1)
    return np.sqrt(2/(n+1))*np.sin(np.pi*np.outer(i,j)/(n+1))

def problem(n,modes):
    h,lam=spectrum(n); grid=h*np.arange(1,n+1)
    u=sum(a*np.sin(m*np.pi*grid) for m,a in modes)
    b=sum(a*(m*np.pi)**2*np.sin(m*np.pi*grid) for m,a in modes)
    V=basis(n); beta=V.T@b; xcoef=beta/lam; x=V@xcoef
    return h,grid,u,b,lam,V,beta,xcoef,x

def fidelity(x,y):
    return float(abs(np.vdot(x,y))**2/(np.vdot(x,x).real*np.vdot(y,y).real))

def angle(F):
    return float(np.arccos(np.sqrt(np.clip(F,0,1))))

def write_csv(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    rows=[]
    for name,modes in PROFILES.items():
        for n in MESHES:
            h,grid,u,b,lam,V,beta,xcoef,x=problem(n,modes)
            for tau in HARD_PARAMS:
                keep=lam>=tau; coeff=np.where(keep,xcoef,0); y=V@coeff
                if np.linalg.norm(y)==0:
                    continue
                F=fidelity(x,y); frac=float(np.sum(abs(coeff)**2)/np.sum(abs(xcoef)**2))
                rows.append(dict(profile=name,N=n,filter="hard_cutoff",parameter=tau,
                    fidelity=F,fs_angle_rad=angle(F),
                    observable_error_bound_unit_norm=2*np.sqrt(max(0,1-F)),
                    solution_norm_fraction=frac,retained_solution_weight=frac,
                    spectral_multiplier_ratio=float(lam[-1]/np.min(lam[keep]))))
            for eta in TIK_PARAMS:
                filt=lam/(lam**2+eta**2); coeff=filt*beta; y=V@coeff
                F=fidelity(x,y); positive=filt[filt>0]
                rows.append(dict(profile=name,N=n,filter="tikhonov",parameter=eta,
                    fidelity=F,fs_angle_rad=angle(F),
                    observable_error_bound_unit_norm=2*np.sqrt(max(0,1-F)),
                    solution_norm_fraction=float(np.sum(abs(coeff)**2)/np.sum(abs(xcoef)**2)),
                    retained_solution_weight="",
                    spectral_multiplier_ratio=float(np.max(positive)/np.min(positive))))
    write_csv(OUT/"filter_comparison_parameter_study.csv",rows)

    fig,ax=plt.subplots(figsize=(6.7,4.4))
    for filt,marker in [("hard_cutoff","o"),("tikhonov","s")]:
        rr=[r for r in rows if r["profile"]=="mixed" and r["N"]==255 and r["filter"]==filt]
        ax.semilogx([r["spectral_multiplier_ratio"] for r in rr],[r["fidelity"] for r in rr],
                    marker+"-",label=filt.replace("_"," "))
    ax.set(xlabel="Spectral-multiplier ratio on the stated support",
           ylabel="Fidelity with unregularized discrete solution",
           ylim=(-.02,1.04),
           title="Hard cutoff and Tikhonov regularization: same Poisson input")
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(OUT/"hard_vs_tikhonov_tradeoff_R1.pdf");plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.7,4.4))
    for name,marker in [("low","o"),("mixed","s"),("high","^")]:
        rr=[r for r in rows if r["profile"]==name and r["N"]==255 and r["filter"]=="hard_cutoff"]
        ax.plot([r["parameter"] for r in rr],[r["fidelity"] for r in rr],marker+"-",label=name)
    ax.set(xlabel=r"Hard cutoff $\tau$",ylabel="Fidelity with unfiltered discrete solution",
           ylim=(-.02,1.04),title="Sensitivity to input spectral composition ($N=255$)")
    ax.grid(alpha=.25);ax.legend(title="input profile");fig.tight_layout()
    fig.savefig(OUT/"input_profile_sensitivity_R1.pdf");plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.7,4.4))
    for filt,par,marker in [("hard_cutoff",50.0,"o"),("tikhonov",50.0,"s")]:
        rr=[r for r in rows if r["profile"]=="mixed" and r["filter"]==filt and r["parameter"]==par]
        ax.plot([r["N"] for r in rr],[r["fidelity"] for r in rr],marker+"-",
                label=f"{filt.replace('_',' ')}; parameter=50")
    ax.set(xlabel="Interior grid points $N$",ylabel="Fidelity with unregularized discrete solution",
           ylim=(-.02,1.04),title="Mesh sensitivity for fixed regularization parameters")
    ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(OUT/"mesh_filter_sensitivity_R1.pdf");plt.close(fig)

    bench=[]
    for n in MESHES:
        h=1/(n+1)
        A=(np.diag(2*np.ones(n))-np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/h**2
        np.linalg.eigh(A); ts=[]
        for _ in range(3):
            t0=time.perf_counter(); vals,vecs=np.linalg.eigh(A); ts.append(time.perf_counter()-t0)
        bench.append(dict(N=n,dense_eigh_median_seconds=float(np.median(ts)),
            dense_eigh_min_seconds=float(np.min(ts)),
            array_footprint_MiB=float((A.nbytes+vecs.nbytes+vals.nbytes)/(1024**2)),
            analytic_spectrum_storage_MiB=float(vals.nbytes/(1024**2))))
    write_csv(OUT/"classical_runtime_memory.csv",bench)
    print("Wrote",len(rows),"filter-comparison rows and",len(bench),"runtime rows")

if __name__=="__main__":
    main()
