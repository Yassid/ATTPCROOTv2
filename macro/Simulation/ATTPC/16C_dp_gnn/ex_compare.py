#!/usr/bin/env python
"""Compare the 17C excitation-energy spectrum with vs without noise cleaning.
Input: ex_A.csv, ex_B.csv (from ex_extract: event,ke,theta,chi2ndf,vz).
Computes Ex via two-body kinematics, applies proton cuts, compares yield + peak.
  ~/gnn_env/bin/python ex_compare.py <ex_A.csv> <ex_B.csv> [chi2cut]
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

u = 931.49401
mC16, md, mp, mC17 = 16.0147013*u, 2.01410178*u, 1.00782503*u, 17.0225864*u
Ebeam = 192.0
CHI2 = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0

def omega2(x, y, z): return np.sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z)
def kine_ex(Kp, thlab, Kej, m1=mC16, m2=md, m3=mp, m4=mC17):
    Et1 = Kp+m1; Et3 = Kej+m3
    s = m1*m1+m2*m2+2*m2*Et1; uu = m2*m2+m3*m3-2*m2*Et3
    inside = (np.cos(thlab)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)
              - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2
    m4ex = np.sqrt(np.clip(inside, 0, None))
    return m4ex - m4

def load(f):
    d = pd.read_csv(f); thd = np.degrees(d.theta)
    sel = (d.ke > 0) & (d.ke < 50) & (thd > 10) & (thd < 170) & (d.chi2ndf < CHI2)
    d = d[sel].copy()
    d["ex"] = np.asarray(kine_ex(Ebeam, d.theta.to_numpy(dtype=float), d.ke.to_numpy(dtype=float)), dtype=float)
    d = d[np.isfinite(d.ex.to_numpy(dtype=float))]
    return d

A, B = load(sys.argv[1]), load(sys.argv[2])
def peak(d, lo=-3, hi=3):
    x = d.ex[(d.ex > lo) & (d.ex < hi)]
    return len(x), (x.mean() if len(x) else np.nan), (x.std() if len(x) else np.nan)
nA, mA, sA = peak(A); nB, mB, sB = peak(B)
print(f"chi2/ndf < {CHI2}")
print(f"{'':<12}{'proton cand':>12}{'gs-peak N':>11}{'gs mean':>10}{'gs sigma':>10}")
print(f"{'NO-CLEAN':<12}{len(A):>12}{nA:>11}{mA:>10.2f}{sA:>10.2f}")
print(f"{'CLEAN':<12}{len(B):>12}{nB:>11}{mB:>10.2f}{sB:>10.2f}")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
bins = np.linspace(-8, 15, 92)
ax[0].hist(A.ex, bins=bins, histtype='step', lw=2, label=f'no-clean (n={len(A)})', color='gray')
ax[0].hist(B.ex, bins=bins, histtype='step', lw=2, label=f'clean (n={len(B)})', color='tab:blue')
ax[0].axvline(0, ls=':', c='r'); ax[0].set_xlabel('Ex(17C) [MeV]'); ax[0].set_ylabel('proton candidates')
ax[0].set_title(f'17C excitation energy (chi2/ndf<{CHI2})'); ax[0].legend()
ax[1].scatter(np.degrees(A.theta), A.ke, s=4, c='gray', label='no-clean', alpha=.5)
ax[1].scatter(np.degrees(B.theta), B.ke, s=4, c='tab:blue', label='clean', alpha=.5)
ax[1].set_xlabel('theta_lab [deg]'); ax[1].set_ylabel('proton KE [MeV]'); ax[1].set_title('kinematics'); ax[1].legend()
plt.tight_layout(); plt.savefig("diagnostics/ex_ab.png", dpi=95); print("wrote diagnostics/ex_ab.png")
