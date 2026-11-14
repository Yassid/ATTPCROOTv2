# ln(sqrt(dE/dx)*Brho) vs Brho: each species is a horizontal band (ln m ladder p 0, d .693, t 1.099, 4He 1.386 up to one constant).
import sys, glob, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
fs = sorted(f for f in glob.glob('/mnt/f/a1975/reco_d2_tau/pid_width_00[0-9][0-9].root') if not f.endswith('0016.root') and __import__('os').path.exists(f+'.done'))
cols=['va','br','sqa','sqw']; A={c:np.concatenate([uproot.open(f)['w'].arrays(cols,library='np')[c] for f in fs]) for c in cols}
v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0)
fig,ax=plt.subplots(1,2,figsize=(15,6))
for j,(lab,q) in enumerate([('peak amplitude (a1975 default)','sqa'),('trace integral (C15d option)','sqw')]):
    y=np.log(A[q][v]*A['br'][v]); a=ax[j]
    a.hist2d(A['br'][v],y,bins=[280,300],range=[[0.05,1.45],np.percentile(y,[0.5,99.7])],norm=LogNorm(),cmap='viridis')
    a.set_xlabel('Brho (Tm)'); a.set_ylabel('ln(sqrt(dE/dx)*Brho)'); a.set_title(f'{len(fs)} D2 runs: {lab}')
    a.grid(alpha=.3)
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_d2_tau/pid_lnm_vs_brho_d2.png',dpi=110); print(len(fs),'runs')
