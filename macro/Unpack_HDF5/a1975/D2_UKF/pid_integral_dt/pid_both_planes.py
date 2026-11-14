import glob, os, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys; DIR=sys.argv[1]; TAG=sys.argv[2]
fs=sorted(f for f in glob.glob(DIR+'/pid_width_0[0-9][0-9][0-9].root') if not f.endswith('0016.root') and os.path.exists(f+'.done'))
c=['va','br','sqa','sqw']; A={k:np.concatenate([uproot.open(f)['w'].arrays(c,library='np')[k] for f in fs]) for k in c}
v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0); br=A['br'][v]
fig,ax=plt.subplots(2,2,figsize=(16,12))
for j,(lab,q) in enumerate([('peak amplitude (a1975 default)','sqa'),('trace integral (C15d option)','sqw')]):
    y=A[q][v]
    ax[0,j].hist2d(y,br,bins=[400,400],range=[[0,np.percentile(y,99.5)],[0,2.5]],norm=LogNorm(),cmap='viridis')
    ax[0,j].set_xlabel('sqrt(dE/dx)'); ax[0,j].set_ylabel('Brho (Tm)'); ax[0,j].set_title(f'{len(fs)} {TAG} runs, {v.sum()} tracks: {lab}')
    l=np.log(y*br); ax[1,j].hist2d(br,l,bins=[400,400],range=[[0.05,2.5],[0.5,3.6]],norm=LogNorm(),cmap='viridis')
    ax[1,j].set_xlabel('Brho (Tm)'); ax[1,j].set_ylabel('ln(sqrt(dE/dx)*Brho)'); ax[1,j].grid(alpha=.3)
plt.tight_layout(); plt.savefig(DIR+f'/pid_both_planes_{TAG}.png',dpi=90)
