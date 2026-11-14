import glob, os, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
def load(D):
    fs=sorted(f for f in glob.glob(D+'/pid_width_0[0-9][0-9][0-9].root') if not f.endswith('0016.root') and os.path.exists(f+'.done'))
    c=['va','br','sqa','sqw']; A={k:np.concatenate([uproot.open(f)['w'].arrays(c,library='np')[k] for f in fs]) for k in c}
    v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0); return len(fs),{k:A[k][v] for k in c}
fig,ax=plt.subplots(2,2,figsize=(16,11))
for i,(lab,D) in enumerate([('D2 (deuterium target)','/mnt/f/a1975/reco_d2_tau'),('H2 (proton target)','/mnt/f/a1975/reco_h2_tau')]):
    n,A=load(D)
    for j,(pl,q) in enumerate([('peak amplitude','sqa'),('trace integral','sqw')]):
        y=np.log(A[q]*A['br']); a=ax[i,j]; a.hist2d(A['br'],y,bins=[300,300],range=[[0.05,2.0],[0.5,3.6]],norm=LogNorm(),cmap='viridis')
        a.set_title(f'{lab}, {n} runs, {len(y)} tracks: {pl}'); a.set_xlabel('Brho (Tm)'); a.set_ylabel('ln(sqrt(dE/dx)*Brho)'); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_h2_tau/pid_ln_D2_vs_H2.png',dpi=90)
