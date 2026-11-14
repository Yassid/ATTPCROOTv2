import glob, os, json, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.path import Path; from matplotlib.colors import LogNorm
fs=sorted(f for f in glob.glob('/mnt/f/a1975/reco_d2_tau/pid_width_00[0-9][0-9].root') if not f.endswith('0016.root') and os.path.exists(f+'.done'))
c=['va','br','sqa','sqw']; A={k:np.concatenate([uproot.open(f)['w'].arrays(c,library='np')[k] for f in fs]) for k in c}
v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0); br,sa,sw=A['br'][v],A['sqa'][v],A['sqw'][v]
G=np.array(json.load(open(os.path.expanduser('~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/triton_d2_dv1104.json')))['vertices'])
ing=Path(G).contains_points(np.c_[sa,br]); print('tracks',v.sum(),'in gate',ing.sum())
ya,yw=np.log(sa*br),np.log(sw*br)
fig,ax=plt.subplots(2,3,figsize=(18,10))
a=ax[0,0]; a.hist2d(br,sa,bins=[300,300],range=[[0.05,2.5],[0,40]],norm=LogNorm()); a.plot(G[:,1],G[:,0],'r-'); a.set_xlabel('Brho'); a.set_ylabel('sqrt dE/dx (peak)'); a.set_title('adopted triton gate on peak plane (test reco)')
a=ax[0,1]; a.hist2d(br,ya,bins=[300,300],range=[[0.05,2.5],[0.5,3.5]],norm=LogNorm()); a.set_title('peak plane ln(sqrt*Brho), gate in red')
a.scatter(br[ing],ya[ing],s=0.3,c='r',alpha=.3)
a=ax[0,2]; a.hist2d(br,yw,bins=[300,300],range=[[0.05,2.5],[0.5,3.5]],norm=LogNorm()); a.scatter(br[ing],yw[ing],s=0.3,c='r',alpha=.3); a.set_title('integral plane, gated tracks in red')
# d anchor per plane at Brho .3-.5 (as in pid_integral_dt_d2.py)
def mode(y,c,w):
    h,e=np.histogram(y,200,range=(c-w,c+w)); from scipy.ndimage import uniform_filter1d; return .5*(e[1:]+e[:-1])[np.argmax(uniform_filter1d(h.astype(float),7))]
s=(br>.3)&(br<.5); dw=mode(yw[s],np.median(yw[s]),.5); print('d anchor integral',round(dw,3))
W=[(0.3,0.6),(0.6,0.9),(0.9,1.2),(1.2,1.5),(1.5,2.0),(2.0,2.5)]
print(f"{'Brho':9s} {'gated':>6s} {'d-like':>7s} {'t-like':>7s} {'a-like':>7s}   (integral y-d: d <0.2, t 0.2-0.6, a >0.6)")
for i,(lo,hi) in enumerate(W):
    w=ing&(br>=lo)&(br<hi); x=yw[w]-dw; n=w.sum()
    fd,ft,fa=np.mean(x<0.2),np.mean((x>=0.2)&(x<0.6)),np.mean(x>=0.6)
    print(f"{lo:.1f}-{hi:.1f}  {n:6d} {fd*100:6.1f}% {ft*100:6.1f}% {fa*100:6.1f}%")
    if i<3: b=ax[1,i]
    else: b=ax[1,i-3]
    b.hist(yw[(br>=lo)&(br<hi)]-dw,200,range=(-0.8,1.2),histtype='step',color='k',label=f'all {lo}-{hi}' if True else '')
    b.hist(x,200,range=(-0.8,1.2),histtype='step',lw=1.4,label=f'in triton gate {lo}-{hi}')
for b in ax[1]: b.set_yscale('log'); b.legend(fontsize=7); b.set_xlabel('integral: ln(sqrt*Brho) - d   (t +0.41, 4He +0.75)'); [b.axvline(z,ls=':') for z in (0.2,0.6)]
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_d2_tau/dt_gate_leak_integral.png',dpi=90)
