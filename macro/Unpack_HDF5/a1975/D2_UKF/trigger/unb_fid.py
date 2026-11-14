import numpy as np, io, contextlib, runpy, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()): G=runpy.run_path('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py')
d,dx,dc,dr=G['d'],G['dx'],G['dc'],G['dr']
base=(d[:,3]<=13.5)&(d[:,1]>=89)&(d[:,1]<=161)&(d[:,0]>=1.2)&np.isfinite(dx)&(d[:,2]>=200)&(d[:,2]<900)
bins=np.arange(-3,8.01,0.2); fig,ax=plt.subplots(1,3,figsize=(18,5))
for i,(a,b) in enumerate([(20,30),(30,40),(40,50)]):
    for R,c in [(0,'0.6'),(20,'C0'),(28,'C3')]:
        w=base&(dc>=a)&(dc<b)&(dr>=R); h,_=np.histogram(dx[w],bins); ax[i].stairs(h,bins,color=c,label=f'r_res>={R}: N={w.sum()}')
        if R==28: print(f'theta_cm {a}-{b} R28: bound(-1.5,.73) {np.sum(w&(dx>-1.5)&(dx<.73))}  2-3.2 {np.sum(w&(dx>2)&(dx<3.2))}  3.2-5 {np.sum(w&(dx>3.2)&(dx<5))}  5-8 {np.sum(w&(dx>5)&(dx<8))}')
    ax[i].set_title(f'theta_cm {a}-{b} (g.s. kinematics), vz 200-900'); ax[i].set_xlabel('Ex(17C) [MeV]'); ax[i].legend(fontsize=8); [ax[i].axvline(z,ls=':',color='k',alpha=.4) for z in (0,2.67,3.8)]
plt.tight_layout(); plt.savefig('/mnt/f/a1975/sim_rebuild/trig/ex_spectra_fiducial.png',dpi=90)
