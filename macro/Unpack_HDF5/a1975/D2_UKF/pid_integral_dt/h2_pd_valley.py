# H2: p/d separation where the (p,d) deuterons are, both planes. Species located per Brho slice from the bands (p = dominant low band, d searched at p+0.55..0.85).
import glob, os, uproot, numpy as np
from scipy.ndimage import uniform_filter1d
fs=sorted(f for f in glob.glob('/mnt/f/a1975/reco_h2_tau/pid_width_0[0-9][0-9][0-9].root') if os.path.exists(f+'.done'))
c=['va','br','sqa','sqw']; A={k:np.concatenate([uproot.open(f)['w'].arrays(c,library='np')[k] for f in fs]) for k in c}
v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0); br=A['br'][v]
J=np.load('/mnt/f/a1975/sim_rebuild/trig/pd/joined_pd.npz'); ke=J['T'][:,2]; bpd=np.sqrt(ke**2+2*1875.613*ke)/299.792458
print(len(fs),'runs; adopted (p,d) deuteron Brho pct 5/50/95:',np.percentile(bpd,[5,50,95]).round(2))
print(f"{'Brho':9s} {'plane':6s} {'p pos':>6s} {'d-p':>6s} {'valley/min peak':>15s} {'N':>6s}")
for lo,hi in [(0.4,0.55),(0.55,0.7),(0.7,0.85),(0.85,1.0),(1.0,1.2)]:
    for q in ['sqa','sqw']:
        y=np.log(A[q][v]*br)[(br>=lo)&(br<hi)]; h,e=np.histogram(y,240,range=(0.6,3.0)); xc=.5*(e[1:]+e[:-1]); hs=uniform_filter1d(h.astype(float),7)
        w=(xc>1.0)&(xc<1.7); mp=xc[w][np.argmax(hs[w])]; w2=(xc>mp+.55)&(xc<mp+.9); md=xc[w2][np.argmax(hs[w2])]
        wv=(xc>mp)&(xc<md); val=hs[wv].min()/min(hs[np.abs(xc-mp)<.02].max(),hs[np.abs(xc-md)<.02].max())
        print(f"{lo:.2f}-{hi:.2f} {q:6s} {mp:6.3f} {md-mp:+6.3f} {val:15.2f} {len(y):6d}")
