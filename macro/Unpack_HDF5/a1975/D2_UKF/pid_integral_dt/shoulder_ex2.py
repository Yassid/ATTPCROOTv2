import glob, os, json, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.path import Path; from scipy.optimize import curve_fit; from scipy.ndimage import uniform_filter1d
fs=sorted(f for f in glob.glob('/mnt/f/a1975/reco_d2_tau/pid_width_00[0-9][0-9].root') if not f.endswith('0016.root') and os.path.exists(f+'.done'))
P=[]
for f in fs:
    a=uproot.open(f)['w'].arrays(['ev','va','br','sqa','sqw'],library='np'); r=int(os.path.basename(f)[10:14])
    P.append(np.c_[np.full(len(a['ev']),r),a['ev'],a['va'],a['br'],a['sqa'],a['sqw']])
P=np.vstack(P); P=P[(P[:,2]>0)&(P[:,3]>0.05)&(P[:,4]>0)&(P[:,5]>0)]
run,ev,br,sa,sw=P[:,0].astype(int),P[:,1].astype(int),P[:,3],P[:,4],P[:,5]; yw=np.log(sw*br)
G=np.array(json.load(open(os.path.expanduser('~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/triton_d2_dv1104.json')))['vertices'])
ing=Path(G).contains_points(np.c_[sa,br])
# per 0.1-Brho slice triton position from the all-track 3-gaussian fit (as gate_leak2.py)
def g(x,a,m,s): return a*np.exp(-.5*((x-m)/s)**2)
def g3(x,*p): return g(x,*p[0:3])+g(x,*p[3:6])+g(x,*p[6:9])+p[9]
edges=np.round(np.arange(0.3,2.01,0.1),2); MT={}
for lo,hi in zip(edges[:-1],edges[1:]):
    y=yw[(br>=lo)&(br<hi)]; h,e=np.histogram(y,160,range=(1.3,3.3)); xc=.5*(e[1:]+e[:-1]); hs=uniform_filter1d(h.astype(float),5)
    w=(xc>1.7)&(xc<2.25); md=xc[w][np.argmax(hs[w])]
    p0=[hs[w].max(),md,.07,max(hs[np.abs(xc-md-.42)<.03].max(),1),md+.42,.07,max(hs[np.abs(xc-md-.75)<.05].max(),1),md+.75,.07,1]
    try: p,_=curve_fit(g3,xc,h,p0=p0,bounds=([0,md-.05,.02,0,md+.30,.03,0,md+.6,.03,0],[np.inf,md+.05,.3,np.inf,md+.52,.15,np.inf,md+.95,.15,np.inf]),sigma=np.sqrt(np.maximum(h,1)),maxfev=40000)
    except Exception: p=p0
    MT[lo]=(p[4],p[5] if p[5]<0.149 else 0.09)
# adopted (d,t) tritons
d=np.load('/mnt/f/a1975/sim_rebuild/trig/dt/joined_dt.npz'); T=d['T']; s=d['sel'].astype(bool)
R,E,KE,EX,TCM=T[s,0].astype(int),T[s,1].astype(int),T[s,3],d['ex'][s],d['tcm'][s]
bgf=np.sqrt(KE**2+2*2808.921*KE)/299.792458
key={}
for i,(r,e) in enumerate(zip(run,ev)): key.setdefault((r,e),[]).append(i)
cls=np.full(len(R),-9.0); nm=0
for j in range(len(R)):
    c=key.get((R[j],E[j]));
    if not c: continue
    c=np.array(c); k=c[np.argmin(np.abs(br[c]-bgf[j]))]
    if abs(br[k]-bgf[j])/bgf[j]>0.15: continue
    nm+=1; lo=np.floor(br[k]*10)/10
    if lo not in MT: cls[j]=-8; continue
    mt,st=MT[lo]; cls[j]=(yw[k]-mt)/st          # pull from the triton peak, integral plane
inrun=np.isin(R,np.unique(run)); print('adopted tritons in these runs',inrun.sum(),'matched',nm, ' (with slice) ', (cls>-8).sum())
m=cls>-8
core=m&(cls<1.5)&(cls>=-1.5); sh=m&(cls>=1.5)&(cls<4.5); far=m&(cls>=4.5); low=m&(cls<-1.5)
print('pull pct 5/25/50/75/95',np.percentile(cls[m],[5,25,50,75,95]).round(2)); print('pull classes: core |pull|<1.5 & low', core.sum(), ' shoulder 1.5-4.5', sh.sum(), ' >4.5', far.sum(), ' below -1.5',low.sum())
fig,ax=plt.subplots(1,3,figsize=(18,5))
ax[0].hist(cls[m],130,range=(-10,8),histtype='step',color='k'); ax[0].set_xlabel('integral pull from triton peak (sigma_t)'); ax[0].set_title('adopted (d,t) tritons, matched'); [ax[0].axvline(z,ls=':') for z in (1.5,4.5)]
bins=np.arange(-4,16,0.25)
for lab,w,col in [('core |pull|<1.5',core,'C0'),('low pull<-1.5 (d side)',low,'C2'),('very low pull<-3',m&(cls<-3),'C1'),('shoulder 1.5-4.5',sh,'C3')]:
    hh,_=np.histogram(EX[w],bins); ax[1].stairs(hh/hh.sum(),bins,color=col,label=f'{lab} N={w.sum()}')
    print(lab, 'Ex<-1 %.1f%%  -1..1.5 %.1f%%  1.5..4.5 %.1f%%  4.5..8 %.1f%%  >8 %.1f%%'%tuple(100*np.mean(c) for c in [EX[w]<-1,(EX[w]>=-1)&(EX[w]<1.5),(EX[w]>=1.5)&(EX[w]<4.5),(EX[w]>=4.5)&(EX[w]<8),EX[w]>=8]))
ax[1].legend(); ax[1].set_xlabel('Ex(15C) [MeV]'); ax[1].set_title('Ex, normalised'); [ax[1].axvline(z,color='0.7',ls=':') for z in (0,0.74,3.1,4.66,6.7)]
ax[2].scatter(bgf[core],TCM[core],s=2,label='core'); ax[2].scatter(bgf[low],TCM[low],s=2,c='C2',label='low'); ax[2].set_xlabel('Brho (genfit)'); ax[2].set_ylabel('theta_cm'); ax[2].legend()
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_d2_tau/dt_shoulder_ex.png',dpi=90)
print('\n--- theta_cm windows ---')
for a,b in [(0,40),(40,80),(80,180)]:
    w=m&(TCM>=a)&(TCM<b); n=w.sum()
    print(f'theta_cm {a}-{b}: matched {n}; core {np.mean(core[w])*100:.1f}%  low(d side) {np.mean(low[w])*100:.1f}%  shoulder {np.mean(sh[w|False])*0+np.sum(sh&w)/max(n,1)*100:.1f}%')
    for lab,c in [('core',core),('low',low)]:
        ww=w&c; e=EX[ww]
        print(f'    {lab:5s} N={ww.sum():4d}  Ex: <-1 {np.mean(e<-1)*100:4.1f}%  g.s.+0.74 (-1..1.5) {np.mean((e>=-1)&(e<1.5))*100:4.1f}%  3.1 (2.3..3.9) {np.mean((e>=2.3)&(e<3.9))*100:4.1f}%  >8 {np.mean(e>=8)*100:4.1f}%')
