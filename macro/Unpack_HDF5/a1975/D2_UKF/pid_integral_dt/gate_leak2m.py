import glob, os, json, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.path import Path; from scipy.optimize import curve_fit; from scipy.ndimage import uniform_filter1d
fs=sorted(f for f in glob.glob('/mnt/f/a1975/reco_d2_tau/pid_width_00[0-9][0-9].root') if not f.endswith('0016.root') and os.path.exists(f+'.done'))
c=['va','br','sqa','sqw']; A={k:np.concatenate([uproot.open(f)['w'].arrays(c,library='np')[k] for f in fs]) for k in c}
v=(A['va']>0)&(A['br']>0.05)&(A['sqa']>0)&(A['sqw']>0); br,sa,sw=A['br'][v],A['sqa'][v],A['sqw'][v]
G=np.array(json.load(open(os.path.expanduser('~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/triton_d2_dv1104.json')))['vertices'])
ing=Path(G).contains_points(np.c_[sa,br]); yw=np.log(sw*br)
def g(x,a,m,s): return a*np.exp(-.5*((x-m)/s)**2)
def g3(x,*p): return g(x,*p[0:3])+g(x,*p[3:6])+g(x,*p[6:9])+p[9]
dt=np.load('/mnt/f/a1975/sim_rebuild/trig/dt/joined_dt.npz'); T=dt['T']; s=dt['sel'].astype(bool); ke=T[s,3]; tb=np.sqrt(ke**2+2*2808.921*ke)/299.792458
edges=np.round(np.arange(0.3,1.41,0.1),2); rows=[]
fig,ax=plt.subplots(4,5,figsize=(22,15)); ax=ax.ravel()
for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
    sl=(br>=lo)&(br<hi); y=yw[sl]; h,e=np.histogram(y,160,range=(1.3,3.3)); xc=.5*(e[1:]+e[:-1]); hs=uniform_filter1d(h.astype(float),5)
    # all tracks in slice: fit d, t, 4He + const; seeds from the ladder around the d mode in [1.7,2.2]
    w=(xc>1.7)&(xc<2.25); md=xc[w][np.argmax(hs[w])]
    p0=[hs[w].max(),md,.07, max(hs[np.abs(xc-md-.42)<.03].max(),1),md+.42,.07, max(hs[np.abs(xc-md-.75)<.05].max(),1),md+.75,.07, 1]
    lb=[0,md-.05,.02, 0,md+.30,.03, 0,md+.6,.03, 0]; ub=[np.inf,md+.05,.3, np.inf,md+.52,.15, np.inf,md+.95,.15, np.inf]
    try: p,_=curve_fit(g3,xc,h,p0=p0,bounds=(lb,ub),sigma=np.sqrt(np.maximum(h,1)),maxfev=40000); ok=1
    except Exception: p=np.array(p0); ok=0
    mt,st,ma,sa_=p[4],p[5],p[7],p[8]
    # gated tracks: t and 4He shapes FIXED from the all-track fit, only amplitudes free
    hg,_=np.histogram(yw[sl&ing],160,range=(1.3,3.3))
    f2=lambda x,at,aa,c0: g(x,at,mt,st)+g(x,aa,ma,sa_)+c0
    try: q,_=curve_fit(f2,xc,hg,p0=[hg.max(),hg.max()/3,0.1],bounds=([0,0,0],[np.inf,np.inf,np.inf]),sigma=np.sqrt(np.maximum(hg,1)))
    except Exception: q=[np.nan]*3
    Nt,Na=q[0]*st*np.sqrt(2*np.pi)/(e[1]-e[0]),q[1]*sa_*np.sqrt(2*np.pi)/(e[1]-e[0]); fa=Na/(Nt+Na)
    yg=yw[sl&ing]; Lc=np.sum((yg<mt)&(yg>mt-0.3)); Rc=np.sum((yg>=mt)&(yg<mt+0.3)); fa=(Rc-Lc)/len(yg)
    wdt=np.mean((tb>=lo)&(tb<hi)); rows.append((lo,hi,(sl&ing).sum(),mt-p[1],ma-p[1],st,sa_,fa,wdt,ok))
    b=ax[i]; b.stairs(h,e,color='0.6',label='all'); b.stairs(hg,e,color='C0',label='in triton gate'); xx=np.linspace(1.3,3.3,500)
    b.plot(xx,g(xx,q[0],mt,st),'C2',lw=1,label='t'); b.plot(xx,g(xx,q[1],ma,sa_),'C3',lw=1,label='4He'); b.set_yscale('log'); b.set_ylim(.7,None)
    b.axvline(mt,color='C2',ls=':'); b.set_title(f'Brho {lo:.1f}-{hi:.1f}: excess above t {fa*100:.0f}%',fontsize=9); b.legend(fontsize=6)
[a.axis('off') for a in ax[len(edges)-1:]]
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_d2_tau/dt_gate_leak_integral_slices_mirror.png',dpi=80)
print(f"{'Brho':8s} {'gated':>6s} {'t-d':>6s} {'a-d':>6s} {'sig_t':>6s} {'sig_a':>6s} {'mirror excess':>8s} {'adopted t share':>15s}")
for r in rows: print(f"{r[0]:.1f}-{r[1]:.1f}  {r[2]:6d} {r[3]:+6.3f} {r[4]:+6.3f} {r[5]:6.3f} {r[6]:6.3f} {r[7]*100:7.1f}% {r[8]*100:14.1f}%"+('' if r[9] else ' FITFAIL'))
R=np.array([(r[7],r[8]) for r in rows]); ok=np.isfinite(R[:,0])
print('weighted by adopted-triton Brho (0.3-1.4 covers %.0f%% of them): 4He fraction %.1f%%'%(R[ok,1].sum()*100, (R[ok,0]*R[ok,1]).sum()/R[ok,1].sum()*100))
