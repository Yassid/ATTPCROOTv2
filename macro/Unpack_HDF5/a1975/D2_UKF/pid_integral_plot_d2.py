import uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
t = uproot.open('/mnt/f/a1975/reco_d2_tau/pid_width_0016_axis.root')['w'].arrays(library='np')
v = (t['va']>0)&(t['br']>0.05)&(t['sqa']>0)&(t['sqw']>0)
br = t['br'][v]
def g(x,a,m,s): return a*np.exp(-.5*((x-m)/s)**2)
def fitpk(x,c,w):
    s=(x>c-w)&(x<c+w); h,e=np.histogram(x[s],40); xc=.5*(e[1:]+e[:-1])
    p,cv=curve_fit(g,xc,h,p0=[h.max(),c,w/2],sigma=np.sqrt(np.maximum(h,1))); return p
fig,ax=plt.subplots(2,2,figsize=(13,10))
for j,(lab,q) in enumerate([('peak amplitude (a1975 default)','sqa'),('trace integral (C15d option)','sqw')]):
    y=t[q][v]
    a=ax[0,j]; a.hist2d(br,y,bins=[300,300],range=[[0,np.percentile(br,99.5)],[0,np.percentile(y,99.5)]],norm=LogNorm(),cmap='viridis')
    a.set_xlabel('Brho (Tm)'); a.set_ylabel('sqrt(dE/dx) [arb]'); a.set_title(f'run 0016: {lab}')
    x=np.log(y*br); x=x-np.median(x[(x>np.percentile(x,5))&(x<np.percentile(x,95))])
    h,e=np.histogram(x,300,range=(np.percentile(x,.5),np.percentile(x,99.5))); xc=.5*(e[1:]+e[:-1])
    d0=xc[np.argmax(h)]; x=x-d0; e=e-d0; xc=xc-d0
    pd=fitpk(x,0,0.12); md,sd=pd[1],abs(pd[2]); x-=md; xc-=md; e-=md
    res={}
    for nm,off in [('p',-0.578),('d',0.),('t',0.37),('4He',0.693)]:
        w=(xc>off-.15)&(xc<off+.15); c=xc[w][np.argmax(h[w])]
        try: p=fitpk(x,c,0.10 if nm!='d' else 0.12); res[nm]=(p[1],abs(p[2]))
        except Exception: res[nm]=(np.nan,np.nan)
    b=ax[1,j]; b.stairs(h,e,color='k'); b.set_yscale('log'); b.set_xlim(-1.2,1.4)
    for nm,(m,s) in res.items(): b.axvline(m,color='C3',ls='--'); b.text(m,h.max()*1.5,nm,ha='center',color='C3')
    for xx,lb in [(-0.693+0.0,''),(0.405,''),(0.693,'')]: b.axvline(xx,color='C0',ls=':',alpha=.7)
    Rpd=(0-res['p'][0])/(res['d'][1]+res['p'][1])
    b.set_title(f"ln(sqrt(dE/dx)*Brho), d at 0   |   4He-d = {res['4He'][0]:+.3f} (expect 0.693),  p-d R = {Rpd:.2f}",fontsize=10)
    b.set_xlabel('ln m (shifted); blue dotted = ideal ln(m/m_d) for p, t, 4He'); b.set_ylim(0.5,h.max()*4)
    print(lab, {k:(round(m,3),round(s,3)) for k,(m,s) in res.items()}, 'Rpd',round(Rpd,2))
plt.tight_layout(); plt.savefig('/mnt/f/a1975/reco_d2_tau/pid_peak_vs_integral_a1975_0016.png',dpi=110)
