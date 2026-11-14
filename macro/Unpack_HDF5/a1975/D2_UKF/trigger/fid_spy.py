# Fiducial (d,p) bound dsigma with the Spyral-validity term (data are PID-gated => need Spyral success; sim acceptance was ungated)
import numpy as np, uproot, io, contextlib, runpy
with contextlib.redirect_stdout(io.StringIO()): G=runpy.run_path('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py')
run,L2,L0,TB=G['run'],G['L2'],G['L0'],G['TB']
X=uproot.open('/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root')['c'].arrays(['thl','vz','ke','ex','tcm','sv'],library='np')
M1=16.0147013*931.49401
def rres(ke,th,vz):
    t=np.radians(th); p=np.sqrt(ke**2+2*ke*938.272); E=183.871-0.011071*vz; pb=np.sqrt(E**2+2*E*M1); return (1000-vz)*np.tan(np.arctan2(p*np.sin(t),pb-p*np.cos(t)))
def fit1(x,er,m,u): S=np.sum((x*m/er**2)[u])/np.sum((m*m/er**2)[u]); return S,np.sqrt(1/np.sum((m*m/er**2)[u])),np.sum((((x-S*m)/er)**2)[u])
for R in (26,28):
    x,er,Y,B,e,Le=run(R)
    b=(X['ex']>-1.5)&(X['ex']<1.0)&(X['vz']>200)&(X['vz']<900)&(X['thl']>=89)&(X['thl']<=161)&(X['ke']>=1.2)&(rres(X['ke'],X['thl'],X['vz'])>=R)
    sv=np.array([np.mean(X['sv'][b&(X['tcm']>=a)&(X['tcm']<c)]>0) if (b&(X['tcm']>=a)&(X['tcm']<c)).sum()>20 else np.nan for a,c in zip(TB[:-1],TB[1:])])
    xc=x/np.where(np.isfinite(sv),sv,1); ec=er/np.where(np.isfinite(sv),sv,1)
    print(f'\nR>={R}  Spyral validity per bin', np.round(sv,2))
    print('  dsdo old ', ' '.join(f'{a:5.2f}({b:4.2f})' if Y[i]>=3 else '     -     ' for i,(a,b) in enumerate(zip(x,er))))
    print('  dsdo new ', ' '.join(f'{a:5.2f}({b:4.2f})' if Y[i]>=3 else '     -     ' for i,(a,b) in enumerate(zip(xc,ec))))
    for lab,u in [('25-40',(Y>=3)&(TB[:-1]>=25)&(TB[1:]<=40)),('25-45',(Y>=3)&(TB[:-1]>=25)&(TB[1:]<=45))]:
        so=fit1(x,er,L2,u); sn=fit1(xc,ec,L2,u); n=u.sum()
        print(f'  l=2 fit {lab}: old S {so[0]:.3f}({so[1]:.3f}) chi2 {so[2]:5.1f}/{n-1}   new S {sn[0]:.3f}({sn[1]:.3f}) chi2 {sn[2]:5.1f}/{n-1}')
