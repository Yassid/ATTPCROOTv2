import numpy as np, uproot, io, contextlib
src=open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("print('theta_cm 10-15")[0]
X=uproot.open('/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root')['c'].arrays(['thl','vz','ke','ex','tcm','sv'],library='np')
M1=16.0147013*931.49401
def rres(ke,th,vz):
    t=np.radians(th); p=np.sqrt(ke**2+2*ke*938.272); E=183.871-0.011071*vz; pb=np.sqrt(E**2+2*E*M1); return (1000-vz)*np.tan(np.arctan2(p*np.sin(t),pb-p*np.cos(t)))
def fit1(x,er,m,u): S=np.sum((x*m/er**2)[u])/np.sum((m*m/er**2)[u]); return S,np.sqrt(1/np.sum((m*m/er**2)[u])),np.sum((((x-S*m)/er)**2)[u])
VAR={'nominal':{},'chi2<5':{'(g[:, 3] <= 13.5)':'(g[:, 3] <= 5)','(d[:, 3] <= 13.5)':'(d[:, 3] <= 5)'},
     'Ex top 1.0':{'EH = -1.5, 0.729':'EH = -1.5, 1.0'},'smear .15':{'SM = 0.238':'SM = 0.15'},'smear .30':{'SM = 0.238':'SM = 0.30'},
     'vz 300-900':{'z0=200':'z0=300'},'vz 200-700':{'z1=900':'z1=700'}}
for name,rep in VAR.items():
    s=src
    for a,b in rep.items():
        assert a in s, (name,a); s=s.replace(a,b)
    G={}; 
    with contextlib.redirect_stdout(io.StringIO()): exec(s,G)
    TB,L2=G['TB'],G['L2']; out=[]
    for R in (24,26,28,30):
        x,er,Y,B,e,Le=G['run'](R)
        b=(X['ex']>-1.5)&(X['ex']<1.0)&(X['vz']>200)&(X['vz']<900)&(X['thl']>=89)&(X['thl']<=161)&(X['ke']>=1.2)&(rres(X['ke'],X['thl'],X['vz'])>=R)
        sv=np.array([np.mean(X['sv'][b&(X['tcm']>=a)&(X['tcm']<c)]>0) if (b&(X['tcm']>=a)&(X['tcm']<c)).sum()>20 else 1.0 for a,c in zip(TB[:-1],TB[1:])])
        xc,ec=x/sv,er/sv; u=(Y>=3)&(e>.3)&(TB[:-1]>=25)&(TB[1:]<=45); S=fit1(xc,ec,L2,u)
        out.append(f'R{R}: S {S[0]:.3f}({S[1]:.3f}) chi2 {S[2]:4.1f}/{u.sum()-1}')
    print(f'{name:12s} '+' | '.join(out))
