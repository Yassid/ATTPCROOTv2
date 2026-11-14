import uproot, numpy as np, glob
c=['thl','vz','ke','ex','tcm','sv','sf','nh']
def load(fs):
    A=[uproot.open(f)['c'].arrays(c,library='np') for f in fs]; return {k:np.concatenate([a[k] for a in A]) for k in c}
S=load(['/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root']); D=load(sorted(glob.glob('/mnt/f/a1975/census_arc/data_*.root')))
def tab(X,lab):
    base=(X['vz']>200)&(X['vz']<900)&(X['thl']>=89)&(X['thl']<=161)&(X['ke']>=1.2)&np.isfinite(X['tcm'])
    win=base&(X['ex']>-1.5)&(X['ex']<1.0); sb=base&(((X['ex']>-4.0)&(X['ex']<-1.5))|((X['ex']>1.5)&(X['ex']<2.0)))   # sidebands below g.s. and in the 1-2 gap
    k=2.5/3.0  # window width / sideband width
    print(f'\n{lab}: window {win.sum()}  sideband {sb.sum()}')
    print('  theta_cm  win(N,valid)  side(N,valid)   net N   net valid   Spyral eff (bkg-subtracted)')
    for a,cc in [(15,20),(20,25),(25,30),(30,35),(35,40),(40,45)]:
        t=(X['tcm']>=a)&(X['tcm']<cc); nw=(win&t).sum(); vw=(win&t&(X['sv']>0)).sum(); ns=(sb&t).sum(); vs=(sb&t&(X['sv']>0)).sum()
        N=nw-k*ns; V=vw-k*vs; e=V/N if N>0 else np.nan; de=np.sqrt(e*(1-e)/max(N,1)+ (k*k*(ns))/max(N,1)**2) if N>0 else np.nan
        print(f'  {a}-{cc}:  {nw:5d} {vw:5d}   {ns:5d} {vs:5d}   {N:7.1f} {V:7.1f}   {e:5.2f} +- {de:.2f}')
tab(S,'SIM g.s.'); tab(D,'DATA (47 runs)')
M1=16.0147013*931.49401
def rres(ke,th,vz):
    t=np.radians(th); p=np.sqrt(ke**2+2*ke*938.272); E=183.871-0.011071*vz; pb=np.sqrt(E**2+2*E*M1); return (1000-vz)*np.tan(np.arctan2(p*np.sin(t),pb-p*np.cos(t)))
for R in (0,20,28):
    for X,lab in [(S,'SIM'),(D,'DATA')]:
        r=rres(X['ke'],X['thl'],X['vz']); keep=r>=R
        Y={k:v[keep] for k,v in X.items()}; print(f'\n===== r_res >= {R}'); tab(Y,lab)
