import numpy as np, uproot, glob
E=[uproot.open(f)['e'].arrays(['mode','found','nhD','tc0','vz0','nhP'],library='np') for f in sorted(glob.glob('/mnt/f/a1975/sim_rebuild/dp_embed/e_r*.root'))]
E={k:np.concatenate([e[k] for e in E]) for k in E[0]}; m=(E['mode']==2)&(E['tc0']>=20)&(E['tc0']<45)
nhD,fd=E['nhD'][m],E['found'][m]>0
edges=np.array([0,50,100,150,200,300,400,600,800,1200,5000])
s=np.array([fd[(nhD>=a)&(nhD<b)].mean() if ((nhD>=a)&(nhD<b)).sum()>30 else np.nan for a,b in zip(edges[:-1],edges[1:])])
nr=np.histogram(nhD,edges)[0]
print('overlay random-event size bins',edges.tolist()); print('  survival s(env) ',np.round(s,3)); print('  random-event freq',np.round(nr/nr.sum(),3), ' mean survival', fd.mean().round(3))
W=np.loadtxt('evsize_win.txt'); B=np.loadtxt('evsize_sb.txt')   # run ev nh_track nh_event
for lab,X in [('bound window',W),('sideband',B)]:
    env=X[:,3]-X[:,2]; print(f'\n{lab}: N {len(X)}  event hits median {np.median(X[:,3]):.0f}, track hits median {np.median(X[:,2]):.0f}, environment median {np.median(env):.0f}  p16/p84 {np.percentile(env,[16,84]).round(0)}')
    if lab=='bound window':
        no=np.histogram(env,edges)[0].astype(float); print('  observed env freq ',np.round(no/no.sum(),3))
        ok=np.isfinite(s)&(s>0); real=no[ok]/s[ok]; print(f'  survival for real (d,p) events = N_obs / sum(N_obs/s) = {no[ok].sum()/real.sum():.3f}   (cf overlay average {fd.mean():.3f})')
# env with sideband subtraction (scale = window width 2.229 / sideband 2.5)
k=2.229/2.5; nw=np.histogram(W[:,3]-W[:,2],edges)[0]; ns=np.histogram(B[:,3]-B[:,2],edges)[0]; net=np.clip(nw-k*ns,0,None); ok=np.isfinite(s)&(s>0)
print(f'\nbackground-subtracted: survival = {net[ok].sum()/(net[ok]/s[ok]).sum():.3f}')
