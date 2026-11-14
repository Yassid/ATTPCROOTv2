import numpy as np, runpy, io, contextlib
src=open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read()
src=src.split("print('theta_cm 10-15")[0]
for TL in (89,85,80,75,70):
    s=src.replace("(g[:, 1] >= 89)",f"(g[:, 1] >= {TL})").replace("(d[:, 1] >= 89)",f"(d[:, 1] >= {TL})")
    G={}; exec(s,G); TB=G['TB']
    for R in (0,28):
        x,er,Y,B,e,Le=G['run'](R)
        print(f'th_lab>={TL} R>={R:2d}', ' '.join(f'{a:5.2f}({b:4.2f})' if Y[i]>=3 else '     -     ' for i,(a,b) in enumerate(zip(x,er))), ' Y',Y.astype(int)[2:], 'eps',np.round(e,2)[2:], 'B',np.round(B,1)[4:])
# lab-angle of bound-state candidates near the edge: data vs sim, theta_cm 38-47
G={}; exec(src,G); d,dx,dc=G['d'],G['dx'],G['dc']; g,sx,sc=G['g'],G['sx'],G['sc']
w=(d[:,2]>=200)&(d[:,2]<900)&(dx>=-1.5)&(dx<=0.729)&(d[:,3]<=13.5)&(d[:,0]>=1.2)
s=(g[:,2]>=200)&(g[:,2]<900)&(g[:,3]<=13.5)&(g[:,0]>=1.2)&np.isfinite(sx)
print('\nbound window, theta_lab histogram 80..100 (2 deg), by theta_cm')
for a,b in [(35,40),(40,45),(45,50)]:
    hd,_=np.histogram(d[w&(dc>=a)&(dc<b),1],np.arange(70,101,2)); hs,_=np.histogram(g[s&(sc>=a)&(sc<b),1],np.arange(70,101,2))
    print(f' {a}-{b} data',hd,'\n        sim ',hs)
print('sim true vs reco th_lab (cols 1 vs 5) near 90: median diff', np.median((g[:,1]-g[:,5])[np.abs(g[:,5]-90)<5]).round(3))
