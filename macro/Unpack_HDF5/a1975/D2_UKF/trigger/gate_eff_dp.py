import numpy as np
T=np.loadtxt('/mnt/f/a1975/sim_rebuild/trig/tracks.txt')  # run ev gate sth svz ske sex stcm zmin zmax nh rzmin rzmax dirB
gate,sth,svz,ske,sex,stcm=T[:,2],T[:,3],T[:,4],T[:,5],T[:,6],T[:,7]
b=(sex>-1.5)&(sex<1.0)&(svz>200)&(svz<900); print('tracks',len(T),'bound-window (Spyral Ex) vz200-900',b.sum())
print('\ngate pass fraction, bound window, by theta_cm (Spyral)')
for a,c in [(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,60)]:
    w=b&(stcm>=a)&(stcm<c); print(f'  {a}-{c}: N {w.sum():5d}  in gate {gate[w].mean()*100 if w.sum() else 0:5.1f}%  median th_lab {np.median(sth[w]) if w.sum() else 0:6.1f}')
print('\ngate pass fraction, bound window, by theta_lab')
for a in range(89,131,3):
    w=b&(sth>=a)&(sth<a+3); print(f'  {a}-{a+3}: N {w.sum():5d}  in gate {gate[w].mean()*100 if w.sum() else 0:5.1f}%')
