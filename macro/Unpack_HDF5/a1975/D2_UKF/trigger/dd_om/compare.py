import numpy as np
def F(n): return np.loadtxt(f'f/outputs/dd_{n}_dsdo.dat')
def P(fn):
    L=open(fn).read().splitlines(); s=next(i for i,l in enumerate(L) if 'ANGLE' in l and 'SIGMA/' in l); d={}
    for l in L[s:]:
        if 'TOTAL REACTION' in l: break
        t=l.split()
        try:
            if len(t)>=5: d.setdefault(float(t[0]),float(t[3]))
        except: pass
    return d
p=P('ptol_D.out'); h=F('D_22p4'); f=F('Dfull_22p4')
print('Daehnick 22.4 MeV, FRESCO/PtolemyCpp with Vso as given and Vso/2')
for a in (20,30,40,50,60,70,90): print(f'  {a:3d} {np.interp(a,f[:,0],f[:,1])/p[a]:.3f} {np.interp(a,h[:,0],h[:,1])/p[a]:.3f}')
print('\ndsigma/dOmega (mb/sr), Ed 22.4 MeV; An-Cai, Han-Shi-Shen, Zhang-Pang-Lou, Daehnick; spread=(max-min)/mean of first 3; dE = ZPL 21.9/23.1')
A,H,Z,D=[F(n+'_22p4') for n in 'AHZD']; a1=F('Z_21p9'); a2=F('Z_23p1')
for a in list(range(10,31,2))+[35,40,45,50]:
    i=a-1; v=np.array([A[i,1],H[i,1],Z[i,1]]); print(f'{a:3d} {A[i,1]:9.2f} {H[i,1]:9.2f} {Z[i,1]:9.2f} {D[i,1]:9.2f}  spread {(v.max()-v.min())/v.mean():.3f}  dE {a1[i,1]/a2[i,1]:.3f}')
