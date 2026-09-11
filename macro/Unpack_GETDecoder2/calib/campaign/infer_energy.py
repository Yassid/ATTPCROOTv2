# Infer the beam energy of a setting from the data alone.
#
# For elastic scattering the two products share the beam's remaining energy, so the sum of
# their ranges falls linearly with how far the beam has already travelled, and extrapolates
# to zero at the beam's TOTAL range. That intercept is a pure length measurement -- no
# assumed beam energy anywhere. Inverting the range-energy relation then gives the energy.
import numpy as np, glob, pickle, warnings, sys
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from cluster_tracks import fit_line, vertex_from_tracks

# Bethe range-energy for He:CO2 90/10, used only to turn a measured range into an energy
fHe,fCO2=0.9,0.1
M=fHe*4.003+fCO2*44.01
wHe,wCO2=fHe*4.003/M, fCO2*44.01/M
ZA=wHe*(2/4.003)+wCO2*(22/44.01)
I=np.exp((wHe*(2/4.003)*np.log(41.8)+wCO2*(22/44.01)*np.log(85.0))/ZA)
Mal,me=3727.379,0.510999
def rng_mm(E,press,EMIN=0.20,step=0.5):
    rho=0.145e-3*press/330.0
    d=0.0
    while E>EMIN and d<20000:
        b2=2*E/Mal
        dd=0.307075*4*ZA/b2*np.log(2*me*1e6*b2/I)*rho*0.1
        if dd<=1e-6: break
        E-=dd*step; d+=step
    return d
def E_from_range(R,press):
    lo,hi=0.5,30.0
    for _ in range(60):
        mid=0.5*(lo+hi)
        if rng_mm(mid,press)<R: lo=mid
        else: hi=mid
    return 0.5*(lo+hi)

T0=-34.0
def measure(setting,press,vD,B=0.0):
    dep,arm=[],[]
    for pth in sorted(glob.glob(f"{setting}/cands_*.pkl")):
        for rec in pickle.load(open(pth,"rb"))["cands"]:
            T=(rec["tb"]-T0)*160.0*1e-3
            P3=np.column_stack([rec["x"],rec["y"],vD*T*10])
            segs={k:fit_line(P3[rec[k]],rec["q"][rec[k]]) for k in ("a1","a2","beam")}
            v=vertex_from_tracks([segs["a1"],segs["a2"]])
            u={};L={}
            for k in ("a1","a2"):
                m=rec[k]; P=P3[m]; c,dd=segs[k]
                t=P@dd; tv=float(v@dd)
                far=P[np.argmax(np.abs(t-tv))]; w=far-v; n=np.linalg.norm(w)
                if n<1e-9: break
                u[k]=w/n; L[k]=float(np.ptp(t))
            if len(u)<2: continue
            op=np.degrees(np.arccos(np.clip(u["a1"]@u["a2"],-1,1)))
            if abs(op-90)>4.0: continue
            cb,db=segs["beam"]; den=db[0]**2+db[1]**2
            if den<1e-12: continue
            tt=-(cb[0]*db[0]+cb[1]*db[1])/den; ent=cb+tt*db
            dep.append(float(np.linalg.norm(v-ent))); arm.append(L["a1"]+L["a2"])
    dep,arm=np.array(dep),np.array(arm)
    if len(dep)<200: return None
    # fit only where the relation is genuinely linear: the well-populated core of the depth
    # distribution. Including the sparse tails flattens the slope and sends the intercept to
    # infinity (that is what gave a nonsense 215660 mm on the first pass).
    lo_d,hi_d=np.percentile(dep,[10,85])
    cs,ys=[],[]
    step=max(40.0,(hi_d-lo_d)/10)
    x=lo_d
    while x<hi_d:
        m=(dep>=x)&(dep<x+step)
        if m.sum()>=25: cs.append(x+step/2); ys.append(np.median(arm[m]))
        x+=step
    if len(cs)<4: return None
    sl,ic=np.polyfit(np.array(cs),np.array(ys),1)
    if sl>=-0.05: return None                      # must fall with depth
    R=-ic/sl
    if not (100<R<3000): return None
    return len(dep),R,E_from_range(R,press)

print(f"  {'setting':>16} {'N':>6} {'range(mm)':>10} {'E_beam':>8} {'MeV/u':>7}  log says")
for s,press,vD,logE in (("P300_B0",299.5,2.142,"1.95"),
                        ("P150_B057",150.0,2.142,"1.95"),
                        ("P300_B102_195",298.7,2.142,"1.95"),
                        ("P300_B102_239",298.7,2.142,"2.39")):
    r=measure(s,press,vD)
    if r is None: print(f"  {s:>16}   (not enough events yet)"); continue
    n,R,E=r
    print(f"  {s:>16} {n:6d} {R:10.0f} {E:8.2f} {E/4:7.2f}  {logE} MeV/u")
