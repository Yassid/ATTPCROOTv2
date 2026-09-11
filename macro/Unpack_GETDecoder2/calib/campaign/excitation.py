# Excitation function for the B=0 set with the CORRECTED beam energy, 11.40 MeV.
import numpy as np, glob, pickle, warnings, sys
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from cluster_tracks import fit_line, vertex_from_tracks
S=299.5/150.0; T0=-34.0; VD=2.142; EB=11.40
d,e=[],[]
for L in open("eloss_alpha_heco2_11p4.txt"):
    if L.startswith("#"): continue
    a,b=L.split(); d.append(float(a)); e.append(float(b))
d,e=np.array(d),np.array(e); RMAX=d[-1]
EofR=lambda R: np.interp(RMAX-np.clip(R*S,0,RMAX),d,e)
Eaft=lambda p: np.interp(p*S,d,e)
Ea,Ep,Tcm=[],[],[]
for pth in sorted(glob.glob("P300_B0/cands_*.pkl")):
    for rec in pickle.load(open(pth,"rb"))["cands"]:
        T=(rec["tb"]-T0)*160.0*1e-3
        P3=np.column_stack([rec["x"],rec["y"],VD*T*10])
        segs={k:fit_line(P3[rec[k]],rec["q"][rec[k]]) for k in ("a1","a2","beam")}
        v=vertex_from_tracks([segs["a1"],segs["a2"]])
        db=segs["beam"][1]
        cb,_=segs["beam"]; den=db[0]**2+db[1]**2
        if den<1e-12: continue
        tt=-(cb[0]*db[0]+cb[1]*db[1])/den; ent=cb+tt*db
        bd=(v-ent); bd=bd/np.linalg.norm(bd)
        u={};L_={}
        for k in ("a1","a2"):
            m=rec[k]; P=P3[m]; c,dd=segs[k]
            t=P@dd; tv=float(v@dd)
            far=P[np.argmax(np.abs(t-tv))]; w=far-v; n=np.linalg.norm(w)
            if n<1e-9: break
            u[k]=w/n; L_[k]=float(np.ptp(t))
        if len(u)<2: continue
        op=np.degrees(np.arccos(np.clip(u["a1"]@u["a2"],-1,1)))
        if abs(op-90)>3.2: continue
        Ea.append(float(EofR(L_["a1"])+EofR(L_["a2"])))
        Ep.append(float(Eaft(np.linalg.norm(v-ent))))
        Tcm.append(2*min(np.degrees(np.arccos(np.clip(u[k]@bd,-1,1))) for k in ("a1","a2")))
Ea,Ep,Tcm=map(np.array,(Ea,Ep,Tcm))
print(f"  {len(Ea)} elastic events")
Ecm=Ea/2.0
print(f"  E_cm from alpha energies: median {np.median(Ecm):.2f}, 16-84% {np.percentile(Ecm,16):.2f}..{np.percentile(Ecm,84):.2f}")
print(f"  E_cm from vertex depth  : median {np.median(Ep/2):.2f}")
print(f"  ceiling {EB/2:.2f} MeV;  above it: {100*(Ecm>EB/2).mean():.0f}%")
dEdx=-np.gradient(e,d)
bins=np.arange(0,6.2,0.12); c=0.5*(bins[1:]+bins[:-1])
y,_=np.histogram(Ecm,bins=bins)
w=np.interp(2*c,e[::-1],dEdx[::-1])
fig,ax=plt.subplots(1,2,figsize=(12.5,4.5))
ax[0].hist(Ecm,bins=bins,color="#4C72B0",label="from alpha energies")
ax[0].hist(Ep/2,bins=bins,color="#DD8452",alpha=.6,label="from vertex depth")
ax[0].axvline(EB/2,color="k",ls="--",lw=1.2,label=f"{EB/2:.2f} MeV = E$_b$/2")
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events / 120 keV"); ax[0].legend(fontsize=8)
ax[0].set_title(f"$E_{{cm}}$ at $E_b$ = {EB:.2f} MeV ({len(Ecm)} events)")
m=(c>0.4)&(c<EB/2)
ax[1].errorbar(c[m],(y*w)[m],yerr=(np.sqrt(np.maximum(y,1))*w)[m],fmt="o-",ms=3,lw=1,color="#C44E52")
ax[1].axvline(2.94,color="k",ls="--",lw=1.2,label="$^8$Be 2$^+$ ($E_x$=3.03)")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel(r"yield $\times$ d$E$/d$x$ [arb.]")
ax[1].set_title("excitation function"); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.tight_layout(); fig.savefig("excitation_11p4.png",dpi=130)
print("  -> excitation_11p4.png")
