# All Dec2014 alpha plots, current configuration:
#   E_beam = 11.40 MeV (paper/data, NOT the log's 7.80), 299.5 torr, v_D = 2.142 cm/us,
#   loose selection (arm >= 25 mm, >= 6 hits), all 68 B=0 runs.
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

TH,EN,THs,ENs,Ecm,Epv,Tcm,Lall=[],[],[],[],[],[],[],[]
for pth in sorted(glob.glob("loose/P300_B0/cands_*.pkl")):
    for rec in pickle.load(open(pth,"rb"))["cands"]:
        T=(rec["tb"]-T0)*160.0*1e-3
        P3=np.column_stack([rec["x"],rec["y"],VD*T*10])
        segs={k:fit_line(P3[rec[k]],rec["q"][rec[k]]) for k in ("a1","a2","beam")}
        v=vertex_from_tracks([segs["a1"],segs["a2"]])
        cb,db=segs["beam"]; den=db[0]**2+db[1]**2
        if den<1e-12: continue
        tt=-(cb[0]*db[0]+cb[1]*db[1])/den; ent=cb+tt*db
        bd=v-ent; n=np.linalg.norm(bd)
        if n<1e-9: continue
        bd=bd/n
        u={};E={};L={}
        ok=True
        for k in ("a1","a2"):
            m=rec[k]; P=P3[m]; c,dd=segs[k]
            t=P@dd; tv=float(v@dd)
            far=P[np.argmax(np.abs(t-tv))]; w=far-v; nn=np.linalg.norm(w)
            if nn<1e-9: ok=False; break
            u[k]=w/nn; L[k]=float(np.ptp(t)); E[k]=float(EofR(L[k]))
        if not ok: continue
        op=np.degrees(np.arccos(np.clip(u["a1"]@u["a2"],-1,1)))
        for k in ("a1","a2"):
            th=np.degrees(np.arccos(np.clip(u[k]@bd,-1,1)))
            TH.append(th); EN.append(E[k]); Lall.append(L[k])
            if abs(op-90)<3.2: THs.append(th); ENs.append(E[k])
        if abs(op-90)<3.2:
            Ecm.append((E["a1"]+E["a2"])/2.0)
            Epv.append(float(Eaft(np.linalg.norm(v-ent)))/2.0)
            Tcm.append(2*min(np.degrees(np.arccos(np.clip(u[k]@bd,-1,1))) for k in ("a1","a2")))
TH,EN,THs,ENs=map(np.array,(TH,EN,THs,ENs))
Ecm,Epv,Tcm,Lall=map(np.array,(Ecm,Epv,Tcm,Lall))
print(f"  {len(TH)} arm tracks, {len(Ecm)} elastic events")

def scat(E,th,m,M):
    s=np.sin(np.radians(th)); c=np.cos(np.radians(th)); disc=M*M-m*m*s*s
    o=np.full(th.shape,np.nan); ok=disc>=0
    o[ok]=E*((m*c[ok]+np.sqrt(disc[ok]))/(m+M))**2; return o
th=np.linspace(0.5,90,700)
TAG=f"$E_b$={EB:.2f} MeV, 299.5 torr, $v_D$={VD} cm/$\\mu$s, arm$\\geq$25 mm, 68 runs"

# 1. kinematics
fig,ax=plt.subplots(1,2,figsize=(15,5.6))
for a,(X,Y,t) in zip(ax,[(TH,EN,f"all arm tracks ({len(TH)})"),(THs,ENs,f"90$^\\circ$ selection ({len(THs)})")]):
    a.hist2d(X,np.clip(Y,0,12),bins=[np.arange(0,90.5,0.4),np.arange(0,12,0.05)],cmap="viridis",cmin=1)
    a.plot(th,scat(EB,th,4,4),"w-",lw=2.2,label=r"$\alpha+\alpha$ elastic")
    for f in (0.7,0.45,0.25): a.plot(th,scat(EB*f,th,4,4),"w-",lw=1,alpha=.45)
    a.plot(th,scat(EB,th,4,12),"-",color="#FF7F3F",lw=1.6,label=r"$\alpha+^{12}$C")
    a.set_xlabel(r"lab angle to the beam $\theta$ [deg]"); a.set_ylabel("energy from range [MeV]")
    a.set_title(t); a.set_xlim(0,90); a.set_ylim(0,12); a.grid(alpha=.2); a.legend(fontsize=8)
fig.suptitle(TAG,fontsize=9,y=0.995); fig.tight_layout(); fig.savefig("plots/01_kinematics.png",dpi=140); plt.close(fig)

# 2. excitation function
dEdx=-np.gradient(e,d)
bins=np.arange(0,6.2,0.12); c=0.5*(bins[1:]+bins[:-1])
y,_=np.histogram(Ecm,bins=bins); w=np.interp(2*c,e[::-1],dEdx[::-1])
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
ax[0].hist(Ecm,bins=bins,color="#4C72B0",label="from alpha ranges")
ax[0].hist(Epv,bins=bins,color="#DD8452",alpha=.6,label="from vertex depth")
ax[0].axvline(EB/2,color="k",ls="--",lw=1.2,label=f"{EB/2:.2f} MeV = $E_b$/2")
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events / 120 keV"); ax[0].legend(fontsize=8)
ax[0].set_title(f"$E_{{cm}}$, two independent determinations ({len(Ecm)} events)")
m=(c>0.4)&(c<EB/2)
ax[1].errorbar(c[m],(y*w)[m],yerr=(np.sqrt(np.maximum(y,1))*w)[m],fmt="o-",ms=3,lw=1,color="#C44E52")
ax[1].axvline(2.94,color="k",ls="--",lw=1.2,label="$^8$Be 2$^+$ ($E_x$=3.03)")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel(r"yield $\times$ d$E$/d$x$ [arb.]")
ax[1].set_title("excitation function"); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.suptitle(TAG,fontsize=9,y=0.995); fig.tight_layout(); fig.savefig("plots/02_excitation.png",dpi=140); plt.close(fig)

# 3. acceptance map
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
h=ax[0].hist2d(Ecm,Tcm,bins=[np.arange(0,6.0,0.1),np.arange(0,95,2)],cmap="viridis")
fig.colorbar(h[3],ax=ax[0]); ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\theta_{cm}$ [deg]")
ax[0].set_title("acceptance map")
for lo,hi,col in ((70,90,"#C44E52"),(50,70,"#4C72B0"),(30,50,"#55A868")):
    mm=(Tcm>=lo)&(Tcm<hi)
    yy,ed=np.histogram(Ecm[mm],bins=np.arange(0,6.0,0.15)); cc=0.5*(ed[1:]+ed[:-1])
    ax[1].step(cc,yy,where="mid",color=col,label=f"$\\theta_{{cm}}$ {lo}-{hi}$^\\circ$ (N={mm.sum()})")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel("events / 150 keV")
ax[1].set_title("yield in $\\theta_{cm}$ bands"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.25)
fig.suptitle(TAG,fontsize=9,y=0.995); fig.tight_layout(); fig.savefig("plots/03_acceptance.png",dpi=140); plt.close(fig)

# 4. range-energy + the arm-length cut
R150=(RMAX-d)[::-1]; Ee=e[::-1]; R300=R150/S
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
ax[0].plot(R300,Ee,lw=2.2,color="#C44E52",label="299.5 torr")
ax[0].plot(R150,Ee,lw=1.4,ls="--",color="#4C72B0",label="150 torr (table density)")
for R,lab,col in ((25,"loose cut","#55A868"),(60,"original cut","#8172B3")):
    Em=float(np.interp(R,R300,Ee)); ax[0].axvline(R,color=col,lw=1.3)
    ax[0].annotate(f" {lab}: {R} mm = {Em:.2f} MeV",(R,Em),color=col,fontsize=8.5,va="bottom")
ax[0].set_xlim(0,600); ax[0].set_ylim(0,12); ax[0].set_xlabel("range [mm]"); ax[0].set_ylabel("alpha energy [MeV]")
ax[0].set_title("range-energy, He:CO$_2$ 90/10"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].hist(Lall,bins=np.arange(0,400,4),color="#4C72B0")
ax[1].axvline(25,color="#55A868",lw=1.5,label="loose cut 25 mm"); ax[1].axvline(60,color="#8172B3",lw=1.5,label="original 60 mm")
ax[1].set_xlabel("arm length [mm]"); ax[1].set_ylabel("arms / 4 mm"); ax[1].legend(fontsize=8)
ax[1].set_title("measured arm lengths"); ax[1].grid(alpha=.25)
fig.suptitle(TAG,fontsize=9,y=0.995); fig.tight_layout(); fig.savefig("plots/04_range_and_cut.png",dpi=140); plt.close(fig)
print("  wrote plots/01_kinematics.png 02_excitation.png 03_acceptance.png 04_range_and_cut.png")
