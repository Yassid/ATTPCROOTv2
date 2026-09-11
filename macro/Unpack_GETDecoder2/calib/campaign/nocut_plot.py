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
EofR=lambda R: float(np.interp(RMAX-np.clip(R*S,0,RMAX),d,e))
TH,EN,LN=[],[],[]
for pth in sorted(glob.glob("nocut/P300_B0/cands_*.pkl")):
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
        u={};L={}
        ok=True
        for k in ("a1","a2"):
            m=rec[k]; P=P3[m]; c,dd=segs[k]
            t=P@dd; tv=float(v@dd)
            far=P[np.argmax(np.abs(t-tv))]; w=far-v; nn=np.linalg.norm(w)
            if nn<1e-9: ok=False; break
            u[k]=w/nn; L[k]=float(np.ptp(t))
        if not ok: continue
        op=np.degrees(np.arccos(np.clip(u["a1"]@u["a2"],-1,1)))
        if abs(op-90)>3.2: continue
        for k in ("a1","a2"):
            TH.append(np.degrees(np.arccos(np.clip(u[k]@bd,-1,1))))
            EN.append(EofR(L[k])); LN.append(L[k])
TH,EN,LN=map(np.array,(TH,EN,LN))
print(f"  {len(TH)} arms in the 90-deg selection, no length cut")
def scat(E,th,m,M):
    s=np.sin(np.radians(th)); c=np.cos(np.radians(th)); disc=M*M-m*m*s*s
    o=np.full(th.shape,np.nan); ok=disc>=0
    o[ok]=E*((m*c[ok]+np.sqrt(disc[ok]))/(m+M))**2; return o
th=np.linspace(0.5,90,700)
fig,ax=plt.subplots(1,2,figsize=(15,5.6))
for a,(sel,t) in zip(ax,[(LN>=0,f"no cut  ({len(TH)} arms)"),
                          (LN<25,f"ONLY the short arms, <25 mm  ({(LN<25).sum()} arms)")]):
    a.hist2d(TH[sel],np.clip(EN[sel],0,12),bins=[np.arange(0,90.5,0.5),np.arange(0,12,0.06)],
             cmap="viridis",cmin=1)
    a.plot(th,scat(EB,th,4,4),"w-",lw=2.2,label=r"$\alpha+\alpha$, $E_b$=11.40")
    for f in (0.7,0.45,0.25,0.12): a.plot(th,scat(EB*f,th,4,4),"w-",lw=1,alpha=.45)
    a.set_xlabel(r"lab angle to the beam $\theta$ [deg]"); a.set_ylabel("energy from range [MeV]")
    a.set_title(t); a.set_xlim(0,90); a.set_ylim(0,12); a.grid(alpha=.2); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig("plots/05_nocut.png",dpi=140)
sh=LN<25
print(f"  short arms (<25 mm): {sh.sum()} = {100*sh.mean():.1f}% of all")
print(f"     their angles: median {np.median(TH[sh]):.0f} deg, 16-84% {np.percentile(TH[sh],16):.0f}..{np.percentile(TH[sh],84):.0f}")
print(f"     long arms   : median {np.median(TH[~sh]):.0f} deg")
lim=np.nan_to_num(scat(EB,TH,4,4),nan=0)
print(f"  above the locus: all {100*(EN>lim).mean():.0f}%, short arms {100*(EN[sh]>lim[sh]).mean():.0f}%, long {100*(EN[~sh]>lim[~sh]).mean():.0f}%")
