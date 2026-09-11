import numpy as np, glob, pickle, warnings, sys
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from cluster_tracks import fit_line, vertex_from_tracks
T0=-34.0; VD=2.142
s1,s2,ss=[],[],[]
for pth in sorted(glob.glob("nocut/P300_B0/cands_*.pkl"))[:15]:
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
        th=[]
        for k in ("a1","a2"):
            m=rec[k]; P=P3[m]; c,dd=segs[k]
            t=P@dd; tv=float(v@dd)
            far=P[np.argmax(np.abs(t-tv))]; w=far-v; nn=np.linalg.norm(w)
            if nn<1e-9: break
            th.append(np.degrees(np.arccos(np.clip((w/nn)@bd,-1,1))))
        if len(th)<2: continue
        th.sort(); s1.append(th[0]); s2.append(th[1]); ss.append(th[0]+th[1])
s1,s2,ss=map(np.array,(s1,s2,ss))
print(f"  {len(ss)} events (15 runs)")
print(f"  theta_1 median {np.median(s1):.1f}, theta_2 median {np.median(s2):.1f}")
print(f"  theta_1 + theta_2: median {np.median(ss):.1f} deg   (must be 90 for elastic equal masses)")
print(f"     16-84%: {np.percentile(ss,16):.1f} .. {np.percentile(ss,84):.1f}")
print(f"     within 90+-5 deg: {100*(np.abs(ss-90)<5).mean():.0f}%")
h,ed=np.histogram(ss,bins=np.arange(40,140,5))
for i in range(len(h)):
    if h[i]>0: print(f"    {ed[i]:3.0f}-{ed[i+1]:3.0f} {'#'*int(50*h[i]/h.max())} {h[i]}")
