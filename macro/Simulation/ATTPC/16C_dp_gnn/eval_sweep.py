import pandas as pd, numpy as np, glob, os
# crossing-event entry indices (global GET # - 34675) that SHOULD separate into >=2 clusters
CROSS=[60,225,309,331,354,407,593]
files=sorted(glob.glob("data/sw_*.csv"))
print(f"{'config':26s} {'clust%':>7s} {'cl/ev':>6s} {'multi%':>7s} {'frag%':>6s} {'sep(cross)':>11s}")
print("-"*72)
for f in files:
    d=pd.read_csv(f)
    clustered=(d.label>=0).mean()*100
    cpe=d[d.label>=0].groupby('event').label.nunique()
    cl_med=cpe.median() if len(cpe) else 0
    multi=(cpe>=2).mean()*100 if len(cpe) else 0
    frag=(cpe>=4).mean()*100 if len(cpe) else 0   # over-segmentation proxy
    sep=tot=0
    for e in CROSS:
        g=d[d.event==e]
        if len(g)<10: continue
        tot+=1
        if g[g.label>=0].label.nunique()>=2: sep+=1
    name=os.path.basename(f).replace('sw_','').replace('.csv','')
    print(f"{name:26s} {clustered:6.0f}% {cl_med:6.0f} {multi:6.0f}% {frag:5.0f}% {sep:>5d}/{tot}")
print("\nGOAL: high clust% (low noise), high sep(cross) (separates overlaps), LOW frag% (not over-segmenting)")
