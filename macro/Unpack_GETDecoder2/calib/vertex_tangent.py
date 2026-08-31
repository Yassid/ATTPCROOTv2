import numpy as np, csv, math, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cluster_tracks import cluster, fit_line, vertex_from_tracks
def load(ev):
    r=[x for x in csv.DictReader(open('/home/yassid/dec2014_calib/events_run0128.csv')) if int(x['event'])==ev]
    return (np.array([int(x['tb']) for x in r]),np.array([float(x['x']) for x in r]),
            np.array([float(x['y']) for x in r]),np.array([float(x['q']) for x in r]))

def analyse(ev,v,R):
    tb,x0,y0,q=load(ev); T=(tb-2.2)*160e-3
    P=np.column_stack([x0,y0,v*T*10]); lab,_,_=cluster(P,q)
    segs=[]
    for L in sorted(set(lab)-{-1}):
        m=lab==L; c,d=fit_line(P[m],q[m]); segs.append((m,c,d,np.ptp(P[m]@d),q[m].mean()))
    if len(segs)<3: return None
    beam=min(segs,key=lambda s:s[4]/(abs(s[2][2])+1e-9))
    arms=sorted([s for s in segs if s is not beam],key=lambda s:-s[3])[:2]
    vtx=vertex_from_tracks([(a[1],a[2]) for a in arms]+[(beam[1],beam[2])])
    ds=[]
    for m,c,d,L,qq in arms:
        Q=P[m]; sel=np.linalg.norm(Q-vtx,axis=1)<R      # tangent: only hits near the vertex
        if sel.sum()<4: return None
        c2,d2=fit_line(Q[sel],q[m][sel])
        ds.append(d2 if np.dot(d2,c2-vtx)>0 else -d2)
    return math.degrees(math.acos(np.clip(np.dot(ds[0],ds[1]),-1,1)))

print('opening angle using only hits within R of the vertex  (full-arm fit = last column)')
for ev in (279,258):
    print(f'\n--- event {ev} ---')
    hdr = '   v_D' + ''.join(('R=%d' % R).rjust(9) for R in (50,80,120,200)) + 'full'.rjust(9)
    print(hdr)
    for v in [1.8,2.0,2.25,2.5,2.75,3.0]:
        row=f'{v:6.2f}'
        for R in (50,80,120,200,10000):
            a=analyse(ev,v,R); row+=f'{a:9.2f}' if a else '      n/a'
        print(row)
