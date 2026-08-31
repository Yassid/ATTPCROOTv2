import numpy as np, csv, math, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cluster_tracks import cluster, fit_line, vertex_from_tracks

def load(ev):
    r=[x for x in csv.DictReader(open('/home/yassid/dec2014_calib/events_run0128.csv')) if int(x['event'])==ev]
    return (np.array([int(x['tb']) for x in r]), np.array([float(x['x']) for x in r]),
            np.array([float(x['y']) for x in r]), np.array([float(x['q']) for x in r]))

TBNS=160.0; TB0=2.2; B=0.5691
def dvec(v,E,tilt):
    t=math.radians(tilt); wt=(B/E)*(v*1e4); f=v/(1+wt*wt)
    return f*wt*math.sin(t), f*wt*wt*math.sin(t)*math.cos(t), f*(1+wt*wt*math.cos(t)**2)

def analyse(ev, v, E=None, tilt=7.0):
    tb,x0,y0,q = load(ev)
    T=(tb-TB0)*TBNS*1e-3
    if E is None:  X,Y,Z = x0, y0, v*T*10
    else:
        vx,vy,vz = dvec(v,E,tilt); X,Y,Z = x0-vx*T*10, y0-vy*T*10, vz*T*10
    P=np.column_stack([X,Y,Z])
    lab,_,_ = cluster(P,q)
    segs=[]
    for L in sorted(set(lab)-{-1}):
        m=lab==L; c,d=fit_line(P[m],q[m]); segs.append((c,d,m.sum(),np.ptp(P[m]@d),q[m].mean()))
    if len(segs)<3: return None
    segs.sort(key=lambda s:-s[3])
    # beam = the segment most aligned with z AND lowest mean charge
    beam=min(segs,key=lambda s: s[4]/ (abs(s[1][2])+1e-9))
    arms=[s for s in segs if s is not beam]
    arms.sort(key=lambda s:-s[3]); arms=arms[:2]
    vtx=vertex_from_tracks([(a[0],a[1]) for a in arms]+[(beam[0],beam[1])])
    ds=[]
    for c,d,n,L,qq in arms:
        dd = d if np.dot(d, c-vtx)>0 else -d
        ds.append(dd)
    op=math.degrees(math.acos(np.clip(np.dot(ds[0],ds[1]),-1,1)))
    bt=math.degrees(math.acos(np.clip(abs(beam[1][2]),-1,1)))
    return op, bt, vtx, [a[3] for a in arms]

print("=== event 279 (both arms long) ===")
print(f"{'v_D':>6} {'opening':>9} {'beamtilt':>9}   vertex(x,y,z)")
for v in [1.8,2.0,2.25,2.5,2.75,3.0,3.5]:
    r=analyse(279,v)
    if r: print(f"{v:6.2f} {r[0]:9.2f} {r[1]:9.2f}   ({r[2][0]:7.1f},{r[2][1]:7.1f},{r[2][2]:7.1f})  armlen={r[3][0]:.0f}/{r[3][1]:.0f}mm")
print()
print("=== event 258 (short recoil, for comparison) ===")
for v in [2.0,2.25,2.5,2.75,3.0]:
    r=analyse(258,v)
    if r: print(f"{v:6.2f} {r[0]:9.2f} {r[1]:9.2f}   ({r[2][0]:7.1f},{r[2][1]:7.1f},{r[2][2]:7.1f})")
