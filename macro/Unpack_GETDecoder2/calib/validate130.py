import math,pickle,sys,warnings,csv
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from batch_scan import load_all, analyse
evs,groups=load_all('/home/yassid/dec2014_calib/all_run0130.csv')
print(f'loaded {len(evs)} events',flush=True)
out=[]
for i,h in enumerate(groups):
    if i%2000==0: print(f'  {i}/{len(evs)} kept={len(out)}',flush=True)
    try: r=analyse(h)
    except Exception: r=None
    if r: out.append(r)
print(f'selected {len(out)} candidates')
pickle.dump(out,open('/home/yassid/dec2014_calib/cands130.pkl','wb'))
