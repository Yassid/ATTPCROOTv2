#!/usr/bin/env python3
"""Parse allang.log into a summary table + the shape parameters vs angle."""
import re, sys
txt=open('allang.log',errors='replace').read()
blocks=re.split(r'###  PEAK SHAPE: ',txt)[1:]
out={}
for bl in blocks:
    shape='GAUSS' if 'PURE GAUSSIAN' in bl.split('\n')[0] else 'GT'
    cur=None
    for ln in bl.split('\n'):
        m=re.match(r'=+ (\S+)\s+theta_CM (\d+)-(\d+) =+',ln.strip())
        if m: cur=(m.group(1),int(m.group(2)),int(m.group(3))); continue
        m=re.match(r'\s+(free|constrained)\s+([\d.]+)\s+([\d.]+)\s+(.*)$',ln)
        if m and cur:
            case,L,c2,rest=m.group(1),float(m.group(2)),float(m.group(3)),m.group(4)
            nums=re.findall(r'-?\d+\.\d+|-?\d+',rest)
            sh=[]
            if case=='free' and len(nums)>=2: sh=[float(nums[0]),float(nums[1])]
            elif case=='constrained' and len(nums)>=1: sh=[None,float(nums[0])]
            out[(shape,cur[0],cur[1],case)]=dict(L=L,c2=c2,sh=sh)
codes=['ATTPCROOT','Spyral']
angs=sorted({k[2] for k in out})
print("\n================ SUMMARY: chi2/ndf, every angle ================")
print("   %-6s %-11s %18s %18s"%("ang","case","ATTPCROOT  G / GT","Spyral  G / GT"))
for a in angs:
    for case in ('free','constrained'):
        row="   %-6s %-11s"%("%d-%d"%(a,a+5),case)
        for c in codes:
            g=out.get(('GAUSS',c,a,case)); t=out.get(('GT',c,a,case))
            row+=" %8s /%8s"%("%.2f"%g['c2'] if g else "-","%.2f"%t['c2'] if t else "-")
        print(row)
print("\n================ shape parameters (FREE case) ================")
print("   %-6s %22s %22s"%("ang","ATTPCROOT  a / k","Spyral  a / k"))
for a in angs:
    row="   %-6s"%("%d-%d"%(a,a+5))
    for c in codes:
        g=out.get(('GT',c,a,'free'))
        row+=" %10s /%10s"%(("%.2f"%g['sh'][0] if g and g['sh'][0] is not None else "-"),
                            ("%.2f"%g['sh'][1] if g and len(g['sh'])>1 else "-"))
    print(row)
print("\n================ peak shape: does GT beat a Gaussian? ================")
print("   %-6s %-11s %14s %14s"%("ang","case","ATTPCROOT dL","Spyral dL"))
for a in angs:
    for case in ('free',):
        row="   %-6s %-11s"%("%d-%d"%(a,a+5),case)
        for c in codes:
            g=out.get(('GAUSS',c,a,case)); t=out.get(('GT',c,a,case))
            row+=" %14s"%("%+.1f"%(g['L']-t['L']) if g and t else "-")
        print(row)
print("\n   positive dL = the GT (tailed) peak fits BETTER than a pure Gaussian\n")
