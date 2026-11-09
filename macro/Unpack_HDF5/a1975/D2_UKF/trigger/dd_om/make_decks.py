# FRESCO type=3 depth is Vso/2 of the global-potential value for spin 1: with it FRESCO = PtolemyCpp to 0.4 %.
# Usage: python3 make_decks.py DIR Ed  (DIR holds the print_potentials binary as p; writes DIR/f/inputs/*.nin and DIR/ptol_*.in)
import sys,re,subprocess
S=sys.argv[1]; Ed=sys.argv[2]
L=subprocess.check_output([S+'/p',Ed],text=True).splitlines()
def deck(n,Ed,p,f=1.0):
 V,r,a,W,rw,aw,Ws,rs,as_,Vso,rso,aso,rc=p
 return f""" d+16C elastic {n} Ed={Ed} surf factor {f}
NAMELIST
 &FRESCO hcm=0.05 rmatch=40.0 jtmin=0.0 jtmax=40.0 absend=-1.0
   thmin=1.0 thmax=180.0 thinc=1.0 iter=0 chans=1 smats=2 xstabl=0 elab(1)={Ed} /
 &PARTITION namep='DEUTERON' massp=2.0141 zp=1 nex=1 pwf=T namet='C-16' masst=16.0147 zt=6 qval=0.0/
 &STATES jp=1.0 ptyp=1 ep=0.0 cpot=1 jt=0.0 ptyt=1 et=0.0/
 &partition /
 &pot kp=1 type=0 p(1:3)=16.0 0.0 {rc}/
 &pot kp=1 type=1 p(1:7)={V} {r} {a} {W} {rw} {aw} 0.0/
 &pot kp=1 type=2 p(1:7)=0.0 0.0 0.0 {f*Ws} {rs} {as_} 0.0/
 &pot kp=1 type=3 p(1:3)={Vso/2} {rso} {aso}/
 &pot /
 &overlap /
 &coupling /
"""
for l in L:
 t=l.split(); n=t[0]
 if t[1]!='ok1': continue
 p=[float(x) for x in t[3::2]]
 tag=Ed.replace('.','p')
 open(f'{S}/f/inputs/dd_{n}_{tag}.nin','w').write(deck(n,Ed,p))
 if n=='D': open(f'{S}/f/inputs/dd_Dx4_{tag}.nin','w').write(deck(n,Ed,p,4.0))
 V,r,a,W,rw,aw,Ws,rs,as_,Vso,rso,aso,rc=p
 open(f'{S}/ptol_{n}_{tag}.in','w').write(f"HEADER d+16C {n}\nCHANNEL: D + 16C\nELAB={Ed}\nPARAMETERSET EL2\nV={V} R0={r} A={a} VI={W} RI0={rw} AI={aw}\nVSI={Ws} RSI0={rs} ASI={as_} VSO={Vso} RSO0={rso} ASO={aso} RC0={rc}\nELASTIC SCATTERING\n;\nRETURN\n")
