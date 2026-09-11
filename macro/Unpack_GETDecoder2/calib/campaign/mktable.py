import numpy as np
fHe,fCO2=0.9,0.1
M=fHe*4.003+fCO2*44.01
wHe,wCO2=fHe*4.003/M, fCO2*44.01/M
ZA=wHe*(2/4.003)+wCO2*(22/44.01)
I=np.exp((wHe*(2/4.003)*np.log(41.8)+wCO2*(22/44.01)*np.log(85.0))/ZA)
Mal,me=3727.379,0.510999
rho150=0.145e-3*150.0/330.0
EMIN=0.20
def dEdx(E):
    b2=2*E/Mal
    v=0.307075*4*ZA/b2*np.log(2*me*1e6*b2/I)*rho150*0.1
    return v if v>1e-6 else 1e-6
d=[0.0]; e=[11.40]; E=11.40; step=0.5
while E>EMIN and len(d)<20000:
    E-=dEdx(E)*step
    if E<=EMIN: break
    d.append(d[-1]+step); e.append(E)
d=np.array(d); e=np.array(e)
with open("eloss_alpha_heco2_11p4.txt","w") as fh:
    fh.write("# distance_mm  energy_MeV   alpha in He:CO2 90/10, 150 torr, E0=11.40 MeV\n")
    for a,b in zip(d,e): fh.write(f"{a:.1f} {b:.5f}\n")
print(f"  {len(d)} points, range {d[-1]:.0f} mm @150 torr = {d[-1]/(299.5/150):.0f} mm @299.5 torr")
i=int(np.argmin(np.abs(e-7.80)))
print(f"  validation: hits 7.80 MeV at {d[i]:.0f} mm, residual range {d[-1]-d[i]:.0f} mm")
print(f"     shipped 7.80 MeV table total = 1100 mm -> agreement {100*(d[-1]-d[i])/1100:.1f}%")
