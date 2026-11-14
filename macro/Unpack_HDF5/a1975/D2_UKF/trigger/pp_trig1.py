# (p,p) elastic, which sets the (p,d) luminosity: does it show the same r_res-organised vz falloff? angdist_gs_pp.C defaults: Ebeam 188, theta-(360/4000)(KE-1.5), chi2<5, IC 1000-1350.
import uproot, numpy as np
U = 931.49401; M1, M2 = 16.0147013*U, 1.00727646688*U; M3, M4 = M2, M1
def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))
def kine(Kp, th, Ke):
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
t = uproot.open('/home/yassid/a1975_analysis/pp/data/pp_kin.root')['pk'].arrays(library='np'); ke, th, vz = t['ke'].astype(float), t['theta'].astype(float), t['vz'].astype(float)
thc = th - (360/4000.)*(ke - 1.5); ex, tcm = kine(188., np.radians(thc), ke)
sel = (t['chi2ndf'] < 5) & (t['ic'] >= 1000) & (t['ic'] <= 1350) & (np.abs(ex) < 1.0) & (ke > 0)
p = np.sqrt(ke**2 + 2*ke*M3); pb = np.sqrt(188.**2 + 2*188.*M1); r = np.radians(thc); rres = (1000 - vz)*np.tan(np.arctan2(p*np.sin(r), pb - p*np.cos(r)))
print('elastic candidates', sel.sum())
ZB = np.arange(0, 1000, 100); TB = [20, 30, 40, 50, 60, 70, 80, 100, 120]
print('elastic density per 100 mm by vz (cols 0..900), relative to vz 0-100 [expected r_res at that vz]')
for a, b in zip(TB[:-1], TB[1:]):
    w = sel & (tcm >= a) & (tcm < b); n = np.histogram(vz[w], np.r_[ZB, 1000])[0]; r0 = np.median(rres[w & (vz < 100)]) if (w & (vz < 100)).sum() else np.nan
    print(f'  tcm {a:3d}-{b:3d} n {w.sum():6d} ref {n[0]:5d}:', ' '.join(f'{v/max(n[0],1):4.2f}[{r0*(950-z0)/950:3.0f}]' for v, z0 in zip(n, ZB)))
