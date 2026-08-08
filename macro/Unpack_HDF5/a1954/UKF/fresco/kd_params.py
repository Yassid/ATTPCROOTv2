#!/usr/bin/env python3
"""Koning-Delaroche (2003) global optical-model parameters for protons, and the FRESCO
inputs for p + 14C at the a1954 energy.

Why compute them here instead of copying the a2091 15C input: that input's real depth
(V_V = 46.89 at 13 MeV) corresponds to the NEUTRON asymmetry term, v1 = 59.30 - 21(N-Z)/A,
and carries no Coulomb correction.  For protons KD has v1 = 59.30 + 21(N-Z)/A plus the
V_C term, which for p+15C gives ~55 MeV, not 46.9.  Reproducing that here would bias the
diffraction pattern we are trying to compare against, so the formulas are implemented
explicitly and the discrepancy is printed.

Energy: the data is 14C at 161 MeV on a proton target.  FRESCO is set up in normal
kinematics (p on 14C), so the lab proton energy is the one giving the same E_cm:
    E_cm  = E_lab(14C) * m_p / (m_p + m_14C)
    E_lab(p) = E_cm * (m_p + m_14C) / m_14C
"""

import math

A, Z = 14, 6
N = A - Z
A13 = A ** (1.0 / 3.0)

M_P = 1.00727646
M_C14 = 14.0032420


def ecm_from_inverse(elab_beam):
    return elab_beam * M_P / (M_P + M_C14)


def elab_proton(ecm):
    return ecm * (M_P + M_C14) / M_C14


def kd_proton(E):
    """KD03 proton parameters at lab energy E (MeV). Returns a dict."""
    Ef = -8.4075 + 0.01378 * A
    dE = E - Ef

    rC = 1.198 + 0.697 * A ** (-2.0 / 3.0) + 12.994 * A ** (-5.0 / 3.0)
    VC = 1.73 / rC * Z * A ** (-1.0 / 3.0)

    v1 = 59.30 + 21.0 * (N - Z) / A - 0.024 * A
    v2 = 0.007067 + 4.23e-6 * A
    v3 = 1.729e-5 + 1.136e-8 * A
    v4 = 7.0e-9

    VV = v1 * (1 - v2 * dE + v3 * dE**2 - v4 * dE**3) \
        + VC * v1 * (v2 - 2 * v3 * dE + 3 * v4 * dE**2)
    rV = 1.3039 - 0.4054 * A ** (-1.0 / 3.0)
    aV = 0.6778 - 1.487e-4 * A

    w1 = 14.667 + 0.009629 * A
    w2 = 73.55 + 0.0795 * A
    WV = w1 * dE**2 / (dE**2 + w2**2)

    d1 = 16.0 + 16.0 * (N - Z) / A
    d2 = 0.0180 + 0.003802 / (1 + math.exp((A - 156.0) / 8.0))
    d3 = 11.5
    WD = d1 * dE**2 / (dE**2 + d3**2) * math.exp(-d2 * dE)
    rD = 1.3424 - 0.01585 * A13
    aD = 0.5187 + 5.205e-4 * A

    vso1, vso2 = 5.922 + 0.0030 * A, 0.0040
    VSO = vso1 * math.exp(-vso2 * dE)
    rSO = 1.1854 - 0.647 * A ** (-1.0 / 3.0)
    aSO = 0.59

    wso1, wso2 = -3.1, 160.0
    WSO = wso1 * dE**2 / (dE**2 + wso2**2)

    return dict(Ef=Ef, dE=dE, rC=rC, VC=VC, VV=VV, rV=rV, aV=aV, WV=WV,
                WD=WD, rD=rD, aD=aD, VSO=VSO, rSO=rSO, aSO=aSO, WSO=WSO)


HEADER = """{title}
NAMELIST
 &FRESCO  hcm= 0.05 rmatch= 30.000
     jtmin=  0.0 jtmax=   60.0 absend=  0.00001
     thmin=  1.00 thmax=180.00 thinc=  1.00
     iblock= {iblock}
     chans= 1 smats= 2 xstabl= 1
     elab(1)= {elab:.4f} pel=1 exl=1 lab=1 lin=1 lex=1 /

 &PARTITION namep='Proton  ' massp= 1.00728 zp= 1 nex= {nex}
            namet=' 14C    ' masst=14.00324 zt= 6 qval= 0.0000 /
"""

POTS_ELASTIC = """ &partition /

 &pot kp= 1 type= 0 p(1:3)= {A:6.3f}  0.0000  {rC:6.4f} /
 &pot kp= 1 type= 1 p(1:7)= {VV:6.3f}  {rV:6.4f}  {aV:6.4f}  {WV:6.3f}  {rV:6.4f}  {aV:6.4f}  0.0000 /
 &pot kp= 1 type= 2 p(1:7)=  0.000  0.0000  0.0000  {WD:6.3f}  {rD:6.4f}  {aD:6.4f}  0.0000 /
 &pot kp= 1 type= 3 p(1:3)= {VSO:6.3f}  {rSO:6.4f}  {aSO:6.4f} /
 &pot /

 &overlap /
 &coupling /
"""

# deformed version: a type=11 deformation follows each deformed potential shape.
POTS_INEL = """ &partition /

 &pot kp= 1 type= 0 p(1:3)= {A:6.3f}  0.0000  {rC:6.4f} /
 &pot kp= 1 type= 1 p(1:7)= {VV:6.3f}  {rV:6.4f}  {aV:6.4f}  {WV:6.3f}  {rV:6.4f}  {aV:6.4f}  0.0000 /
 &pot kp= 1 type=11 p(1:5)= {d1:6.4f}  {d2:6.4f}  {d3:6.4f}  {d4:6.4f}  0.0000 /
 &pot kp= 1 type= 2 p(1:7)=  0.000  0.0000  0.0000  {WD:6.3f}  {rD:6.4f}  {aD:6.4f}  0.0000 /
 &pot kp= 1 type=11 p(1:5)= {d1:6.4f}  {d2:6.4f}  {d3:6.4f}  {d4:6.4f}  0.0000 /
 &pot kp= 1 type= 3 p(1:3)= {VSO:6.3f}  {rSO:6.4f}  {aSO:6.4f} /
 &pot /

 &overlap /
 &coupling /
"""


def write_elastic(path, elab, p, title):
    with open(path, "w") as f:
        f.write(HEADER.format(title=title, iblock=1, elab=elab, nex=1))
        f.write(" &STATES jp= 0.5 ptyp= 1 ep= 0.0000 cpot= 1 jt= 0.0 ptyt= 1 et= 0.0000 /\n")
        f.write(POTS_ELASTIC.format(A=float(A), **p))


def write_inelastic(path, elab, p, title, L, jt, ptyt, ex, betaR=0.281):
    """One excited state reached by a single multipole L (deformation length betaR, fm)."""
    d = [0.0, 0.0, 0.0, 0.0]
    d[L - 1] = betaR
    with open(path, "w") as f:
        f.write(HEADER.format(title=title, iblock=2, elab=elab, nex=2))
        f.write(" &STATES jp= 0.5 ptyp= 1 ep= 0.0000 cpot= 1 jt= 0.0 ptyt= 1 et= 0.0000 /\n")
        f.write(" &STATES                             cpot= 1 jt= %3.1f ptyt=%2d et= %6.4f /\n"
                % (jt, ptyt, ex))
        f.write(POTS_INEL.format(A=float(A), d1=d[0], d2=d[1], d3=d[2], d4=d[3], **p))


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    ind = os.path.join(here, "inputs")
    os.makedirs(ind, exist_ok=True)

    for ebeam, tag in ((161.0, "161"), (155.0, "155")):
        ecm = ecm_from_inverse(ebeam)
        ep = elab_proton(ecm)
        p = kd_proton(ep)
        print("\n=== 14C beam %.1f MeV -> E_cm %.3f MeV -> E_lab(p) %.3f MeV ===" % (ebeam, ecm, ep))
        for k in ("Ef", "rC", "VC", "VV", "rV", "aV", "WV", "WD", "rD", "aD", "VSO", "rSO", "aSO", "WSO"):
            print("   %-4s %9.4f" % (k, p[k]))

        write_elastic(os.path.join(ind, "p14C_el_%s.nin" % tag), ep, p,
                      "p + 14C elastic, KD03 proton OMP, Elab(p)=%.3f MeV" % ep)
        # 14C levels in the analysis window, each by its natural multipole
        for L, jt, ptyt, ex, nm in ((1, 1.0, -1, 6.0940, "6094_L1"),
                                    (3, 3.0, -1, 6.7280, "6728_L3"),
                                    (2, 2.0, +1, 7.0120, "7012_L2")):
            write_inelastic(os.path.join(ind, "p14C_inel_%s_%s.nin" % (tag, nm)), ep, p,
                            "p + 14C -> %s MeV (L=%d), KD03, Elab(p)=%.3f MeV" % (nm.split("_")[0], L, ep),
                            L, jt, ptyt, ex)
    print("\ninputs written to", ind)
