#!/usr/bin/env python3
"""Step 5 -- the geometric acceptance.  This is the heart of the analysis.

What fraction of reactions at a given (E_cm, theta_cm) do we actually detect?  For this
measurement that is almost purely GEOMETRY, so it can be computed directly -- no Geant4.

The chain is fully determined once (E_cm, theta_cm) is fixed:

    E_beam at the vertex = 2 * E_cm            (equal masses)
    depth of the vertex  = how far the beam must travel to be left with that energy
    theta_lab            = theta_cm / 2        (and 90 - theta_cm/2 for the partner)
    E1, E2               = E_beam cos^2(theta_lab1), E_beam cos^2(theta_lab2)
    R1, R2               = ranges of those energies at this pressure

The only freedom left is the azimuth of the reaction plane.  So: scan the azimuth, and
for each one ask whether BOTH arms stay inside the active volume and are long enough to
be found by the clustering.  The fraction that pass is the acceptance.

Makes  out/05_acceptance.png  and  out/acceptance.npz

    python3 05_acceptance.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config
from lib.eloss import RangeTable

tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

# --- beam geometry ---------------------------------------------------------------
# The tilt pivots at the entrance window, so the beam starts ON the detector axis at
# z = DRIFT_LENGTH and diverges as it goes deeper (towards the pad plane at z = 0).
TILT   = np.radians(config.TILT)
PHI_B  = np.radians(13.0)              # azimuth of the tilt, measured from the data
BEAM_DIR = np.array([np.sin(TILT)*np.cos(PHI_B), np.sin(TILT)*np.sin(PHI_B), -np.cos(TILT)])
ENTRANCE = np.array([0.0, 0.0, config.DRIFT_LENGTH])

MIN_TRACK_MM = 5.0        # clustering floor: shorter than this and no cluster forms
N_AZIMUTH    = 180        # azimuthal samples per cell

# an orthonormal frame around the beam
_u = np.array([0.0, 0.0, 1.0])
E1_ = np.cross(BEAM_DIR, _u); E1_ /= np.linalg.norm(E1_)
E2_ = np.cross(BEAM_DIR, E1_)


def arm_direction(theta_deg, psi):
    """Unit vector at polar angle theta from the beam, azimuth psi around it."""
    t = np.radians(theta_deg)
    return (np.cos(t) * BEAM_DIR
            + np.sin(t) * (np.cos(psi) * E1_ + np.sin(psi) * E2_))


def contained(point):
    """Is this point inside the active volume?  (cylinder: pad radius and drift length)"""
    r = np.hypot(point[..., 0], point[..., 1])
    z = point[..., 2]
    return (r < config.R_PAD) & (z > 0) & (z < config.DRIFT_LENGTH)


def acceptance(E_cm, theta_cm):
    """Fraction of azimuths for which both arms are contained and long enough."""
    E_beam = 2.0 * E_cm
    if E_beam > config.E_BEAM:              # cannot react with more energy than it has
        return 0.0
    depth = float(tab.depth_of_energy(E_beam))
    vertex = ENTRANCE + depth * BEAM_DIR
    if not contained(vertex):
        return 0.0

    th1 = theta_cm / 2.0                    # equal masses: theta_lab = theta_cm / 2
    th2 = 90.0 - th1                        # the partner, since theta1 + theta2 = 90
    e1 = E_beam * np.cos(np.radians(th1)) ** 2
    e2 = E_beam * np.cos(np.radians(th2)) ** 2
    r1 = float(tab.range_of_energy(e1))
    r2 = float(tab.range_of_energy(e2))
    if r1 < MIN_TRACK_MM or r2 < MIN_TRACK_MM:
        return 0.0

    psi = np.linspace(0, 2 * np.pi, N_AZIMUTH, endpoint=False)
    # the two arms are back-to-back in azimuth about the beam
    d1 = np.stack([arm_direction(th1, p) for p in psi])
    d2 = np.stack([arm_direction(th2, p + np.pi) for p in psi])
    ok = contained(vertex + r1 * d1) & contained(vertex + r2 * d2)
    return float(ok.mean())


E_EDGES = np.arange(0.2, config.E_BEAM / 2 + 0.2, 0.1)
T_EDGES = np.arange(2.0, 180.0, 2.0)
Ec = 0.5 * (E_EDGES[1:] + E_EDGES[:-1])
Tc = 0.5 * (T_EDGES[1:] + T_EDGES[:-1])
EPS = np.zeros((len(Ec), len(Tc)))
for i, E in enumerate(Ec):
    for j, T in enumerate(Tc):
        EPS[i, j] = acceptance(E, T)
print(f"  acceptance grid {EPS.shape}, max {EPS.max():.2f}, "
      f"non-zero in {100*(EPS>0).mean():.0f}% of cells")

np.savez_compressed(os.path.join(config.OUT_DIR, "acceptance.npz"),
                    E=Ec, T=Tc, eps=EPS)

fig, ax = plt.subplots(1, 2, figsize=(14, 5.0))
im = ax[0].pcolormesh(Ec, Tc, EPS.T, cmap="magma", vmin=0, vmax=1, shading="auto")
fig.colorbar(im, ax=ax[0], label=r"$\varepsilon$")
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\theta_{cm}$ [deg]")
ax[0].set_title("geometric acceptance")
for T in (90, 70, 50, 30):
    j = np.argmin(np.abs(Tc - T))
    ax[1].plot(Ec, EPS[:, j], lw=1.8, label=rf"$\theta_{{cm}}$ = {T}$^\circ$")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel(r"$\varepsilon$")
ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=9); ax[1].grid(alpha=.25)
ax[1].set_title("acceptance vs energy at fixed angle")
fig.tight_layout()
out = os.path.join(config.OUT_DIR, "05_acceptance.png"); fig.savefig(out, dpi=140)

j90 = np.argmin(np.abs(Tc - 90))
print("  at theta_cm = 90 deg:")
for E in (1.0, 2.0, 3.0, 4.0, 5.0):
    i = np.argmin(np.abs(Ec - E))
    print(f"    E_cm {E:.1f} MeV -> eps = {EPS[i, j90]:.2f}")
print(f"  wrote {out} and acceptance.npz")
