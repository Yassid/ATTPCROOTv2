#!/usr/bin/env python3
"""
Colored remake of Figure 3 from:

  S. Julier & J. Uhlmann, "A New Extension of the Kalman Filter to Nonlinear
  Systems", 1997.

It compares:
  - "True" mean/covariance via Monte Carlo
  - Linearization / EKF-style transform
  - Unscented Transform (UT)

The example is the polar → Cartesian conversion of a noisy range/bearing
measurement to (x, y).
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


# ---------------------------------------------------------------------------
# Problem setup: range–bearing measurement of a target at (0, 1)
# ---------------------------------------------------------------------------

# True Cartesian location
x_true = np.array([0.0, 1.0])

# Mean measurement in polar coordinates (r, theta)
mu_polar = np.array([1.0, math.pi / 2.0])   # r = 1 m, theta = 90 deg

# Sensor noise: good range, poor bearing
sigma_r = 0.02                              # 2 cm
sigma_theta = math.radians(15.0)            # 15 degrees
P_polar = np.diag([sigma_r**2, sigma_theta**2])


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def polar_to_cartesian(z: np.ndarray) -> np.ndarray:
    """Nonlinear transform f: (r, theta) -> (x, y)."""
    r, theta = z
    return np.array([r * np.cos(theta), r * np.sin(theta)])


def jacobian_polar_to_cartesian(r: float, theta: float) -> np.ndarray:
    """Jacobian of f at (r, theta) for the EKF-style linearization."""
    return np.array([
        [np.cos(theta), -r * np.sin(theta)],
        [np.sin(theta),  r * np.cos(theta)],
    ])


def unscented_transform(mu: np.ndarray,
                        P: np.ndarray,
                        f,
                        kappa: float | None = None):
    """
    Basic Julier & Uhlmann (1997) unscented transform for mean/covariance.

    mu : (n,) mean
    P  : (n,n) covariance
    f  : nonlinear mapping R^n -> R^m
    kappa : scaling parameter; default enforces n + kappa = 3.
    """
    mu = np.atleast_1d(mu)
    n = mu.size

    if kappa is None:
        # Heuristic from the paper: choose n + kappa = 3
        kappa = 3.0 - n

    lam = n + kappa
    if lam <= 0:
        raise ValueError("Require n + kappa > 0")

    # Sigma points (2n+1) with scaling
    sqrtP = np.linalg.cholesky(lam * P)
    sigma = np.zeros((2 * n + 1, n))
    sigma[0] = mu
    for i in range(n):
        sigma[1 + i]     = mu + sqrtP[:, i]
        sigma[1 + n + i] = mu - sqrtP[:, i]

    # Weights
    Wm = np.full(2 * n + 1, 1.0 / (2.0 * lam))
    Wc = np.full(2 * n + 1, 1.0 / (2.0 * lam))
    Wm[0] = kappa / lam
    Wc[0] = kappa / lam

    # Propagate through nonlinear map
    Y = np.array([f(sp) for sp in sigma])   # (2n+1, m)
    y_mean = np.sum(Wm[:, None] * Y, axis=0)

    diff = Y - y_mean
    m = y_mean.size
    y_cov = np.zeros((m, m))
    for i in range(2 * n + 1):
        y_cov += Wc[i] * np.outer(diff[i], diff[i])

    return y_mean, y_cov, sigma, Y


def plot_cov_ellipse(mean: np.ndarray,
                     cov: np.ndarray,
                     ax: plt.Axes,
                     n_std: float = 1.0,
                     **kwargs) -> Ellipse:
    """
    Plot an n_std covariance ellipse for a 2×2 covariance matrix.
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    # Ellipse geometry
    width, height = 2.0 * n_std * np.sqrt(vals)
    angle = math.degrees(math.atan2(vecs[1, 0], vecs[0, 0]))

    ell = Ellipse(xy=mean, width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ell)
    return ell


# ---------------------------------------------------------------------------
# Compute mean/covariance with three different methods
# ---------------------------------------------------------------------------

print("Starting MC")
# 1. "True" via Monte Carlo
N_SAMPLES = 500_000  # increase if you want even smoother estimates
samples_polar = np.random.multivariate_normal(mu_polar, P_polar, size=N_SAMPLES)
r = samples_polar[:, 0]
theta = samples_polar[:, 1]
print("Mean theta:", np.degrees(theta.mean()))
print("Mean r:", r.mean())
samples_cart = np.column_stack((r * np.cos(theta), r * np.sin(theta)))

mu_true = samples_cart.mean(axis=0)
print("True mean:", mu_true)
P_true = np.cov(samples_cart, rowvar=False)

print("Running linear")

# 2. Linearization / EKF
mu_ekf = polar_to_cartesian(mu_polar)
J = jacobian_polar_to_cartesian(*mu_polar)
P_ekf = J @ P_polar @ J.T

print("Running UT")
# 3. Unscented transform
mu_ut, P_ut, sigma_pts, sigma_cart = unscented_transform(mu_polar, P_polar,
                                                        polar_to_cartesian)
print("Finished computations")

# ---------------------------------------------------------------------------
# Make the figure
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6, 6))

# Means
ax.scatter(mu_true[0], mu_true[1], color="C0", marker="x", s=80,
           label="True mean (Monte Carlo)")
ax.scatter(mu_ekf[0], mu_ekf[1], color="C1", marker="o", s=60,
           label="Linearisation / EKF mean")
ax.scatter(mu_ut[0], mu_ut[1], color="C2", marker="^", s=70,
           label="Unscented transform mean")

# 1σ covariance ellipses
plot_cov_ellipse(mu_true, P_true, ax,
                 edgecolor="C0", facecolor="none", linewidth=2)
plot_cov_ellipse(mu_ekf, P_ekf, ax,
                 edgecolor="C1", linestyle="--", facecolor="none", linewidth=2)
plot_cov_ellipse(mu_ut, P_ut, ax,
                 edgecolor="C2", linestyle="-.", facecolor="none", linewidth=2)

# Optionally: show UT sigma points in Cartesian space
print(sigma_cart.shape, "sigma points in Cartesian")
i = 2
#ax.scatter(sigma_cart[:, 0], sigma_cart[:, 1],
#           color="b", alpha=1, s=20, label="Sigma points")
print("Sigma point", i, "in polar:", sigma_pts[i])
print("Sigma point", i, "in Cartesian:", sigma_cart[i])


# Formatting
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax.set_aspect("auto", adjustable="box")

# Match the visual window of the original figures
#ax.set_xlim(-0.4, 0.4)
#x.set_ylim(0.9, 1.04)

# De-duplicate legend entries (sigma points share label color)
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc="lower center")

plt.tight_layout()
plt.show()

# Optionally save straight to file:
fig.savefig("julier97_fig3_color.png", dpi=300)
