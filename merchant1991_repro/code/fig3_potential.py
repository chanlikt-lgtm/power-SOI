"""Fig. 3 -- equipotential contour map for the SOI device of Fig. 2
(t_s=0.1 um, t_ox=2 um, L=70 um), at the breakdown voltage.

Unlike Nakagawa 1990 (which needs a numerical 2-D Poisson solve because its
doping profile is uniform, not tailored), this paper's whole point is a closed
form potential (eqs. 5, 7, 9, 15): the lateral doping N(x) is *designed* so that

    phi(x,y) = (V x/W) [1 - y^2/(2 t t_s)]                 0<=x<=W, 0<=y<=t_s
    phi(x,y) = phi(x,t_s) - (3 V x /(t W)) (y - t_s)        0<=x<=W, t_s<y<=t_s+t_ox

for the depleted region 0<=x<=W. Beyond the depletion edge (W<x<=L) the SOI is
assumed quasi-neutral at the applied bias V, and the buried oxide there just
carries a plain 1-D vertical drop V -> 0 into the grounded substrate. We plot
this closed-form phi directly -- it is what the paper's optimized doping
profile is designed to produce, and is the same field distribution used to
derive breakdown_voltage() in model.py.

The design width W is found by inverting the model so that the breakdown
voltage at (t_s, t_ox, W) equals the paper's stated Fig.3 breakdown of 750 V
(a specific L=70um device); this pins down the one free design parameter.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from model import UM, SI, EPS_OX, t_eff, breakdown_voltage

OUTDIR = "../figs"


def find_design_W(t_s, t_ox, V_target, W_lo=5 * UM, W_hi=69 * UM):
    """Invert breakdown_voltage(t_s, t_ox, W) = V_target for W."""
    f = lambda W: breakdown_voltage(t_s, t_ox, W) - V_target
    return brentq(f, W_lo, W_hi, xtol=1e-9)


def phi(x, y, V, W, t_s, t_ox):
    t = t_eff(t_s, t_ox)
    ratio = SI.eps_r / EPS_OX   # = 3
    x = np.minimum(x, W)                      # clamp: flat beyond the depletion edge
    in_si = y <= t_s
    phi_si = (V * x / W) * (1.0 - y**2 / (2.0 * t * t_s))
    phi_ts = (V * x / W) * (1.0 - t_s / (2.0 * t))
    Ey_ox = ratio * V * x / (t * W)
    phi_ox = phi_ts - Ey_ox * (y - t_s)
    return np.where(in_si, phi_si, phi_ox)


def main():
    t_s, t_ox, L, V = 0.1 * UM, 2.0 * UM, 70.0 * UM, 750.0
    W = find_design_W(t_s, t_ox, V)
    print(f"Fig.3  design width W = {W/UM:.1f} um  (must be < L=70 um)")

    nx, ny = 700, 400
    x = np.linspace(0, L, nx)
    y = np.linspace(0, t_s + t_ox, ny)
    X, Y = np.meshgrid(x, y)

    # For x > W the SOI is quasi-neutral at V; the oxide there drops linearly
    # from V to 0 with no lateral variation.
    Phi = np.where(
        X <= W,
        phi(X, Y, V, W, t_s, t_ox),
        np.where(Y <= t_s, V, V * (1.0 - (Y - t_s) / t_ox)),
    )

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    levels = np.array([0, 79.5, 159, 239, 319, 399, 479, 559, 639, 700, 750])
    cs = ax.contour(X / UM, Y / UM, Phi, levels=levels, colors="k", linewidths=0.9)
    ax.clabel(cs, fmt="%.0f", fontsize=7)
    ax.axhline(t_s / UM, color="gray", lw=0.6, ls=":")
    ax.axvline(W / UM, color="gray", lw=0.6, ls="--")
    ax.invert_yaxis()
    ax.set_xlabel("lateral position $x$  ($\\mu$m)")
    ax.set_ylabel("depth $y$  ($\\mu$m)")
    ax.text(2, (t_s + t_ox) / UM - 0.25, "Buried Oxide", fontsize=8)
    ax.set_title(f"Fig.3  Equipotential contours at V={V:.0f} V "
                 f"($t_s$=0.1 $\\mu$m, $t_{{ox}}$=2 $\\mu$m, $L$=70 $\\mu$m)")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig3_potential.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
