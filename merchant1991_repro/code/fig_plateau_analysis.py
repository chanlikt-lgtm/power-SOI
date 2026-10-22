"""The "plateau problem": why the minimum of BV(t_s) (Fig. 4's "U" curve) is a
broad, flat plateau for Si but a sharp "V" for 4H-SiC/GaN.

Since V_BR(t_s,t_ox,mat) = critical_field(t_s,mat) * t_eff(t_s,t_ox,mat) exactly
(model.critical_field, model.t_eff -- verified numerically in model.py's
docstring), and critical_field(t_s) does NOT depend on t_ox, the location and
shape of the minimum is controlled by how fast E_c(t_s) falls as t_s grows,
relative to how fast the geometric factor t_eff grows. Taking logs,

    d ln(V_BR)/d ln(t_s) = d ln(E_c)/d ln(t_s) + d ln(t_eff)/d ln(t_s)
                          = d ln(E_c)/d ln(t_s) + (t_s/2)/(t_s/2 + r*t_ox)

so the minimum (d V_BR/d t_s = 0) occurs where

    -d ln(E_c)/d ln(t_s)  =  (t_s/2) / (t_s/2 + r*t_ox),     r = eps_r/eps_ox   (*)

The left-hand side depends only on t_s and the material's ionization prefactor
a_ion (not on t_ox, b_ion, or eps_r -- see model.g_star); the right-hand side
depends on t_s, t_ox and the material's permittivity ratio r (not on a_ion or
b_ion). Where the two intersect is (approximately) the minimum of BV(t_s).

This script produces:
  fig_dbv_dts.png       raw derivative dV_BR/dt_s vs t_s, Si/4H-SiC/GaN, one
                        panel per t_ox in {2,4,6} um -- the direct answer to
                        "plot the derivative vs t_s for different material
                        and t_ox".
  fig_log_slope.png     the (*) LHS/RHS reproduction: LHS is material-only
                        (3 solid curves), RHS is material+t_ox (dashed
                        curves); marked intersections are cross-checked
                        against a direct numerical minimization of BV(t_s).
  fig_plateau_width.png "plateau width" (the t_s-interval around the minimum
                        over which |dV_BR/dt_s / BV| stays below 1%/um) vs
                        t_ox, per material -- quantifies why Si's bottom is
                        much broader than 4H-SiC's/GaN's.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, minimize_scalar

from model import UM, EPS_OX, SI, SIC_4H, GAN, breakdown_voltage, critical_field

OUTDIR = "../figs"
MATERIALS = (SI, SIC_4H, GAN)
COLORS = {"Si": "tab:blue", "4H-SiC": "tab:orange", "GaN": "tab:green"}
TOX_LIST = (2.0, 4.0, 6.0)   # um


def dbv_dts(ts_um, tox_um, mat, rel_step=1e-3):
    """Central finite difference of breakdown_voltage w.r.t. t_s (V per um)."""
    ts = ts_um * UM
    h = max(ts_um * rel_step, 1e-6) * UM
    v_plus = breakdown_voltage(ts + h, tox_um * UM, W=None, mat=mat)
    v_minus = breakdown_voltage(ts - h, tox_um * UM, W=None, mat=mat)
    return (v_plus - v_minus) / (2 * h) * UM   # V/um


def dlnEc_dlnts(ts_um, mat, rel_step=1e-3):
    """Central finite difference of ln(E_c) w.r.t. ln(t_s) (dimensionless)."""
    ts = ts_um * UM
    h = max(ts_um * rel_step, 1e-6) * UM
    ec_plus = critical_field(ts + h, mat)
    ec_minus = critical_field(ts - h, mat)
    return (np.log(ec_plus) - np.log(ec_minus)) / (np.log(ts + h) - np.log(ts - h))


def rhs_slope(ts_um, tox_um, mat):
    """(t_s/2) / (t_s/2 + r*t_ox), r = eps_r/eps_ox."""
    ts, tox = ts_um * UM, tox_um * UM
    r = mat.eps_r / EPS_OX
    return (ts / 2.0) / (ts / 2.0 + r * tox)


# ----------------------------------------------------------------------------
def fig_dbv_dts():
    ts_um = np.logspace(np.log10(0.02), np.log10(50), 160)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4), sharey=False)
    for ax, tox_um in zip(axes, TOX_LIST):
        for mat in MATERIALS:
            d = np.array([dbv_dts(t, tox_um, mat) for t in ts_um])
            ax.plot(ts_um, d, color=COLORS[mat.name], lw=1.8, label=mat.name)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=10)
        ax.set_xlabel("$t_s$  ($\\mu$m)")
        ax.set_title(f"$t_{{ox}}$={tox_um:.0f} $\\mu$m")
        ax.grid(which="both", alpha=0.25)
    axes[0].set_ylabel("$dV_{BR}/dt_s$  (V/$\\mu$m)")
    axes[0].legend(fontsize=9)
    fig.suptitle("Slope of the breakdown-voltage \"U\" curve vs film thickness")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig_dbv_dts.png", dpi=140)
    plt.close(fig)


def fig_log_slope():
    ts_um = np.logspace(np.log10(0.02), np.log10(50), 200)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    for mat in MATERIALS:
        c = COLORS[mat.name]
        lhs = np.array([-dlnEc_dlnts(t, mat) for t in ts_um])
        ax.plot(ts_um, lhs, color=c, lw=2.2, label=f"{mat.name}: $-d\\ln E_c/d\\ln t_s$")
        for tox_um, ls in zip(TOX_LIST, ("--", "-.", ":")):
            rhs = np.array([rhs_slope(t, tox_um, mat) for t in ts_um])
            ax.plot(ts_um, rhs, color=c, lw=1.1, ls=ls, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel("$t_s$  ($\\mu$m)")
    ax.set_ylabel("slope (dimensionless)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Log-slope minimum condition: solid = LHS (material only,\n"
                 "$t_{ox}$-independent), dashed/dash-dot/dotted = RHS at "
                 "$t_{ox}$=2/4/6 $\\mu$m")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig_log_slope.png", dpi=150)
    plt.close(fig)

    print("Log-slope intersection (approx. U-curve minimum) vs direct numerical minimization:")
    for mat in MATERIALS:
        for tox_um in TOX_LIST:
            f = lambda logt: -dlnEc_dlnts(np.exp(logt), mat) - rhs_slope(np.exp(logt), tox_um, mat)
            lo, hi = np.log(0.02), np.log(50)
            # bracket a sign change
            grid = np.linspace(lo, hi, 400)
            vals = np.array([f(g) for g in grid])
            sign_changes = np.where(np.diff(np.sign(vals)) != 0)[0]
            if len(sign_changes) == 0:
                print(f"  {mat.name:8s} tox={tox_um:.0f}um: no intersection found")
                continue
            i = sign_changes[0]
            ts_cross = np.exp(brentq(f, grid[i], grid[i + 1]))
            res = minimize_scalar(lambda t: breakdown_voltage(t * UM, tox_um * UM, W=None, mat=mat),
                                   bounds=(0.02, 50), method="bounded",
                                   options={"xatol": 1e-6})
            print(f"  {mat.name:8s} tox={tox_um:.0f}um: log-slope crossing t_s="
                  f"{ts_cross:6.3f}um   direct numerical minimum t_s={res.x:6.3f}um  "
                  f"(BV_min={res.fun:7.1f} V)")


def fig_plateau_width():
    """Width (in um) of the t_s-interval around the minimum over which the
    relative slope |dV_BR/dt_s| / BV_min stays below a small threshold --
    a direct, quantitative measure of how "broad" the plateau is."""
    threshold = 0.01   # 1%/um relative slope
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    tox_um = np.linspace(1.0, 10.0, 19)

    for mat in MATERIALS:
        widths = []
        for tox in tox_um:
            res = minimize_scalar(lambda t: breakdown_voltage(t * UM, tox * UM, W=None, mat=mat),
                                   bounds=(0.02, 50), method="bounded",
                                   options={"xatol": 1e-7})
            ts_min, bv_min = res.x, res.fun

            def rel_slope(ts_um):
                return abs(dbv_dts(ts_um, tox, mat)) / bv_min - threshold

            # search outward from the minimum until the relative slope exceeds threshold
            lo = ts_min
            while lo > 0.021 and rel_slope(lo) < 0:
                lo /= 1.05
            hi = ts_min
            while hi < 49 and rel_slope(hi) < 0:
                hi *= 1.05
            widths.append(hi - lo)

        ax.plot(tox_um, widths, color=COLORS[mat.name], lw=2.0, marker="o", ms=4,
                label=mat.name)

    ax.set_xlabel("$t_{ox}$  ($\\mu$m)")
    ax.set_ylabel(f"plateau width  ($\\mu$m film thickness, |slope|$/BV_{{min}}$<{threshold*100:.0f}\\%/$\\mu$m)")
    ax.set_yscale("log")
    ax.set_title("Plateau width around the BV(t_s) minimum")
    ax.legend(fontsize=9)
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig_plateau_width.png", dpi=150)
    plt.close(fig)


def main():
    fig_dbv_dts()
    fig_log_slope()
    fig_plateau_width()


if __name__ == "__main__":
    main()
