"""Extension beyond the paper: re-evaluate the same closed-form Merchant 1991
breakdown model (model.py) with 4H-SiC and GaN Chynoweth ionization
coefficients in place of silicon's, holding the device geometry and the
buried SiO2 oxide fixed. See derivation.md sec.6 for the physics and full
literature citations for the 4H-SiC/GaN coefficients.

Produces three figures:
  fig4_materials.png   BV vs SOI/film thickness t_s, t_ox=2um, Si/4H-SiC/GaN,
                        W->inf limit (materials analogue of Fig.4), with each
                        material's own oxide-breakdown cap.
  fig5_materials.png   BV vs lateral depletion width W, t_s=0.1um, t_ox=2um,
                        Si/4H-SiC/GaN (materials analogue of Fig.5).
  fig_oxide_crossover.png   For a range of buried-oxide thicknesses t_ox, the
                        film thickness t_s below which the device is
                        oxide-limited rather than semiconductor-limited (the
                        crossing point between the flat oxide cap and the
                        falling left branch of each material's BV(t_s) "U"
                        curve) -- i.e. how thick a 4H-SiC/GaN film must stay
                        before the buried SiO2, not the semiconductor, sets
                        the achievable breakdown voltage.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from model import (UM, SI, SIC_4H, GAN, breakdown_voltage,
                    oxide_breakdown_voltage, pin_diode_breakdown_voltage,
                    crossover_film_thickness)

OUTDIR = "../figs"
MATERIALS = (SI, SIC_4H, GAN)
COLORS = {"Si": "tab:blue", "4H-SiC": "tab:orange", "GaN": "tab:green"}


def fig4_materials():
    ts_um = np.logspace(np.log10(0.01), np.log10(100), 200)
    t_ox = 2.0 * UM
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mat in MATERIALS:
        bv = np.array([breakdown_voltage(t * UM, t_ox, W=None, mat=mat) for t in ts_um])
        vox = np.array([oxide_breakdown_voltage(t * UM, t_ox, mat) for t in ts_um])
        capped = np.minimum(bv, vox)
        c = COLORS[mat.name]
        ax.plot(ts_um, capped, lw=2.0, color=c, label=f"{mat.name} (oxide-capped)")
        ax.plot(ts_um, bv, lw=1.0, color=c, ls=":", alpha=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("film thickness $t_s$  ($\\mu$m)")
    ax.set_ylabel("breakdown voltage (V)")
    ax.set_xlim(0.01, 100)
    ax.set_title("BV vs film thickness, $t_{ox}$=2 $\\mu$m ($W\\to\\infty$ limit)\n"
                 "solid = oxide-capped, dotted = semiconductor-limited (uncapped)")
    ax.legend(fontsize=9)
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig4_materials.png", dpi=150)
    plt.close(fig)


def fig5_materials():
    W_um = np.linspace(15, 100, 90)
    t_s, t_ox = 0.1 * UM, 2.0 * UM
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mat in MATERIALS:
        bv = np.array([breakdown_voltage(t_s, t_ox, w * UM, mat) for w in W_um])
        vox = oxide_breakdown_voltage(t_s, t_ox, mat)
        c = COLORS[mat.name]
        ax.plot(W_um, np.minimum(bv, vox), lw=2.0, color=c, label=f"{mat.name}")
        ax.axhline(vox, color=c, lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel("lateral depletion width  ($\\mu$m)")
    ax.set_ylabel("breakdown voltage (V)")
    ax.set_xlim(0, 100)
    ax.set_title("BV vs lateral depletion width $W$\n"
                 "($t_s$=0.1 $\\mu$m, $t_{ox}$=2 $\\mu$m; dashed = each material's oxide cap)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig5_materials.png", dpi=150)
    plt.close(fig)


def fig_oxide_crossover():
    """Verify crossover_film_thickness() (closed form, exactly t_ox-independent)
    against a brute-force numerical scan over t_ox, and plot both."""
    tox_um = np.logspace(np.log10(0.5), np.log10(20), 15)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mat in MATERIALS:
        c = COLORS[mat.name]
        ts_closed = crossover_film_thickness(mat) / UM
        ts_numeric = []
        for t in tox_um:
            f = lambda ts: (breakdown_voltage(ts * UM, t * UM, W=None, mat=mat)
                             - oxide_breakdown_voltage(ts * UM, t * UM, mat))
            ts_numeric.append(brentq(f, 1e-4, 90.0, xtol=1e-10))
        ax.plot(tox_um, ts_numeric, color=c, lw=2.0, marker="o", ms=4, label=f"{mat.name} (numeric scan)")
        ax.axhline(ts_closed, color=c, lw=1.0, ls="--", alpha=0.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("buried-oxide thickness $t_{ox}$  ($\\mu$m)")
    ax.set_ylabel("crossover film thickness $t_s$  ($\\mu$m)")
    ax.set_title("Oxide-limited vs semiconductor-limited crossover\n"
                 "below the curve: oxide limits BV; above: film limits BV\n"
                 "(dashed = closed-form, t_ox-independent; markers = numeric scan)")
    ax.legend(fontsize=8.5)
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig_oxide_crossover.png", dpi=150)
    plt.close(fig)

    print(f"\nOxide-limited/semiconductor-limited crossover t_s "
          f"(closed form, exactly independent of t_ox):")
    for mat in MATERIALS:
        tc = crossover_film_thickness(mat)
        print(f"  {mat.name:8s}  t_s < {tc/UM:6.3f} um -> oxide-limited; "
              f"t_s > {tc/UM:6.3f} um -> {mat.name}-film-limited")


def main():
    fig4_materials()
    fig5_materials()
    fig_oxide_crossover()

    t_s, t_ox = 0.1 * UM, 2.0 * UM
    print(f"\nAt t_s=0.1um, t_ox=2um, W->inf:")
    for mat in MATERIALS:
        bv = breakdown_voltage(t_s, t_ox, W=None, mat=mat)
        vox = oxide_breakdown_voltage(t_s, t_ox, mat)
        print(f"  {mat.name:8s}  semiconductor-limit BV={bv:8.0f} V   "
              f"oxide cap={vox:6.0f} V   actual (capped)={min(bv,vox):6.0f} V   "
              f"[{mat.citation}]")


if __name__ == "__main__":
    main()
