"""Fig.5 reproduction — electric-field distribution along symmetry axis AA'.

Structure B, t_ox = 3 um, V = 500 V, for four shallow-layer doses:
    A: 2.5e12   B: 2.0e12   C: 1.5e12   D: 0.0   (cm^-2)
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import model as m

UM = m.UM
V = 500.0
CURVES = [("A", 2.5e12), ("B", 2.0e12), ("C", 1.5e12), ("D", 0.0)]
COLORS = {"A": "#c0392b", "B": "#e67e22", "C": "#27ae60", "D": "#2c3e50"}


def main() -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for lab, SD in CURVES:
        g = m.Geom(t_ox=3.0 * UM, S_D=SD)
        y, phi, E = m.solve_potential(g, V)
        ax.plot(y / UM, np.abs(E), color=COLORS[lab], lw=1.8,
                label=f"{lab}: $S_D$={SD/1e12:.1f}$\\times10^{{12}}$")

    L_si = 20.0
    ax.axvspan(L_si, L_si + 3.0, color="#d0d0d0", alpha=0.5, zorder=0)
    ax.text(L_si + 1.5, 0.05e6, "Oxide", ha="center", fontsize=9)
    ax.text(10, 0.05e6, "Silicon", ha="center", fontsize=9)
    ax.axvline(L_si, color="k", lw=0.7, ls="--")

    ax.set_xlabel("Depth  y  ($\\mu$m)")
    ax.set_ylabel("Electric field  |E|  (V/cm)")
    ax.set_title("Fig.5  Field along symmetry axis AA'  "
                 "(structure B, $t_{ox}$=3 $\\mu$m, 500 V)")
    ax.set_xlim(0, 23)
    ax.set_ylim(0, None)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = "../figs/fig5_field.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)

    # numeric summary
    print("\nBulk-Si field @10um and oxide field @21.5um, V=500:")
    for lab, SD in CURVES:
        g = m.Geom(t_ox=3.0 * UM, S_D=SD)
        y, phi, E = m.solve_potential(g, V)
        Ea = np.abs(E)
        i10 = np.searchsorted(y, 10 * UM)
        iox = np.searchsorted(y, 21.5 * UM)
        print(f"  {lab}: E_bulk={Ea[i10]:.2e}  E_ox={Ea[iox]:.2e}  "
              f"Vox/V={m.oxide_voltage_fraction(g, V):.2f}")


if __name__ == "__main__":
    main()
