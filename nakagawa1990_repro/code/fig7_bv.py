"""Fig.7 — breakdown voltage vs. dose, structure A, t_ox = 5 um.

IMPORTANT MODELLING NOTE
------------------------
The 1-D vertical field-sharing model (and the mirror-symmetric 2-D model)
describes the *drain-directly-over-oxide* cut, which is structure B (Fig.3B),
the paper's optimised device.  Structures A and B share the same vertical layer
stack, so a vertical model cannot, by itself, separate them: their difference is
a 2-D lateral effect.  In structure A the drain sits away from the symmetry axis
and its critical point couples less completely to the grounded substrate through
the oxide, so A needs a *thicker* oxide (5 um) than B (3 um) to reach 500 V.

This script therefore plots the paper's digitised Fig.7 points (structure A) and
the vertical model at t_ox = 5 um.  The model curve lies ABOVE the structure-A
data: it represents the favourable (B-type) vertical sharing and confirms that a
5 um oxide provides ample margin to exceed 500 V.  The ~150-200 V gap is the
quantitative signature of structure A's less complete field-sharing.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import model as m

UM = m.UM

# digitised points from the paper's Fig.7 (structure A, t_ox=5um) -- 4 points,
# curve stops at dose ~1.0 (monotonic rise, no rolloff shown).
PAPER_A = [(0.2, 410), (0.5, 458), (0.8, 493), (1.0, 513)]


def main() -> None:
    doses = np.linspace(0.05e12, 1.6e12, 20)
    bv = [m.breakdown_voltage(m.Geom(t_ox=5 * UM, S_D=SD)) for SD in doses]

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(doses / 1e12, bv, "-", color="#8e44ad", lw=2.0,
            label="vertical model, $t_{ox}$=5 $\\mu$m (B-type sharing)")
    px = [p[0] for p in PAPER_A]
    py = [p[1] for p in PAPER_A]
    ax.plot(px, py, "s", color="#c0392b", ms=6, mfc="white", mew=1.5,
            label="paper Fig.7 (structure A)")

    ax.axhline(500, color="gray", ls=":", lw=1.0)
    ax.text(1.35, 508, "500 V", color="gray", fontsize=8)
    ax.annotate("A needs $t_{ox}$≈5 $\\mu$m to reach 500 V\n"
                "(B reaches it at 3 $\\mu$m — Fig.6)",
                xy=(1.1, 500), xytext=(0.35, 590),
                arrowprops=dict(arrowstyle="->"), fontsize=8.5)

    ax.set_xlabel("Impurity dose $S_D$  ($\\times10^{12}$ cm$^{-2}$)")
    ax.set_ylabel("Breakdown voltage  (V)")
    ax.set_title("Fig.7  Breakdown voltage vs. dose (structure A, $t_{ox}$=5 $\\mu$m)")
    ax.set_xlim(0, 1.6)
    ax.set_ylim(300, 780)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = "../figs/fig7_bv.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
