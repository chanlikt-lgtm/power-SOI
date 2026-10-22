"""Fig.6 reproduction — breakdown voltage vs. shallow-layer dose, structure B.

Curves for t_ox = 1, 2, 3 um.  Paper reference points (digitized from Fig.6)
overlaid as markers for comparison.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

import model as m

UM = m.UM

# digitized approximate reference points from the paper's Fig.6 (structure B)
PAPER = {
    3: [(0.25, 535), (0.5, 560), (0.75, 578), (1.0, 575), (1.25, 568),
        (1.5, 558), (2.0, 500), (2.5, 430)],
    2: [(0.25, 448), (0.5, 470), (0.75, 480), (1.0, 478), (1.25, 470),
        (1.5, 455), (2.0, 405), (2.5, 350)],
    1: [(0.25, 355), (0.5, 375), (0.75, 382), (1.0, 378), (1.25, 368),
        (1.5, 352), (2.0, 305), (2.5, 232)],
}
COLORS = {3: "#c0392b", 2: "#2980b9", 1: "#27ae60"}


def main() -> None:
    doses = np.linspace(0.05e12, 2.6e12, 26)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    for tox in (3, 2, 1):
        bv = [m.breakdown_voltage(m.Geom(t_ox=tox * UM, S_D=SD)) for SD in doses]
        ax.plot(doses / 1e12, bv, color=COLORS[tox], lw=2.0,
                label=f"model  $t_{{ox}}$={tox} $\\mu$m")
        px = [p[0] for p in PAPER[tox]]
        py = [p[1] for p in PAPER[tox]]
        ax.plot(px, py, "o", color=COLORS[tox], ms=5, mfc="white", mew=1.4,
                label=f"paper  $t_{{ox}}$={tox} $\\mu$m")

    ax.axhline(500, color="gray", ls=":", lw=1.0)
    ax.text(2.35, 508, "500 V", color="gray", fontsize=8)
    ax.set_xlabel("Impurity dose $S_D$  ($\\times10^{12}$ cm$^{-2}$)")
    ax.set_ylabel("Breakdown voltage  (V)")
    ax.set_title("Fig.6  Breakdown voltage vs. dose (structure B)")
    ax.set_xlim(0, 2.6)
    ax.set_ylim(0, 650)
    ax.legend(loc="lower left", fontsize=8, ncol=3, framealpha=0.95)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = "../figs/fig6_bv.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)

    print("\nModel peak BV per oxide thickness:")
    for tox in (1, 2, 3):
        bv = [m.breakdown_voltage(m.Geom(t_ox=tox * UM, S_D=SD)) for SD in doses]
        i = int(np.argmax(bv))
        print(f"  t_ox={tox}um: peak {bv[i]:.0f} V @ dose {doses[i]/1e12:.2f}e12")


if __name__ == "__main__":
    main()
