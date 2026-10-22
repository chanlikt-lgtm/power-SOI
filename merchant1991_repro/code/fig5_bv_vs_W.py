"""Fig. 5 -- breakdown voltage vs lateral depletion width W at breakdown, for the
SOI device of Fig. 2 (t_s=0.1um, t_ox=2um, full eq. 23: max(I1, I2'+I2'')=1)
compared against the reference 1-D p-i-n diode of the same width W (I1=1 alone).
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from model import UM, breakdown_voltage, pin_diode_breakdown_voltage

OUTDIR = "../figs"


def main():
    W_um = np.linspace(15, 100, 90)
    t_s, t_ox = 0.1 * UM, 2.0 * UM

    bv_soi = np.array([breakdown_voltage(t_s, t_ox, w * UM) for w in W_um])
    bv_pin = np.array([pin_diode_breakdown_voltage(w * UM) for w in W_um])

    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.plot(W_um, bv_pin, "k--", lw=1.6, label="1-D $p$-$i$-$n$ diode")
    ax.plot(W_um, bv_soi, "b-", lw=2.0, label="SOI device (Fig.2)")
    ax.set_xlabel("lateral depletion width  ($\\mu$m)")
    ax.set_ylabel("breakdown voltage (V)")
    ax.set_xlim(0, 100)
    ax.set_title("Fig.5  BV vs lateral depletion width at breakdown\n"
                  "($t_s$=0.1 $\\mu$m, $t_{ox}$=2 $\\mu$m)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig5_bv_vs_W.png", dpi=150)
    plt.close(fig)

    for w in (20, 40, 55, 70, 100):
        print(f"W={w:3d} um:  BV(SOI)={breakdown_voltage(t_s, t_ox, w*UM):7.1f} V   "
              f"BV(1-D p-i-n)={pin_diode_breakdown_voltage(w*UM):7.1f} V")


if __name__ == "__main__":
    main()
