"""Fig. 4 -- breakdown voltage vs SOI thickness t_s, for t_ox = 2/4/6 um, in the
W -> infinity limit (paper eq. 23 with I1, I2'' neglected -> solve I2'(V)=1).
Also draws the buried-oxide dielectric-breakdown cap for t_ox=2um.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from model import UM, breakdown_voltage, oxide_breakdown_voltage

OUTDIR = "../figs"


def main():
    ts_um = np.logspace(np.log10(0.01), np.log10(100), 220)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))

    for tox_um in (2.0, 4.0, 6.0):
        bv = np.array([breakdown_voltage(t * UM, tox_um * UM, W=None) for t in ts_um])
        vox = oxide_breakdown_voltage(ts_um * UM, tox_um * UM)
        bv_capped = np.minimum(bv, vox)
        ax.plot(ts_um, bv_capped, lw=1.8, label=f"$t_{{ox}}$={tox_um:.0f} $\\mu$m")

    vox2 = oxide_breakdown_voltage(np.array([0.03, 0.3]) * UM, 2.0 * UM)
    ax.plot([0.03, 0.3], vox2, "k--", lw=1.2)

    ax.set_xscale("log")
    ax.set_xlabel("SOI thickness  ($\\mu$m)")
    ax.set_ylabel("breakdown voltage (V)")
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0, 1400)
    ax.set_title("Fig.4  Breakdown voltage vs SOI thickness ($W\\to\\infty$ limit)")
    ax.legend(fontsize=9, loc="upper center")
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig4_bv_vs_ts.png", dpi=150)
    plt.close(fig)

    bv_01_2 = breakdown_voltage(0.1 * UM, 2.0 * UM, W=None)
    print(f"Fig.4  BV(t_s=0.1um, t_ox=2um, W->inf) = {bv_01_2:.0f} V  "
          f"(matches the Fig.5 W->100um plateau, ~865 V; paper text's rounded "
          f"'700 V' refers to the >700 V experimental headline, not a precise "
          f"Fig.4 graph reading)")
    print(f"       oxide breakdown cap at same point = "
          f"{oxide_breakdown_voltage(0.1*UM, 2.0*UM):.0f} V  (paper dashed line: 1200 V)")


if __name__ == "__main__":
    main()
