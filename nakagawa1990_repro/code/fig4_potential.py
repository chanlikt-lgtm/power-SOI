"""Fig.4 reproduction — potential distribution, structure B, 500 V.

The paper draws Fig.4 as a 3-D wireframe of potential (height) over the
cross-section: a high plateau on the drain/n-layer side sloping down to the p+
side, with a steep face on the left across the *oxide film* -- that face is the
>half-voltage drop across the buried oxide.

This script reproduces that 3-D view (to match the paper) plus a 2-D contour map
and a vertical-cut panel that quantifies the oxide share.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d proj)

import solver2d as s2
from model import UM


def main() -> None:
    d = s2.Device2D(structure="B", S_D=1.5e12, t_ox=3 * UM, Lx=50 * UM)
    phi = d.solve(500.0)
    Xum, Yum = d.X / UM, d.Y / UM
    i_if = np.searchsorted(d.y, d.L_si)

    # ================= 3-D wireframe view (matches the paper) =================
    v_ox = phi[i_if, 0]                       # potential at Si/oxide interface, x=0
    ox_share = v_ox / 500.0                   # oxide's share of the applied 500 V
    si_share = 1.0 - ox_share

    fig = plt.figure(figsize=(9.4, 6.8))
    ax = fig.add_subplot(111, projection="3d")

    # orient like the paper: drain/n-layer high plateau on the left, p+ low-right.
    # depth axis kept coarse so the oxide "cliff" reads clearly.
    step_x = max(1, d.nx // 60)
    step_y = max(1, d.ny // 45)
    Xp = Xum[::step_y, ::step_x]
    Yp = Yum[::step_y, ::step_x]
    Pp = phi[::step_y, ::step_x]

    surf = ax.plot_surface(Xp, Yp, Pp, rcount=45, ccount=60, cmap="turbo",
                           edgecolor="k", linewidth=0.22, alpha=0.9,
                           antialiased=True)

    # Si/buried-oxide interface traced across the surface (y = L_si = 20 um)
    yif = d.L_si / UM
    ax.plot(Xum[i_if, :], np.full(d.nx, yif), phi[i_if, :],
            color="crimson", lw=2.2, label="Si / buried-oxide interface")

    # the drain-edge (x=0) depth profile -- the "voltage cliff": gentle drop
    # through the silicon, then a STEEP drop across the thin buried oxide
    yline = d.y / UM
    zline = phi[:, 0]
    ax.plot(np.zeros_like(yline), yline, zline, color="k", lw=1.4)
    ax.plot(np.zeros(2), [yif, (d.L_si + d.t_ox) / UM], [v_ox, 0.0],
            color="magenta", lw=3.2, label="oxide voltage drop (drain edge)")

    # ---- region + electrode labels (placed in free space, at true depths) ----
    ax.set_xlim(0, 50)
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9)
    ax.text(1.0, 0.0, 525, "n layer / drain  (500 V)", color="darkred",
            fontsize=9, ha="left", va="bottom", fontweight="bold")
    ax.text(50, 1.0, 55, "p$^+$ layer / source  (0 V)", color="navy",
            fontsize=9, ha="right", va="bottom", fontweight="bold")
    ax.text(6, 9.0, 360, f"Silicon film 0-20 $\\mu$m\n($\\approx${100*si_share:.0f}% of $V$)",
            color="black", fontsize=8, ha="center", va="center", bbox=box)
    ax.text(3, 21.5, 150,
            f"Buried oxide 20-23 $\\mu$m\n$>$half of $V$ ($\\approx${100*ox_share:.0f}%)",
            color="magenta", fontsize=8, ha="center", va="center", bbox=box)

    ax.set_xlabel("x  ($\\mu$m)")
    ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_zlabel("potential (V)")
    ax.set_zlim(0, 500)
    ax.set_title("Nakagawa 1990 Fig.4 in 3-D -- potential distribution "
                 "(structure B, 500 V)\n"
                 f"more than half the applied voltage falls across the buried oxide "
                 f"($\\approx${100*ox_share:.0f}%)")
    ax.view_init(elev=28, azim=-60)
    ax.invert_yaxis()
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.02,
                 label="potential (V)")
    fig.tight_layout()
    out3d = "../figs/fig4_potential_3d.png"
    fig.savefig(out3d, dpi=150)
    plt.close(fig)
    print("wrote", out3d)

    # ================= 2-D contour + vertical cut ============================
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                   gridspec_kw={"height_ratios": [2.3, 1]})
    levels = np.arange(0, 501, 25)
    cf = ax1.contourf(Xum, Yum, phi, levels=levels, cmap="turbo")
    ax1.contour(Xum, Yum, phi, levels=levels, colors="k", linewidths=0.4)
    ax1.axhline(20.0, color="white", lw=1.5, ls="--")
    ax1.text(25, 21.6, "buried oxide", color="white", ha="center", fontsize=9)
    ax1.text(4, 2.0, "drain 500 V", color="white", fontsize=9, weight="bold")
    ax1.text(40, 2.0, "source 0 V", color="white", fontsize=9, weight="bold")
    ax1.invert_yaxis()
    ax1.set_xlabel("x  ($\\mu$m)")
    ax1.set_ylabel("depth y  ($\\mu$m)")
    ax1.set_title("Fig.4  Potential distribution — 2-D map (structure B, 500 V)")
    fig.colorbar(cf, ax=ax1, label="potential (V)", shrink=0.9)

    col = phi[:, 0]
    ax2.plot(d.y / UM, col, "b-", lw=2)
    ax2.axvspan(20, 23, color="#d0d0d0", alpha=0.6)
    ax2.axvline(20, color="k", ls="--", lw=0.7)
    v_ox = col[i_if]
    ax2.annotate(f"{v_ox/500*100:.0f}% of 500 V\nacross the oxide",
                 xy=(21.5, v_ox / 2), xytext=(11, 120),
                 arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax2.text(9, 380, "Silicon", fontsize=9)
    ax2.text(21.5, 380, "Oxide", fontsize=9, ha="center")
    ax2.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain, x=0)")
    ax2.set_ylabel("potential (V)")
    ax2.set_xlim(0, 23)
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    out2d = "../figs/fig4_potential.png"
    fig.savefig(out2d, dpi=150)
    plt.close(fig)
    print("wrote", out2d)
    print(f"oxide voltage share under drain = {v_ox/500:.2f}")


if __name__ == "__main__":
    main()
