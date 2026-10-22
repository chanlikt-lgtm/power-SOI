"""Fig. 3 (Merchant 1991) -- 3-D reproduction of the equipotential landscape.

The paper's Fig. 3 is a 2-D *numerical* Poisson solution: an equipotential
contour map of the optimized SOI device (t_s=0.1 um, t_ox=2 um, L=70 um). The
existing fig3_potential.py plots the closed-form phi(x,y) as 2-D contours. This
companion does three things the paper's one figure implies -- derive, simulate,
and plot in 3-D:

  DERIVE   -- restate the closed form (paper eqs. 5-9) that the optimized
              lateral doping N(x) is *designed* to produce, in the depleted
              window 0<=x<=W:

                phi(x,y) = (V x / W) [1 - y^2 / (2 t t_s)]        0<=y<=t_s   (Si)
                phi(x,y) = phi(x,t_s) - (K_s/K_ox)(V x)/(t W)(y-t_s)  t_s<y   (oxide)
                t        = t_s/2 + (K_s/K_ox) t_ox                            (eq 8)

              and the linear donor sheet that forces it (from eq. 3, E_y=Vxy/(tWt_s)):

                N(x) = K_s eps0 V x / (q t W t_s)      (linear in x)

  SIMULATE -- independently solve the 2-D Poisson equation

                div( K(y) grad phi ) = -q N(x)/eps0

              on 0<=x<=W, 0<=y<=t_s+t_ox with a variable-permittivity
              finite-volume scheme (K_s in the film, K_ox in the buried oxide,
              the Si/oxide displacement-continuity x3 field jump emerging from
              the face permittivities -- not imposed). BCs: E_y=0 at the top
              (y=0, fully-depleted surface), grounded substrate phi=0 at the
              bottom (y=t_s+t_ox), phi=0 at the P+ edge (x=0), and the closed-
              form vertical profile at the depletion edge (x=W). The interior
              linearity phi ~ x is NOT imposed -- it is what the linear N(x)
              produces, so recovering it validates the derivation.

  PLOT 3D  -- a hero 3-D surface of phi over the full device cross-section
              0<=x<=L (closed form for x<W, quasi-neutral plateau for x>W),
              and a second 3-D figure overlaying the numerical solve (surface)
              on the closed form (wireframe) over the depleted window, with the
              max/RMS discrepancy printed.

Units: CGS-practical (cm, V, cm^-2), UM=1e-4 cm, as elsewhere in this repo.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

from model import UM, Q, EPS0, EPS_OX, SI, t_eff
from fig3_potential import find_design_W, phi as phi_closed

OUTDIR = "../figs"
K_S = SI.eps_r          # 11.7
K_OX = EPS_OX           # 3.9


# ---------------------------------------------------------------------------
# SIMULATE -- variable-permittivity 2-D finite-volume Poisson solve
# ---------------------------------------------------------------------------
def solve_poisson_2d(t_s, t_ox, W, V, n_si=16, n_ox=90, nx=161):
    """Solve div(K grad phi) = -q N(x)/eps0 on the depleted window [0,W]x[0,t_s+t_ox].

    Returns (X, Y, PHI) meshes (cm, cm, V). Node-based finite volume with the
    sharp Si/oxide interface placed *on a cell face* at y=t_s -- no node
    straddles the permittivity jump. The interface face uses the distance-
    weighted harmonic permittivity (exact series combination of the Si and
    oxide half-gaps), so the x3 vertical-field jump is captured to second order
    and every silicon node carries a full control-volume of donor charge.
    """
    # --- grid: NO node at y=t_s; the interface falls on the face between the
    # last silicon node (y<t_s) and the first oxide node (y>t_s) ---
    dsi = t_s / n_si
    y_si = np.arange(n_si) * dsi                       # 0 .. t_s-dsi   (all < t_s)
    y_ox = np.linspace(t_s, t_s + t_ox, n_ox)[1:]      # t_s+.. .. t_s+t_ox (>t_s)
    y = np.concatenate([y_si, y_ox])
    x = np.linspace(0.0, W, nx)
    ny = y.size
    j0 = n_si - 1                                      # last Si node; face j0 = interface

    # lateral permittivity per row (unambiguous -- no node sits on t_s)
    Kx = np.where(y < t_s, K_S, K_OX)

    # vertical face permittivity between rows; interface face gets the distance-
    # weighted harmonic mean with the material boundary pinned at y=t_s
    Kyf = np.where(0.5 * (y[:-1] + y[1:]) < t_s, K_S, K_OX)
    Kyf[j0] = (y[j0 + 1] - y[j0]) / ((t_s - y[j0]) / K_S + (y[j0 + 1] - t_s) / K_OX)

    # face positions (midpoints, except the interface face pinned at t_s), from
    # which the control-volume heights hy are built so the last Si / first oxide
    # cells meet exactly at the interface
    yf = 0.5 * (y[:-1] + y[1:])
    yf[j0] = t_s
    hy = np.empty(ny)
    hy[0] = yf[0] - y[0]
    hy[1:-1] = yf[1:] - yf[:-1]
    hy[-1] = y[-1] - yf[-1]

    dx = W / (nx - 1)
    hx = np.full(nx, dx); hx[0] = hx[-1] = 0.5 * dx

    # derived linear donor sheet N(x) (cm^-3), silicon rows only
    t = t_eff(t_s, t_ox)
    Nx = K_S * EPS0 * V * x / (Q * t * W * t_s)                  # N(x) ~ x

    # Dirichlet edge value at x=W (closed-form vertical profile), and phi=0 walls
    phiW = phi_closed(np.full_like(y, W), y, V, W, t_s, t_ox)

    # unknowns: i in 1..nx-2 (x=0,W Dirichlet), j in 0..ny-2 (bottom Dirichlet;
    # top j=0 is a no-flux Neumann unknown)
    ii = np.arange(1, nx - 1)
    jj = np.arange(0, ny - 1)
    nI, nJ = ii.size, jj.size
    kof = -np.ones((nx, ny), dtype=int)
    kof[1:nx - 1, 0:ny - 1] = np.arange(nI * nJ).reshape(nI, nJ)

    rows, cols, vals = [], [], []
    b = np.zeros(nI * nJ)
    for a, i in enumerate(ii):
        for c, j in enumerate(jj):
            k = kof[i, j]
            aP = 0.0
            # east / west (lateral, K uniform along the row)
            aE = Kx[j] * hy[j] / (x[i + 1] - x[i]); aP += aE
            aW = Kx[j] * hy[j] / (x[i] - x[i - 1]); aP += aW
            if i + 1 == nx - 1:
                b[k] += aE * phiW[j]                 # Dirichlet neighbour @ x=W
            else:
                rows.append(k); cols.append(kof[i + 1, j]); vals.append(-aE)
            if i - 1 == 0:
                b[k] += aW * 0.0                     # Dirichlet neighbour @ x=0
            else:
                rows.append(k); cols.append(kof[i - 1, j]); vals.append(-aW)
            # north (toward larger y / substrate)
            aN = Kyf[j] * hx[i] / (y[j + 1] - y[j]); aP += aN
            if j + 1 == ny - 1:
                b[k] += aN * 0.0                     # grounded substrate
            else:
                rows.append(k); cols.append(kof[i, j + 1]); vals.append(-aN)
            # south (toward smaller y / surface): omitted at j=0 -> E_y=0 Neumann
            if j - 1 >= 0:
                aS = Kyf[j - 1] * hx[i] / (y[j] - y[j - 1]); aP += aS
                rows.append(k); cols.append(kof[i, j - 1]); vals.append(-aS)
            # source: ionized donors, silicon rows only (full control volume)
            if y[j] < t_s:
                b[k] += Q * Nx[i] / EPS0 * hx[i] * hy[j]
            rows.append(k); cols.append(k); vals.append(aP)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(nI * nJ, nI * nJ))
    sol = spsolve(A, b)

    PHI = np.zeros((nx, ny))
    PHI[0, :] = 0.0
    PHI[nx - 1, :] = phiW
    PHI[:, ny - 1] = 0.0
    for a, i in enumerate(ii):
        for c, j in enumerate(jj):
            PHI[i, j] = sol[kof[i, j]]
    X, Yg = np.meshgrid(x, y, indexing="ij")
    return X, Yg, PHI


# ---------------------------------------------------------------------------
# PLOT 3D
# ---------------------------------------------------------------------------
def hero_surface(t_s, t_ox, L, W, V):
    """Closed-form phi over the full device cross-section, as a 3-D surface."""
    nx, ny = 260, 150
    x = np.linspace(0, L, nx)
    y = np.linspace(0, t_s + t_ox, ny)
    X, Yg = np.meshgrid(x, y, indexing="ij")
    # depleted window: closed form; beyond W: quasi-neutral film @V, oxide V->0
    Phi = np.where(
        X <= W,
        phi_closed(X, Yg, V, W, t_s, t_ox),
        np.where(Yg <= t_s, V, V * (1.0 - (Yg - t_s) / t_ox)),
    )
    fig = plt.figure(figsize=(9.6, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X / UM, Yg / UM, Phi, cmap="viridis", rstride=2, cstride=2,
        linewidth=0, antialiased=True, alpha=0.88,
    )
    # project the paper's own equipotential contours onto the floor (z=0):
    # this is literally the paper's Fig.3 contour map sitting under the surface.
    paper_levels = np.array([79.5, 159, 239, 319, 399, 479, 559, 639, 700])
    cset = ax.contour(X / UM, Yg / UM, Phi, levels=paper_levels,
                      zdir="z", offset=0.0, cmap="viridis", linewidths=1.1)
    ax.clabel(cset, fmt="%.0f", fontsize=6)
    # and trace the same contours on the surface itself
    ax.contour(X / UM, Yg / UM, Phi, levels=paper_levels,
               colors="k", linewidths=0.5, alpha=0.5)
    # mark the Si/oxide interface and the depletion edge on the surface
    yi = t_s / UM
    ax.plot(x / UM, np.full_like(x, yi),
            np.where(x <= W, phi_closed(x, np.full_like(x, t_s), V, W, t_s, t_ox),
                     V), color="crimson", lw=2.0, label="Si / buried-oxide interface")
    yv = np.linspace(0, t_s + t_ox, 60)
    ax.plot(np.full_like(yv, W / UM), yv / UM,
            phi_closed(np.full_like(yv, W), yv, V, W, t_s, t_ox),
            color="k", lw=1.8, ls="--", label="depletion edge  $x=W$")

    # --- extended region + electrode labels (Si film / oxide / substrate) ----
    # Region labels sit in free space just past the drain edge, each at the TRUE
    # depth y of the layer it names, so they line up with the drain-side curtain.
    ax.set_xlim(0, L / UM + 15)                     # right-hand room for labels
    xr = L / UM + 1.5
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9)
    ax.text(xr, 0.0, V - 5, "Silicon film  $t_s$=0.1 $\\mu$m",
            color="crimson", fontsize=8, ha="left", va="center", bbox=box)
    ax.text(xr, 1.0, V * (1 - 0.9 / 2), "Buried oxide  $t_{ox}$=2 $\\mu$m",
            color="saddlebrown", fontsize=8, ha="left", va="center", bbox=box)
    ax.text(xr, (t_s + t_ox) / UM, 20.0, "N$^+$ substrate (0 V)",
            color="dimgray", fontsize=8, ha="left", va="center", bbox=box)
    # electrodes on the top (y=0) surface: P+ at x=0 (0 V), N+ contact at x=L (+V)
    ax.text(0.0, 0.0, 70, "P$^+$", color="navy", fontsize=10,
            ha="center", va="bottom", fontweight="bold")
    ax.text(L / UM, 0.0, V + 55, "N$^+$ contact (+$V$)", color="darkgreen",
            fontsize=8, ha="right", va="bottom", fontweight="bold")
    ax.set_xlabel("lateral  $x$  ($\\mu$m)")
    ax.set_ylabel("depth  $y$  ($\\mu$m)")
    ax.set_zlabel("potential  $\\phi$  (V)")
    ax.set_title("Merchant 1991 Fig. 3 in 3-D -- optimized SOI potential surface\n"
                 "(paper's equipotential contours projected on the floor)\n"
                 f"$t_s$=0.1 $\\mu$m, $t_{{ox}}$=2 $\\mu$m, $L$=70 $\\mu$m, "
                 f"$V$={V:.0f} V, $W$={W/UM:.1f} $\\mu$m")
    ax.set_zlim(0, V * 1.12)               # headroom for the film/contact labels
    ax.invert_yaxis()                      # depth increases into the page
    ax.view_init(elev=26, azim=-125)
    ax.set_box_aspect((2.0, 1.0, 1.1))
    ax.legend(loc="upper left", fontsize=8)
    fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.02, label="$\\phi$ (V)")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig3_potential_3d.png", dpi=150)
    plt.close(fig)


def validation_surface(t_s, t_ox, W, V):
    """3-D overlay: numerical Poisson solve (surface) vs closed form (wireframe)."""
    X, Yg, PHI = solve_poisson_2d(t_s, t_ox, W, V)
    PHIc = phi_closed(X, Yg, V, W, t_s, t_ox)
    err = PHI - PHIc
    emax = np.max(np.abs(err))
    erms = np.sqrt(np.mean(err**2))
    print(f"  numerical vs closed-form phi:  max |err| = {emax:.3f} V "
          f"({100*emax/V:.2f}% of V),  RMS = {erms:.3f} V")

    fig = plt.figure(figsize=(9.2, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X / UM, Yg / UM, PHI, cmap="plasma", rstride=3, cstride=6,
                    linewidth=0, antialiased=True, alpha=0.9)
    ax.plot_wireframe(X / UM, Yg / UM, PHIc, rstride=12, cstride=24,
                      color="k", linewidth=0.6, label="closed form (eqs. 5-9)")
    ax.set_xlabel("lateral  $x$  ($\\mu$m)")
    ax.set_ylabel("depth  $y$  ($\\mu$m)")
    ax.set_zlabel("potential  $\\phi$  (V)")
    ax.set_title("Validation: 2-D variable-$\\varepsilon$ Poisson solve (surface)\n"
                 f"vs closed form (wireframe) -- max err {emax:.2f} V "
                 f"({100*emax/V:.2f}% of $V$)")
    ax.invert_yaxis()
    ax.view_init(elev=24, azim=-122)
    ax.set_box_aspect((2.0, 1.0, 1.1))
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/fig3_potential_3d_check.png", dpi=150)
    plt.close(fig)
    return emax, erms


def main():
    t_s, t_ox, L, V = 0.1 * UM, 2.0 * UM, 70.0 * UM, 750.0
    W = find_design_W(t_s, t_ox, V)
    t = t_eff(t_s, t_ox)
    print(f"Fig.3 (3-D)  design width W = {W/UM:.2f} um  (< L=70 um),  "
          f"t_eff = {t/UM:.3f} um,  K_s/K_ox = {K_S/K_OX:.2f}")
    print(f"  surface top-edge ramp: phi(0,0)=0 V  ->  phi(W,0)={V:.0f} V "
          f"(phi ~ V x/W, linear along the surface)")
    hero_surface(t_s, t_ox, L, W, V)
    validation_surface(t_s, t_ox, W, V)
    print(f"  wrote {OUTDIR}/fig3_potential_3d.png  and  "
          f"{OUTDIR}/fig3_potential_3d_check.png")


if __name__ == "__main__":
    main()
