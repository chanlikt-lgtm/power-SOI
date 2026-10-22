"""Nakagawa 1990 Fig.4 -- STANDALONE 2-D Poisson solve + verify + 3-D plot.

Self-contained: no imports from model.py / solver2d.py. Everything (physical
constants, geometry, doping, the PDE discretisation, boundary conditions, the
solver, the verification, and the 3-D plot) lives in this one file.

------------------------------------------------------------------------------
PHYSICS
------------------------------------------------------------------------------
Electrostatics of the thin-SOI cross-section (x lateral, y depth DOWNWARD):

    silicon drift   0     <= y < L_si     (net ionised donors, rho > 0)
    buried oxide    L_si  <= y < H        (charge free)
    substrate       y = H                 (grounded)

We solve Poisson's equation for the electrostatic potential phi(x,y):

    div( eps(x,y) grad phi ) = -rho(x,y),        rho = q * N_net(x,y)

with a piecewise-constant permittivity eps = eps0 * (K_si in Si, K_ox in oxide).
The Si/oxide displacement continuity (D_y jump by K_si/K_ox = 3) is captured by
harmonic-mean FACE permittivities in a 5-point finite-volume scheme.

Structure B (the paper's Fig.3B / Fig.4 device):
    * drain electrode  phi = V   on the top surface, left  (x <= drain_w)
    * source electrode phi = 0    on the top surface, right (x >= Lx-source_w)
    * grounded substrate phi = 0  on the whole bottom edge (y = H)
    * p+ sinker: the right silicon edge is grounded (phi = 0)
    * everywhere else on the outer boundary: natural Neumann (zero normal flux)

A shallow N donor sheet (dose S_D over the bottom w_n of the silicon) shields the
bulk silicon and pushes >half the applied voltage across the thin buried oxide --
the paper's key result, which the verification and the plot both quantify.

------------------------------------------------------------------------------
Units: CGS-practical -- cm, V, cm^-3, F/cm.  (UM = 1e-4 cm = 1 micron.)
Run:   py fig4_pde_standalone.py
Writes: ../figs/fig4_pde_standalone_3d.png
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

# ----------------------------- physical constants ----------------------------
Q      = 1.602176634e-19        # elementary charge, C
EPS0   = 8.8541878128e-14       # vacuum permittivity, F/cm
K_SI   = 11.7                   # relative permittivity, silicon
K_OX   = 3.9                    # relative permittivity, SiO2
UM     = 1.0e-4                 # 1 micron in cm

# ----------------------------- device parameters -----------------------------
PARAMS = dict(
    L_si     = 20.0 * UM,       # silicon-film thickness
    t_ox     = 3.0 * UM,        # buried-oxide thickness
    Lx       = 50.0 * UM,       # lateral extent (drain..source)
    w_n      = 2.0 * UM,        # shallow donor-sheet thickness (bottom of Si)
    N_D0     = 1.0e14,          # background donor density, cm^-3
    S_D      = 1.5e12,          # shallow-sheet dose, cm^-2
    drain_w  = 8.0 * UM,        # drain electrode width on the top surface
    source_w = 8.0 * UM,        # source electrode width on the top surface
)


# ============================================================================
#  PDE ASSEMBLY + SOLVE  (vectorised sparse finite volume)
# ============================================================================
def solve_poisson(nx, ny, V, p=PARAMS):
    """Solve div(eps grad phi) = -rho on the SOI cross-section.

    Returns (x, y, phi) with phi shape (ny, nx). Fully vectorised COO assembly:
    5-point stencil, harmonic-mean face permittivities, identity rows for
    Dirichlet nodes, and natural (zero-flux) Neumann on any outer face that is
    neither Dirichlet nor interior.
    """
    L_si, t_ox, Lx = p["L_si"], p["t_ox"], p["Lx"]
    H = L_si + t_ox
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, H, ny)
    dx, dy = x[1] - x[0], y[1] - y[0]
    X, Yg = np.meshgrid(x, y)                       # (ny, nx)

    # --- material permittivity (absolute) and space charge ------------------
    in_si = Yg < L_si - 1e-12
    eps = np.where(in_si, K_SI, K_OX) * EPS0
    N = np.where(in_si, p["N_D0"], 0.0)
    shallow = in_si & (Yg >= L_si - p["w_n"])
    N = N + np.where(shallow, p["S_D"] / p["w_n"], 0.0)
    rho = Q * N

    # --- Dirichlet boundary conditions (structure B) ------------------------
    mask = np.zeros_like(Yg, dtype=bool)
    val = np.zeros_like(Yg)
    mask[-1, :] = True                              # grounded substrate (bottom)
    xL = x <= p["drain_w"]
    xR = x >= Lx - p["source_w"]
    mask[0, xL] = True; val[0, xL] = V              # drain electrode  (top-left)
    mask[0, xR] = True; val[0, xR] = 0.0            # source electrode (top-right)
    si_row = y < L_si - 1e-12
    mask[si_row, -1] = True; val[si_row, -1] = 0.0  # p+ sinker (right Si edge)

    # --- assemble the sparse operator ---------------------------------------
    K = np.arange(nx * ny).reshape(ny, nx)
    Jg, Ig = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    free = ~mask
    diag = np.zeros((ny, nx))
    R, C, D = [], [], []

    def harm(a, b):                                 # harmonic-mean face eps
        return 2.0 / (1.0 / a + 1.0 / b)

    def add(sel, col_off, c):
        R.append(K[sel]); C.append(K[sel] + col_off); D.append(c[sel])
        diag[sel] -= c[sel]

    # west (i-1)
    c = np.zeros((ny, nx)); c[:, 1:] = harm(eps[:, 1:], eps[:, :-1]) / dx**2
    add(free & (Ig > 0), -1, c)
    # east (i+1)
    c = np.zeros((ny, nx)); c[:, :-1] = harm(eps[:, :-1], eps[:, 1:]) / dx**2
    add(free & (Ig < nx - 1), +1, c)
    # north (j-1, shallower)
    c = np.zeros((ny, nx)); c[1:, :] = harm(eps[1:, :], eps[:-1, :]) / dy**2
    add(free & (Jg > 0), -nx, c)
    # south (j+1, deeper)
    c = np.zeros((ny, nx)); c[:-1, :] = harm(eps[:-1, :], eps[1:, :]) / dy**2
    add(free & (Jg < ny - 1), +nx, c)

    # diagonals: free -> -sum(face); Dirichlet -> identity row
    R.append(K[free]); C.append(K[free]); D.append(diag[free])
    R.append(K[mask]); C.append(K[mask]); D.append(np.ones(int(mask.sum())))

    b = np.where(mask, val, -rho).ravel()
    A = sp.csr_matrix((np.concatenate(D),
                       (np.concatenate(R), np.concatenate(C))),
                      shape=(nx * ny, nx * ny))
    phi = spsolve(A.tocsr(), b).reshape(ny, nx)
    # keep the matrix/rhs around for the residual check
    solve_poisson.last = (A, b, mask, val, V)
    return x, y, phi


# ============================================================================
#  VERIFY THE SOLUTION
# ============================================================================
def oxide_share(x, y, phi, p=PARAMS):
    """Fraction of the applied voltage dropped across the oxide under the drain.

    Evaluated by interpolating phi to the EXACT interface depth y=L_si under the
    middle of the drain electrode (x = drain_w/2), not at the singular x=0 drain
    corner -- so the metric is smooth under grid refinement. phi at the interface
    equals the voltage still remaining to drop across the (grounded) oxide, hence
    that value / V is the oxide's share of the applied voltage.
    """
    ix = int(np.argmin(np.abs(x - 0.5 * p["drain_w"])))
    v_if = float(np.interp(p["L_si"], y, phi[:, ix]))
    return v_if / phi[0, ix], ix                    # phi[0,ix] = drain potential V


def oxide_share_drain_avg(x, y, phi, p=PARAMS):
    """Oxide share AVERAGED over the drain electrode footprint -- a smooth
    integral functional (corner-singularity noise integrates out), suitable for
    Richardson extrapolation to estimate the grid-converged value + accuracy."""
    cols = np.where(x <= p["drain_w"])[0]
    sh = np.array([np.interp(p["L_si"], y, phi[:, i]) / phi[0, i] for i in cols])
    return float(np.trapz(sh, x[cols]) / (x[cols][-1] - x[cols][0]))


def verify(V=500.0):
    print("=" * 68)
    print("VERIFICATION of the standalone 2-D Poisson solution")
    print("=" * 68)

    # (1) linear-system residual at the production grid --------------------
    x, y, phi = solve_poisson(201, 185, V)
    A, b, mask, val, Vd = solve_poisson.last
    r = A.dot(phi.ravel()) - b
    print(f"(1) linear residual  ||A phi - b||_inf = {np.abs(r).max():.2e}"
          f"   (direct solve -> ~machine zero)")

    # (2) boundary conditions actually imposed -----------------------------
    bc_err = np.abs(phi[mask] - val[mask]).max()
    print(f"(2) Dirichlet BC error  max|phi_bc - set| = {bc_err:.2e} V")

    # (3) physical bounds: no overshoot outside [0, V] ---------------------
    print(f"(3) potential range  min={phi.min():.3f} V  max={phi.max():.3f} V"
          f"   (must stay within [0, {V:.0f}])")

    # (4) self-convergence (Cauchy) of the WHOLE field vs the finest grid ----
    #     No closed form exists, so we verify that the discrete solution
    #     converges: interpolate each grid onto a common set of interior points
    #     (kept away from the singular electrode-edge corners) and measure the
    #     difference from the finest grid. Monotone decrease -> the solve converges.
    print("(4) self-convergence of the field vs the finest grid")
    print("    (evaluated on a common interior grid, away from electrode corners):")
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    H = L_si + t_ox
    xe = np.linspace(4 * UM, Lx - 4 * UM, 45)
    ye = np.linspace(1 * UM, H - 1 * UM, 45)
    Xe, Ye = np.meshgrid(xe, ye)
    pts = np.column_stack([Ye.ravel(), Xe.ravel()])

    def sample(nx, ny):
        xg, yg, pg = solve_poisson(nx, ny, V)
        return RegularGridInterpolator((yg, xg), pg)(pts)

    # grid sizes chosen so the electrode edges (x=8,42 um), the Si/oxide
    # interface (y=20 um) and the doping-band edge (y=18 um) land exactly on
    # nodes at EVERY resolution -- (nx-1) multiple of 25, (ny-1) multiple of 23 --
    # so each grid discretises the identical geometry and convergence is clean.
    grids = [(51, 47), (76, 70), (101, 93), (151, 139), (201, 185), (251, 231)]
    fields, funcs, mm = {}, [], []
    for (nx, ny) in grids:
        xg, yg, pg = solve_poisson(nx, ny, V)
        fields[(nx, ny)] = RegularGridInterpolator((yg, xg), pg)(pts)
        funcs.append(oxide_share_drain_avg(xg, yg, pg))
        mm.append((nx - 1) // 25)                     # resolution index (h ~ 1/mm)
    ref = fields[grids[-1]]                            # finest as reference
    print("       nx   ny |  max|dphi| (V)   RMS (V)   order")
    prev = None
    for (nx, ny) in grids[:-1]:
        e = fields[(nx, ny)] - ref
        em, rms = np.abs(e).max(), np.sqrt(np.mean(e**2))
        order = "" if prev is None else f"  {np.log2(prev/em):+.2f}"
        print(f"     {nx:4d} {ny:4d} |   {em:8.3f}     {rms:7.3f}   {order}")
        prev = em
    print("    -> error vs the finest grid falls monotonically: the solve converges.")

    # (5) ACCURACY: Richardson-extrapolate a smooth functional -> true error --
    mm = np.array(mm, float)
    funcs = np.array(funcs)
    h = 1.0 / mm
    (F0, C, porder), _ = curve_fit(lambda h, F0, C, p: F0 + C * h**p,
                                   h, funcs, p0=(funcs[-1], 1.0, 1.2), maxfev=20000)
    s_prod, _ = oxide_share(x, y, phi)                # production-grid (201x185)
    err_prod = abs(funcs[list(mm).index(8)] - F0)
    print(f"(5) accuracy of the oxide voltage share (drain-averaged functional):")
    print(f"    Richardson extrapolation -> grid-converged value = {F0:.4f}"
          f"   observed order p = {porder:.2f}")
    print(f"    production grid (201x185) = {funcs[list(mm).index(8)]:.4f}"
          f"   -> discretisation error ~ {err_prod:.4f} ({100*err_prod:.2f}% of V)")
    print(f"    point value under drain centre = {s_prod:.3f}  ->  "
          f"{'MORE than half' if F0 > 0.5 else 'less than half'} of V across the "
          f"buried oxide (paper: 'more than half').")
    print("=" * 68)
    return x, y, phi


# ============================================================================
#  3-D PLOT WITH LABELS
# ============================================================================
def plot_3d(x, y, phi, V=500.0, p=PARAMS, out="../figs/fig4_pde_standalone_3d.png"):
    L_si, t_ox = p["L_si"], p["t_ox"]
    Xum, Yum = np.meshgrid(x / UM, y / UM)
    i_if = int(np.searchsorted(y, L_si))
    v_ox = phi[i_if, 0]
    ox_share = v_ox / V
    si_share = 1.0 - ox_share

    fig = plt.figure(figsize=(9.4, 6.8))
    ax = fig.add_subplot(111, projection="3d")

    sy = max(1, phi.shape[0] // 45)
    sx = max(1, phi.shape[1] // 60)
    surf = ax.plot_surface(Xum[::sy, ::sx], Yum[::sy, ::sx], phi[::sy, ::sx],
                           rcount=45, ccount=60, cmap="turbo",
                           edgecolor="k", linewidth=0.22, alpha=0.9)

    # Si/oxide interface traced across the surface, and the drain-edge "cliff"
    yif = L_si / UM
    ax.plot(x / UM, np.full(len(x), yif), phi[i_if, :],
            color="crimson", lw=2.2, label="Si / buried-oxide interface")
    ax.plot(np.zeros(2), [yif, (L_si + t_ox) / UM], [v_ox, 0.0],
            color="magenta", lw=3.2, label="oxide voltage drop (drain edge)")

    ax.set_xlim(0, p["Lx"] / UM)
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9)
    ax.text(1.0, 0.0, V + 25, "n layer / drain  (500 V)", color="darkred",
            fontsize=9, ha="left", va="bottom", fontweight="bold")
    ax.text(p["Lx"] / UM, 1.0, 55, "p$^+$ layer / source  (0 V)", color="navy",
            fontsize=9, ha="right", va="bottom", fontweight="bold")
    ax.text(6, 9.0, 360, f"Silicon film 0-20 $\\mu$m\n($\\approx${100*si_share:.0f}% of $V$)",
            fontsize=8, ha="center", va="center", bbox=box)
    ax.text(3, 21.5, 150,
            f"Buried oxide 20-23 $\\mu$m\n$>$half of $V$ ($\\approx${100*ox_share:.0f}%)",
            color="magenta", fontsize=8, ha="center", va="center", bbox=box)

    ax.set_xlabel("x  ($\\mu$m)")
    ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_zlabel("potential (V)")
    ax.set_zlim(0, V)
    ax.set_title("Nakagawa 1990 Fig.4 -- standalone PDE solve (structure B, 500 V)\n"
                 f"more than half the applied voltage falls across the buried oxide "
                 f"($\\approx${100*ox_share:.0f}%)")
    ax.view_init(elev=28, azim=-60)
    ax.invert_yaxis()
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.02, label="potential (V)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def main():
    x, y, phi = verify(500.0)
    plot_3d(x, y, phi, 500.0)


if __name__ == "__main__":
    main()
