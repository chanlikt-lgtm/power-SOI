"""Nakagawa 1990 Fig.4 -- a THIRD method: GREEN'S FUNCTION (eigenfunction series).

Companion to `fig4_pde_standalone.py` (finite-volume Poisson) and
`fig4_conformal.py` (Schwarz-Christoffel conformal map).  This file solves the
SAME structure-B cross-section with the Dirichlet GREEN'S FUNCTION of the
rectangle, expanded in its Laplacian eigenfunctions.  The point is a mutual
check: the conformal map and the Green's function solve the *identical* idealised
boundary-value problem by completely different mathematics, so they must agree to
series-truncation precision -- which then pins the residual ~5% gap to the
finite-volume solve firmly on the shared boundary-condition idealisation, not on
either analytic method.

------------------------------------------------------------------------------
THE GREEN'S-FUNCTION SOLUTION
------------------------------------------------------------------------------
Work in the effective single-dielectric coordinate  Y(y)  (silicon: Y=y; oxide
stretched by K_si/K_ox), exactly as in fig4_conformal.py, so the two-dielectric
source-free operator becomes the plain Laplacian on the rectangle
[0,a]x[0,b],  a=Lx, b=H_eff.  Poisson's equation there is

    lap phi = -f(x,Y),     f = rho(Y)/eps_si   in the silicon band (Y<t_si), else 0.

The Dirichlet Green's function of the rectangle gives, in one shot,

    phi(r) = INT_V G(r,r') f(r') dV'   -   OINT_S phi_bc dG/dn dS'
           =  u (space-charge response, =0 on all walls)
           +  B + v'  (harmonic extension of the wall data).

Concretely, everything is a sine/sinh eigenfunction series:

  * u   -- double sine series  sum_{n,m} u_nm sin(n pi x/a) sin(m pi Y/b),
           u_nm = f_nm / ((n pi/a)^2 + (m pi/b)^2)      [the Green's-function kernel]
  * B   -- the (harmonic) bilinear  V (1-x/a)(1-Y/b)  that carries the one non-zero
           corner value (drain touches the top-left corner), so the remaining wall
           data vanishes at every corner and its sine series converge cleanly.
  * v'  -- single sine/sinh series for the reduced top and left wall data.

No PDE grid, no conformal map -- just the rectangle's eigenfunctions.  The space
charge is handled DIRECTLY by the volume Green's-function integral (no phi_p
superposition trick is needed for it here; B+v' carry only the electrode data).

------------------------------------------------------------------------------
Units: CGS-practical -- cm, V, cm^-3, F/cm.  (UM = 1e-4 cm = 1 micron.)
Run:   py fig4_green.py
Writes: ../figs/fig4_green_3d.png , ../figs/fig4_green.png ,
        ../figs/fig4_three_methods.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

# reuse the shared model pieces so the three methods are strictly comparable
from fig4_conformal import (Q, EPS0, K_SI, K_OX, UM, PARAMS,
                            particular_1d, Y_of_y, solve as conformal_solve,
                            _oxide_share)


def _y_of_Y(Yv, p=PARAMS):
    """Inverse of Y_of_y: effective depth -> true depth."""
    L_si = p["L_si"]
    Yv = np.asarray(Yv, float)
    return np.where(Yv <= L_si, Yv, L_si + (K_OX / K_SI) * (Yv - L_si))


def _sinh_ratio(coef, arg, span):
    """sinh(coef*arg)/sinh(coef*span), evaluated without overflow for large coef.
    (both arg,span >= 0, arg <= span)."""
    # = exp(-coef*(span-arg)) * (1-exp(-2 coef arg)) / (1-exp(-2 coef span))
    return (np.exp(-coef * (span - arg)) * (1.0 - np.exp(-2.0 * coef * arg))
            / (1.0 - np.exp(-2.0 * coef * span)))


def solve(V, p=PARAMS, nx=221, ny=171, N=400, charge_free=False):
    """Green's-function (eigenfunction-series) solution.  Returns (x, y, phi, meta),
    phi on the TRUE (x, y) grid for direct comparison with the other two methods."""
    L_si, t_ox, Lx = p["L_si"], p["t_ox"], p["Lx"]
    H = L_si + t_ox
    a = Lx
    b = L_si + (K_SI / K_OX) * t_ox                    # effective height H_eff
    eps_si = K_SI * EPS0

    # output grid (true coords) and its effective depths
    x = np.linspace(0.0, a, nx)
    y = np.linspace(0.0, H, ny)
    Ye = Y_of_y(y, p)                                  # effective depth of each true y

    n = np.arange(1, N + 1)                            # mode indices
    kx = n * np.pi / a
    kY = n * np.pi / b

    # basis matrices on the output grid
    SXx = np.sin(np.outer(x, kx))                      # (nx, N)  sin(n pi x/a)
    SYe = np.sin(np.outer(Ye, kY))                     # (ny, N)  sin(m pi Y/b)

    # ---- 1-D charged profile phi_p(y) (for the left-wall closure & metadata) ----
    y1d, phi_p = particular_1d(V, p, charge_free=charge_free)
    phip_if = float(np.interp(L_si, y1d, phi_p))

    # ---- (u) space-charge response: double sine series ----------------------
    # f depends on Y only -> x-part is the sine expansion of 1 (odd n only)
    bn = np.where(n % 2 == 1, 4.0 / (n * np.pi), 0.0)  # (2/a) INT sin(n pi x/a) dx
    if charge_free:
        em = np.zeros(N)
    else:
        Yq = np.linspace(0.0, L_si, 4001)              # silicon band (Y=y here)
        NY = np.where(Yq >= L_si - p["w_n"], p["N_D0"] + p["S_D"] / p["w_n"],
                      p["N_D0"])
        fY = Q * NY / eps_si                            # source density (1/cm^2 * ...)
        # e_m = (2/b) INT_0^{t_si} fY sin(m pi Y/b) dY
        em = (2.0 / b) * np.trapz(fY[None, :] * np.sin(np.outer(kY, Yq)), Yq, axis=1)
    denom = kx[:, None] ** 2 + kY[None, :] ** 2        # (N,N)  eigenvalues
    U = (bn[:, None] * em[None, :]) / denom            # u_nm  (N,N)
    u = SYe @ U.T @ SXx.T                               # (ny, nx)

    # ---- (B) harmonic bilinear carrying the top-left (drain) corner value ----
    Xg, Yg = np.meshgrid(x, Ye)                         # (ny, nx) effective
    B = V * (1.0 - Xg / a) * (1.0 - Yg / b)

    # ---- (v') reduced harmonic wall data (vanishes at every corner) ----------
    # top wall  Y=0:  data = top_actual(x) - V(1 - x/a)
    xf = np.linspace(0.0, a, 6001)
    a_d, b_d = p["drain_w"], Lx - p["source_w"]
    top_actual = np.where(xf <= a_d, V,
                  np.where(xf >= b_d, 0.0,
                           V * 0.5 * (1 + np.cos(np.pi * (xf - a_d) / (b_d - a_d)))))
    dtop = top_actual - V * (1.0 - xf / a)
    An = (2.0 / a) * np.trapz(dtop[None, :] * np.sin(np.outer(kx, xf)), xf, axis=1)
    # v'_top = sum_n A_n sin(n pi x/a) sinh(n pi (b-Y)/a)/sinh(n pi b/a)
    Ttop = _sinh_ratio(kx[None, :], (b - Ye)[:, None], b)     # (ny, N)
    vtop = (Ttop * An[None, :]) @ SXx.T                        # (ny, nx)

    # left wall  x=0:  data = phi_p(y(Y)) - V(1 - Y/b)
    Yf = np.linspace(0.0, b, 6001)
    phip_left = np.interp(_y_of_Y(Yf, p), y1d, phi_p)
    dleft = phip_left - V * (1.0 - Yf / b)
    Cm = (2.0 / b) * np.trapz(dleft[None, :] * np.sin(np.outer(kY, Yf)), Yf, axis=1)
    # v'_left = sum_m C_m sin(m pi Y/b) sinh(m pi (a-x)/b)/sinh(m pi a/b)
    Lft = _sinh_ratio(kY[None, :], (a - x)[:, None], a)        # (nx, N)
    vleft = (SYe * Cm[None, :]) @ Lft.T                        # (ny, nx)

    phi = u + B + vtop + vleft
    meta = dict(N=N, H_eff=b, phip_if=phip_if, eps_si=eps_si,
                em=em, denom=denom, u=u, charge_free=charge_free)
    return x, y, phi, meta


# ============================================================================
#  VERIFY
# ============================================================================
def verify(V=500.0):
    print("=" * 70)
    print("VERIFICATION of the GREEN'S-FUNCTION solution (Nakagawa Fig.4)")
    print("=" * 70)
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    H = L_si + t_ox
    x, y, phi, meta = solve(V, nx=221, ny=171, N=400)

    print(f"(1) eigenfunction series: {meta['N']} modes per direction"
          f"   (Green's function G = sum sin.sin / ((n pi/a)^2+(m pi/b)^2))")

    # (2) Poisson residual in EFFECTIVE coords: lap phi + f == 0 (validates series)
    b = meta["H_eff"]
    nxg = nYg = 161
    xs = np.linspace(0, Lx, nxg); Ys = np.linspace(0, b, nYg)
    # evaluate the series directly at effective depths Ys (uniform effective grid)
    phiE = _eval_effective(V, xs, Ys, meta["N"])
    dx = xs[1] - xs[0]; dY = Ys[1] - Ys[0]
    lap = ((phiE[1:-1, 2:] - 2*phiE[1:-1, 1:-1] + phiE[1:-1, :-2]) / dx**2
           + (phiE[2:, 1:-1] - 2*phiE[1:-1, 1:-1] + phiE[:-2, 1:-1]) / dY**2)
    # source f on the interior effective grid (silicon band Y<t_si)
    Xi, Yi = np.meshgrid(xs[1:-1], Ys[1:-1])
    NY = np.where((Yi < L_si) & (Yi >= L_si - PARAMS["w_n"]),
                  PARAMS["N_D0"] + PARAMS["S_D"] / PARAMS["w_n"],
                  np.where(Yi < L_si, PARAMS["N_D0"], 0.0))
    f = Q * NY / meta["eps_si"]
    res = np.abs(lap + f)
    scale = V / b**2
    # exclude the electrode-edge singularities on the top wall
    near = (((np.abs(Xi - PARAMS["drain_w"]) < 3*UM) |
             (np.abs(Xi - (Lx-PARAMS["source_w"])) < 3*UM)) & (Yi < 5*UM))
    clean = (res / scale)[~near]
    print(f"(2) Poisson residual (effective coords, away from electrode edges):"
          f"\n    max|lap phi + f|/(V/H_eff^2) = {clean.max():.2e},"
          f"  RMS = {np.sqrt(np.mean(clean**2)):.2e}  (-> 0 validates the series)")

    # (3) boundary conditions
    a_d, b_d = PARAMS["drain_w"], Lx - PARAMS["source_w"]
    top = phi[0, :]
    print(f"(3) top-surface BC:  drain max|phi-V| = {np.abs(top[x<=a_d]-V).max():.2f} V,"
          f"   source max|phi| = {np.abs(top[x>=b_d]).max():.2f} V")
    print(f"    right p+ edge max|phi| = {np.abs(phi[y<L_si, -1]).max():.2f} V,"
          f"   substrate max|phi| = {np.abs(phi[-1, :]).max():.2f} V")

    # (4) physical range
    print(f"(4) potential range  min={phi.min():.2f} V  max={phi.max():.2f} V")

    # (5) oxide voltage share, and agreement with the conformal method
    sc, ix = _oxide_share(x, y, phi)
    print(f"(5) oxide voltage share under drain centre = {sc:.3f}"
          f"   ({'MORE' if sc>0.5 else 'less'} than half -- paper: 'more than half')")
    xc, yc, phic, _ = conformal_solve(V, nx=221, ny=171)
    diff = phi - phic
    print(f"(6) Green's function vs conformal map (same idealised BVP):"
          f"\n    max|Green - conformal| = {np.abs(diff).max():.2f} V,"
          f"  RMS = {np.sqrt(np.mean(diff**2)):.2f} V"
          f"   (-> 0: two independent analytic methods agree)")
    print("=" * 70)
    return x, y, phi, meta


def _eval_effective(V, xs, Ys, N, p=PARAMS, charge_free=False):
    """Evaluate the Green's-function series directly at effective coords (xs, Ys)
    -- used only by the Poisson-residual check (harmonic + charge parts)."""
    L_si, t_ox, Lx = p["L_si"], p["t_ox"], p["Lx"]
    a = Lx; b = L_si + (K_SI / K_OX) * t_ox
    eps_si = K_SI * EPS0
    n = np.arange(1, N + 1); kx = n * np.pi / a; kY = n * np.pi / b
    SXx = np.sin(np.outer(xs, kx)); SYe = np.sin(np.outer(Ys, kY))
    y1d, phi_p = particular_1d(V, p, charge_free=charge_free)

    bn = np.where(n % 2 == 1, 4.0 / (n * np.pi), 0.0)
    Yq = np.linspace(0.0, L_si, 4001)
    NY = np.where(Yq >= L_si - p["w_n"], p["N_D0"] + p["S_D"] / p["w_n"], p["N_D0"])
    fY = Q * NY / eps_si
    em = (2.0 / b) * np.trapz(fY[None, :] * np.sin(np.outer(kY, Yq)), Yq, axis=1)
    U = (bn[:, None] * em[None, :]) / (kx[:, None]**2 + kY[None, :]**2)
    u = SYe @ U.T @ SXx.T
    Xg, Yg = np.meshgrid(xs, Ys)
    B = V * (1 - Xg/a) * (1 - Yg/b)
    xf = np.linspace(0, a, 6001); a_d, b_d = p["drain_w"], Lx - p["source_w"]
    top_actual = np.where(xf <= a_d, V, np.where(xf >= b_d, 0.0,
                          V*0.5*(1+np.cos(np.pi*(xf-a_d)/(b_d-a_d)))))
    An = (2/a)*np.trapz((top_actual - V*(1-xf/a))[None, :]*np.sin(np.outer(kx, xf)),
                        xf, axis=1)
    Ttop = _sinh_ratio(kx[None, :], (b - Ys)[:, None], b)
    vtop = (Ttop * An[None, :]) @ SXx.T
    Yf = np.linspace(0, b, 6001)
    dleft = np.interp(_y_of_Y(Yf, p), y1d, phi_p) - V*(1-Yf/b)
    Cm = (2/b)*np.trapz(dleft[None, :]*np.sin(np.outer(kY, Yf)), Yf, axis=1)
    Lft = _sinh_ratio(kY[None, :], (a - xs)[:, None], a)
    vleft = (SYe * Cm[None, :]) @ Lft.T
    return u + B + vtop + vleft


# ============================================================================
#  PLOTS
# ============================================================================
def make_plots(V=500.0):
    x, y, phi, _ = solve(V, nx=221, ny=171, N=400)
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    Xum, Yum = np.meshgrid(x / UM, y / UM)
    i_if = int(np.argmin(np.abs(y - L_si)))
    ox_share, ix = _oxide_share(x, y, phi)

    # 3-D surface
    fig = plt.figure(figsize=(9.4, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    sx = max(1, (len(x) - 1) // 60); sy = max(1, (len(y) - 1) // 45)
    surf = ax.plot_surface(Xum[::sy, ::sx], Yum[::sy, ::sx], phi[::sy, ::sx],
                           rcount=45, ccount=60, cmap="turbo",
                           edgecolor="k", linewidth=0.22, alpha=0.9)
    ax.plot(x / UM, np.full_like(x, L_si / UM), phi[i_if, :], color="crimson",
            lw=2.2, label="Si / buried-oxide interface")
    v_ox = phi[i_if, 0]
    ax.plot(np.zeros(2), [L_si/UM, (L_si+t_ox)/UM], [v_ox, 0.0], color="magenta",
            lw=3.2, label="oxide voltage drop (drain edge)")
    ax.set_xlabel("x  ($\\mu$m)"); ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_zlabel("potential (V)"); ax.set_zlim(0, V)
    ax.set_title("Nakagawa 1990 Fig.4 by GREEN'S FUNCTION (structure B, 500 V)\n"
                 f"rectangle eigenfunction series -- oxide carries "
                 f"$\\approx${100*ox_share:.0f}% of $V$")
    ax.view_init(elev=28, azim=-60); ax.invert_yaxis()
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.02, label="potential (V)")
    fig.tight_layout()
    out3d = "../figs/fig4_green_3d.png"
    fig.savefig(out3d, dpi=150); plt.close(fig); print("wrote", out3d)

    # 2-D contour + vertical cut
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                   gridspec_kw={"height_ratios": [2.3, 1]})
    levels = np.arange(0, V + 1, 25)
    cf = ax1.contourf(Xum, Yum, phi, levels=levels, cmap="turbo")
    ax1.contour(Xum, Yum, phi, levels=levels, colors="k", linewidths=0.4)
    ax1.axhline(20.0, color="white", lw=1.5, ls="--")
    ax1.text(25, 21.6, "buried oxide", color="white", ha="center", fontsize=9)
    ax1.text(4, 2.0, "drain 500 V", color="white", fontsize=9, weight="bold")
    ax1.text(40, 2.0, "source 0 V", color="white", fontsize=9, weight="bold")
    ax1.invert_yaxis()
    ax1.set_xlabel("x  ($\\mu$m)"); ax1.set_ylabel("depth y  ($\\mu$m)")
    ax1.set_title("Fig.4 by Green's function — 2-D map (structure B, 500 V)")
    fig.colorbar(cf, ax=ax1, label="potential (V)", shrink=0.9)
    col = phi[:, ix]
    ax2.plot(y / UM, col, "b-", lw=2)
    ax2.axvspan(20, 23, color="#d0d0d0", alpha=0.6); ax2.axvline(20, color="k", ls="--", lw=0.7)
    v_ox = float(np.interp(L_si, y, col))
    ax2.annotate(f"{ox_share*100:.0f}% of {V:.0f} V\nacross the oxide",
                 xy=(21.5, v_ox/2), xytext=(11, 120),
                 arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax2.text(9, 380, "Silicon", fontsize=9); ax2.text(21.5, 380, "Oxide", fontsize=9, ha="center")
    ax2.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain centre)")
    ax2.set_ylabel("potential (V)"); ax2.set_xlim(0, 23); ax2.grid(alpha=0.25)
    fig.tight_layout()
    out2d = "../figs/fig4_green.png"
    fig.savefig(out2d, dpi=150); plt.close(fig); print("wrote", out2d)


def compare_three_methods(V=500.0):
    """Overlay all three methods: 1-D-free finite-volume Poisson, conformal map,
    and Green's function -- vertical cut under the drain + share table."""
    import fig4_pde_standalone as fv
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    xg, yg, phig, _ = solve(V, nx=221, ny=171, N=400)             # Green
    xc, yc, phic, _ = conformal_solve(V, nx=221, ny=171)          # conformal
    xf, yf, phif = fv.solve_poisson(201, 185, V)                  # finite volume
    sg, ig = _oxide_share(xg, yg, phig)
    sc, ic = _oxide_share(xc, yc, phic)
    sf, iff = _oxide_share(xf, yf, phif, fv.PARAMS)

    yq = np.linspace(0, L_si + t_ox, 400)
    cg = np.interp(yq, yg, phig[:, ig])
    cc = np.interp(yq, yc, phic[:, ic])
    cf_ = np.interp(yq, yf, phif[:, iff])

    print("-" * 70)
    print("THREE-METHOD comparison -- oxide voltage share under the drain centre:")
    print(f"    Green's function      = {sg:.3f}")
    print(f"    conformal map         = {sc:.3f}")
    print(f"    finite-volume Poisson = {sf:.3f}   (paper: 'more than half')")
    print(f"    |Green - conformal| RMS on the cut = "
          f"{np.sqrt(np.mean((cg-cc)**2)):.1f} V   (two analytic methods)")
    print(f"    |analytic - FV|      RMS on the cut ~ "
          f"{np.sqrt(np.mean((0.5*(cg+cc)-cf_)**2)):.1f} V   (BC idealisation)")

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.plot(yq/UM, cf_, color="0.35", lw=3.2, label=f"finite-volume Poisson  ({100*sf:.0f}%)")
    ax.plot(yq/UM, cc, "b-", lw=2.0, label=f"conformal map  ({100*sc:.0f}%)")
    ax.plot(yq/UM, cg, "r--", lw=2.0, label=f"Green's function  ({100*sg:.0f}%)")
    ax.axvspan(20, 23, color="#d0d0d0", alpha=0.6); ax.axvline(20, color="k", ls=":", lw=0.8)
    ax.text(9, 460, "Silicon", fontsize=10); ax.text(21.4, 460, "Oxide", fontsize=10, ha="center")
    ax.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain centre)")
    ax.set_ylabel("potential (V)")
    ax.set_title("Nakagawa Fig.4 — three independent methods\n"
                 "conformal map and Green's function coincide; both track the FV solve")
    ax.set_xlim(0, 23); ax.grid(alpha=0.25); ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out = "../figs/fig4_three_methods.png"
    fig.savefig(out, dpi=150); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    verify(500.0)
    make_plots(500.0)
    compare_three_methods(500.0)
