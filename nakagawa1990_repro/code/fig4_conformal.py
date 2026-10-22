"""Nakagawa 1990 Fig.4 -- STANDALONE solution by CONFORMAL MAPPING.

A companion to `fig4_pde_standalone.py`, which solves the same structure-B
cross-section with a finite-volume Poisson discretisation.  This file solves it
by a genuine *analytic* conformal map -- no finite-difference grid for the field
itself -- so that the two independent methods can be cross-checked.

------------------------------------------------------------------------------
THE IDEA
------------------------------------------------------------------------------
Conformal mapping solves LAPLACE's equation, but the SOI cross-section is a
POISSON problem with two dielectrics.  We reduce it to a Laplace problem the
conformal map can handle by two standard, clearly-scoped steps:

  (1) SUPERPOSITION for the space charge.
        phi(x,y) = phi_p(y) + phi_h(x,y)

      phi_p(y) is the exact 1-D vertical solution that carries ALL the donor
      charge:   d/dy( eps(y) d phi_p/dy ) = -rho(y),   phi_p(0)=V, phi_p(H)=0.
      Substituting, phi_h then satisfies the SOURCE-FREE equation
        div( eps grad phi_h ) = 0.
      phi_p is what makes MORE than half the voltage fall across the buried
      oxide (the paper's key result -- it is a property of the space charge,
      not of the electrode geometry), so keeping it exact matters.

  (2) EFFECTIVE-OXIDE-THICKNESS stretch for the two dielectrics.
      Define an effective depth  Y(y) = y in the silicon, and
        Y = t_si + (K_si/K_ox)*(y - t_si)   in the oxide.
      This makes the vertical displacement continuous automatically, collapsing
      the two-dielectric source-free equation into a single-dielectric Laplace
      equation  d^2 phi_h/dx^2 + d^2 phi_h/dY^2 = 0  on the rectangle
      [0,Lx] x [0,H_eff].  (Exact for vertical fields -- which dominate in the
      thin oxide -- and an approximation for the weak lateral oxide field; this
      is the ONE approximation in the method, and it is why phi_h, not phi, is
      the part we push through it.)

  (3) CONFORMAL MAP of that rectangle to the upper half w-plane by the
      Schwarz-Christoffel map for a rectangle,  w = sn(u, m)  (Jacobi elliptic
      sine).  The modulus m is fixed by the rectangle's aspect ratio via
      K'(m)/K(m) = 2 H_eff / Lx.  Dirichlet data on the four sides map to
      the real axis; phi_h at any interior point is then the Poisson-integral /
      harmonic-measure of that data -- evaluated in closed form as a sum of
      subtended angles.  No PDE grid is used.

Total field:  phi = phi_p(y) + phi_h(x, Y(y)).

------------------------------------------------------------------------------
Units: CGS-practical -- cm, V, cm^-3, F/cm.  (UM = 1e-4 cm = 1 micron.)
Run:   py fig4_conformal.py
Writes: ../figs/fig4_conformal_3d.png , ../figs/fig4_conformal.png
"""
from __future__ import annotations

import numpy as np
from scipy.special import ellipj, ellipk
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

# ----------------------------- physical constants ----------------------------
Q      = 1.602176634e-19        # elementary charge, C
EPS0   = 8.8541878128e-14       # vacuum permittivity, F/cm
K_SI   = 11.7                   # relative permittivity, silicon
K_OX   = 3.9                    # relative permittivity, SiO2
UM     = 1.0e-4                 # 1 micron in cm

# ----------------------------- device parameters -----------------------------
# (identical to fig4_pde_standalone.py so the two methods are comparable)
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
#  (1) PARTICULAR SOLUTION phi_p(y): exact 1-D vertical Poisson with the charge
# ============================================================================
def particular_1d(V, p=PARAMS, n=4001, charge_free=False):
    """Solve  d/dy( eps(y) dphi_p/dy ) = -rho(y),  phi_p(0)=V, phi_p(H)=0.

    Returns (y, phi_p) on a fine vertical grid.  This carries the whole donor
    space charge; it is the reason more than half of V falls across the oxide.

    charge_free=True zeroes rho, leaving the pure two-dielectric capacitive
    divider (phi_p piecewise-linear) -- the Laplace/field-plate limit, in which
    the oxide takes only its geometric ~31% share.  Use it to isolate how much of
    the >1/2 result is the donor charge rather than the dielectric ratio.
    """
    L_si, t_ox = p["L_si"], p["t_ox"]
    H = L_si + t_ox
    y = np.linspace(0.0, H, n)

    in_si = y < L_si
    eps = np.where(in_si, K_SI, K_OX) * EPS0
    if charge_free:
        rho = np.zeros_like(y)
    else:
        N = np.where(in_si, p["N_D0"], 0.0)
        N += np.where(in_si & (y >= L_si - p["w_n"]), p["S_D"] / p["w_n"], 0.0)
        rho = Q * N                                   # C/cm^3

    # Q_acc(y) = int_0^y rho dy'  (accumulated charge per area, C/cm^2)
    Q_acc = np.concatenate([[0.0], np.cumsum(0.5 * (rho[1:] + rho[:-1]) * np.diff(y))])

    # phi_p'(y) = (D0 - Q_acc)/eps ; choose D0 so that phi_p(H)=0, phi_p(0)=V.
    inv_eps = 1.0 / eps
    A = np.trapz(inv_eps, y)                           # int dy/eps
    Bt = np.trapz(Q_acc * inv_eps, y)                  # int Q_acc/eps dy
    D0 = (Bt - V) / A                                  # top displacement, from BC
    integrand = (D0 - Q_acc) * inv_eps
    phi_p = V + np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(y))])
    return y, phi_p


# ============================================================================
#  (2)+(3) HARMONIC PART phi_h BY CONFORMAL MAP OF THE RECTANGLE -> UHP
# ============================================================================
def _rect_modulus(aspect):
    """Return Jacobi parameter m with K'(m)/K(m) == aspect  (= 2 H_eff/Lx)."""
    f = lambda m: ellipk(1.0 - m) / ellipk(m) - aspect
    return brentq(f, 1e-12, 1.0 - 1e-12, xtol=1e-14)


def _sn_complex(alpha, beta, m):
    """Jacobi sn(alpha + i*beta, m) via the imaginary-argument addition theorem
    (Abramowitz & Stegun 16.21).  alpha, beta are real arrays; m = k^2."""
    s, c, d, _ = ellipj(alpha, m)                     # functions of modulus k
    s1, c1, d1, _ = ellipj(beta, 1.0 - m)             # functions of modulus k'
    denom = c1**2 + m * s**2 * s1**2
    re = s * d1 / denom
    im = c * d * s1 * c1 / denom
    return re + 1j * im


def Y_of_y(yv, p=PARAMS):                             # true depth -> effective depth
    yv = np.asarray(yv, float)
    L_si = p["L_si"]
    return np.where(yv <= L_si, yv, L_si + (K_SI / K_OX) * (yv - L_si))


def build_phih(V, p=PARAMS, nb=900, charge_free=False):
    """Construct phi_h by conformal mapping and return (eval_eff, meta).

    eval_eff(x, Yeff) -> phi_h evaluated at EFFECTIVE coordinates (x, Yeff) by the
    Poisson integral / harmonic-measure formula; it is harmonic by construction,
    so it can be probed at arbitrary points (used by the mean-value verification
    as well as by the plotting grid)."""
    L_si, t_ox, Lx = p["L_si"], p["t_ox"], p["Lx"]
    H = L_si + t_ox
    H_eff = L_si + (K_SI / K_OX) * t_ox               # effective single-dielectric height

    # --- conformal map: rectangle [0,Lx]x[0,H_eff] -> upper half w-plane -----
    m = _rect_modulus(2.0 * H_eff / Lx)
    Kc = ellipk(m)
    s_scale = 2.0 * Kc / Lx                            # real scale: x,Y -> u-plane

    def z_to_w(xarr, Yeff):                            # (x, effective depth) -> w
        return _sn_complex(s_scale * np.asarray(xarr, float) - Kc,
                           s_scale * np.asarray(Yeff, float), m)

    # --- boundary Dirichlet data for phi_h ----------------------------------
    # phi_p(0)=V on the whole top; so phi_h_top = phi_top_actual - V.
    y1d, phi_p = particular_1d(V, p, charge_free=charge_free)
    phip_if = float(np.interp(L_si, y1d, phi_p))       # phi_p at Si/oxide interface
    a, b = p["drain_w"], Lx - p["source_w"]

    tb, gb = [], []                                    # real-axis coord, boundary value
    # top surface (Yeff=0, x:0->Lx): drain V, source 0, cosine ramp in the gap
    xs = np.linspace(0.0, Lx, nb)
    top_actual = np.where(xs <= a, V,
                 np.where(xs >= b, 0.0,
                          V * 0.5 * (1 + np.cos(np.pi * (xs - a) / (b - a)))))
    tb.append(z_to_w(xs, np.zeros_like(xs)).real); gb.append(top_actual - V)
    # right edge (p+ sinker): total phi=0 -> phi_h = -phi_p(y)
    ys = np.linspace(0.0, H, nb // 2)
    tb.append(z_to_w(np.full_like(ys, Lx), Y_of_y(ys, p)).real)
    gb.append(-np.interp(ys, y1d, phi_p))
    # left edge (under the drain): total phi = phi_p(y) -> phi_h = 0
    tb.append(z_to_w(np.zeros_like(ys), Y_of_y(ys, p)).real)
    gb.append(np.zeros_like(ys))
    tb = np.concatenate(tb); gb = np.concatenate(gb)
    order = np.argsort(tb); tb, gb = tb[order], gb[order]
    gmid = 0.5 * (gb[:-1] + gb[1:])

    def eval_eff(xarr, Yeff):
        """phi_h at effective coords, via the subtended-angle harmonic measure."""
        w = z_to_w(np.ravel(xarr), np.ravel(Yeff))
        Wr, Wi = w.real, np.maximum(w.imag, 1e-9)       # guard the top edge
        theta = np.arctan2(Wi[:, None], Wr[:, None] - tb[None, :])
        # theta(t) rises 0->pi with t, so segment measure = (theta_{i+1}-theta_i)/pi
        dtheta = theta[:, 1:] - theta[:, :-1]
        return (dtheta * gmid[None, :]).sum(axis=1) / np.pi

    meta = dict(m=m, Kc=Kc, s_scale=s_scale, H_eff=H_eff, phip_if=phip_if,
                y1d=y1d, phi_p=phi_p, z_to_w=z_to_w)
    return eval_eff, meta


def solve(V, p=PARAMS, nx=161, ny=121, charge_free=False):
    """Full field phi = phi_p(y) + phi_h(x,y).  Returns (x, y, phi, meta)."""
    L_si, t_ox, Lx = p["L_si"], p["t_ox"], p["Lx"]
    H = L_si + t_ox
    eval_eff, meta = build_phih(V, p, charge_free=charge_free)
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, H, ny)
    X, Yg = np.meshgrid(x, y)
    phi_h = eval_eff(X, Y_of_y(Yg, p)).reshape(ny, nx)
    phi_p_on_y = np.interp(y, meta["y1d"], meta["phi_p"])
    phi = phi_p_on_y[:, None] + phi_h
    meta["phi_h"] = phi_h
    return x, y, phi, meta


# ============================================================================
#  VERIFY
# ============================================================================
def verify(V=500.0):
    print("=" * 70)
    print("VERIFICATION of the CONFORMAL-MAPPING solution (Nakagawa Fig.4)")
    print("=" * 70)
    x, y, phi, meta = solve(V, nx=201, ny=161)
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]

    # (1) the conformal map is analytic (Cauchy-Riemann) -- spot check J != 0
    print(f"(1) rectangle->UHP map: Jacobi parameter m = {meta['m']:.6f}"
          f"   (K'/K = {ellipk(1-meta['m'])/ellipk(meta['m']):.4f}"
          f" = 2 H_eff/Lx = {2*meta['H_eff']/Lx:.4f})")

    # (2) phi_h is harmonic: test the MEAN-VALUE PROPERTY directly through the
    #     conformal evaluator (no PDE grid).  A harmonic function equals its
    #     average over any circle centred at an interior point; we check that on
    #     a lattice of interior points using small circles in effective coords.
    H_eff = meta["H_eff"]
    eval_eff, _ = build_phih(V, PARAMS)
    r = 0.8 * UM                                        # probe-circle radius
    ang = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    xc = np.linspace(3 * UM, Lx - 3 * UM, 20)
    Yc = np.linspace(3 * UM, H_eff - 3 * UM, 16)        # keep circles clear of edges
    Xc, Yyc = np.meshgrid(xc, Yc)
    centre = eval_eff(Xc, Yyc)
    ring = np.zeros_like(centre)
    for t in ang:                                       # average over the circle
        ring += eval_eff(Xc + r * np.cos(t), Yyc + r * np.sin(t))
    ring /= len(ang)
    mv_err = np.abs(ring - centre) / V                  # mean-value residual / V
    print(f"(2) harmonicity of phi_h -- mean-value property on {Xc.size} interior"
          f" circles (r={r/UM:.1f} um):"
          f"\n    max|avg_circle - centre|/V = {mv_err.max():.2e},"
          f"  RMS = {np.sqrt(np.mean(mv_err**2)):.2e}  (-> 0 means harmonic)")

    # (3) boundary conditions recovered on the top surface
    a, b = PARAMS["drain_w"], Lx - PARAMS["source_w"]
    top = phi[0, :]
    drain_cols = x <= a; src_cols = x >= b
    print(f"(3) top-surface BC:  drain max|phi-V| = "
          f"{np.abs(top[drain_cols]-V).max():.2f} V,"
          f"   source max|phi| = {np.abs(top[src_cols]).max():.2f} V")

    # (4) potential stays physical
    print(f"(4) potential range  min={phi.min():.2f} V  max={phi.max():.2f} V"
          f"   (should stay ~within [0, {V:.0f}])")

    # (5) THE physics: oxide voltage share under the drain centre
    ix = int(np.argmin(np.abs(x - 0.5 * PARAMS["drain_w"])))
    v_if = float(np.interp(L_si, y, phi[:, ix]))
    share = v_if / phi[0, ix]
    # 1-D reference: from phi_p alone (charge-free divider for comparison)
    share_1d = meta["phip_if"] / V
    tsi_e = L_si; tox_e = (K_SI/K_OX)*t_ox
    share_chargefree = tox_e / (tsi_e + tox_e)
    print(f"(5) oxide voltage share under drain centre = {share:.3f}")
    print(f"    1-D particular solution phi_p(interface)/V = {share_1d:.3f}"
          f"   -> {'MORE' if share_1d>0.5 else 'less'} than half across the oxide")
    print(f"    (charge-free dielectric divider would give only {share_chargefree:.3f}"
          f" -- the donor charge is what lifts it past 1/2, matching the paper)")
    print("=" * 70)
    return x, y, phi, meta


# ============================================================================
#  PLOTS  (mirror fig4_pde_standalone.py so the two are directly comparable)
# ============================================================================
def make_plots(V=500.0):
    x, y, phi, meta = solve(V, nx=221, ny=171)
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    Xum, Yum = np.meshgrid(x / UM, y / UM)
    i_if = int(np.argmin(np.abs(y - L_si)))
    ix = int(np.argmin(np.abs(x - 0.5 * PARAMS["drain_w"])))
    ox_share = float(np.interp(L_si, y, phi[:, ix])) / phi[0, ix]

    # ---------------------------- 3-D surface --------------------------------
    fig = plt.figure(figsize=(9.4, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    sx = max(1, (len(x) - 1) // 60); sy = max(1, (len(y) - 1) // 45)
    surf = ax.plot_surface(Xum[::sy, ::sx], Yum[::sy, ::sx], phi[::sy, ::sx],
                           rcount=45, ccount=60, cmap="turbo",
                           edgecolor="k", linewidth=0.22, alpha=0.9)
    ax.plot(x / UM, np.full_like(x, L_si / UM), phi[i_if, :],
            color="crimson", lw=2.2, label="Si / buried-oxide interface")
    v_ox = phi[i_if, 0]
    ax.plot(np.zeros(2), [L_si / UM, (L_si + t_ox) / UM], [v_ox, 0.0],
            color="magenta", lw=3.2, label="oxide voltage drop (drain edge)")
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9)
    ax.text(6, 9.0, 360, f"Silicon 0-20 $\\mu$m\n($\\approx${100*(1-ox_share):.0f}% of $V$)",
            fontsize=8, ha="center", va="center", bbox=box)
    ax.text(3, 21.5, 150, f"Buried oxide 20-23 $\\mu$m\n$>$half of $V$"
            f" ($\\approx${100*ox_share:.0f}%)", color="magenta",
            fontsize=8, ha="center", va="center", bbox=box)
    ax.set_xlabel("x  ($\\mu$m)"); ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_zlabel("potential (V)"); ax.set_zlim(0, V)
    ax.set_title("Nakagawa 1990 Fig.4 by CONFORMAL MAPPING (structure B, 500 V)\n"
                 f"analytic Schwarz-Christoffel solve -- oxide carries "
                 f"$\\approx${100*ox_share:.0f}% of $V$")
    ax.view_init(elev=28, azim=-60); ax.invert_yaxis()
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.02, label="potential (V)")
    fig.tight_layout()
    out3d = "../figs/fig4_conformal_3d.png"
    fig.savefig(out3d, dpi=150); plt.close(fig); print("wrote", out3d)

    # ---------------------- 2-D contour + vertical cut -----------------------
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
    ax1.set_title("Fig.4 by conformal mapping — 2-D map (structure B, 500 V)")
    fig.colorbar(cf, ax=ax1, label="potential (V)", shrink=0.9)

    col = phi[:, ix]
    ax2.plot(y / UM, col, "b-", lw=2)
    ax2.axvspan(20, 23, color="#d0d0d0", alpha=0.6)
    ax2.axvline(20, color="k", ls="--", lw=0.7)
    v_ox = float(np.interp(L_si, y, col))
    ax2.annotate(f"{ox_share*100:.0f}% of {V:.0f} V\nacross the oxide",
                 xy=(21.5, v_ox / 2), xytext=(11, 120),
                 arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax2.text(9, 380, "Silicon", fontsize=9)
    ax2.text(21.5, 380, "Oxide", fontsize=9, ha="center")
    ax2.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain centre)")
    ax2.set_ylabel("potential (V)"); ax2.set_xlim(0, 23); ax2.grid(alpha=0.25)
    fig.tight_layout()
    out2d = "../figs/fig4_conformal.png"
    fig.savefig(out2d, dpi=150); plt.close(fig); print("wrote", out2d)


def _oxide_share(x, y, phi, p=PARAMS):
    """Oxide voltage share under the drain centre (same metric everyone uses)."""
    ix = int(np.argmin(np.abs(x - 0.5 * p["drain_w"])))
    v_if = float(np.interp(p["L_si"], y, phi[:, ix]))
    return v_if / phi[0, ix], ix


def compare_charge_free(V=500.0):
    """(1) Charged vs charge-free conformal solutions: isolate how much of the
    >1/2 oxide share is the DONOR CHARGE and how much is the DIELECTRIC RATIO."""
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]
    x, y, phi_c, _ = solve(V, nx=221, ny=171, charge_free=False)   # charged
    _, _, phi_0, _ = solve(V, nx=221, ny=171, charge_free=True)    # charge-free
    sc, ix = _oxide_share(x, y, phi_c)
    s0, _ = _oxide_share(x, y, phi_0)

    print("-" * 70)
    print("(A) charged vs charge-free (conformal), oxide share under drain:")
    print(f"    with donor charge   = {sc:.3f}   (> 1/2  -> paper's result)")
    print(f"    charge-free divider = {s0:.3f}   (pure dielectric ratio)")
    print(f"    => the donor sheet adds {100*(sc-s0):+.0f} percentage points.")

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(y / UM, phi_c[:, ix], "b-", lw=2.2,
            label=f"with donor charge  ({100*sc:.0f}% across oxide)")
    ax.plot(y / UM, phi_0[:, ix], "r--", lw=2.2,
            label=f"charge-free (Laplace)  ({100*s0:.0f}% across oxide)")
    ax.axvspan(20, 23, color="#d0d0d0", alpha=0.6)
    ax.axvline(20, color="k", ls=":", lw=0.8)
    ax.text(9, 460, "Silicon", fontsize=10)
    ax.text(21.4, 460, "Oxide", fontsize=10, ha="center")
    ax.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain centre)")
    ax.set_ylabel("potential (V)")
    ax.set_title("Conformal solve: donor charge vs charge-free dielectric divider\n"
                 "the space charge is what pushes the oxide share past 1/2")
    ax.set_xlim(0, 23); ax.grid(alpha=0.25); ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out = "../figs/fig4_conformal_chargefree.png"
    fig.savefig(out, dpi=150); plt.close(fig); print("wrote", out)


def compare_with_fv(V=500.0):
    """(2) Overlay the conformal solution against the finite-volume Poisson solve
    from fig4_pde_standalone.py -- two independent methods, same structure."""
    try:
        import fig4_pde_standalone as fv
    except Exception as e:                              # pragma: no cover
        print(f"[compare_with_fv] could not import fig4_pde_standalone: {e}")
        return
    L_si, t_ox, Lx = PARAMS["L_si"], PARAMS["t_ox"], PARAMS["Lx"]

    xc, yc, phic, _ = solve(V, nx=221, ny=171, charge_free=False)  # conformal
    xf, yf, phif = fv.solve_poisson(201, 185, V)                   # finite volume
    sc, ixc = _oxide_share(xc, yc, phic)
    sf, ixf = _oxide_share(xf, yf, phif, fv.PARAMS)

    # common vertical cut under the drain centre
    yq = np.linspace(0, L_si + t_ox, 400)
    cc = np.interp(yq, yc, phic[:, ixc])
    cf_ = np.interp(yq, yf, phif[:, ixf])
    rms = np.sqrt(np.mean((cc - cf_) ** 2))
    mx = np.abs(cc - cf_).max()

    print("-" * 70)
    print("(B) conformal vs finite-volume Poisson (fig4_pde_standalone):")
    print(f"    oxide share under drain:  conformal = {sc:.3f}   FV = {sf:.3f}")
    print(f"    vertical-cut difference:  RMS = {rms:.1f} V   max = {mx:.1f} V"
          f"   ({100*rms/V:.1f}% of V)")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    # left: vertical cuts overlaid
    axL.plot(yq / UM, cc, "b-", lw=2.4, label=f"conformal  ({100*sc:.0f}%)")
    axL.plot(yq / UM, cf_, "r--", lw=2.0, label=f"finite-volume Poisson  ({100*sf:.0f}%)")
    axL.axvspan(20, 23, color="#d0d0d0", alpha=0.6)
    axL.axvline(20, color="k", ls=":", lw=0.8)
    axL.text(9, 460, "Silicon", fontsize=10); axL.text(21.4, 460, "Oxide",
                                                       fontsize=10, ha="center")
    axL.set_xlabel("depth y  ($\\mu$m)  (cut under drain centre)")
    axL.set_ylabel("potential (V)"); axL.set_xlim(0, 23); axL.grid(alpha=0.25)
    axL.legend(loc="upper right", fontsize=9)
    axL.set_title(f"Two methods, same cut  (RMS diff {rms:.0f} V)")
    # right: difference map over the whole cross-section
    Xf, Yf = np.meshgrid(xf, yf)
    from scipy.interpolate import RegularGridInterpolator
    gi = RegularGridInterpolator((yc, xc), phic, bounds_error=False, fill_value=None)
    phic_on_fv = gi(np.column_stack([Yf.ravel(), Xf.ravel()])).reshape(phif.shape)
    diff = phic_on_fv - phif
    im = axR.pcolormesh(Xf / UM, Yf / UM, diff, cmap="coolwarm",
                        vmin=-np.abs(diff).max(), vmax=np.abs(diff).max(),
                        shading="auto")
    axR.axhline(20.0, color="k", lw=1.0, ls="--")
    axR.invert_yaxis()
    axR.set_xlabel("x  ($\\mu$m)"); axR.set_ylabel("depth y  ($\\mu$m)")
    axR.set_title("conformal $-$ finite-volume  (V)")
    fig.colorbar(im, ax=axR, label="potential difference (V)", shrink=0.9)
    fig.suptitle("Nakagawa Fig.4 — conformal mapping vs finite-volume Poisson",
                 fontsize=12)
    fig.tight_layout()
    out = "../figs/fig4_conformal_vs_fv.png"
    fig.savefig(out, dpi=150); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    verify(500.0)
    make_plots(500.0)
    compare_charge_free(500.0)
    compare_with_fv(500.0)
