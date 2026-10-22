"""
================================================================================
 Nakagawa, Yasuhara & Baba, ISPSD 1990 -- thin-SOI 500 V lateral device
 FULL self-contained reproduction: from the PDE to the 2-D/3-D plots.
================================================================================

Paper: A. Nakagawa, N. Yasuhara, Y. Baba, "New 500 V Output Device Structures
for Thin Silicon Layer on Silicon Dioxide Film," Proc. ISPSD 1990, Tokyo,
pp. 97-101.

--------------------------------------------------------------------------------
THE PDE
--------------------------------------------------------------------------------
Electrostatics under full depletion (no mobile carriers in the blocking state):

        div( eps(r) grad phi(r) ) = -rho(r),      rho = q * N_net    [Poisson]
        E = -grad phi

Materials / charge:
  * silicon drift layer (0 <= y < L_si): net donor charge q*N_D0, plus a shallow
    N donor layer of areal dose S_D near the bottom (the "shield").
  * buried oxide (L_si <= y < L_si+t_ox): charge free.
  * grounded silicon substrate at the bottom.

Interface condition (built into the finite-volume discretization via face
permittivities): the normal displacement D = eps*E is continuous, so the field
jumps up into the oxide by eps_si/eps_ox ~ 3.

We solve this PDE in two forms:
  (1) 1-D vertical cut along the drain symmetry axis AA'  -> Figs. 5, 6, 7
  (2) 2-D cross-section                                    -> Fig. 4 (2-D + 3-D)

Avalanche breakdown (Figs. 6, 7): breakdown when the ionization integral in the
silicon reaches unity,   I(V) = integral alpha(|E|) dy = 1,  with the Chynoweth
coefficient  alpha(E) = A exp(-B/E).

Run:  python nakagawa1990_full.py            # writes all figures to ../figs
Units: cm, V, cm^-3, cm^-2.
================================================================================
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D            # noqa: F401 (registers 3d)
from scipy.linalg import solve_banded
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# numpy>=2 renamed trapz -> trapezoid
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ------------------------------------------------------------------ constants --
Q      = 1.602176634e-19      # elementary charge, C
EPS0   = 8.8541878128e-14     # vacuum permittivity, F/cm
EPS_SI = 11.7                 # relative permittivity, silicon
EPS_OX = 3.9                  # relative permittivity, SiO2
UM     = 1.0e-4               # 1 micron in cm

# Chynoweth ionization coefficient for silicon (high-field branch)
A_ION = 7.03e5               # cm^-1
B_ION = 1.231e6              # V/cm

OUTDIR = "../figs"


# ============================================================================ #
#  PART 1 -- 1-D VERTICAL POISSON  (drain symmetry axis AA', structure B)       #
# ============================================================================ #
class Geom:
    """Vertical device geometry + doping along the drain axis."""
    def __init__(self, L_si=20 * UM, t_ox=3 * UM, w_n=2 * UM,
                 N_D0=1.0e14, S_D=1.5e12):
        self.L_si = L_si          # silicon drift thickness
        self.t_ox = t_ox          # buried-oxide thickness
        self.w_n  = w_n           # shallow N-layer thickness
        self.N_D0 = N_D0          # background n- doping (cm^-3)
        self.S_D  = S_D           # shallow-layer areal dose (cm^-2)

    @property
    def L_tot(self):
        return self.L_si + self.t_ox


def build_grid_1d(g: Geom, n_si=2000, n_ox=400):
    """Node coordinates y (cm), permittivity eps(y), net donor N(y) (cm^-3)."""
    y = np.unique(np.concatenate([
        np.linspace(0.0, g.L_si, n_si, endpoint=False),
        np.linspace(g.L_si, g.L_tot, n_ox + 1),
    ]))
    eps = np.where(y < g.L_si - 1e-12, EPS_SI, EPS_OX) * EPS0
    N = np.zeros_like(y)
    in_si = y < g.L_si - 1e-12
    N[in_si] = g.N_D0
    N[in_si & (y >= g.L_si - g.w_n)] += g.S_D / g.w_n     # shallow layer
    return y, eps, N


def solve_1d(g: Geom, V_applied, n_si=2000, n_ox=400):
    """Solve the 1-D Poisson PDE with Dirichlet BCs phi(0)=V, phi(L_tot)=0.

    Finite-volume discretization with harmonic-mean face permittivities;
    tridiagonal system solved in banded form.  Returns (y, phi, E), E=-dphi/dy.
    """
    y, eps, N = build_grid_1d(g, n_si, n_ox)
    n = y.size
    dy = np.diff(y)
    eps_face = 2.0 / (1.0 / eps[:-1] + 1.0 / eps[1:])     # harmonic mean
    cond = eps_face / dy

    ab = np.zeros((3, n))     # banded matrix: rows = [upper, diag, lower]
    b = np.zeros(n)
    rho = Q * N

    cl, cr = cond[:-1], cond[1:]
    h = 0.5 * (y[2:] - y[:-2])          # control-volume widths, interior nodes
    idx = np.arange(1, n - 1)
    ab[0, idx + 1] = cr                 # upper diagonal
    ab[1, idx]     = -(cl + cr)         # main diagonal
    ab[2, idx - 1] = cl                 # lower diagonal
    b[idx] = -rho[idx] * h

    ab[1, 0]  = 1.0; b[0]  = V_applied  # Dirichlet: drain
    ab[1, -1] = 1.0; b[-1] = 0.0        # Dirichlet: substrate

    phi = solve_banded((1, 1), ab, b)
    E = -np.gradient(phi, y)
    return y, phi, E


def oxide_fraction(g: Geom, V_applied):
    """Fraction of the applied voltage dropped across the buried oxide."""
    y, phi, _ = solve_1d(g, V_applied)
    i_if = np.searchsorted(y, g.L_si)
    return (phi[i_if] - phi[-1]) / V_applied


# ============================================================================ #
#  PART 2 -- AVALANCHE IONIZATION + BREAKDOWN VOLTAGE                            #
# ============================================================================ #
def alpha_ion(E):
    """Chynoweth ionization coefficient alpha(E) for silicon (cm^-1)."""
    E = np.asarray(E, float)
    a = np.zeros_like(E)
    m = E > 1.0e4
    a[m] = A_ION * np.exp(-B_ION / E[m])
    return a


def ionization_integral(g: Geom, V_applied):
    """I = integral alpha(|E|) dy through the silicon along AA'."""
    y, _, E = solve_1d(g, V_applied)
    in_si = y < g.L_si - 1e-12
    return _trapz(alpha_ion(np.abs(E[in_si])), y[in_si])


def breakdown_voltage(g: Geom, v_lo=20.0, v_hi=2000.0, tol=0.5):
    """Bisection on V so that the ionization integral equals 1."""
    f = lambda V: ionization_integral(g, V) - 1.0
    if f(v_lo) > 0:
        return v_lo
    if f(v_hi) < 0:
        return v_hi
    for _ in range(80):
        vm = 0.5 * (v_lo + v_hi)
        fm = f(vm)
        if abs(fm) < 1e-3 or (v_hi - v_lo) < tol:
            return vm
        if fm > 0:
            v_hi = vm
        else:
            v_lo = vm
    return 0.5 * (v_lo + v_hi)


# ============================================================================ #
#  PART 3 -- 2-D POISSON  (cross-section, for the Fig. 4 potential map)          #
# ============================================================================ #
class Device2D:
    """2-D cross-section (x lateral, y depth). Structure B: drain top-left over
    the oxide; source ground top-right, tied to the substrate by a right-edge
    p+ sinker; substrate grounded at the bottom."""
    def __init__(self, L_si=20 * UM, t_ox=3 * UM, Lx=50 * UM, w_n=2 * UM,
                 N_D0=1.0e14, S_D=1.5e12, nx=201, ny=185,
                 drain_w=8 * UM, source_w=8 * UM, structure="B"):
        self.L_si, self.t_ox, self.Lx = L_si, t_ox, Lx
        self.H = L_si + t_ox
        self.w_n, self.N_D0, self.S_D = w_n, N_D0, S_D
        self.nx, self.ny = nx, ny
        self.drain_w, self.source_w = drain_w, source_w
        self.structure = structure
        self.x = np.linspace(0, Lx, nx)
        self.y = np.linspace(0, self.H, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)      # (ny, nx)

    def eps_map(self):
        return np.where(self.Y < self.L_si - 1e-12, EPS_SI, EPS_OX) * EPS0

    def charge_map(self):
        N = np.zeros_like(self.Y)
        in_si = self.Y < self.L_si - 1e-12
        N[in_si] = self.N_D0
        N[in_si & (self.Y >= self.L_si - self.w_n)] += self.S_D / self.w_n
        return Q * N

    def dirichlet(self, V):
        mask = np.zeros_like(self.Y, bool)
        val  = np.zeros_like(self.Y)
        mask[-1, :] = True                                # substrate grounded
        xL = self.x <= self.drain_w
        xR = self.x >= self.Lx - self.source_w
        in_si_row = self.y < self.L_si - 1e-12
        if self.structure == "B":
            mask[0, xL] = True; val[0, xL] = V            # drain top-left
            mask[0, xR] = True                            # source top-right (0)
            mask[in_si_row, -1] = True                    # right p+ sinker (0)
        else:                                             # structure A (mirror)
            mask[0, xR] = True; val[0, xR] = V
            mask[0, xL] = True
            mask[in_si_row, 0] = True
        return mask, val

    def solve(self, V):
        """Assemble the 5-point finite-volume Laplacian (harmonic-mean face eps)
        and solve the sparse linear system.  Returns phi (ny, nx)."""
        ny, nx = self.ny, self.nx
        eps, rho = self.eps_map(), self.charge_map()
        mask, val = self.dirichlet(V)
        dx = self.x[1] - self.x[0]
        dy = self.y[1] - self.y[0]
        idx = lambda j, i: j * nx + i
        rows, cols, data = [], [], []
        b = np.zeros(nx * ny)
        for j in range(ny):
            for i in range(nx):
                k = idx(j, i)
                if mask[j, i]:
                    rows.append(k); cols.append(k); data.append(1.0)
                    b[k] = val[j, i]
                    continue
                diag = 0.0
                for jj, ii, hh in ((j, i - 1, dx), (j, i + 1, dx),
                                   (j - 1, i, dy), (j + 1, i, dy)):
                    if 0 <= jj < ny and 0 <= ii < nx:
                        ef = 2.0 / (1.0 / eps[j, i] + 1.0 / eps[jj, ii])
                        c = ef / (hh * hh)
                        rows.append(k); cols.append(idx(jj, ii)); data.append(c)
                        diag -= c
                    # else: Neumann (zero-flux) free surface -> no term
                rows.append(k); cols.append(k); data.append(diag)
                b[k] = -rho[j, i]
        A = sp.csr_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny))
        return spla.spsolve(A, b).reshape(ny, nx)


# ============================================================================ #
#  PART 4 -- PLOTS  (Fig. 4 2-D/3-D, Fig. 5, Fig. 6, Fig. 7)                     #
# ============================================================================ #
def plot_fig4(outdir=OUTDIR):
    """Fig. 4 -- potential distribution (structure B, 500 V): 3-D surface + 2-D."""
    d = Device2D(structure="B", S_D=1.5e12, t_ox=3 * UM, Lx=50 * UM)
    phi = d.solve(500.0)
    Xum, Yum = d.X / UM, d.Y / UM
    i_if = np.searchsorted(d.y, d.L_si)

    # ---- 3-D surface (matches the paper's wireframe view) ----
    fig = plt.figure(figsize=(8.0, 6.2))
    ax = fig.add_subplot(111, projection="3d")
    sy, sx = max(1, d.ny // 45), max(1, d.nx // 60)
    ax.plot_surface(Xum[::sy, ::sx], Yum[::sy, ::sx], phi[::sy, ::sx],
                    rcount=45, ccount=60, cmap="turbo",
                    edgecolor="k", linewidth=0.25, alpha=0.95)
    ax.text(3, 2, 505, "n layer (500 V)", fontsize=9)
    ax.text(2, 21.5, 40, "oxide film", fontsize=9)
    ax.text(46, 2, 30, "p$^+$ layer (0 V)", fontsize=9)
    ax.set_xlabel("x  ($\\mu$m)"); ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_zlabel("potential (V)")
    ax.set_title("Fig.4  Potential distribution - 3-D view (structure B, 500 V)")
    ax.view_init(elev=28, azim=-60); ax.invert_yaxis()
    fig.tight_layout(); fig.savefig(f"{outdir}/fig4_potential_3d.png", dpi=150)
    plt.close(fig)

    # ---- 2-D contour + vertical cut ----
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                 gridspec_kw={"height_ratios": [2.3, 1]})
    lv = np.arange(0, 501, 25)
    cf = a1.contourf(Xum, Yum, phi, levels=lv, cmap="turbo")
    a1.contour(Xum, Yum, phi, levels=lv, colors="k", linewidths=0.4)
    a1.axhline(20, color="white", ls="--", lw=1.5)
    a1.text(25, 21.6, "buried oxide", color="white", ha="center", fontsize=9)
    a1.text(4, 2, "drain 500 V", color="white", fontsize=9, weight="bold")
    a1.text(40, 2, "source 0 V", color="white", fontsize=9, weight="bold")
    a1.invert_yaxis(); a1.set_xlabel("x  ($\\mu$m)"); a1.set_ylabel("depth y  ($\\mu$m)")
    a1.set_title("Fig.4  Potential distribution - 2-D map (structure B, 500 V)")
    fig.colorbar(cf, ax=a1, label="potential (V)", shrink=0.9)

    col = phi[:, 0]
    a2.plot(d.y / UM, col, "b-", lw=2)
    a2.axvspan(20, 23, color="#d0d0d0", alpha=0.6); a2.axvline(20, color="k", ls="--", lw=0.7)
    v_ox = col[i_if]
    a2.annotate(f"{v_ox/500*100:.0f}% of 500 V\nacross the oxide",
                xy=(21.5, v_ox / 2), xytext=(11, 120),
                arrowprops=dict(arrowstyle="->"), fontsize=9)
    a2.text(9, 380, "Silicon", fontsize=9); a2.text(21.5, 380, "Oxide", ha="center", fontsize=9)
    a2.set_xlabel("depth y  ($\\mu$m)   (vertical cut under the drain, x=0)")
    a2.set_ylabel("potential (V)"); a2.set_xlim(0, 23); a2.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig4_potential.png", dpi=150)
    plt.close(fig)
    print(f"Fig.4  oxide voltage share under drain = {v_ox/500:.2f}")


def plot_fig5(outdir=OUTDIR):
    """Fig. 5 -- field vs depth along AA' for four shallow-layer doses."""
    curves = [("A", 2.5e12, "#c0392b"), ("B", 2.0e12, "#e67e22"),
              ("C", 1.5e12, "#27ae60"), ("D", 0.0, "#2c3e50")]
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for lab, SD, c in curves:
        y, _, E = solve_1d(Geom(t_ox=3 * UM, S_D=SD), 500.0)
        ax.plot(y / UM, np.abs(E), color=c, lw=1.8,
                label=f"{lab}: $S_D$={SD/1e12:.1f}$\\times10^{{12}}$")
    ax.axvspan(20, 23, color="#d0d0d0", alpha=0.5, zorder=0)
    ax.axvline(20, color="k", lw=0.7, ls="--")
    ax.text(10, 0.05e6, "Silicon", ha="center", fontsize=9)
    ax.text(21.5, 0.05e6, "Oxide", ha="center", fontsize=9)
    ax.set_xlabel("Depth  y  ($\\mu$m)"); ax.set_ylabel("Electric field  |E|  (V/cm)")
    ax.set_title("Fig.5  Field along symmetry axis AA'  (structure B, $t_{ox}$=3 $\\mu$m, 500 V)")
    ax.set_xlim(0, 23); ax.set_ylim(0, None)
    ax.legend(loc="upper left", fontsize=9); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig5_field.png", dpi=150)
    plt.close(fig)
    print("Fig.5  written")


def plot_fig6(outdir=OUTDIR):
    """Fig. 6 -- breakdown voltage vs dose, structure B, tox = 1,2,3 um."""
    paper = {
        3: [(0.25, 535), (0.5, 560), (0.75, 578), (1.0, 575), (1.25, 568),
            (1.5, 558), (2.0, 500), (2.5, 430)],
        2: [(0.25, 448), (0.5, 470), (0.75, 480), (1.0, 478), (1.25, 470),
            (1.5, 455), (2.0, 405), (2.5, 350)],
        1: [(0.25, 355), (0.5, 375), (0.75, 382), (1.0, 378), (1.25, 368),
            (1.5, 352), (2.0, 305), (2.5, 232)],
    }
    colors = {3: "#c0392b", 2: "#2980b9", 1: "#27ae60"}
    doses = np.linspace(0.05e12, 2.6e12, 26)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for tox in (3, 2, 1):
        bv = [breakdown_voltage(Geom(t_ox=tox * UM, S_D=SD)) for SD in doses]
        ax.plot(doses / 1e12, bv, color=colors[tox], lw=2.0,
                label=f"model  $t_{{ox}}$={tox} $\\mu$m")
        ax.plot([p[0] for p in paper[tox]], [p[1] for p in paper[tox]], "o",
                color=colors[tox], ms=5, mfc="white", mew=1.4,
                label=f"paper  $t_{{ox}}$={tox} $\\mu$m")
        i = int(np.argmax(bv))
        print(f"Fig.6  tox={tox}um peak {bv[i]:.0f} V @ dose {doses[i]/1e12:.2f}e12")
    ax.axhline(500, color="gray", ls=":", lw=1.0); ax.text(2.35, 508, "500 V", color="gray", fontsize=8)
    ax.set_xlabel("Impurity dose $S_D$  ($\\times10^{12}$ cm$^{-2}$)")
    ax.set_ylabel("Breakdown voltage  (V)")
    ax.set_title("Fig.6  Breakdown voltage vs. dose (structure B)")
    ax.set_xlim(0, 2.6); ax.set_ylim(0, 650)
    ax.legend(loc="lower left", fontsize=8, ncol=3); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig6_bv.png", dpi=150)
    plt.close(fig)


def plot_fig7(outdir=OUTDIR):
    """Fig. 7 -- breakdown voltage vs dose, structure A, tox = 5 um.

    The vertical model represents structure B's favourable (drain-over-oxide)
    sharing; A and B share the same vertical stack, so the model lies above the
    structure-A data -- their difference is a 2-D effect (A needs 5um vs B's 3um).
    """
    paper_A = [(0.2, 410), (0.5, 458), (0.8, 493), (1.0, 513)]
    doses = np.linspace(0.05e12, 1.6e12, 20)
    bv = [breakdown_voltage(Geom(t_ox=5 * UM, S_D=SD)) for SD in doses]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(doses / 1e12, bv, "-", color="#8e44ad", lw=2.0,
            label="vertical model, $t_{ox}$=5 $\\mu$m (B-type sharing)")
    ax.plot([p[0] for p in paper_A], [p[1] for p in paper_A], "s",
            color="#c0392b", ms=6, mfc="white", mew=1.5,
            label="paper Fig.7 (structure A)")
    ax.axhline(500, color="gray", ls=":", lw=1.0); ax.text(1.35, 508, "500 V", color="gray", fontsize=8)
    ax.annotate("A needs $t_{ox}$~5 $\\mu$m to reach 500 V\n(B reaches it at 3 $\\mu$m - Fig.6)",
                xy=(1.0, 513), xytext=(0.30, 600), arrowprops=dict(arrowstyle="->"), fontsize=8.5)
    ax.set_xlabel("Impurity dose $S_D$  ($\\times10^{12}$ cm$^{-2}$)")
    ax.set_ylabel("Breakdown voltage  (V)")
    ax.set_title("Fig.7  Breakdown voltage vs. dose (structure A, $t_{ox}$=5 $\\mu$m)")
    ax.set_xlim(0, 1.6); ax.set_ylim(300, 780)
    ax.legend(loc="lower right", fontsize=8.5); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig7_bv.png", dpi=150)
    plt.close(fig)
    print("Fig.7  written")


# ============================================================================ #
#  MAIN                                                                          #
# ============================================================================ #
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("=" * 68)
    print("Nakagawa ISPSD 1990 -- full reproduction (PDE -> plots)")
    print("=" * 68)

    # quick sanity numbers straight from the PDE
    g = Geom(t_ox=3 * UM, S_D=1.5e12)
    print(f"oxide voltage share (tox=3um, 500 V) = {oxide_fraction(g, 500.0):.3f}")
    print(f"BV (tox=3um, S_D=1.5e12)             = {breakdown_voltage(g):.1f} V")
    print("-" * 68)

    plot_fig4()      # 2-D + 3-D potential
    plot_fig5()      # field vs depth
    plot_fig6()      # BV vs dose, structure B
    plot_fig7()      # BV vs dose, structure A
    print("-" * 68)
    print(f"All figures written to {os.path.abspath(OUTDIR)}")


if __name__ == "__main__":
    main()
