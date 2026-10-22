"""Fig. 6 reproduction -- step-thickness drift-region SOI device (2-D Poisson).

Companion to the analytical core: this is one of the paper's 2-D DEVICE figures
(Fig. 6), reproduced qualitatively with a 2-D Poisson solve (not MEDICI).

Geometry (y downward, y=0 at the flat top of the structure):
  * buried oxide fixed at depth [TS_MAX, TS_MAX+tI]     (tI = 1 um)
  * top silicon film, stepped thickness tS(x): 0.1 / 0.3 / 0.5 um in three
    lateral segments (thin at the source, thick at the drain)
  * field oxide fills the etched region above the thin steps (up to y=0)
  * grounded silicon substrate below the buried oxide
  * drain over the tS3=0.5um step (ES,C=67, EI=201 V/um analytically);
    source over tS1=0.1um; linear (RESURF) drift doping in between.

Produces:
  fig6b_field3d.png  -- 3-D |E|(x,y) surface  (the paper's Fig. 6(b))
  fig6a_potential.png -- equipotential contours (the paper's Fig. 6(a), bonus)

The key features to match: silicon field limited to ~ES,C, the buried-oxide field
enhanced to ~3x that (EI ~ eps_si/eps_ox * ES,C), and NEW FIELD PEAKS at each
step interface.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import scipy.sparse as sp
import scipy.sparse.linalg as spla

Q, EPS0, EPS_S, EPS_OX = 1.602176634e-19, 8.8541878128e-14, 11.7, 3.9
UM = 1e-4

# geometry (cm)
LD   = 30 * UM             # lateral drift length shown
TS   = [0.1, 0.3, 0.5]     # step silicon thicknesses (um)
TS_MAX = 0.5 * UM
TI   = 1.0 * UM            # buried oxide
Y_OX0, Y_OX1 = TS_MAX, TS_MAX + TI
Y_BOT = Y_OX1 + 2.0 * UM   # substrate bottom (grounded)

# doping: linear RESURF ramp, source -> drain (tuned so the vertical field at
# the Si/oxide interface under the drain reaches ~ES,C = 61 V/um as in the paper)
ND0  = 8.5e15             # cm^-3 at source
ND_SLOPE = 12.0            # drain doping = (1+SLOPE)*ND0
V_DRAIN = 206.0            # applied drain voltage (paper BV=206 V)

NX, NY = 241, 176


def ts_of_x(xc):
    """silicon thickness (cm) vs lateral position, 3 equal steps."""
    seg = min(2, int(xc / LD * 3))
    return TS[seg] * UM


def build():
    x = np.linspace(0, LD, NX)
    y = np.linspace(0, Y_BOT, NY)
    X, Y = np.meshgrid(x, y)                      # (ny,nx)
    eps = np.full_like(X, EPS_OX)                 # default: oxide (field + buried)
    Nnet = np.zeros_like(X)                       # net donor (cm^-3)
    mat = np.zeros_like(X, dtype=int)             # 0 oxide,1 si-film,2 substrate

    ysurf = np.array([TS_MAX - ts_of_x(xx) for xx in x])   # silicon top per column
    for j in range(NY):
        for i in range(NX):
            yy = y[j]
            if yy >= Y_OX1:                       # substrate
                eps[j, i] = EPS_S; mat[j, i] = 2
            elif yy >= Y_OX0:                     # buried oxide
                eps[j, i] = EPS_OX
            elif yy >= ysurf[i]:                  # silicon film
                eps[j, i] = EPS_S; mat[j, i] = 1
                nd = ND0 * (1.0 + ND_SLOPE * x[i] / LD)
                Nnet[j, i] = nd
            else:                                 # field oxide (etched region)
                eps[j, i] = EPS_OX
    return x, y, X, Y, eps * EPS0, Q * Nnet, mat, ysurf


def solve():
    x, y, X, Y, eps, rho, mat, ysurf = build()
    dx, dy = x[1] - x[0], y[1] - y[0]
    idx = lambda j, i: j * NX + i
    rows, cols, data, b = [], [], [], np.zeros(NX * NY)

    # Dirichlet masks: substrate bottom = 0; drain (right, in silicon) = V; source (left) = 0
    for j in range(NY):
        for i in range(NX):
            k = idx(j, i)
            dirich = None
            if j == NY - 1:
                dirich = 0.0                              # substrate grounded
            elif i == NX - 1 and mat[j, i] == 1:
                dirich = V_DRAIN                          # drain contact (thick-Si end)
            elif i == 0 and mat[j, i] == 1:
                dirich = 0.0                              # source contact
            if dirich is not None:
                rows.append(k); cols.append(k); data.append(1.0); b[k] = dirich
                continue
            diag = 0.0
            for jj, ii, h in ((j, i-1, dx), (j, i+1, dx), (j-1, i, dy), (j+1, i, dy)):
                if 0 <= jj < NY and 0 <= ii < NX:
                    ef = 2.0 / (1.0/eps[j, i] + 1.0/eps[jj, ii])
                    c = ef / (h*h)
                    rows.append(k); cols.append(idx(jj, ii)); data.append(c); diag -= c
            rows.append(k); cols.append(k); data.append(diag); b[k] = -rho[j, i]
    A = sp.csr_matrix((data, (rows, cols)), shape=(NX*NY, NX*NY))
    phi = spla.spsolve(A, b).reshape(NY, NX)

    Ey, Ex = np.gradient(-phi, y, x)
    Emag = np.hypot(Ex, Ey)
    return x, y, X, Y, phi, Emag, np.abs(Ey), mat, ysurf


def main():
    x, y, X, Y, phi, Emag, Eyabs, mat, ysurf = solve()
    Xum, Yum = X/UM, Y/UM

    # report the VERTICAL field Ey across the buried-oxide interface under the
    # drain -- the quantity the continuity relation EI = (eps_si/eps_ox) ES,C
    # governs. Sampled mid-step (away from the contact corner).
    icut = np.searchsorted(x, 0.80 * LD)          # inside the tS3 (drain) step
    j_if = np.searchsorted(y, Y_OX0)
    esc = Eyabs[j_if - 2, icut] / 1e4             # silicon side of interface
    ei = Eyabs[j_if + 2, icut] / 1e4              # oxide side of interface
    print(f"vertical field across Si/oxide interface under drain (x={x[icut]/UM:.0f}um):")
    print(f"  ES,C (silicon side) = {esc:6.1f} V/um   (paper ~61, analytical 67)")
    print(f"  EI   (oxide side)   = {ei:6.1f} V/um   (paper ~190, analytical 201)")
    print(f"  ratio EI/ES,C       = {ei/esc:6.2f}     (continuity eps_si/eps_ox = 3.0)")

    # ---- Fig 6(b): 3-D field surface ----
    # clip the immediate contact-corner singularities for display (paper shows a
    # smoothed device field); cap at 2.0e6 V/cm = 200 V/um like the paper axis.
    Edisp = np.clip(Emag, 0, 2.0e6) / 1e6
    fig = plt.figure(figsize=(8.2, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    sy, sx = max(1, NY//70), max(1, NX//90)
    ax.plot_surface(Xum[::sy, ::sx], Yum[::sy, ::sx], Edisp[::sy, ::sx],
                    cmap="viridis", edgecolor="k", linewidth=0.2, alpha=0.96,
                    rcount=70, ccount=90)
    ax.set_xlabel("x  ($\\mu$m)"); ax.set_ylabel("y  ($\\mu$m)")
    ax.set_zlabel("E  ($\\times10^{6}$ V/cm)")
    ax.set_zlim(0, 2.0)
    ax.set_title("Fig.6(b)  3-D electric-field distribution (step-thickness SOI)")
    ax.text2D(0.60, 0.82, "$E_I$ peak (buried oxide\nunder drain, ~2x10$^6$ V/cm)",
              transform=ax.transAxes, fontsize=8, color="#083")
    ax.text2D(0.10, 0.60, "new peaks at\nstep interfaces", transform=ax.transAxes,
              fontsize=8, color="#a00")
    ax.view_init(elev=24, azim=-58)
    fig.tight_layout(); fig.savefig("../figs/fig6b_field3d.png", dpi=150)
    plt.close(fig)
    print("wrote fig6b_field3d.png")

    # ---- Fig 6(a): equipotential contours (bonus) ----
    fig, ax = plt.subplots(figsize=(8.0, 3.6))
    lv = np.linspace(0, V_DRAIN, 22)
    ax.contour(Xum, Yum, phi, levels=lv, colors="k", linewidths=0.5)
    ax.plot(x/UM, ysurf/UM, color="C0", lw=1.2)
    ax.axhline(Y_OX0/UM, color="gray", lw=1); ax.axhline(Y_OX1/UM, color="gray", lw=1)
    ax.text(1, 0.75, "$t_{S1}$"); ax.text(11, 0.55, "$t_{S2}$"); ax.text(21, 0.35, "$t_{S3}$")
    ax.text(13, 1.0, "buried oxide", color="gray", fontsize=8)
    ax.text(13, 2.2, "substrate", color="gray", fontsize=8)
    ax.set_xlabel("x  ($\\mu$m)"); ax.set_ylabel("depth y  ($\\mu$m)")
    ax.set_title("Fig.6(a)  Equipotential contours (step-thickness SOI)")
    ax.invert_yaxis(); ax.set_ylim(Y_BOT/UM, 0)
    fig.tight_layout(); fig.savefig("../figs/fig6a_potential.png", dpi=150)
    plt.close(fig)
    print("wrote fig6a_potential.png")


if __name__ == "__main__":
    main()
