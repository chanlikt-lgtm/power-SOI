"""2-D Poisson solver for the SOI cross-section (Nakagawa ISPSD 1990, Fig.4).

Cross-section (x lateral, y depth downward):
    silicon drift  0 <= y < L_si   (net donor charge: background + shallow layer)
    buried oxide   L_si <= y < H    (charge free)
    substrate      y = H            (grounded)

Structure B (Fig.3B): drain (V) at the top-left, directly over the buried oxide;
source/body ground at the top-right, tied to the substrate by a p+ sinker on the
right edge.  This is the paper's optimized device and the one Fig.4 depicts.

5-point finite-volume Laplacian with face permittivities; sparse direct solve.
Units: cm, V, cm^-3.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from model import Q, EPS0, EPS_SI, EPS_OX, UM

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


class Device2D:
    def __init__(self, L_si=20 * UM, t_ox=3 * UM, Lx=50 * UM,
                 w_n=2 * UM, N_D0=1.0e14, S_D=1.5e12,
                 nx=201, ny=185, drain_w=8 * UM, source_w=8 * UM,
                 structure="B"):
        self.L_si, self.t_ox, self.Lx = L_si, t_ox, Lx
        self.H = L_si + t_ox
        self.w_n, self.N_D0, self.S_D = w_n, N_D0, S_D
        self.nx, self.ny = nx, ny
        self.drain_w, self.source_w = drain_w, source_w
        self.structure = structure
        self.x = np.linspace(0, Lx, nx)
        self.y = np.linspace(0, self.H, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)   # shape (ny, nx)

    # ---- material / charge maps -------------------------------------------
    def eps_map(self):
        e = np.where(self.Y < self.L_si - 1e-12, EPS_SI, EPS_OX) * EPS0
        return e

    def charge_map(self):
        """rho = q * N_net (C/cm^3), positive donors."""
        N = np.zeros_like(self.Y)
        in_si = self.Y < self.L_si - 1e-12
        N[in_si] = self.N_D0
        in_shallow = in_si & (self.Y >= self.L_si - self.w_n)
        N[in_shallow] += self.S_D / self.w_n
        return Q * N

    # ---- boundary mask -----------------------------------------------------
    def dirichlet(self, V):
        """Return (mask, values) arrays for Dirichlet nodes."""
        mask = np.zeros_like(self.Y, dtype=bool)
        val = np.zeros_like(self.Y)

        # substrate bottom grounded
        mask[-1, :] = True
        val[-1, :] = 0.0

        xL = self.x <= self.drain_w
        xR = self.x >= self.Lx - self.source_w

        if self.structure == "B":
            # drain (V) top-left, source ground top-right
            mask[0, xL] = True; val[0, xL] = V
            mask[0, xR] = True; val[0, xR] = 0.0
            # p+ sinker: right edge silicon grounded
            in_si_row = self.y < self.L_si - 1e-12
            mask[in_si_row, -1] = True; val[in_si_row, -1] = 0.0
        else:  # structure A — mirror: drain top-right, ground top-left + left sinker
            mask[0, xR] = True; val[0, xR] = V
            mask[0, xL] = True; val[0, xL] = 0.0
            in_si_row = self.y < self.L_si - 1e-12
            mask[in_si_row, 0] = True; val[in_si_row, 0] = 0.0
        return mask, val

    # ---- assemble & solve --------------------------------------------------
    def solve(self, V):
        ny, nx = self.ny, self.nx
        eps = self.eps_map()
        rho = self.charge_map()
        mask, val = self.dirichlet(V)

        dx = self.x[1] - self.x[0]
        dy = self.y[1] - self.y[0]

        def idx(j, i):
            return j * nx + i

        N = nx * ny
        rows, cols, data = [], [], []
        b = np.zeros(N)

        for j in range(ny):
            for i in range(nx):
                k = idx(j, i)
                if mask[j, i]:
                    rows.append(k); cols.append(k); data.append(1.0)
                    b[k] = val[j, i]
                    continue
                diag = 0.0
                # neighbours with harmonic-mean face permittivity
                for (jj, ii, h) in ((j, i - 1, dx), (j, i + 1, dx),
                                    (j - 1, i, dy), (j + 1, i, dy)):
                    if 0 <= jj < ny and 0 <= ii < nx:
                        ef = 2.0 / (1.0 / eps[j, i] + 1.0 / eps[jj, ii])
                        c = ef / (h * h)
                        rows.append(k); cols.append(idx(jj, ii)); data.append(c)
                        diag -= c
                    else:
                        # Neumann (zero-flux) surface: reflect -> no term
                        pass
                rows.append(k); cols.append(k); data.append(diag)
                b[k] = -rho[j, i]

        A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
        phi = spla.spsolve(A, b).reshape(ny, nx)
        return phi

    def fields(self, phi):
        Ey, Ex = np.gradient(-phi, self.y, self.x)
        Emag = np.hypot(Ex, Ey)
        return Ex, Ey, Emag


if __name__ == "__main__":
    d = Device2D(structure="B", S_D=1.5e12, t_ox=3 * UM)
    phi = d.solve(500.0)
    # oxide voltage fraction under the drain (x=0 column)
    col = phi[:, 0]
    i_if = np.searchsorted(d.y, d.L_si)
    print("structure B: phi at drain-bottom (Si/ox) =", round(col[i_if], 1),
          "V  -> oxide share", round(col[i_if] / 500.0, 2))
    Ex, Ey, Em = d.fields(phi)
    print("peak |E| in silicon:", f"{Em[d.y < d.L_si-1e-12].max():.2e} V/cm")
