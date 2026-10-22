"""1-D vertical electrostatics + avalanche model for Nakagawa et al. ISPSD 1990.

Paper: A. Nakagawa, N. Yasuhara, Y. Baba, "New 500V Output Device Structures for
Thin Silicon Layer on Silicon Dioxide Film," Proc. ISPSD 1990, Tokyo, pp. 97-101.

Physical picture (symmetry axis AA' of structure B, Fig.3B)
----------------------------------------------------------
The drain (top, n+) sits at the applied voltage V; the silicon substrate below
the buried oxide (BOX) is grounded.  Along the vertical cut AA' the full voltage
therefore drops across  [ depleted Si drift layer  +  buried oxide ].

Vertical stack (y measured downward from the top of the Si):
    0 .......... L_si         n- silicon drift layer  (background donors N_D0)
                              + shallow N donor layer of dose S_D near the bottom
    L_si ....... L_si+t_ox    buried oxide (charge free)
    substrate (y = L_si+t_ox) grounded (phi = 0)

Governing equation (electrostatics, full depletion assumed):
    d/dy( eps(y) dphi/dy ) = -rho(y),      rho = q * N_net(y)   (donors > 0)
    E(y) = -dphi/dy
At the Si/oxide interface the normal displacement D = eps*E is continuous, so
    E_ox = (eps_si/eps_ox) * E_si   (~3x jump upward into the oxide).

The shallow donor sheet supplies areal charge q*S_D.  Crossing it upward from the
oxide, the displacement drops by q*S_D, so the field in the *bulk* silicon is
    E_si,bulk = (D_ox - q*S_D)/eps_si,
which decreases as S_D increases -> the donor layer "shields" the bulk silicon
while a high field is sustained in the oxide (the paper's key mechanism).

Breakdown: avalanche is evaluated in the silicon via the ionization integral
    I = \int alpha(E) dy  along AA';   breakdown when I = 1.

Everything is in CGS-practical units: lengths in cm, fields in V/cm, dose in cm^-2.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

# numpy>=2 renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ---- physical constants -------------------------------------------------------
Q = 1.602176634e-19            # C
EPS0 = 8.8541878128e-14        # F/cm
EPS_SI = 11.7
EPS_OX = 3.9
UM = 1.0e-4                    # 1 micron in cm


# ---- default geometry / doping (structure B, Fig.5 conditions) ----------------
class Geom:
    def __init__(
        self,
        L_si: float = 20.0 * UM,     # silicon drift thickness
        t_ox: float = 3.0 * UM,      # buried-oxide thickness
        w_n: float = 2.0 * UM,       # shallow N diffusion-layer thickness
        N_D0: float = 1.0e14,        # background n- doping (cm^-3)
        S_D: float = 1.5e12,         # shallow-layer areal dose (cm^-2)
    ):
        self.L_si = L_si
        self.t_ox = t_ox
        self.w_n = w_n
        self.N_D0 = N_D0
        self.S_D = S_D

    @property
    def L_tot(self) -> float:
        return self.L_si + self.t_ox


def build_grid(g: Geom, n_si: int = 2000, n_ox: int = 400):
    """Uniform-ish grid; returns cell-centered y (cm), permittivity, net doping."""
    y_si = np.linspace(0.0, g.L_si, n_si, endpoint=False)
    y_ox = np.linspace(g.L_si, g.L_tot, n_ox + 1)
    y = np.concatenate([y_si, y_ox])
    y = np.unique(y)

    eps = np.where(y < g.L_si - 1e-12, EPS_SI, EPS_OX) * EPS0

    # net donor concentration (cm^-3): background everywhere in Si, plus the
    # shallow N layer distributed over w_n at the bottom of the silicon.
    N = np.zeros_like(y)
    in_si = y < g.L_si - 1e-12
    N[in_si] = g.N_D0
    n_conc = g.S_D / g.w_n                       # cm^-3 in the shallow layer
    in_shallow = (y >= g.L_si - g.w_n) & in_si
    N[in_shallow] += n_conc
    return y, eps, N


def solve_potential(g: Geom, V_applied: float, n_si: int = 2000, n_ox: int = 400):
    """Solve 1-D Poisson with Dirichlet BCs phi(0)=V_applied, phi(L_tot)=0.

    Returns (y, phi, E) with E = -dphi/dy (>0 pointing downward).
    Finite-volume discretization with face permittivities (harmonic mean).
    """
    y, eps, N = build_grid(g, n_si, n_ox)
    n = y.size

    # face midpoints and conductances eps/dy
    dy = np.diff(y)
    eps_face = 2.0 / (1.0 / eps[:-1] + 1.0 / eps[1:])   # harmonic mean
    cond = eps_face / dy                                  # F/cm per (V) -> flux coeff

    # assemble tridiagonal system in banded form (rows: upper, diag, lower)
    ab = np.zeros((3, n))
    b = np.zeros(n)
    rho = Q * N   # C/cm^3

    cl = cond[:-1]                       # cond[i-1] for i=1..n-2
    cr = cond[1:]                        # cond[i]   for i=1..n-2
    h = 0.5 * (y[2:] - y[:-2])           # control-volume width, i=1..n-2
    idx = np.arange(1, n - 1)
    ab[0, idx + 1] = cr                  # upper diagonal
    ab[1, idx] = -(cl + cr)              # main diagonal
    ab[2, idx - 1] = cl                  # lower diagonal
    b[idx] = -rho[idx] * h

    # Dirichlet BCs
    ab[1, 0] = 1.0
    b[0] = V_applied
    ab[1, -1] = 1.0
    b[-1] = 0.0

    phi = solve_banded((1, 1), ab, b)
    E = -np.gradient(phi, y)     # V/cm
    return y, phi, E


def oxide_voltage_fraction(g: Geom, V_applied: float) -> float:
    """Fraction of the applied voltage dropped across the buried oxide."""
    y, phi, E = solve_potential(g, V_applied)
    i_iface = np.searchsorted(y, g.L_si)
    V_ox = phi[i_iface] - phi[-1]          # potential across the oxide
    return V_ox / V_applied


# ---- avalanche / ionization ---------------------------------------------------
# van Overstraeten-de Man (Si, high-field branch), electrons.
_A_ION = 7.03e5     # cm^-1
_B_ION = 1.231e6    # V/cm

def alpha_ion(E: np.ndarray) -> np.ndarray:
    """Chynoweth ionization coefficient alpha(E) for silicon (cm^-1)."""
    E = np.asarray(E, dtype=float)
    a = np.zeros_like(E)
    m = E > 1.0e4
    a[m] = _A_ION * np.exp(-_B_ION / E[m])
    return a


def ionization_integral(g: Geom, V_applied: float) -> float:
    """I = \int alpha(E) dy through the silicon along AA'."""
    y, phi, E = solve_potential(g, V_applied)
    in_si = y < g.L_si - 1e-12
    ys = y[in_si]
    a = alpha_ion(np.abs(E[in_si]))
    return _trapz(a, ys)


def breakdown_voltage(g: Geom, v_lo: float = 20.0, v_hi: float = 2000.0,
                      tol: float = 0.5) -> float:
    """Bisection on V so that ionization_integral == 1 (avalanche breakdown)."""
    def f(V):
        return ionization_integral(g, V) - 1.0

    flo, fhi = f(v_lo), f(v_hi)
    if flo > 0:            # already breaks at v_lo
        return v_lo
    if fhi < 0:            # never reaches avalanche in range
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


if __name__ == "__main__":
    g = Geom(t_ox=3.0 * UM, S_D=1.5e12)
    for V in (500.0,):
        frac = oxide_voltage_fraction(g, V)
        I = ionization_integral(g, V)
        print(f"V={V:.0f}  oxide-fraction={frac:.3f}  ionization-integral={I:.3f}")
    print("BV(S_D=1.5e12, tox=3um) =", round(breakdown_voltage(g), 1), "V")
