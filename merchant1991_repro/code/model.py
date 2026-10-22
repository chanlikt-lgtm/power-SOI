"""Closed-form breakdown-voltage model for Merchant, Arnold, Baumgart, Mukherjee,
Pein & Pinker, "Realization of High Breakdown Voltage (>700 V) in Thin SOI
Devices," Proc. ISPSD 1991, pp. 31-35.

Physical picture (Fig. 2 of the paper)
---------------------------------------
A lateral P+-N diode in a thin SOI film sits over a buried oxide on a grounded
N+ substrate (shorted to the P+ region). x is lateral (0 at the P+/N junction),
y is vertical, downward, 0 at the top of the SOI. t_s = SOI thickness, t_ox =
buried-oxide thickness, W = lateral depletion width at breakdown (W < L, the
P+ to N+ contact spacing).

Key result (paper eqs. 1-9): choosing the *lateral* doping profile N(x) so that
the lateral field E_x is independent of y (and, with E_y(x,0)=0, the SOI layer
fully depletes) forces N(x) to be linear in x and yields

    E_y(x,y) = V x y / (t W t_s)                                        (15)
    E_x(y)   = -(V/W) (1 - y^2 / (2 t_s t))                              (16)
    t        = t_s/2 + (eps_semi/eps_ox) t_ox                            (8)

For silicon eps_si/eps_ox = 11.7/3.9 = 3 exactly (the Si/oxide displacement-
continuity jump), so t is the "effective" series thickness of half-depleted
silicon plus 3x the oxide. E_x is largest (in magnitude) at y=0; E_y is
largest at (x=W, y=t_s) -- the bottom outer corner of the depletion region.
The buried oxide is always SiO2 here regardless of the semiconductor on top
(eps_ox=3.9 is a material constant of the *device*, not swept).

Avalanche (paper eqs. 10-23): ionization only occurs in the semiconductor
film, using a single Chynoweth coefficient alpha(E) = a exp(-b/E). The paper's
own fit for silicon is a=7.03e5 cm^-1, b=1.468e6 V/cm (valid
1.75e5 <= E <= 6.4e5 V/cm -- note this b differs from the "high-field" Van
Overstraeten-de Man branch used in some later papers, e.g. Nakagawa 1990's
b=1.231e6). This module also carries literature Chynoweth fits for 4H-SiC and
GaN (see the MATERIALS registry below) so the same closed-form breakdown model
can be re-evaluated for wide-bandgap SOI-style films -- see derivation.md sec.6.

Two hole-electron generation points dominate:
  path 1 (eq 17): generated at y=0, 0<x<W -> horizontal ionization integral
      I1 = a W exp(-b W / V)                     (identical to a 1-D p-i-n
                                                    diode of width W)
  path 2 (eqs 18-22): generated at the corner (x=W, y=t_s) -> a vertical leg
  (through the silicon depth, closed form via the exponential integral Ei)
  plus a horizontal leg at y=t_s:
      I2' (eq 20)  = a t_s g [Ei(-g) + exp(-g)/g],   g = b t / V   (eq 21)
      I2''(eq 22)  = a W exp(-b W / (V (1 - t_s/(2t))))
Breakdown (eq 23): the lowest V such that (1-I1)(1-I2'-I2'') = 0, i.e. the
lowest V at which either I1 or (I2'+I2'') first reaches 1.

In the W -> infinity limit (Fig. 4) I1 and I2'' vanish (they decay like
W exp(-const*W)) and breakdown reduces to solving I2'(V) = 1 alone -- a
BV(t_s, t_ox) relation independent of W.

Units: cm, V, cm^-2 unless a variable name says micrometres (t_s_um etc).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import expi

# ---- physical constants -------------------------------------------------------
Q = 1.602176634e-19            # C
EPS0 = 8.8541878128e-14        # F/cm
EPS_OX = 3.9                   # buried oxide is always SiO2, regardless of film
UM = 1.0e-4                    # 1 micron in cm

# ---- buried-oxide dielectric strength (paper: "typically 600 V/um") ------------
E_OX_CRIT = 600.0 / UM          # V/cm


@dataclass(frozen=True)
class Material:
    """A semiconductor film: relative permittivity + single-carrier Chynoweth
    ionization fit alpha(E) = a_ion * exp(-b_ion/E)."""
    name: str
    eps_r: float
    a_ion: float     # cm^-1
    b_ion: float     # V/cm
    citation: str


# Si: the paper's own fit (eq. 10), effectively a combined/undifferentiated
# ionization coefficient. 4H-SiC and GaN: single-carrier (hole) Chynoweth fits
# from the literature, chosen because (a) 4H-SiC avalanche along the field
# directions relevant here is hole-dominated, and (b) using the same carrier
# type for both wide-bandgap materials keeps the three-material comparison
# apples-to-apples. See derivation.md sec.6 for discussion and full citations.
SI = Material("Si", 11.7, 7.03e5, 1.468e6,
              "Merchant et al., Proc. ISPSD 1991, eq.10 (paper's own fit)")
SIC_4H = Material("4H-SiC", 9.7, 3.25e6, 1.71e7,
                   "Raghunathan & Baliga, Appl. Phys. Lett. 72, 3503 (1998) "
                   "-- holes, room temperature")
GAN = Material("GaN", 8.9, 6.4e5, 1.454e7,
               "Cao et al., Appl. Phys. Lett. 112, 262103 (2018) "
               "-- holes, bulk-GaN substrate")

MATERIALS = {"Si": SI, "4H-SiC": SIC_4H, "GaN": GAN}


def t_eff(t_s: float, t_ox: float, mat: Material = SI) -> float:
    """Effective series thickness t = t_s/2 + (eps_semi/eps_ox) t_ox   (eq. 8)."""
    return t_s / 2.0 + (mat.eps_r / EPS_OX) * t_ox


def Ex_top(V: float, W: float) -> float:
    """|E_x| at y=0 (uniform in x by construction) = V/W."""
    return V / W


def Ey_corner(V: float, t_s: float, t_ox: float, mat: Material = SI) -> float:
    """|E_y| at the (x=W, y=t_s) corner = V/t."""
    return V / t_eff(t_s, t_ox, mat)


def oxide_field(V: float, t_s: float, t_ox: float, mat: Material = SI) -> float:
    """Vertical field in the buried oxide (displacement continuity jump)."""
    return (mat.eps_r / EPS_OX) * Ey_corner(V, t_s, t_ox, mat)


def oxide_breakdown_voltage(t_s: float, t_ox: float, mat: Material = SI) -> float:
    """V at which the buried-oxide field reaches E_OX_CRIT."""
    return E_OX_CRIT * t_eff(t_s, t_ox, mat) / (mat.eps_r / EPS_OX)


# ---- ionization integrals (eqs. 17, 20, 22) ------------------------------------
def I1(V: float, W: float, mat: Material = SI) -> float:
    """Path 1: horizontal, generated at y=0 -- same as a 1-D p-i-n diode."""
    return mat.a_ion * W * np.exp(-mat.b_ion * W / V)


def I2_prime(V: float, t_s: float, t_ox: float, mat: Material = SI) -> float:
    """Path 2, vertical leg (through the film depth at x=W), closed form (eq. 20)."""
    t = t_eff(t_s, t_ox, mat)
    g = mat.b_ion * t / V
    return mat.a_ion * t_s * g * (expi(-g) + np.exp(-g) / g)


def I2_dprime(V: float, W: float, t_s: float, t_ox: float, mat: Material = SI) -> float:
    """Path 2, horizontal leg (along y=t_s from x=W to 0), eq. 22."""
    t = t_eff(t_s, t_ox, mat)
    return mat.a_ion * W * np.exp(-mat.b_ion * W / (V * (1.0 - t_s / (2.0 * t))))


def ionization_max(V: float, t_s: float, t_ox: float, W: float | None,
                    mat: Material = SI) -> float:
    """max(I1, I2'+I2'') -- breakdown is where this first equals 1 as V rises.

    W=None selects the W -> infinity limit (Fig. 4): I1, I2'' -> 0, leaving
    only I2'(V, t_s, t_ox).
    """
    i2 = I2_prime(V, t_s, t_ox, mat)
    if W is None:
        return i2
    i1 = I1(V, W, mat)
    i2dp = I2_dprime(V, W, t_s, t_ox, mat)
    return max(i1, i2 + i2dp)


def _bisect_bv(f, v_lo=1.0, v_hi=200000.0):
    """Smallest V in (v_lo, v_hi) with f(V)=ionization-1 == 0, f increasing."""
    flo, fhi = f(v_lo), f(v_hi)
    if flo >= 0:
        return v_lo
    if fhi <= 0:
        return v_hi
    return brentq(f, v_lo, v_hi, xtol=1e-2)


def breakdown_voltage(t_s: float, t_ox: float, W: float | None = None,
                       mat: Material = SI) -> float:
    """Lowest V solving eq. 23 (or its W->infinity limit if W is None)."""
    return _bisect_bv(lambda V: ionization_max(V, t_s, t_ox, W, mat) - 1.0)


def pin_diode_breakdown_voltage(W: float, mat: Material = SI) -> float:
    """Reference 1-D p-i-n diode of intrinsic width W: solve I1(V,W)=1 alone."""
    return _bisect_bv(lambda V: I1(V, W, mat) - 1.0)


def g_star(t_s: float, mat: Material = SI) -> float:
    """Root of a_ion * t_s * g * [Ei(-g) + exp(-g)/g] = 1 for g.

    This is exactly the g = b*t/V that appears inside I2_prime's breakdown
    condition (W->inf limit), but the defining equation involves only t_s and
    mat.a_ion -- t, t_ox, mat.b_ion and mat.eps_r never enter it. So g_star (and
    hence the critical field below) is a function of t_s and the material's
    ionization prefactor alone, NOT of t_ox. See critical_field() and
    derivation.md sec.7 ("the plateau problem").
    """
    f = lambda g: mat.a_ion * t_s * g * (expi(-g) + np.exp(-g) / g) - 1.0
    return brentq(f, 1e-8, 1.0e5, xtol=1e-14, rtol=1e-13)


def critical_field(t_s: float, mat: Material = SI) -> float:
    """Avalanche critical field E_c(t_s) = b_ion / g_star(t_s, mat).

    breakdown_voltage(t_s, t_ox, W=None, mat) == critical_field(t_s, mat) *
    t_eff(t_s, t_ox, mat) exactly (verified numerically) -- i.e. V_BR factors
    into a t_ox-independent "critical field vs t_s" part and a purely
    geometric "effective thickness" part that carries all the t_ox dependence.
    """
    return mat.b_ion / g_star(t_s, mat)


def crossover_film_thickness(mat: Material = SI) -> float:
    """Closed form for the film thickness t_s at which the W->inf semiconductor
    breakdown voltage equals the buried-oxide breakdown voltage.

    Both bv and vox are proportional to t = t_s/2 + (eps_r/eps_ox) t_ox (bv via
    V=b*t/g, vox via E_OX_CRIT*t/ratio), so t cancels out of bv(t_s,t_ox) ==
    vox(t_s,t_ox), leaving a condition on t_s ALONE:

        g_cross = b_ion * (eps_r/eps_ox) / E_OX_CRIT
        a_ion * t_s * g_cross * [Ei(-g_cross) + exp(-g_cross)/g_cross] = 1

    i.e. this crossover is exactly independent of t_ox -- a fact of the model,
    not a numerical coincidence (see derivation.md sec.6). Below this t_s the
    device is oxide-limited; above it, film-limited (at least until the right
    branch of the BV(t_s) "U" climbs back above the oxide cap at large t_s).
    """
    ratio = mat.eps_r / EPS_OX
    g_cross = mat.b_ion * ratio / E_OX_CRIT
    bracket = expi(-g_cross) + np.exp(-g_cross) / g_cross
    return 1.0 / (mat.a_ion * g_cross * bracket)


if __name__ == "__main__":
    t_s, t_ox = 0.1 * UM, 2.0 * UM
    bv_inf = breakdown_voltage(t_s, t_ox)
    print(f"t_s=0.1um, t_ox=2um, W->inf:  BV = {bv_inf:.0f} V  "
          f"(paper Fig.5 plateau as W->100um: ~865 V)")
    print(f"oxide breakdown voltage at this (t_s,t_ox):  {oxide_breakdown_voltage(t_s, t_ox):.0f} V "
          f"(paper: 1200 V dashed line, Fig.4)")
    for W_um in (20, 40, 55, 70):
        bv = breakdown_voltage(t_s, t_ox, W_um * UM)
        pin = pin_diode_breakdown_voltage(W_um * UM)
        print(f"W={W_um:3.0f} um:  BV(SOI)={bv:7.1f} V   BV(1-D p-i-n)={pin:7.1f} V")

    print("\n4H-SiC / GaN extension (t_s=0.1um, t_ox=2um, W->inf):")
    for mat in (SI, SIC_4H, GAN):
        bv = breakdown_voltage(t_s, t_ox, W=None, mat=mat)
        vox = oxide_breakdown_voltage(t_s, t_ox, mat)
        capped = min(bv, vox)
        print(f"  {mat.name:8s}  BV={bv:8.0f} V   oxide_cap={vox:8.0f} V   "
              f"BV_capped={capped:8.0f} V")
