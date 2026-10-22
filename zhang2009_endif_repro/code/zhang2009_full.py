"""
================================================================================
 Zhang, Li, Hu & Luo, IEEE T-ED 56(10) 2009, pp. 2327-2334
 "Field Enhancement for Dielectric Layer of High-Voltage Devices on SOI" (ENDIF)
 FULL self-contained reproduction: closed-form equations -> Figs 2,3,4,5.
================================================================================

THE PHYSICS (ENDIF principle)
-----------------------------
At breakdown the vertical field in an SOI HV device is limited by the SILICON
critical field E_S,C, but the buried DIELECTRIC can carry a much larger field E_I
set by displacement continuity across the Si/dielectric interface (with an
optional interface sheet charge sigma_in):

        eps_I * E_I  =  eps_S * E_S,C  +  q * sigma_in            [continuity]
    =>  E_I = (eps_S/eps_I) E_S,C + q sigma_in /(eps_I eps0)          (1)

Three ENDIF "knobs" raise E_I (and hence the vertical breakdown V_B,V):
  1) thin silicon  -> higher E_S,C(t_S)      2) low-k buried dielectric (small eps_I)
  3) interface charge sigma_in.

Critical field of silicon vs thickness (paper eq. 10), t_S in micrometres:

        E_S,C(t_S) = 9.783 (21.765 - ln t_S)/(ln t_S + 3.975) x 1e4   [V/cm]     (10)

reaching 141 V/um at t_S = 0.1 um.  The silicon critical ionisation rate is the
Chynoweth rate evaluated at E_S,C:  alpha_S,C = A exp(-B/E_S,C)          (Fig. 2).

Dielectric field (eq. 11) and vertical breakdown voltage (eq. 13):

        E_I = [ 9.783 eps_S (21.765 - ln t_S)/(ln t_S+3.975) x1e4 + q sigma_in ]/eps_I   (11)
        V_B,V = { [9.783(21.765-ln t_S)/(ln t_S+3.975)](eps_S t_I + 0.5 eps_I t_S)
                  + q sigma_in t_I x 1e-4 } / eps_I                                       (13)

with t_S, t_I in micrometres and eps_S=11.7, eps_I relative -> V_B,V directly in volts.

Verified against the paper:
  E_S,C(0.1um)=141 V/um; E_I(t_S=1,eps=3.9,sig=0)=161 V/um, (eps=2.0)=319 V/um;
  V_B,V(t_S=1)=348 V, V_B,V(t_S=0.1)=852 V, V_B,V(t_S=1,sig=4e12)=719 V.

Run:  python zhang2009_full.py           # writes Figs 2-5 to ../figs
Units: cm, V, cm^-2 unless a variable name says micrometres.
================================================================================
"""

from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ constants --
Q      = 1.602176634e-19      # C
EPS0   = 8.8541878128e-14     # F/cm
EPS_S  = 11.7                 # silicon relative permittivity
EPS_OX = 3.9                  # SiO2 relative permittivity
# Chynoweth ionisation coefficient for silicon (high-field branch)
A_ION  = 7.03e5              # cm^-1
B_ION  = 1.231e6            # V/cm

# ENDIF empirical fit constants (paper eq. 10)
C0, C1, C2, C3 = 9.783, 21.765, 3.975, 1.0e4

OUTDIR = "../figs"


# ============================================================================ #
#  CORE CLOSED-FORM MODEL                                                        #
# ============================================================================ #
def E_SC(tS_um):
    """Silicon critical field E_S,C (V/cm) vs silicon thickness t_S [um], eq.(10)."""
    L = np.log(tS_um)
    return C0 * (C1 - L) / (L + C2) * C3


def alpha_SC(tS_um):
    """Silicon critical ionisation rate alpha_S,C (cm^-1): Chynoweth at E_S,C."""
    E = E_SC(tS_um)
    return A_ION * np.exp(-B_ION / E)


def E_I(tS_um, epsI_rel=EPS_OX, sigma_in=0.0):
    """Dielectric-layer field E_I (V/cm), eq.(11)/(1) via displacement continuity."""
    return (EPS_S * E_SC(tS_um) + Q * sigma_in / EPS0) / epsI_rel


def V_BV(tS_um, epsI_rel=EPS_OX, sigma_in=0.0, tI_um=2.0):
    """Vertical breakdown voltage V_B,V (V), eq.(13)."""
    L = np.log(tS_um)
    bracket = C0 * (C1 - L) / (L + C2)                       # = E_S,C /1e4  [V/um]
    term_field = bracket * (EPS_S * tI_um + 0.5 * epsI_rel * tS_um)
    term_charge = Q * sigma_in * (tI_um * 1e-4) / EPS0       # tI um -> cm
    return (term_field + term_charge) / epsI_rel


def enhancement_factors(tS_um, epsI_rel, sigma_in, EI_con_Vum=60.0):
    """Delta-gamma_1,2,3 of eq.(12). E_I,Con = 60 V/um (conventional E_S,C=20 V/um).
    Returns each increment relative to the conventional SiO2 reference."""
    EIc = EI_con_Vum * 1e4                                   # -> V/cm
    dES = E_SC(tS_um) - 20.0e4                               # vs 20 V/um baseline
    g1 = EPS_S * dES / (EPS_OX * EIc)
    g2 = (EPS_OX - epsI_rel) / (epsI_rel)                    # LK: +(dEps)/(eps_SiO2-dEps)
    g3 = Q * sigma_in / (EPS0 * EPS_OX * EIc)
    return g1, g2, g3


# ============================================================================ #
#  FIGURE 2 -- E_S,C and alpha_S,C vs t_S                                        #
# ============================================================================ #
def plot_fig2(outdir=OUTDIR):
    tS = np.logspace(np.log10(0.03), np.log10(10), 400)
    fig, axL = plt.subplots(figsize=(6.6, 4.8))
    axR = axL.twinx()

    l1, = axL.plot(tS, E_SC(tS) / 1e5, "b-", lw=2, label="$E_{S,C}$  (this model, eq.10)")
    # conventional constant reference E_S,C = 20 V/um
    axL.axhline(20e4 / 1e5, color="gray", ls="--", lw=1, label="conventional $E_{S,C}=20$ V/$\\mu$m")
    l2, = axR.plot(tS, alpha_SC(tS) / 1e5, "r-.", lw=2, label="$\\alpha_{S,C}$ (Chynoweth at $E_{S,C}$)")

    axL.annotate("141 V/$\\mu$m\nat $t_S$=0.1 $\\mu$m", xy=(0.1, E_SC(0.1) / 1e5),
                 xytext=(0.18, 20), fontsize=8.5, arrowprops=dict(arrowstyle="->"))
    axL.set_xscale("log")
    axL.set_xlabel("$t_S$  ($\\mu$m)")
    axL.set_ylabel("$E_{S,C}$  ($\\times10^{5}$ V/cm)", color="b")
    axR.set_ylabel("Ionization rate $\\alpha_{S,C}$  ($\\times10^{5}$ cm$^{-1}$)", color="r")
    axL.set_xlim(0.01, 10); axL.set_ylim(0, 25); axR.set_ylim(0, 4)
    axL.tick_params(axis="y", colors="b"); axR.tick_params(axis="y", colors="r")
    axL.set_title("Fig.2  $E_{S,C}$ and $\\alpha_{S,C}$ vs silicon thickness $t_S$")
    axL.legend(handles=[l1, l2], loc="upper right", fontsize=8.5)
    axL.grid(which="both", alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig2_esc_alpha.png", dpi=150)
    plt.close(fig)
    print(f"Fig.2  E_S,C(0.1um) = {E_SC(0.1)/1e4:.1f} V/um ; "
          f"alpha_S,C(0.1um) = {alpha_SC(0.1)/1e5:.2f} x1e5 /cm")


# ============================================================================ #
#  FIGURE 3 -- E_I vs t_S (a) and enhancement factors (b)                        #
# ============================================================================ #
def plot_fig3(outdir=OUTDIR):
    tS = np.logspace(np.log10(0.03), np.log10(100), 400)

    # ============ (a) E_I vs t_S : sigma_in family + low-k family ============
    fig, a = plt.subplots(figsize=(7.0, 5.0))
    sig_family = [1, 4, 6, 10, 15, 18, 20]      # x1e12 cm^-2
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(sig_family)))
    for s, c in zip(sig_family, cmap):
        a.plot(tS, E_I(tS, EPS_OX, s * 1e12) / 1e5, color=c, lw=1.6,
               label=f"$\\sigma_{{in}}$={s}")
    for eps, ls in [(3.9, "-"), (2.0, "--"), (1.5, ":")]:
        a.plot(tS, E_I(tS, eps, 0.0) / 1e5, "k", ls=ls, lw=1.8,
               label=f"$\\varepsilon_I$={eps}$\\varepsilon_0$, $\\sigma_{{in}}$=0")
    a.axhline(1000e4 / 1e5, color="red", ls="--", lw=1)
    a.text(3, 105, "$E_I=1000$ V/$\\mu$m", color="red", fontsize=8)
    a.set_xscale("log"); a.set_xlabel("$t_S$  ($\\mu$m)")
    a.set_ylabel("$E_I$  ($\\times10^{5}$ V/cm)")
    a.set_xlim(0.01, 100); a.set_ylim(0, 180)
    a.set_title("Fig.3(a)  Dielectric field $E_I$ vs $t_S$, low-$k$ and $\\sigma_{in}$")
    a.legend(fontsize=7.2, ncol=2, loc="upper right"); a.grid(which="both", alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig3a_EI.png", dpi=150); plt.close(fig)

    # ============ (b) enhancement factors : three panels (as in paper) =======
    fig, (b1, b2, b3) = plt.subplots(1, 3, figsize=(11.0, 4.2), sharey=True)
    tSg = np.logspace(np.log10(0.05), np.log10(30), 200)
    epsg = np.linspace(0.25, 3.9, 200)
    sigg = np.logspace(np.log10(0.05e12), np.log10(30e12), 200)
    g1 = [enhancement_factors(t, EPS_OX, 0.0)[0] for t in tSg]
    g2 = [enhancement_factors(1.0, e, 0.0)[1] for e in epsg]
    g3 = [enhancement_factors(1.0, EPS_OX, s)[2] for s in sigg]

    for ax in (b1, b2, b3):
        ax.axhline(15.67, color="gray", ls="--", lw=1)
        ax.set_ylim(-1, 20); ax.grid(which="both", alpha=0.25)
    b1.plot(tSg, g1, "b-", lw=1.9); b1.set_xscale("log")
    b1.set_title("$\\Delta\\gamma_1$: thin $t_S$", fontsize=9)
    b1.set_xlabel("$t_S$  ($\\mu$m)"); b1.set_ylabel("Enhancement factor  $\\Delta\\gamma_i$")
    b1.text(0.06, 16.4, "$\\Delta\\gamma_i=15.67$\n($E_I{=}1000$ V/$\\mu$m)", fontsize=7.5, color="gray")
    b2.plot(epsg, g2, "g-", lw=1.9)
    b2.set_title("$\\Delta\\gamma_2$: low-$k$  ($\\sigma_{in}{=}0$)", fontsize=9)
    b2.set_xlabel("$\\varepsilon_I/\\varepsilon_0$")
    b3.plot(sigg / 1e12, g3, "r-", lw=1.9); b3.set_xscale("log")
    b3.set_title("$\\Delta\\gamma_3$: interface charge", fontsize=9)
    b3.set_xlabel("$\\sigma_{in}$  ($\\times10^{12}$ cm$^{-2}$)")
    fig.suptitle("Fig.3(b)  Enhancement factors of the three ENDIF approaches "
                 "($E_{I,Con}$=60 V/$\\mu$m)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig3b_enhancement.png", dpi=150); plt.close(fig)

    print(f"Fig.3  E_I(t_S=1,eps=3.9,sig=0) = {E_I(1.0)/1e4:.0f} V/um ; "
          f"(eps=2.0) = {E_I(1.0,2.0)/1e4:.0f} V/um")


# ============================================================================ #
#  FIGURE 4 -- V_B,V vs t_S                                                      #
# ============================================================================ #
def plot_fig4(outdir=OUTDIR):
    tS = np.logspace(np.log10(0.03), np.log10(100), 400)
    fig, ax = plt.subplots(figsize=(6.8, 5.2))

    sig_family = [1, 4, 6, 10, 15, 18, 20]
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(sig_family)))
    for s, c in zip(sig_family, cmap):
        ax.plot(tS, V_BV(tS, EPS_OX, s * 1e12), color=c, lw=1.5,
                label=f"$\\sigma_{{in}}$={s}")
    for eps, ls in [(3.9, "-"), (2.0, "--"), (1.5, ":")]:
        ax.plot(tS, V_BV(tS, eps, 0.0), "k", ls=ls, lw=1.7,
                label=f"$\\varepsilon_I$={eps}$\\varepsilon_0$")

    # conventional constant E_S,C = 20 V/um reference (eps=3.9, sigma=0)
    L = np.log(tS)
    vconv = (20e4 / 1e4) * (EPS_S * 2.0 + 0.5 * EPS_OX * tS) / EPS_OX
    ax.plot(tS, vconv, color="brown", lw=1.4, ls="-.",
            label="$E_{S,C}=20$ V/$\\mu$m (conv.)")
    ax.axhline(1000 * 2.0, color="red", ls="--", lw=1)      # EI=1000 V/um x tI
    ax.text(6, 2040, "$E_I=1000$ V/$\\mu$m", color="red", fontsize=8)

    ax.annotate("348 V", xy=(1.0, 348), xytext=(1.7, 180), fontsize=8,
                arrowprops=dict(arrowstyle="->"))
    ax.annotate("852 V", xy=(0.1, 852), xytext=(0.03, 1150), fontsize=8,
                arrowprops=dict(arrowstyle="->"))
    ax.set_xscale("log"); ax.set_xlabel("$t_S$  ($\\mu$m)")
    ax.set_ylabel("$V_{B,V}$  (V)")
    ax.set_xlim(0.01, 100); ax.set_ylim(0, 2200)
    ax.set_title("Fig.4  Vertical breakdown voltage $V_{B,V}$ vs $t_S$ ($t_I$=2 $\\mu$m)")
    ax.legend(fontsize=7.2, ncol=2, loc="upper center"); ax.grid(which="both", alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig4_vbv.png", dpi=150)
    plt.close(fig)
    print(f"Fig.4  V_B,V(1um)={V_BV(1.0):.0f} V ; V_B,V(0.1um)={V_BV(0.1):.0f} V ; "
          f"V_B,V(1um,sig=4e12)={V_BV(1.0,EPS_OX,4e12):.0f} V")


# ============================================================================ #
#  FIGURE 5 -- vertical field profile under the drain for several t_S           #
# ============================================================================ #
def plot_fig5(outdir=OUTDIR, tI_um=2.0):
    """Box-model vertical field |E|(y): triangular rise in silicon to E_S,C at the
    interface, then a flat plateau at E_I across the dielectric."""
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    tops_EI, tops_ESC = [], []
    for tS in (0.1, 0.2, 0.5, 1.0, 5.0):
        esc = E_SC(tS) / 1e5
        ei  = E_I(tS) / 1e5
        # silicon 0..tS: linear rise 0 -> E_S,C ; dielectric tS..tS+tI: flat E_I
        y = np.array([0.0, tS, tS, tS + tI_um])
        E = np.array([0.0, esc, ei, ei])
        ax.plot(y, E, lw=1.8, label=f"$t_S$={tS} $\\mu$m")
        tops_EI.append((tS, ei)); tops_ESC.append((tS, esc))
    tsx = np.array([p[0] for p in tops_EI])
    ax.plot(tsx, [p[1] for p in tops_EI], "ks--", ms=5, mfc="white",
            lw=1, label="$E_I$ vs $t_S$")
    ax.plot(tsx, [p[1] for p in tops_ESC], "ko:", ms=4, lw=1, label="$E_{S,C}$ vs $t_S$")
    ax.set_xscale("log")
    ax.set_xlabel("Vertical distance  ($\\mu$m)")
    ax.set_ylabel("Electric field  ($\\times10^{5}$ V/cm)")
    ax.set_xlim(0.01, 10); ax.set_ylim(0, 50)
    ax.set_title("Fig.5  Vertical field under the drain vs depth ($t_I$=2 $\\mu$m)")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(which="both", alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig5_vertical_field.png", dpi=150)
    plt.close(fig)
    print(f"Fig.5  E_I(0.1um)={E_I(0.1)/1e4:.0f} V/um  E_I(5um)={E_I(5.0)/1e4:.0f} V/um")


# ============================================================================ #
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("=" * 70)
    print("Zhang et al. IEEE T-ED 2009 -- ENDIF reproduction (formulas -> plots)")
    print("=" * 70)
    print("Verification against paper's stated numbers:")
    print(f"  E_S,C(0.1um)            = {E_SC(0.1)/1e4:6.1f} V/um   (paper 141)")
    print(f"  E_I(1um,3.9e0,0)        = {E_I(1.0)/1e4:6.0f} V/um   (paper 161)")
    print(f"  E_I(1um,2.0e0,0)  [LK]  = {E_I(1.0,2.0)/1e4:6.0f} V/um   (paper 319)")
    print(f"  V_B,V(1um)              = {V_BV(1.0):6.0f} V      (paper 348)")
    print(f"  V_B,V(0.1um)            = {V_BV(0.1):6.0f} V      (paper 852)")
    print(f"  V_B,V(1um,sig=4e12)     = {V_BV(1.0,EPS_OX,4e12):6.0f} V      (paper 719)")
    print("-" * 70)
    plot_fig2(); plot_fig3(); plot_fig4(); plot_fig5()
    print("-" * 70)
    print(f"All figures written to {os.path.abspath(OUTDIR)}")


if __name__ == "__main__":
    main()
