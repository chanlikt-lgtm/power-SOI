"""Proof plots for the 681 V/um erratum.

(1) E_I vs sigma_in at t_S=1um: the correct eps_I=3.9eps0 line passes through
    (4e12, 347); the claimed 681 lies exactly on the eps_I=2.0eps0 (low-k) line
    -> visual proof the wrong permittivity was used.
(2) Decomposition bars: 347 = 161(field)+186(charge) at 3.9eps0 vs the paper's
    681 = 319(field)+362(charge) at 2.0eps0.
"""
import numpy as np
import matplotlib.pyplot as plt
import zhang2009_full as z

Q, EPS0, EPS_S = z.Q, z.EPS0, z.EPS_S
OUT = "../figs_paper"
tS = 1.0
esc = z.E_SC(tS) / 1e4                      # 53.57 V/um

def EI(epsI, sig):                          # V/um
    return EPS_S * esc / epsI + Q * sig / (epsI * EPS0) * 1e-4

sig = np.linspace(0, 8e12, 400)

# ---------- Plot 1: E_I vs sigma_in ------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 5.4))
ax.plot(sig / 1e12, EI(3.9, sig), color="#0a7d0a", lw=2.6,
        label=r"correct: $\varepsilon_I=3.9\varepsilon_0$  (eq. 11)")
ax.plot(sig / 1e12, EI(2.0, sig), color="#8a8a8a", lw=2.2, ls="--",
        label=r"low-$k$: $\varepsilon_I=2.0\varepsilon_0$")

# correct point
ax.plot(4, EI(3.9, 4e12), "o", color="#0a7d0a", ms=13, zorder=5)
ax.annotate(f"CORRECT\n$E_I$ = 161 + 186 = {EI(3.9,4e12):.0f} V/$\\mu$m",
            xy=(4, EI(3.9, 4e12)), xytext=(4.5, 250),
            fontsize=11, color="#0a7d0a", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#0a7d0a", lw=1.8))
# claimed wrong point
ax.plot(4, 681, "X", color="#c00000", ms=15, zorder=5)
ax.annotate("PAPER TEXT says 681 V/$\\mu$m\n(2x too large)",
            xy=(4, 681), xytext=(0.2, 720),
            fontsize=11, color="#c00000", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#c00000", lw=1.8))
# show 681 sits on the 2.0eps0 line
ax.annotate("681 lands on the $\\varepsilon_I=2.0\\varepsilon_0$ line\n"
            "-> the low-$k$ value was used by mistake",
            xy=(4, EI(2.0, 4e12)), xytext=(4.6, 560),
            fontsize=10, color="#555555",
            arrowprops=dict(arrowstyle="->", color="#555555", lw=1.5))

ax.axvline(4, color="k", ls=":", lw=0.8, alpha=0.5)
ax.set_xlabel("Interface charge  $\\sigma_{in}$  ($\\times10^{12}$ cm$^{-2}$)", fontsize=12)
ax.set_ylabel("Dielectric field  $E_I$  (V/$\\mu$m)", fontsize=12)
ax.set_title("Proof: $E_I$ vs $\\sigma_{in}$ at $t_S$=1 $\\mu$m  "
             "($E_I=\\varepsilon_S E_{S,C}/\\varepsilon_I + q\\sigma_{in}/\\varepsilon_I$)",
             fontsize=12)
ax.set_xlim(0, 8); ax.set_ylim(0, 820)
ax.legend(loc="lower right", fontsize=11); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(f"{OUT}/proof_EI_vs_sigma.png", dpi=160)
plt.close(fig)
print("wrote proof_EI_vs_sigma.png")

# ---------- Plot 2: decomposition bars ---------------------------------------
fig, ax = plt.subplots(figsize=(6.6, 5.2))
labels = ["CORRECT\n$\\varepsilon_I=3.9\\varepsilon_0$", "PAPER TEXT\n(used $2.0\\varepsilon_0$)"]
field = [EPS_S * esc / 3.9, EPS_S * esc / 2.0]      # 161, 319
charge = [Q * 4e12 / (3.9 * EPS0) * 1e-4, Q * 4e12 / (2.0 * EPS0) * 1e-4]  # 186, 362
x = [0, 1]
b1 = ax.bar(x, field, width=0.55, color="#4c78a8", label="field term  $\\varepsilon_S E_{S,C}/\\varepsilon_I$")
b2 = ax.bar(x, charge, width=0.55, bottom=field, color="#e45756",
            label="charge term  $q\\sigma_{in}/\\varepsilon_I$")
totals = ["= 347", "= 675  (paper prints 681)"]
for xi, f, c, tot in zip(x, field, charge, totals):
    ax.text(xi, f / 2, f"{f:.0f}", ha="center", va="center", color="white", weight="bold")
    ax.text(xi, f + c / 2, f"{c:.0f}", ha="center", va="center", color="white", weight="bold")
    ax.text(xi, f + c + 14, tot, ha="center", weight="bold", fontsize=12.5,
            color=("#0a7d0a" if xi == 0 else "#c00000"))
ax.axhline(347, color="#0a7d0a", ls=":", lw=1.2)
ax.text(-0.42, 353, "347 correct", color="#0a7d0a", fontsize=9, ha="left")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("$E_I$  (V/$\\mu$m)", fontsize=12)
ax.set_title("$E_I$ decomposition at $t_S$=1 $\\mu$m, $\\sigma_{in}$=4$\\times10^{12}$ cm$^{-2}$",
             fontsize=12)
ax.set_ylim(0, 780); ax.legend(loc="upper left", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/proof_decomposition.png", dpi=160)
plt.close(fig)
print("wrote proof_decomposition.png")
