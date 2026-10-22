# Derivation — optimum-RESURF SOI lateral diode (Merchant et al., ISPSD 1991)

**Paper:** S. Merchant, E. Arnold, H. Baumgart, S. Mukherjee, H. Pein, R. Pinker,
*"Realization of High Breakdown Voltage (>700 V) in Thin SOI Devices,"*
Proc. ISPSD 1991, pp. 31–35.

## 1. The device and the design goal

A lateral P+–N diode is built in a thin SOI film over a buried oxide (BOX) on a
grounded N+ substrate, with the substrate externally shorted to the P+ region
(Fig. 2). Coordinates: `x` lateral (0 at the P+/N junction), `y` vertical,
downward, 0 at the top of the SOI. `t_s` = SOI thickness, `t_ox` = buried-oxide
thickness, `L` = P+-to-N+ spacing, `W` = lateral depletion width at breakdown
(`W < L`).

Unlike a conventional RESURF design (start from a doping profile, solve for
the field), this paper works **backward**: it picks the field distribution
that gives the highest breakdown voltage for a given `W`, then derives the
lateral doping profile `N(x)` that produces it via Poisson's equation. The
target field distribution is a **laterally uniform** `E_x` — independent of
`y` in the depletion region — so that, looking along `x`, the device behaves
like an ideal 1-D p-i-n diode of width `W`.

## 2. From "uniform lateral field" to the doping profile

Gauss's law relates the doping to the field (paper eq. 2):

$$N(x) = \frac{K_s\epsilon_0}{q}\,\frac{\partial E_x(x,y)}{\partial x}
       = -\frac{K_s\epsilon_0}{q}\,\frac{\partial^2\phi}{\partial x^2}$$

Requiring `∂E_x/∂x = 0` (laterally uniform `E_x`) turns this into (eq. 3–4):

$$N(x) = \frac{K_s\epsilon_0}{q}\frac{\partial E_y(x,y)}{\partial y}
\;\Longrightarrow\;
Q(x) \equiv \int_0^{t_s} N(x)\,dy = \frac{K_s\epsilon_0}{q}\big[E_y(x,t_s)-E_y(x,0)\big]$$

Choosing `E_y(x,0) = 0` (full depletion, minimum surface field) and combining
with the boundary condition `φ(x,0) = Vx/W` (linear top-surface potential, the
signature of a laterally uniform field `E_x(x,0) = -V/W`), the paper obtains
closed forms for everything in the depleted region `0 ≤ x ≤ W`, `0 ≤ y ≤ t_s`:

$$E_y(x,y) = \frac{Vxy}{tWt_s} \qquad(15)$$
$$E_x(y) = -\frac{V}{W}\left(1-\frac{y^2}{2t_st}\right) \qquad(16)$$
$$t = \frac{t_s}{2} + \frac{K_s}{K_{ox}}t_{ox}, \qquad \frac{K_s}{K_{ox}}=\frac{11.7}{3.9}=3 \qquad(8)$$

`t` is the effective series thickness of half-depleted silicon plus 3× the
oxide (the 3× is the Si/oxide displacement-continuity jump). `E_x` is largest
in magnitude at `y=0`; `E_y` is largest at the corner `(x=W, y=t_s)`.
Integrating `Q(x)` (eq. 9) shows the **required lateral doping is linear in
`x`** — the SOI analogue of the Stengl–Gösele "variation of lateral doping"
concept cited in the paper (ref. 3), except derived here directly from the
uniform-field requirement rather than fit to a target field taper.

The potential itself follows directly (used for [fig3_potential.py](code/fig3_potential.py)):

$$\phi(x,y) = \frac{Vx}{W}\left(1-\frac{y^2}{2tt_s}\right), \quad 0\le y\le t_s$$

and, since `E_y` is vertically uniform and `E_x` laterally uniform inside the
charge-free oxide, `φ(x,y) ∝ xy` there — the field jumps by ×3 at the
interface and continues linearly down to the grounded substrate.

## 3. Avalanche breakdown — two dominant ionization paths

Ionization only happens in the silicon, via the Chynoweth coefficient
(paper's own fit, eq. 10):

$$\alpha(E) = a\,e^{-b/E}, \qquad a = 7.03\times10^5\,\text{cm}^{-1},\;\;
b = 1.468\times10^{6}\,\text{V/cm}$$

(valid for `1.75e5 ≤ E ≤ 6.4e5 V/cm`). Breakdown is the lowest `V` at which
*any* generation-to-collection path's ionization integral reaches 1 (eq. 13).
Two points dominate:

**Path 1** — generated at `y=0` (where `|E_x|` is largest and, by
construction, constant across `0<x<W`): hole and electron travel horizontally
to `x=0` and `x=W`. This is exactly the ionization integral of a 1-D p-i-n
diode of width `W` (eq. 17):

$$I_1 = aW\exp\!\left(-\frac{bW}{V}\right)$$

**Path 2** — generated at the corner `(x=W, y=t_s)` (where `|E_y|` is
largest): the electron travels vertically up through the silicon (`I_2'`,
closed-form via the exponential integral `Ei`, eqs. 19–21), the hole travels
horizontally along `y=t_s` back to `x=0` (`I_2''`, eq. 22):

$$I_2' = at_s\gamma\left[\mathrm{Ei}(-\gamma)+\frac{e^{-\gamma}}{\gamma}\right],
\qquad \gamma = \frac{bt}{V} \qquad(20,21)$$
$$I_2'' = aW\exp\!\left(-\frac{bW}{V(1-t_s/2t)}\right) \qquad(22)$$

Breakdown voltage (eq. 23) is the lowest `V` solving
`(1-I_1)(1-I_2'-I_2'') = 0`, i.e. the lowest `V` at which
`max(I_1, I_2'+I_2'')` first reaches 1 — implemented as
[`model.breakdown_voltage`](code/model.py).

## 4. Why thinner SOI can give *higher* breakdown (the paper's headline result)

In the `W → ∞` limit, `I_1` and `I_2''` both vanish (they decay like
`W e^{-\text{const}\cdot W}`), leaving only `I_2'(V, t_s, t_ox) = 1` to solve —
a `BV(t_s, t_ox)` relation independent of `W` ([fig4_bv_vs_ts.py](code/fig4_bv_vs_ts.py)).
Writing `BV ≈ (t_s/2 + 3t_ox)·E_c(t_s)` (using eq. 8 and 16 at the corner),
for fixed critical field `BV` would grow with `t_s`. But the *critical field*
`E_c` itself is set by the ionization integral over a path of length ~`t_s`
(eq. 19): a **shorter** silicon path needs a **higher** field to reach
ionization integral 1, exactly as in a one-sided step junction where the
critical field rises as the depletion width shrinks. Below `t_s ≈ 1 µm` this
"critical-field-rises-fast" effect overtakes the "shorter path" effect, so
`BV` *increases* again as `t_s` decreases — the characteristic "U" shape of
Fig. 4, and the reason `t_s = 0.1 µm` gives *higher* breakdown than `t_s = 1 µm`
even though the silicon is thinner. This trades off against the buried
oxide's own dielectric strength (typically 600 V/µm): for very thin `t_s` the
oxide field (3× the silicon corner field, roughly `3V/t ≈ V/t_ox` as `t_s→0`)
saturates the model's prediction — the dashed cap in Fig. 4
([`model.oxide_breakdown_voltage`](code/model.py)).

## 5. Parameters used

| Quantity | Value | Note |
|---|---|---|
| `K_s, K_ox` | 11.7, 3.9 | silicon, SiO2 relative permittivity |
| Chynoweth `a, b` | 7.03e5 cm⁻¹, 1.468e6 V/cm | paper's own fit (eq. 10); *not* the same `b` used by Nakagawa 1990 or later ENDIF-style papers, which use a different branch/reference |
| `E_ox,crit` | 600 V/µm | paper: "typically 600 V µm⁻¹" |
| Fig. 3/6 device | `t_s`=0.1 µm, `t_ox`=2 µm, `L`=70 µm | paper's fabricated/simulated device |
| Fig. 4 sweep | `t_ox` = 2/4/6 µm, `t_s` = 0.01–100 µm, `W→∞` | |
| Fig. 5 sweep | `t_s`=0.1 µm, `t_ox`=2 µm, `W` = 15–100 µm | |

## 6. Extension — 4H-SiC and GaN films

Nothing in secs. 1–3 is silicon-specific except the single Chynoweth
ionization fit `alpha(E) = a·e^{-b/E}` and the relative permittivity `eps_r`
that enters `t = t_s/2 + (eps_r/eps_ox) t_ox` (eq. 8) — the buried oxide is
still SiO2 (`eps_ox=3.9`) regardless of the semiconductor film on top. Swapping
in literature Chynoweth fits for 4H-SiC and GaN re-evaluates the exact same
closed-form breakdown model ([`model.py`](code/model.py)) for a
"4H-SiC-on-insulator" / "GaN-on-insulator" version of the paper's device.

**Material parameters used** (single-carrier fits, chosen for
comparability — see the module docstring in `model.py`):

| Material | `eps_r` | `a` (cm⁻¹) | `b` (V/cm) | Source |
|---|---|---|---|---|
| Si | 11.7 | 7.03×10⁵ | 1.468×10⁶ | Merchant et al. 1991, eq. 10 (paper's own fit) |
| 4H-SiC | 9.7 | 3.25×10⁶ | 1.71×10⁷ | Raghunathan & Baliga, *Appl. Phys. Lett.* **72**, 3503 (1998) — holes, room temp |
| GaN | 8.9 | 6.4×10⁵ | 1.454×10⁷ | Cao et al., *Appl. Phys. Lett.* **112**, 262103 (2018) — holes, bulk-GaN substrate |

4H-SiC's and GaN's much larger `b` (higher threshold field for ionization)
means their `W→∞` breakdown voltage at the paper's own design point
(`t_s`=0.1 µm, `t_ox`=2 µm) is enormous if only the semiconductor's own
avalanche limit is considered — **3990 V (4H-SiC) and 6947 V (GaN)** vs.
**863 V (Si)** — see `code/fig_materials_comparison.py`. But the buried SiO2
itself breaks down first: its dielectric strength (600 V/µm, unchanged) caps
all three materials at essentially the same **~1210–1213 V**, since the oxide
cap is set almost entirely by `t_ox` (`V_ox,crit ≈ E_ox,crit·t_ox` once
`t_s ≪ t_ox`). **In this architecture, at this film thickness, 4H-SiC and GaN
gain nothing over silicon — the oxide, not the semiconductor, is the
bottleneck.**

**When does the film, not the oxide, become the limit?** Both `bv(t_s,t_ox)`
(the `W→∞` semiconductor limit) and `V_ox(t_s,t_ox)` are proportional to the
same `t = t_s/2 + (eps_r/eps_ox)t_ox` — `bv = b·t/g*(t_s)` where `g*(t_s)` is
the root of eq. 20/21's `I2'=1` condition, and `V_ox = E_ox,crit·t/(eps_r/eps_ox)`.
Setting `bv = V_ox`, the common factor `t` **cancels**, leaving a condition on
`t_s` alone:

$$g_{cross} = \frac{b\,(eps_r/eps_{ox})}{E_{ox,crit}}, \qquad
a\,t_s\,g_{cross}\Big[\mathrm{Ei}(-g_{cross})+\frac{e^{-g_{cross}}}{g_{cross}}\Big] = 1$$

i.e. **the oxide/film crossover thickness is exactly independent of
`t_ox`** — a closed-form fact of the model
([`model.crossover_film_thickness`](code/model.py)), not a numerical
coincidence (confirmed by brute-force root-finding across `t_ox` in
`fig_materials_comparison.py`). Solving it:

| Material | Crossover `t_s` |
|---|---|
| Si | 0.064 µm |
| 4H-SiC | 32.8 µm |
| GaN | 28.8 µm |

Below its crossover thickness, a film is oxide-limited; above it,
semiconductor-limited. Silicon's crossover (0.064 µm) sits right at the edge
of the paper's own thin-film regime, which is why Si's Fig. 4 "U" is
genuinely film-limited over most of its interesting range. 4H-SiC and GaN's
crossover thicknesses (~29–33 µm) are two orders of magnitude larger: **the
whole point of this paper — pushing the film to sub-micron thickness to raise
BV — does not transfer to wide-bandgap films with a plain 2 µm SiO2 BOX**,
because at any sub-micron thickness they are deep in the oxide-limited
regime. Exploiting 4H-SiC's/GaN's intrinsic advantage in this architecture
would require either a much thicker buried oxide or a higher-dielectric-
strength insulator, not a thinner film.

## 7. The plateau problem — why Si's minimum is broad and 4H-SiC's/GaN's is sharp

Fig. 4's "U" curve has a minimum somewhere in `t_s`. Visually, silicon's
minimum is a wide, flat plateau; 4H-SiC's and GaN's are narrow, sharp "V"s.
This falls out of the exact factorization already used in sec. 6:

$$V_{BR}(t_s,t_{ox}) = E_c(t_s)\cdot t_{eff}(t_s,t_{ox}), \qquad
E_c(t_s) := \frac{b}{g^\ast(t_s)}$$

where `g*(t_s)` is the root of `a·t_s·g·[Ei(-g)+e^{-g}/g] = 1` — an equation
in `t_s` and the material's ionization prefactor `a` **alone**
([`model.g_star`](code/model.py), [`model.critical_field`](code/model.py)).
Neither `t_ox` nor `eps_r` enters `E_c(t_s)` at all; they only enter through
the second, purely geometric factor `t_eff`. This factorization was verified
numerically (`breakdown_voltage(t_s,t_ox,W=None,mat)` equals
`critical_field(t_s,mat)*t_eff(t_s,t_ox,mat)` to solver tolerance for all
three materials).

**Minimum condition.** Taking logarithmic derivatives,

$$\frac{d\ln V_{BR}}{d\ln t_s} = \frac{d\ln E_c}{d\ln t_s}
+ \underbrace{\frac{t_s/2}{t_s/2 + r\,t_{ox}}}_{d\ln t_{eff}/d\ln t_s},
\qquad r=\frac{eps_r}{eps_{ox}}$$

so `dV_BR/dt_s = 0` where

$$\underbrace{-\frac{d\ln E_c}{d\ln t_s}}_{\text{LHS: material \& }t_s\text{ only}}
\;=\;
\underbrace{\frac{t_s/2}{t_s/2 + r\,t_{ox}}}_{\text{RHS: material's }r\text{, }t_s\text{, }t_{ox}}
\qquad (\star)$$

— confirmed against a direct numerical minimization of `BV(t_s)` for all 9
(material, `t_ox`) combinations tested: the crossing of $(\star)$ locates the
minimum to within 0.01 µm every time (see the printed check in
[`fig_plateau_analysis.py`](code/fig_plateau_analysis.py)).

**Why Si's plateau is broader.** The LHS of $(\star)$ — the log-slope of the
critical field — depends on material *only* through `a` (not `b`, not
`eps_r`). Numerically, LHS rises far more gently with `t_s` for Si than for
4H-SiC/GaN, because Si's `b` (and hence the whole ionization law's field
sensitivity) is an order of magnitude smaller, so a given fractional change in
`t_s` moves Si's critical field only mildly while it swings 4H-SiC's/GaN's
sharply (the ionization integral itself is far more field-sensitive for the
wide-bandgap fits — see `model.py`'s `b_ion` values, ~10× larger for 4H-SiC and
GaN). Since $(\star)$'s RHS is a slowly-varying, bounded (0 to 1) function of
`t_s` for all materials, a gently-varying LHS crosses it over a *wide* range
of `t_s` where the two sides stay nearly equal (broad plateau); a
steeply-varying LHS crosses it *sharply*, at one well-defined point (narrow
"V").

**Quantified.** Defining "plateau width" as the `t_s`-interval, centered on
the numerically-found minimum, over which `|dV_BR/dt_s|/BV_min` stays below
1%/µm ([`fig_plateau_analysis.fig_plateau_width`](code/fig_plateau_analysis.py)):

| `t_ox` | Si plateau width | 4H-SiC plateau width | GaN plateau width |
|---|---|---|---|
| 2 µm | 1.19 µm | 0.58 µm | 0.73 µm |
| 4 µm | 4.47 µm | 2.27 µm | 2.78 µm |
| 6 µm | 13.03 µm | 5.24 µm | 6.34 µm |

Silicon's plateau is consistently **~2× wider** than 4H-SiC's or GaN's at the
same `t_ox` (and all three widen with `t_ox`, since a thicker oxide dilutes
the `t_s`-dependence of `t_eff` — see the RHS of $(\star)$ shrinking towards 0
as `t_ox` grows for fixed `t_s`, which flattens the whole curve near its
minimum). Permittivity (`eps_r/eps_ox`) only enters through the RHS and is a
secondary effect here; the dominant cause of Si's broader bottom is its much
weaker field-sensitivity of impact ionization (`b_ion`), not permittivity.

Scripts: [`fig_plateau_analysis.py`](code/fig_plateau_analysis.py) produces
`fig_dbv_dts.png` (raw `dV_BR/dt_s` vs `t_s`, one panel per `t_ox`, three
materials each), `fig_log_slope.png` (the $(\star)$ LHS/RHS reproduction with
marked crossings), and `fig_plateau_width.png` (the table above, plotted vs
`t_ox`).
