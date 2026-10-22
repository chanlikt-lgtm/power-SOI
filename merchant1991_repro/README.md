# Merchant et al. (ISPSD 1991) — thin-SOI >700 V RESURF diode reproduction

Reproduction of **S. Merchant, E. Arnold, H. Baumgart, S. Mukherjee, H. Pein,
R. Pinker, "Realization of High Breakdown Voltage (>700 V) in Thin SOI
Devices," Proc. ISPSD 1991, pp. 31–35** (source PDF:
`../Realization of High Breakdown Voltage 700 V in Thin SOI Devices.pdf`).

## What this is

A lateral P+–N SOI diode whose **lateral doping profile is derived**, not
assumed: the paper picks the field distribution that maximizes breakdown
voltage for a given lateral depletion width `W` — a laterally uniform `E_x`
— and works backward through Poisson's equation to the (linear) doping
profile `N(x)` that produces it. The resulting closed-form field/potential
gives, via two dominant avalanche-ionization paths, a breakdown-voltage model
with a striking result: **below ~1 µm SOI thickness, breakdown voltage rises
as the silicon gets thinner**, because the silicon's own critical field rises
faster than the shrinking path length loses. This predicts (and the paper
experimentally confirms) >700 V breakdown from a 0.1 µm-thick SOI film over a
2 µm buried oxide.

## Results (model vs paper)

| Quantity | Model | Paper |
|---|---|---|
| Oxide breakdown voltage, `t_s`=0.1 µm `t_ox`=2 µm | 1210 V | 1200 V (Fig. 4 dashed line) |
| `BV`, `t_s`=0.1 µm `t_ox`=2 µm, `W`→∞ | 863 V | ~865 V (Fig. 5 plateau) |
| `BV` vs `W` shape (Fig. 5, `t_s`=0.1, `t_ox`=2) | matches closely | — |
| `BV`, `t_s`=0.1 µm `t_ox`=2 µm, `W`=70 µm | 849 V | 750 V (Fig. 3, `L`=70 µm device) |
| Design `W` for `BV`=750 V (Fig. 3 device) | 49.4 µm | < `L`=70 µm ✓ self-consistent |
| Fig. 4 "U" shape (BV vs `t_s`, 3 values of `t_ox`) | reproduced | — |

**Fig. 5** (BV vs lateral depletion width) matches the paper closely, both in
shape and in the ~865 V saturation value. **Fig. 4**'s "U" shape and oxide-cap
mechanism reproduce correctly, but the absolute BV level runs ~15% higher than
the paper's text ("700 V may be achieved" at `t_s`=0.1 µm, `t_ox`=2 µm); since
Fig. 5's plateau and the oxide-breakdown cap both match closely, this is most
likely the paper's rounded reference to its >700 V experimental headline
rather than a precise Fig. 4 graph reading — see `derivation.md` §4.
**Fig. 3**'s equipotential map is reproduced by plotting the model's
closed-form potential directly (no numerical Poisson solve needed, since the
whole point of the paper's design is a closed-form field). **Fig. 6**
(measured I-V oscilloscope trace) is experimental data and out of scope.

## Extension — 4H-SiC and GaN films

Nothing in the derivation is silicon-specific except the single-carrier
Chynoweth ionization fit and the film's relative permittivity (the buried
oxide stays SiO2). Swapping in literature fits for 4H-SiC (Raghunathan &
Baliga 1998) and GaN (Cao et al. 2018) re-evaluates the same closed-form model
for a "wide-bandgap-on-insulator" version of the device
(`code/model.py`'s `Material`/`MATERIALS` registry,
`code/fig_materials_comparison.py`). See `derivation.md` §6 for the full
derivation, citations, and discussion.

| Material | Semiconductor-limit BV (`t_s`=0.1 µm, `t_ox`=2 µm, `W→∞`) | Oxide cap | Actual (capped) BV | Oxide/film crossover `t_s` |
|---|---|---|---|---|
| Si | 863 V | 1210 V | 863 V | 0.064 µm |
| 4H-SiC | 3990 V | 1212 V | 1212 V | 32.8 µm |
| GaN | 6947 V | 1213 V | 1213 V | 28.8 µm |

**Key finding:** at the paper's own thin-film design point, 4H-SiC and GaN's
much higher intrinsic critical field is entirely wasted — the buried SiO2
breaks down first, capping all three materials at essentially the same
~1210–1213 V. The oxide/film crossover thickness turns out to be **exactly
independent of `t_ox`** (a closed-form fact of the model, not a numerical
coincidence — the common factor `t` cancels out of `bv=V_ox`), and for
4H-SiC/GaN it sits around 29–33 µm — two orders of magnitude above the
sub-micron regime that makes this paper's Si result interesting. Thinning a
4H-SiC or GaN film the way this paper thins silicon does not raise its
breakdown voltage in this architecture; only a thicker (or higher-strength)
buried dielectric would let the wide-bandgap film's advantage show through.

## The plateau problem — why Si's BV(t_s) minimum is broad and 4H-SiC's/GaN's is sharp

`V_BR(t_s,t_ox) = E_c(t_s)·t_eff(t_s,t_ox)` factors *exactly* into a
`t_ox`-independent critical-field term `E_c(t_s)` (`code/model.py`'s
`critical_field`/`g_star`) and a purely geometric term. The minimum of the "U"
curve occurs where `-d ln E_c/d ln t_s = (t_s/2)/(t_s/2+r·t_ox)` — confirmed
against direct numerical minimization to within 0.01 µm for all 9
(material, `t_ox`) combinations tested (`code/fig_plateau_analysis.py`).
Because Si's Chynoweth `b` is ~10× smaller than 4H-SiC's/GaN's, its critical
field's log-slope varies far more gently with `t_s`, so it stays "nearly
matched" to the geometric term's slope over a much wider `t_s`-range —
producing a broad plateau instead of a sharp "V". Quantified (plateau width =
`t_s`-range around the minimum with `|dV_BR/dt_s|/BV_min` < 1%/µm):

| `t_ox` | Si | 4H-SiC | GaN |
|---|---|---|---|
| 2 µm | 1.19 µm | 0.58 µm | 0.73 µm |
| 4 µm | 4.47 µm | 2.27 µm | 2.78 µm |
| 6 µm | 13.03 µm | 5.24 µm | 6.34 µm |

Si's plateau is consistently ~2× wider than 4H-SiC's/GaN's at the same
`t_ox`. See `derivation.md` §7 for the full derivation.

## Layout

```
code/model.py               core: closed-form fields/potential (eqs 1-9), ionization
                             integrals I1/I2'/I2'' (eqs 17-22), breakdown_voltage();
                             Material/MATERIALS registry (Si, 4H-SiC, GaN) and
                             crossover_film_thickness() for the materials extension
code/fig3_potential.py       Fig. 3 — equipotential contours (closed-form phi(x,y))
code/fig4_bv_vs_ts.py        Fig. 4 — BV vs SOI thickness, tox=2/4/6um (W->inf limit)
code/fig5_bv_vs_W.py         Fig. 5 — BV vs lateral depletion width W, vs 1-D p-i-n
code/fig_materials_comparison.py   4H-SiC/GaN extension: BV vs t_s and vs W per
                             material, plus the oxide/film crossover-thickness check
code/build_report.py         assembles report/merchant1991_report_pdfpages.pdf (no LaTeX)
derivation.md                full analytic derivation (eqs 1-23) + materials extension (sec.6)
figs/                        generated figures
figs_paper/                  rendered paper pages + cropped paper figures
report/merchant1991_report.tex/.pdf        LaTeX side-by-side comparison report (theory + figs + extension)
report/merchant1991_report_pdfpages.pdf    same report via matplotlib-PdfPages (no LaTeX needed)
```

## Run

```bash
cd code
python fig3_potential.py
python fig4_bv_vs_ts.py
python fig5_bv_vs_W.py
python fig_materials_comparison.py   # 4H-SiC / GaN extension
python build_report.py               # matplotlib-PdfPages report (no LaTeX needed)
```

Needs numpy, scipy, matplotlib, pymupdf (used only for rendering `figs_paper/`
crops from the source PDF).

To rebuild the LaTeX report (WSL TeX Live):
```bash
wsl bash -c "cd /mnt/e/MATERIAL/SOI/merchant1991_repro/report && pdflatex merchant1991_report.tex"
```
(run twice to resolve cross-references).

## Key modelling assumptions

- Ionization coefficients are the paper's own stated fit (eq. 10):
  `a = 7.03×10⁵ cm⁻¹`, `b = 1.468×10⁶ V/cm`, valid `1.75×10⁵ ≤ E ≤ 6.4×10⁵ V/cm`
  — a different `b` than Nakagawa (1990)'s `1.231×10⁶ V/cm`.
- Buried-oxide dielectric strength: `600 V/µm` (paper: "typically").
- Only the two paper-identified dominant ionization paths (generated at
  `y=0` and at the corner `x=W, y=t_s`) are evaluated; eq. 13's product
  formally has infinitely many possible paths.
- Fig. 3's device is assumed to follow the derived linear `N(x)` out to a
  single design width `W` (found by inverting `BV(W)=750 V`); beyond `W` the
  SOI is treated as quasi-neutral at the applied bias, with a plain 1-D
  vertical drop across the oxide there.
