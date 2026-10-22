# Nakagawa (ISPSD 1990) — thin-SOI 500 V device reproduction

Reproduction of **A. Nakagawa, N. Yasuhara, Y. Baba, "New 500 V Output Device
Structures for Thin Silicon Layer on Silicon Dioxide Film," Proc. ISPSD 1990,
Tokyo, pp. 97–101** (source PDF: `../1990ispsd.pdf`).

## What this is

A thin (20 µm) silicon layer on a buried oxide (BOX) on a grounded substrate.
When the silicon fully depletes, **more than half of the applied 500 V is forced
across the buried oxide**, and a shallow N donor layer shields the bulk silicon
so avalanche is delayed. The reproduction derives the physics, simulates it, and
compares figure-by-figure with the paper.

## Results (model vs paper)

| Quantity | Model | Paper |
|---|---|---|
| Oxide voltage share, tox=3 µm | 54–60 % | > 50 % (Fig. 4) |
| Peak BV, tox = 1 / 2 / 3 µm | 401 / 479 / 565 V | ~380 / 480 / 580 V (Fig. 6) |
| RESURF optimum dose | ~1×10¹² cm⁻² | ~1×10¹² cm⁻² |
| Structure A vs B oxide need | 2-D effect (over-predicted) | 5 µm vs 3 µm (Fig. 7) |

**Figs. 4, 5, 6** (structure B, the optimized device) reproduce quantitatively.
**Fig. 7** (structure A) is a 2-D effect a vertical model captures only in trend —
see `derivation.md` §4 and the Fig. 7 script header.

## Layout

```
code/nakagawa1990_full.py              SINGLE-FILE: PDE -> 1D/2D solvers -> all 2D/3D plots
derivation.md                          full analytic derivation
report/nakagawa1990_full_report.pdf    MAIN: theory + per-figure side-by-side
report/nakagawa1990_reproduction_report.pdf   earlier (model-only) report
code/model.py                          1-D vertical electrostatics + Chynoweth avalanche
code/solver2d.py                       2-D Poisson (potential map)
code/fig4_potential.py                 Fig. 4  — 2-D + 3-D potential distribution
code/fig5_field.py                     Fig. 5  — field vs depth, 4 doses
code/fig6_bv.py                        Fig. 6  — BV vs dose, tox 1/2/3 µm (structure B)
code/fig7_bv.py                        Fig. 7  — BV vs dose, structure A (with caveat)
code/build_full_report.py              assembles the MAIN side-by-side report
code/build_report.py                   assembles the earlier model-only report
figs/                                  generated figures
figs_paper/                            rendered paper pages + cropped paper figures
```

## Run

```bash
# uses the RAG repo venv (numpy, scipy, matplotlib, pymupdf)
cd code
python fig4_potential.py
python fig5_field.py
python fig6_bv.py
python fig7_bv.py
python build_report.py
```

## Key modelling assumptions

- Full depletion of the 20 µm silicon at blocking voltage.
- Background n⁻ doping `N_D0 = 1×10¹⁴ cm⁻³` (not given in the paper; chosen to
  centre the Fig. 6 RESURF peak; standard for a 500 V-class 20 µm drift layer).
- Shallow N layer 2 µm thick ("a few microns", paper).
- Silicon ionization: `α = 7.03×10⁵ · exp(−1.231×10⁶ / E)` (Chynoweth).
