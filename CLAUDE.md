# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is **not a software project** — it is a physics reproduction workspace for papers on
SOI (silicon-on-insulator) high-voltage device breakdown. Each source paper (PDF, at the
repo root) has a matching `<author><year>_repro/` directory containing a from-scratch
numerical/analytical reproduction of that paper's key results: derive the governing
equations, implement a solver, regenerate the paper's figures, and produce a side-by-side
comparison report (model vs. paper, numbers and plots).

Root-level PDFs are the source papers:
- `1990ispsd.pdf` / `nakagawa1990---*.pdf` — same paper, two copies (5-page conference
  version and a 12-page longer version) → reproduced in `nakagawa1990_repro/`.
- `Field Enhancement for Dielectric Layer of  HVIC  for SOI.pdf` — Zhang, Li, Hu & Luo,
  IEEE T-ED 2009, the ENDIF paper → reproduced in `zhang2009_endif_repro/`.
- `Realization of High Breakdown Voltage 700 V in Thin SOI Devices.pdf` — Merchant, Arnold,
  Baumgart, Mukherjee, Pein & Pinker, ISPSD 1991 → reproduced in `merchant1991_repro/`.
- `Deriving Fig4 Breakdown Voltage.pdf` — a hand-written derivation note for Merchant's
  Fig. 4 BV formula (supporting material, not itself a source paper).

Some results are cross-checked in **Mathematica notebooks** (`.nb`) that live outside the
`code/` dirs and are not part of the Python pipeline: `VBd vs SOI group by BOX.nb` and
`VBd vs SOI group by BOX  for Si.nb` at the repo root (BV vs SOI thickness grouped by BOX
thickness), and `nakagawa1990_repro/2D   Potential plot.nb` (2-D potential map). Treat these
as reference/verification artifacts — the Python scripts are the primary reproduction.

## Running the code

There's no project-wide build/lint/test tooling — each `<repro>/code/*.py` script is a
standalone entry point run with a plain Python 3 (numpy, scipy, matplotlib, and pymupdf
must be importable; confirmed present under
`C:\Users\User\AppData\Local\Programs\Python\Python311\python.exe`).

```bash
# Nakagawa 1990 (thin-SOI 500V device)
cd nakagawa1990_repro/code
python fig4_potential.py   # Fig.4 — 2-D/3-D potential
python fig4_pde_standalone.py  # Fig.4 — self-contained 2-D Poisson solve + 3-D plot
python fig4_conformal.py   # Fig.4 — same structure by conformal mapping (Schwarz-Christoffel)
python fig4_green.py       # Fig.4 — same structure by Green's function (eigenfunction series) + 3-method compare
python fig5_field.py       # Fig.5 — field vs depth, 4 doses
python fig6_bv.py          # Fig.6 — BV vs dose, tox 1/2/3 um (structure B)
python fig7_bv.py          # Fig.7 — BV vs dose, structure A (2-D caveat)
python build_report.py     # assembles the (earlier) model-only PDF report
python build_full_report.py  # assembles the MAIN theory+figures side-by-side PDF report

# Zhang 2009 (ENDIF)
cd zhang2009_endif_repro/code
python zhang2009_full.py         # closed-form model -> Figs 2-5, prints verification table
python verify_681.py             # three-way check of the paper's "681 V/um" erratum
python make_erratum_assets.py    # annotated paper crops for the erratum
python proof_plot.py             # E_I-vs-sigma proof plots backing the 681 V/um erratum
python build_erratum_ppt.py      # erratum slide deck (pptx)
python build_pdfpages_report.py  # matplotlib-PdfPages side-by-side report (no LaTeX needed)
python fig6_step_device.py       # Fig.6 — step-thickness 2-D Poisson device (qualitative)

# Merchant 1991 (>700V thin-SOI RESURF diode)
cd merchant1991_repro/code
python fig3_potential.py   # Fig.3 — equipotential contours (closed-form phi(x,y))
python fig4_bv_vs_ts.py    # Fig.4 — BV vs SOI thickness, tox 2/4/6 um (W->inf limit)
python fig5_bv_vs_W.py     # Fig.5 — BV vs lateral depletion width W, vs 1-D p-i-n diode
python fig_materials_comparison.py  # extension: same model in 4H-SiC / GaN films
python fig_plateau_analysis.py  # dV_BR/dt_s, log-slope & plateau-width (why Si's U-curve minimum is broad)
python build_report.py     # matplotlib-PdfPages side-by-side report (no LaTeX needed)
```

The Merchant LaTeX report additionally needs a TeX toolchain (this environment uses WSL):
```bash
wsl bash -c "cd /mnt/e/MATERIAL/SOI/merchant1991_repro/report && pdflatex merchant1991_report.tex"
```

The Zhang LaTeX report additionally needs a TeX toolchain (this environment uses WSL):
```bash
wsl bash -c "cd /mnt/e/MATERIAL/SOI/zhang2009_endif_repro/report && pdflatex zhang2009_report.tex"
```

There is no automated test suite. Verification is done by printing model numbers next to
the paper's stated numbers (see each script's `print(...)` calls and the tables in the
READMEs) and by visually comparing generated figures in `figs/` against the paper crops in
`figs_paper/`.

## Repository layout (per reproduction)

Each `*_repro/` directory follows the same shape:
```
README.md            summary, results table (model vs paper), file map, run instructions
derivation.md         (Nakagawa, Merchant) full analytic derivation of the physics
code/                 all Python — models, per-figure plotters, report builders
figs/                 figures generated by the model (output of code/*.py)
figs_paper/            paper pages/figures rendered to PNG for side-by-side comparison
report/               assembled PDF/PPTX deliverables (+ .tex source where applicable)
```

## Architecture notes

**Nakagawa 1990** (`nakagawa1990_repro/`) models a thin (20 µm) Si layer fully depleted
over a buried oxide (BOX) on a grounded substrate, where a shallow N donor sheet shields
the bulk silicon and forces most of the voltage across the oxide (RESURF-like behavior).
- `code/model.py` — 1-D vertical electrostatics: builds a graded grid over
  silicon+oxide, solves Poisson via a tridiagonal (banded) finite-volume discretization
  with Dirichlet BCs, and finds breakdown voltage by bisection on the Chynoweth
  ionization integral (`I = ∫α(E)dy == 1`). This module (`Geom`, `solve_potential`,
  `oxide_voltage_fraction`, `breakdown_voltage`, `alpha_ion`) is imported by nearly
  everything else.
- `code/solver2d.py` — 2-D Poisson solve (sparse direct) over the (x, y) cross-section
  for structure B, used only for the Fig.4 potential map; imports constants from `model.py`.
- `code/fig4_potential.py` … `code/fig7_bv.py` — one script per paper figure; each
  produces a PNG in `figs/`.
- `code/build_report.py` / `code/build_full_report.py` — assemble the PNGs plus
  hand-written theory text into a single PDF via `matplotlib.backends.backend_pdf.PdfPages`
  (no LaTeX dependency). `build_full_report.py` is the maintained, complete report;
  `build_report.py` is an earlier model-only version — prefer the "full" one.
- `code/nakagawa1990_full.py` is a later, single-file consolidation of the PDE setup,
  1D/2D solvers, and all figure generation in one script — check this file first when
  making changes, since it may supersede the split `model.py`/`solver2d.py`/`fig*.py` files.
- `code/fig4_pde_standalone.py` — a fully self-contained (no imports from
  `model.py`/`solver2d.py`) 2-D Poisson solve of the Fig.4 structure-B cross-section:
  a variable-permittivity 5-point finite-volume scheme with harmonic-mean face
  permittivities across the Si/oxide interface (the ×3 D_y jump), drain/source/substrate
  Dirichlet BCs plus a grounded p+ sinker on the right Si edge, followed by a verification
  of the donor-sheet shielding and a 3-D potential plot (`figs/fig4_pde_standalone_3d.png`).
  It is the Nakagawa analogue of Merchant's `fig3_potential_3d.py` — an independent
  corroboration of the electrostatics, kept standalone on purpose.
- `code/fig4_conformal.py` — a THIRD, independent Fig.4 solution by **conformal
  mapping** (no finite-difference grid for the field). It splits the problem by
  superposition `phi = phi_p(y) + phi_h(x,y)`: an exact 1-D vertical particular
  solution `phi_p` carries the donor space charge (this is what makes >half the
  voltage fall across the oxide — a charge effect, not a geometry effect), leaving
  `phi_h` source-free. The two dielectrics are folded in by the effective-oxide-
  thickness stretch `Y(y)` (exact for vertical fields — the one documented
  approximation), and `phi_h` is solved by a genuine Schwarz–Christoffel map of the
  rectangle to the upper half-plane (`w = sn(u,m)`, Jacobi elliptic sine with
  complex argument via the A&S 16.21 addition theorem; modulus `m` fixed by
  `K'(m)/K(m) = 2 H_eff/Lx`) followed by a closed-form Poisson/harmonic-measure
  (subtended-angle) evaluation of the Dirichlet data. `build_phih` returns a point
  evaluator so `verify()` can confirm `phi_h` is harmonic to machine precision via
  the mean-value property (no PDE grid). Reproduces the ~58% oxide voltage share
  (paper: "more than half") and writes `figs/fig4_conformal_3d.png` /
  `figs/fig4_conformal.png`, mirroring `fig4_pde_standalone.py` for side-by-side
  comparison of the two methods. A `charge_free` switch threads through
  `solve`/`build_phih`/`particular_1d` to drop the space charge and recover the
  pure ~29% dielectric divider (the Laplace/field-plate limit). Two extra
  self-checks run from `__main__`: `compare_charge_free()` isolates the donor
  charge's contribution (charged 58% vs charge-free 29% -> the sheet adds ~29
  points; `figs/fig4_conformal_chargefree.png`), and `compare_with_fv()` imports
  `fig4_pde_standalone` and overlays the two independent methods (oxide share
  0.582 vs 0.539; vertical-cut RMS ~24 V ≈ 4.8% of V, with the difference map
  showing the two agree in the bulk/oxide and diverge only in near-surface
  dipoles at the electrode edges, where their top-surface BC idealizations differ;
  `figs/fig4_conformal_vs_fv.png`).
- `code/fig4_green.py` — a THIRD Fig.4 method: the **Green's-function / eigenfunction
  series** of the rectangle. In the same effective-thickness single-dielectric
  coordinate as the conformal file, it writes the Dirichlet Green's function of the
  rectangle as a Laplacian-eigenfunction (sine/sinh) series and gets the whole field
  in one shot: a double-sine sum `u_nm = f_nm/((n pi/a)^2+(m pi/b)^2)` integrates the
  donor space charge directly (no `phi_p` superposition needed for it), a harmonic
  bilinear `B = V(1-x/a)(1-Y/b)` absorbs the single non-zero (drain) corner so the
  wall data vanishes at every corner, and single sine/sinh series carry the reduced
  top/left electrode data. It imports the shared pieces from `fig4_conformal`
  (`PARAMS`, `particular_1d`, `Y_of_y`, `solve`, `_oxide_share`) so the methods are
  strictly comparable. `verify()` checks the Poisson residual in effective coords,
  exact BC recovery, and — the key cross-check — that the Green's function and the
  conformal map agree to ~0.09 V RMS (0.02% of V; both give the ~58% oxide share),
  which pins the residual ~5% gap to the finite-volume solve on the shared BC
  idealization rather than on either analytic method. Writes `figs/fig4_green_3d.png`,
  `figs/fig4_green.png`, and `figs/fig4_three_methods.png` (the vertical-cut overlay
  of all three methods, via `compare_three_methods()`, which imports
  `fig4_pde_standalone`). Together the four `fig4_*` scripts solve the identical
  device by the 1-D model, finite volumes, conformal mapping, and Green's functions.
- The physics reasoning (why the donor sheet shields the silicon, why structure A needs a
  thicker oxide than B — a 2-D effect a 1-D model can't fully capture) is written out in
  `derivation.md`; read it before changing the electrostatics.

**Zhang 2009 / ENDIF** (`zhang2009_endif_repro/`) is purely analytical/closed-form (no PDE
solve): displacement continuity across the Si/dielectric interface sets a dielectric field
`E_I` that can exceed the silicon critical field, and three "knobs" (thin silicon, low-k
dielectric, interface charge) raise it and the vertical breakdown voltage `V_B,V`.
- `code/zhang2009_full.py` is the single-file model: closed-form functions `E_SC`,
  `alpha_SC`, `E_I`, `V_BV`, `enhancement_factors` (paper eqs. 10–13) followed by one
  `plot_figN()` per figure (Figs. 2–5) and a `main()` that prints a model-vs-paper
  verification table before generating plots. This is the one file most changes touch.
- Figs. 6–8 of the paper are the paper's specific 2-D device structures (vs MEDICI +
  experiment). Figs. 7–8 remain out of scope for the analytical core — don't reproduce them
  without discussing scope first. Fig. 6 (the step-thickness drift-region device) now has a
  *qualitative* companion, `code/fig6_step_device.py`, which does a genuine 2-D Poisson solve
  (not MEDICI) to show the silicon field capped near `E_SC`, the buried-oxide field enhanced
  ~3x, and new field peaks at each thickness step — writing `fig6a_potential.png` /
  `fig6b_field3d.png` to `figs/`. It is a standalone illustration, separate from the
  closed-form Figs. 2–5 pipeline. Note the README still lists Figs 6–8 as "not yet done" under
  its Scope section — that text predates `fig6_step_device.py`.
- `code/verify_681.py` documents/proves a specific erratum in the paper (a typo'd field
  value on p.2329); `make_erratum_assets.py`, `proof_plot.py`, and `build_erratum_ppt.py`
  build the annotated proof plots and slide deck explaining it (`proof_plot.py` shows the 681
  value falls on the wrong-permittivity `eps_I=2.0 eps0` line, evidence a low-k value was used
  by mistake) — these are one-off investigative scripts, not part of the main figure pipeline.
- `code/build_pdfpages_report.py` builds the side-by-side report without LaTeX (same
  `PdfPages` approach as Nakagawa's `build_full_report.py`); `report/zhang2009_report.tex`
  is the LaTeX alternative.

**Merchant 1991** (`merchant1991_repro/`) is also purely closed-form/analytical: the
paper picks the field distribution that maximizes breakdown voltage for a given lateral
depletion width `W` (a laterally uniform `E_x`) and derives, via Poisson's equation, the
(linear) lateral doping profile `N(x)` that produces it — the reverse of the usual
RESURF analysis. The payoff is a striking, experimentally-confirmed result: below
~1 µm SOI thickness, breakdown voltage *increases* as the silicon gets thinner.
- `code/model.py` — closed-form fields/potential (`t_eff`, `Ex_top`, `Ey_corner`,
  `oxide_field`/`oxide_breakdown_voltage`) and the two dominant avalanche-ionization
  paths (`I1`, `I2_prime` via the exponential integral `Ei`, `I2_dprime`); breakdown
  voltage is found by bisection (`breakdown_voltage`) on `max(I1, I2'+I2'')==1`, with
  `W=None` selecting the `W→∞` limit used for Fig. 4. `pin_diode_breakdown_voltage` is
  the 1-D p-i-n reference curve (`I1==1` alone). Everything else imports this module.
- `code/fig3_potential.py` — reproduces the Fig. 3 equipotential map by plotting the
  model's closed-form `phi(x,y)` directly (no numerical Poisson solve needed, unlike
  Nakagawa's `solver2d.py`); the one free design parameter `W` is found by inverting
  `breakdown_voltage(t_s, t_ox, W) == 750 V` to match the paper's stated Fig. 3 result.
- `code/fig3_potential_3d.py` — 3-D companion to `fig3_potential.py`: derive/simulate/plot.
  Besides plotting the closed-form `phi(x,y)` as a 3-D surface over the full `0..L` cross-
  section (`fig3_potential_3d.png`), it *independently* solves the 2-D Poisson equation
  `div(K(y) grad phi) = -qN(x)/eps0` over the depleted window with a variable-permittivity
  finite-volume scheme (its own solver, sharp Si/oxide interface placed *on a cell face* at
  y=t_s with a distance-weighted harmonic face permittivity, so the ×3 displacement jump
  emerges to 2nd order rather than being imposed and every Si node carries a full donor-charge
  cell; E_y=0 Neumann top, grounded-substrate and P+ Dirichlet walls, closed-form profile at
  `x=W`) and overlays it on the closed form (`fig3_potential_3d_check.png`), printing the
  numeric-vs-closed-form discrepancy. The closed form is an *exact* continuous solution
  (verified: Si Poisson to machine precision, oxide harmonic, exact ×3 interface displacement
  continuity), and the solver converges to it at 2nd order — ~0.02 V max / 0.005 V RMS
  (~0.003% of V) at the default 16·90·161 grid — so this is a genuine corroboration of the
  derived linear `N(x)`, not a new physics result. NB: an earlier interface treatment (node
  straddling y=t_s) sat at ~2.9%; the face-based interface is what gets it to ~0.003%. Imports
  `phi` and `find_design_W` from `fig3_potential.py`.
- `code/fig4_bv_vs_ts.py` / `code/fig5_bv_vs_W.py` — one script per remaining
  reproducible figure; Fig. 5 (BV vs `W`) matches the paper closely including the
  saturation value, while Fig. 4's absolute BV level runs ~15% high relative to the
  paper's rounded text despite the correct "U" shape and an oxide-breakdown cap that
  matches almost exactly (1210 V modeled vs. 1200 V paper) — see `derivation.md` §4
  for the likely explanation (graph-reading/rounding tolerance, not a formula error).
- Fig. 6 (a photographed oscilloscope I-V trace) is experimental data and out of scope,
  same treatment as Zhang's Figs. 6-8.
- `code/build_report.py` builds the side-by-side report without LaTeX (same `PdfPages`
  approach as Zhang's `build_pdfpages_report.py`, writing
  `report/merchant1991_report_pdfpages.pdf`); `report/merchant1991_report.tex` is the
  LaTeX alternative (writes `report/merchant1991_report.pdf` — note the two report
  formats deliberately use different output filenames so rebuilding one doesn't clobber
  the other).
- `code/model.py` also carries a `Material`/`MATERIALS` registry (Si, 4H-SiC, GaN — the
  latter two from literature Chynoweth fits, see the module docstring for citations) and
  `crossover_film_thickness()`, an extension re-evaluating the same closed-form model for
  wide-bandgap SOI-style films. Every function takes an optional `mat: Material = SI`
  kwarg, defaulting to silicon so the original Figs. 3-5 are unaffected.
  `code/fig_materials_comparison.py` produces the Si/4H-SiC/GaN comparison figures. Key
  finding: at the paper's own thin-film design point, the buried SiO2 (not the
  semiconductor) caps all three materials at nearly the same ~1210-1213 V, and the
  oxide/film crossover thickness is *exactly* independent of `t_ox` (a closed-form
  identity, not a coincidence — see `derivation.md` sec.6). Don't add new materials to
  `MATERIALS` without also citing a real literature source for `a_ion`/`b_ion`/`eps_r` in
  the same style, per user preference for defensible, cited physics constants.
- `code/model.py`'s `g_star`/`critical_field` expose another exact factorization:
  `breakdown_voltage(t_s, t_ox, W=None, mat) == critical_field(t_s, mat) * t_eff(t_s, t_ox, mat)`,
  where `critical_field` depends on `t_s` and the material's `a_ion` alone (not `t_ox`,
  `b_ion`, or `eps_r`). `code/fig_plateau_analysis.py` uses this to explain the "plateau
  problem" (Si's `BV(t_s)` minimum is a broad plateau, 4H-SiC's/GaN's a sharp "V"): the
  minimum occurs where `-d ln(E_c)/d ln(t_s) == (t_s/2)/(t_s/2 + (eps_r/eps_ox) t_ox)`,
  confirmed to match direct numerical minimization to within 0.01 µm. Si's plateau is
  ~2x wider than 4H-SiC's/GaN's at the same `t_ox` because Si's Chynoweth `b` is ~10x
  smaller, making its critical field far less field-sensitive — see `derivation.md` §7.

## Conventions to preserve

- All physics/code uses **CGS-practical units**: lengths in cm, fields in V/cm, doses in
  cm⁻², voltages in V — except where a variable name explicitly says micrometres (`t_S`,
  `tS_um`, etc.), which is why unit conversions (`UM = 1e-4`, `×1e4`, `/1e4`) appear
  throughout. Keep new code consistent with whichever convention the surrounding function
  uses, and comment the unit when it's not obvious from the name.
- Every model script prints its key numbers next to the paper's stated values (see the
  tables in each README) — when changing a model, update these printed comparisons rather
  than removing them, since they're the primary correctness check in lieu of a test suite.
- Figures are written to `figs/`, paper reference crops live in `figs_paper/`; reports pull
  from both. Keep filenames stable (`fig4_potential.png`, `fig3a_EI.png`, etc.) since report
  builders hardcode these paths.
