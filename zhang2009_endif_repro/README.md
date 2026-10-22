# ENDIF reproduction — Zhang, Li, Hu & Luo, IEEE T-ED 2009

Reproduction of **"Field Enhancement for Dielectric Layer of High-Voltage Devices
on Silicon on Insulator"**, IEEE Trans. Electron Devices **56**(10), 2009,
pp. 2327–2334 — the **ENDIF** (enhanced dielectric-layer field) principle.

## The idea
In an SOI HV device the vertical breakdown is limited by the *silicon* critical
field `E_S,C`, but the buried *dielectric* holds a much larger field `E_I` set by
displacement continuity: `ε_I·E_I = ε_S·E_S,C + q·σ_in`. Three knobs raise `E_I`
(and the vertical breakdown `V_B,V`): thin silicon (higher `E_S,C`), low-k buried
dielectric, and interface charge `σ_in`.

## Scope
The paper's **analytical core** (Figs. 2–5, closed-form). Figs. 6–8 are specific
2-D device structures (step-thickness, VK/LK, trench SOI) compared to MEDICI +
experiment — out of scope for this analytical reproduction.

## Verification (model vs. paper — all match to rounding)
| Quantity | Model | Paper |
|---|---|---|
| `E_S,C` at t_S=0.1 µm | 140.8 V/µm | 141 |
| `E_I` (t_S=1, ε=3.9, σ=0) | 161 V/µm | 161 |
| `E_I` (t_S=1, ε=2.0, low-k) | 313 V/µm | 319 |
| `E_I` at t_S=0.1 µm | 422 V/µm | 423 |
| `V_B,V` at t_S=1 µm | 348 V | 348 |
| `V_B,V` at t_S=0.1 µm | 852 V | 852 |
| `V_B,V` (t_S=1, σ=4e12) | 719 V | 719 |

## Erratum found in the paper (p.2329)
The text "E_I = 681 V/µm for σ_in=4e12, ε_I=3.9ε₀" is a **typo**: eq.(11) gives
347, the paper's own V_B,V=719 V implies 347, and its own Fig.3(a) reads ~350.
The 681 comes from using the low-k ε_I=2.0ε₀ (319+362=681). See
`code/verify_681.py` (three-way proof) and the annotated slide deck.

## Files
```
code/zhang2009_full.py              SINGLE-FILE: closed-form model -> Figs 2-5 + self-check
code/verify_681.py                  three-way verification of the 681 erratum
code/make_erratum_assets.py         annotated paper crops (681 circled, Fig3a markers)
code/build_erratum_ppt.py           builds the erratum slide deck
code/build_pdfpages_report.py       matplotlib-PdfPages report (no LaTeX)
report/zhang2009_report.tex/.pdf    LaTeX side-by-side comparison report (theory + figs)
report/zhang2009_report_pdfpages.pdf  same report via matplotlib-PdfPages
report/zhang2009_erratum.pptx       annotated erratum slide deck (5 slides)
figs/                               generated figures (fig2, fig3a, fig3b, fig4, fig5)
figs_paper/                         rendered paper pages, cropped figures, erratum crops
```

Figs 6-8 (step-thickness / VK-LK / trench SOI devices) need a 2-D Poisson device
solve and are not yet done -- see "Scope".

## Run
```bash
python code/zhang2009_full.py       # writes figs/ and prints the verification table
# rebuild the report (WSL TeX Live):
wsl bash -c "cd /mnt/e/MATERIAL/SOI/zhang2009_endif_repro/report && pdflatex zhang2009_report.tex"
```
