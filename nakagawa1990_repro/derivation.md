# Derivation — vertical field-sharing model for the thin-SOI 500 V device

**Paper:** A. Nakagawa, N. Yasuhara, Y. Baba, *"New 500 V Output Device Structures
for Thin Silicon Layer on Silicon Dioxide Film,"* Proc. ISPSD 1990, Tokyo,
pp. 97–101.

## 1. The device and the key idea

A thin (20 µm) silicon layer sits on a buried silicon-dioxide film (BOX) on a
grounded silicon substrate. A lateral high-voltage diode is built in the silicon.
When a large reverse voltage is applied, the thin silicon layer **fully depletes**,
and — this is the paper's central point — **a large fraction of the applied
voltage is forced across the buried oxide**, so the device can block ≥ 500 V even
though the silicon is only 20 µm thick.

A **shallow N-type diffusion layer** (areal dose `S_D`) at the bottom of the
silicon, just above the oxide, controls how the field is shared. Its ionized
(positive) donors "shield" the bulk silicon: they let the field be large in the
oxide while keeping it low in the silicon, so avalanche in the silicon is delayed.

## 2. One-dimensional vertical model (symmetry axis AA′ of structure B)

Along the vertical cut at the drain (structure B, Fig. 3B) the drain contact
(top) is at the applied voltage `V` and the substrate (below the oxide) is at 0.
Take `y` downward, `y=0` at the top of the silicon:

```
  0 ........ L_si        silicon: background donors N_D0 (+ shallow layer S_D near bottom)
  L_si ..... L_si+t_ox   buried oxide (charge free)
  y = L_si+t_ox          substrate, grounded
```

**Electrostatics (full depletion).** With `E = −dφ/dy` and charge density
`ρ = qN` (donors positive),

$$\frac{d}{dy}\!\left(\varepsilon\,\frac{d\phi}{dy}\right) = -qN(y),
\qquad \frac{dE}{dy} = \frac{qN(y)}{\varepsilon}.$$

**Interface condition.** The normal displacement `D = εE` is continuous across
the Si/oxide interface, so the field jumps *up* into the oxide:

$$E_{ox} = \frac{\varepsilon_{si}}{\varepsilon_{ox}}\,E_{si}
        \approx \frac{11.7}{3.9}\,E_{si} \approx 3\,E_{si}.$$

**Shielding by the donor sheet.** Model the shallow layer as an areal sheet
`q S_D` at `y = L_si`. Crossing it *upward* from the oxide, the displacement drops
by `q S_D`, so the field in the bulk silicon is

$$\boxed{\,E_{si,\text{bulk}} = \frac{D_{ox} - q S_D}{\varepsilon_{si}}\,}$$

The oxide displacement `D_ox = ε_ox E_ox` is fixed mostly by the applied voltage
and oxide thickness. As `S_D` increases, `E_si,bulk` **decreases** — the donor
sheet absorbs the oxide's displacement and shields the silicon. This is exactly
the paper's Fig. 5 behaviour: the zero-dose curve (D) has the highest silicon
field; the highest-dose curve (A) the lowest.

**Voltage sharing.** Integrating `E` over each region,

$$V = \underbrace{\int_0^{L_{si}}\! E_{si}\,dy}_{V_{si}}
      + \underbrace{E_{ox}\,t_{ox}}_{V_{ox}},\qquad
  \frac{V_{ox}}{V} \;\text{grows with } S_D \text{ and } t_{ox}.$$

With realistic numbers (`t_ox = 3 µm`, `S_D = 1.5×10¹² cm⁻²`, `V = 500 V`) the
model gives `V_ox/V ≈ 0.60`, reproducing Fig. 4's *"more than half of the applied
voltage falls across the silicon dioxide film."*

## 3. Avalanche breakdown (Figs. 6, 7)

Avalanche happens in the **silicon** (the oxide does not ionize). Using the
Chynoweth ionization coefficient for silicon,

$$\alpha(E) = A\,e^{-B/E},\qquad A = 7.03\times10^{5}\,\text{cm}^{-1},\;
  B = 1.231\times10^{6}\,\text{V/cm},$$

breakdown occurs when the ionization integral along AA′ reaches unity:

$$I(V) = \int_{\text{silicon}} \alpha\big(|E(y)|\big)\,dy = 1
   \;\Longrightarrow\; V = BV.$$

**RESURF optimum.** Two competing effects set `BV(S_D)`:

* too little dose → poor shielding → high silicon field → early avalanche (low BV);
* too much dose → the silicon field reverses/over-shields and the *oxide* field
  (and the field at the bottom corner) grows → BV falls again.

So `BV(S_D)` rises, peaks near the RESURF charge (~1×10¹² cm⁻²), then falls —
the shape of Fig. 6. Thicker oxide shares more voltage → higher peak BV
(≈ 400 / 480 / 565 V for `t_ox` = 1 / 2 / 3 µm in the model, vs. the paper's
≈ 380 / 480 / 580 V).

## 4. Why structure A needs a thicker oxide (Fig. 7) — a 2-D effect

Structures A and B have the **same vertical layer stack**, so a purely vertical
(1-D) model cannot separate them. Their difference is lateral: in structure B the
drain sits directly over the oxide at the symmetry axis, giving the most complete
vertical field-sharing (breakdown point #1 right under the drain). In structure A
the drain is displaced from the axis and its critical point couples less
completely to the grounded substrate through the oxide, so **A needs a thicker
oxide (5 µm) than B (3 µm)** to reach 500 V. The vertical model represents the
favourable B-type sharing and therefore over-predicts structure A; the gap is the
quantitative signature of A's incomplete 2-D field-sharing.

## 5. Assumptions and parameters

| Quantity | Value | Note |
|---|---|---|
| Si layer `L_si` | 20 µm | paper |
| Buried oxide `t_ox` | 1–5 µm | paper (swept) |
| Shallow-layer depth `w_n` | 2 µm | "a few microns" (paper); assumed |
| Background n⁻ doping `N_D0` | 1×10¹⁴ cm⁻³ | **assumed** (not given); sets RESURF optimum |
| Dose `S_D` | 0–2.5×10¹² cm⁻² | paper (swept) |
| `ε_si , ε_ox` | 11.7, 3.9 | — |
| `α(E)` | Chynoweth `A e^{−B/E}` | Si, high-field branch |

`N_D0` and `w_n` are the only free parameters; `N_D0 = 1×10¹⁴ cm⁻³` centres the
Fig. 6 RESURF peak near the paper's optimum and is a physically standard light
drift doping for a 500 V-class 20 µm layer.
