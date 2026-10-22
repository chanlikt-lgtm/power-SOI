"""Full reproduction report: theory + per-figure side-by-side (paper vs model).

One self-contained PDF (no LaTeX):
  1  cover + abstract + model-vs-paper results table
  2  theory I   - device idea + 1-D vertical electrostatics
  3  theory II  - avalanche/RESURF + structure-A (2-D) + assumptions
  4  Fig.4  paper | model (3-D), + 2-D map
  5  Fig.5  paper | model
  6  Fig.6  paper | model
  7  Fig.7  paper | model
  8  notes / limitations
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import image as mpimg

import model as m

UM = m.UM
FIG = "../figs"
PAP = "../figs_paper"
OUT = "../report/nakagawa1990_full_report.pdf"

TITLE = "Reproduction of Nakagawa, Yasuhara & Baba (ISPSD 1990)"
SUB = ('"New 500 V Output Device Structures for Thin Silicon Layer on\n'
       'Silicon Dioxide Film", Proc. ISPSD 1990, Tokyo, pp. 97-101')

A4 = (8.27, 11.69)


def text_page(pdf, title, blocks, footer=None):
    """blocks: list of (kind, content). kind in {'h','t','m'}: header/text/math."""
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.95, title, fontsize=15, weight="bold", va="top")
    y = 0.89
    for kind, content in blocks:
        if kind == "h":
            y -= 0.008
            fig.text(0.08, y, content, fontsize=11.5, weight="bold", va="top")
            y -= 0.028
        elif kind == "m":
            fig.text(0.11, y, content, fontsize=12, va="top")
            y -= 0.052
        else:
            n = content.count("\n") + 1
            fig.text(0.08, y, content, fontsize=10.2, va="top", family="monospace",
                     linespacing=1.45)
            y -= 0.020 * n + 0.010
    if footer:
        fig.text(0.08, 0.035, footer, fontsize=8, color="gray")
    pdf.savefig(fig)
    plt.close(fig)


def compare_page(pdf, title, paper_img, model_img, caption, extra_img=None):
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.955, title, fontsize=14, weight="bold", va="top")

    if extra_img:
        top, ph = 0.60, 0.32
    else:
        top, ph = 0.44, 0.46
    axp = fig.add_axes([0.03, top, 0.47, ph])
    axp.imshow(mpimg.imread(paper_img)); axp.axis("off")
    axp.set_title("Paper (original)", fontsize=10)
    axm = fig.add_axes([0.51, top, 0.47, ph])
    axm.imshow(mpimg.imread(model_img)); axm.axis("off")
    axm.set_title("This reproduction", fontsize=10)

    if extra_img:
        axe = fig.add_axes([0.14, 0.20, 0.72, 0.36])
        axe.imshow(mpimg.imread(extra_img)); axe.axis("off")

    cap_y = 0.16 if extra_img else 0.34
    fig.text(0.08, cap_y, caption, fontsize=9.8, va="top", wrap=True,
             family="monospace", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    g3 = m.Geom(t_ox=3 * UM, S_D=1.5e12)
    frac = m.oxide_voltage_fraction(g3, 500.0)
    doses = np.linspace(0.05e12, 2.6e12, 26)
    peaks = {}
    for tox in (1, 2, 3):
        bv = [m.breakdown_voltage(m.Geom(t_ox=tox * UM, S_D=SD)) for SD in doses]
        i = int(np.argmax(bv)); peaks[tox] = (bv[i], doses[i] / 1e12)

    with PdfPages(OUT) as pdf:
        # ---------- 1. cover ----------
        cover = (
            f"{SUB}\n\n"
            "PURPOSE\n"
            "  Independent derivation, simulation and figure-by-figure comparison\n"
            "  of the paper's thin-SOI 500 V lateral device.\n\n"
            "MODEL (one line)\n"
            "  1-D vertical full-depletion electrostatics along the drain axis +\n"
            "  Chynoweth avalanche integral; 2-D Poisson for the potential map.\n\n"
            "KEY RESULTS  (model vs paper)\n"
            "  ------------------------------------------------------------------\n"
            "  Quantity                         Model        Paper\n"
            "  ------------------------------------------------------------------\n"
            f"  Oxide voltage share (tox=3um)    {frac*100:4.0f}%        > 50%  (Fig.4)\n"
            f"  Peak BV, tox=1um                 {peaks[1][0]:4.0f} V      ~380 V (Fig.6)\n"
            f"  Peak BV, tox=2um                 {peaks[2][0]:4.0f} V      ~480 V (Fig.6)\n"
            f"  Peak BV, tox=3um                 {peaks[3][0]:4.0f} V      ~580 V (Fig.6)\n"
            "  RESURF optimum dose              ~1e12        ~1e12  cm^-2\n"
            "  Structure A needs thicker oxide  yes(2-D)     5um vs 3um (Fig.7)\n"
            "  ------------------------------------------------------------------\n\n"
            "VERDICT\n"
            "  Figs. 4, 5, 6 (structure B, the optimized device) reproduce\n"
            "  quantitatively. Fig. 7 (structure A) is a 2-D effect the vertical\n"
            "  model captures only in trend (see that page)."
        )
        text_page(pdf, TITLE, [("t", cover)],
                  footer="Full report - theory + side-by-side comparison. Generated from code/.")

        # ---------- 2. theory I ----------
        text_page(pdf, "Theory I - device idea and vertical electrostatics", [
            ("h", "1.  The device and the key idea"),
            ("t",
             "A thin (20 um) silicon layer sits on a buried silicon-dioxide film\n"
             "(BOX) on a grounded substrate. Under high reverse bias the thin silicon\n"
             "FULLY DEPLETES, and a large fraction of the applied voltage is forced\n"
             "across the buried oxide -- so the device blocks >= 500 V with only 20 um\n"
             "of silicon. A shallow N donor layer (dose S_D) at the Si/BOX interface\n"
             "'shields' the bulk silicon: it lets the field be large in the oxide\n"
             "while keeping it low in the silicon, delaying avalanche."),
            ("h", "2.  1-D vertical model (drain symmetry axis AA', structure B)"),
            ("t",
             "y downward, y=0 at top of Si; drain at V (top), substrate at 0 (bottom).\n"
             "Full-depletion electrostatics with donor charge rho = qN:"),
            ("m", r"$\frac{d}{dy}\left(\varepsilon\,\frac{d\phi}{dy}\right)=-qN(y),"
                  r"\qquad E=-\frac{d\phi}{dy}$"),
            ("t", "Displacement D = eps*E is continuous at the Si/oxide interface, so\n"
                  "the field jumps UP into the oxide:"),
            ("m", r"$E_{ox}=\frac{\varepsilon_{si}}{\varepsilon_{ox}}\,E_{si}"
                  r"\approx\frac{11.7}{3.9}\,E_{si}\approx 3\,E_{si}$"),
            ("t", "Modelling the shallow layer as an areal sheet q*S_D, crossing it\n"
                  "upward drops the displacement by q*S_D, so the BULK-silicon field is"),
            ("m", r"$E_{si,bulk}=\frac{D_{ox}-q\,S_D}{\varepsilon_{si}}$"),
            ("t", "which DECREASES as S_D rises -> the donor sheet shields the silicon\n"
                  "(exactly the Fig.5 ordering D>C>B>A). Integrating E over each region,\n"
                  "V_ox/V grows with S_D and t_ox; at tox=3um, S_D=1.5e12, V=500 V the\n"
                  "model gives V_ox/V ~ 0.60 -> Fig.4's 'more than half across the oxide'."),
        ])

        # ---------- 3. theory II ----------
        text_page(pdf, "Theory II - avalanche, RESURF, and structure A", [
            ("h", "3.  Avalanche breakdown (Figs. 6, 7)"),
            ("t", "Avalanche occurs in the silicon (oxide does not ionize). With the\n"
                  "Chynoweth ionization coefficient for silicon,"),
            ("m", r"$\alpha(E)=A\,e^{-B/E},\;\;A=7.03\times10^{5}\,\mathrm{cm}^{-1},"
                  r"\;B=1.231\times10^{6}\,\mathrm{V/cm}$"),
            ("t", "breakdown occurs when the ionization integral reaches unity:"),
            ("m", r"$I(V)=\int_{Si}\alpha(|E(y)|)\,dy=1"
                  r"\;\Rightarrow\;V=BV$"),
            ("h", "RESURF optimum"),
            ("t", "Too little dose -> poor shielding -> high Si field -> early avalanche.\n"
                  "Too much dose -> field reverses/over-shields, corner field grows.\n"
                  "So BV(S_D) rises, peaks near the RESURF charge (~1e12 cm^-2), then\n"
                  "falls -- the shape of Fig.6. Thicker oxide shares more voltage ->\n"
                  "higher peak BV (model ~400/480/565 V for tox=1/2/3 um)."),
            ("h", "4.  Why structure A needs a thicker oxide (Fig.7) -- a 2-D effect"),
            ("t", "A and B share the SAME vertical stack, so a vertical (or mirror-\n"
                  "symmetric 2-D) model cannot separate them. In B the drain sits over\n"
                  "the oxide at the axis (full vertical sharing, breakdown pt #1 under\n"
                  "the drain); in A the drain is displaced and couples less to the\n"
                  "substrate, so A needs 5 um oxide vs B's 3 um to reach 500 V. The\n"
                  "vertical model represents the favourable B-type sharing and thus\n"
                  "over-predicts A -- the gap is the signature of A's 2-D disadvantage."),
            ("h", "Assumptions / parameters"),
            ("t", "  L_si=20um; t_ox=1-5um (swept); shallow-layer w_n=2um ('a few um').\n"
                  "  Background n- doping N_D0=1e14 cm^-3 (NOT given in paper; centres\n"
                  "  the Fig.6 RESURF peak; standard for a 500V-class 20um drift layer).\n"
                  "  eps_si=11.7, eps_ox=3.9."),
        ])

        # ---------- 4-7. side-by-side figure pages ----------
        compare_page(
            pdf, "Fig. 4 - Potential distribution (structure B, 500 V)",
            f"{PAP}/fig4_paper_crop.png", f"{FIG}/fig4_potential_3d.png",
            "Potential as height over the cross-section. Both show a 500 V plateau at\n"
            "the n-layer/drain sloping to the p+ side, with a steep face over the OXIDE\n"
            "FILM (left) -- the >half-voltage drop across the buried oxide. The 2-D map\n"
            "below quantifies it: 54% of 500 V falls across the oxide under the drain.",
            extra_img=f"{FIG}/fig4_potential.png")

        compare_page(
            pdf, "Fig. 5 - Electric field along the symmetry axis AA'",
            f"{PAP}/fig5_paper_crop.png", f"{FIG}/fig5_field.png",
            "Field vs depth for four doses (tox=3um, 500 V). In the silicon the field\n"
            "is low and ordered D>C>B>A (more dose = more shielding); it jumps ~3x at\n"
            "the Si/oxide interface and is constant across the charge-free oxide. Note\n"
            "the paper's oxide bars (A highest, D lowest) match the model's oxide field.")

        compare_page(
            pdf, "Fig. 6 - Breakdown voltage vs. dose (structure B)",
            f"{PAP}/fig6_paper_crop.png", f"{FIG}/fig6_bv.png",
            "Model curves (lines) vs digitized paper points (circles) for tox=1,2,3um.\n"
            "The RESURF optimum -- BV rising with dose, peaking near 1e12 cm^-2, then\n"
            "falling -- is reproduced, with peak BV ~400/480/565 V vs the paper's\n"
            "~380/480/580 V. This is the paper's headline quantitative result.")

        compare_page(
            pdf, "Fig. 7 - Breakdown voltage vs. dose (structure A, tox=5um)",
            f"{PAP}/fig7_paper_crop.png", f"{FIG}/fig7_bv.png",
            "Structure A. The vertical model (purple, B-type sharing) lies ABOVE the\n"
            "structure-A data because A and B share the same vertical stack -- their\n"
            "difference is 2-D. The gap (~150-200 V) is the signature of A's less\n"
            "complete field-sharing, consistent with A needing 5um oxide where B needs\n"
            "3um. Quantitative A modelling requires the full 2-D asymmetric device.")

        # ---------- 8. notes ----------
        text_page(pdf, "Notes, files and limitations", [
            ("h", "What matches"),
            ("t", "  * Oxide carries > half the voltage (Fig.4): model 54-60%.\n"
                  "  * Field shielding by the donor sheet (Fig.5): correct ordering,\n"
                  "    ~3x interface jump, constant oxide field.\n"
                  "  * BV(dose) RESURF curves and magnitudes (Fig.6): within ~5%."),
            ("h", "Known limitation"),
            ("t", "  * Structure A (Fig.7) is a 2-D effect; the vertical model over-\n"
                  "    predicts its absolute BV. Reproducing it quantitatively needs the\n"
                  "    full 2-D device with the actual asymmetric drain/sinker doping."),
            ("h", "Files"),
            ("t", "  derivation.md          full analytic derivation\n"
                  "  code/model.py          1-D electrostatics + avalanche\n"
                  "  code/solver2d.py       2-D Poisson (Fig.4)\n"
                  "  code/fig4..fig7_*.py   figure scripts\n"
                  "  code/build_full_report.py  this report\n"
                  "  figs/ , figs_paper/    model figures and cropped paper figures"),
        ])

    print("wrote", OUT)


if __name__ == "__main__":
    main()
