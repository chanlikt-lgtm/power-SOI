"""Assemble the reproduction report PDF (self-contained, no LaTeX required)."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import image as mpimg

import model as m

UM = m.UM
FIGDIR = "../figs"
OUT = "../report/nakagawa1990_reproduction_report.pdf"

TITLE = "Reproduction of Nakagawa, Yasuhara & Baba (ISPSD 1990)"
SUB = ("\"New 500 V Output Device Structures for Thin Silicon Layer on\n"
       "Silicon Dioxide Film\", Proc. ISPSD 1990, Tokyo, pp. 97-101")


def _text_page(pdf, title, body, footer=None):
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    fig.text(0.08, 0.94, title, fontsize=16, weight="bold", va="top")
    fig.text(0.08, 0.89, body, fontsize=10.3, va="top", family="monospace",
             linespacing=1.5)
    if footer:
        fig.text(0.08, 0.04, footer, fontsize=8, color="gray")
    pdf.savefig(fig)
    plt.close(fig)


def _image_page(pdf, title, img_path, caption):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.08, 0.95, title, fontsize=14, weight="bold", va="top")
    ax = fig.add_axes([0.06, 0.30, 0.88, 0.60])
    ax.imshow(mpimg.imread(img_path))
    ax.axis("off")
    fig.text(0.08, 0.25, caption, fontsize=10, va="top", wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    # recompute a few numbers for the results table
    g3 = m.Geom(t_ox=3 * UM, S_D=1.5e12)
    frac = m.oxide_voltage_fraction(g3, 500.0)
    peaks = {}
    doses = np.linspace(0.05e12, 2.6e12, 26)
    for tox in (1, 2, 3):
        bv = [m.breakdown_voltage(m.Geom(t_ox=tox * UM, S_D=SD)) for SD in doses]
        i = int(np.argmax(bv))
        peaks[tox] = (bv[i], doses[i] / 1e12)

    with PdfPages(OUT) as pdf:
        # ---- cover / summary ----
        cover = (
            f"{SUB}\n\n"
            "PURPOSE\n"
            "  Independent derivation, simulation and figure-by-figure comparison\n"
            "  of the paper's thin-SOI 500 V lateral device.\n\n"
            "MODEL\n"
            "  1-D vertical full-depletion electrostatics along the drain symmetry\n"
            "  axis (silicon drift + shallow N donor layer + buried oxide over a\n"
            "  grounded substrate); D = eps*E continuity gives the ~3x field jump\n"
            "  into the oxide.  Avalanche via the Chynoweth ionization integral\n"
            "  (breakdown when INT alpha dy = 1).  A 2-D Poisson solve is used for\n"
            "  the potential map (Fig.4).\n\n"
            "KEY RESULTS  (model vs paper)\n"
            "  ------------------------------------------------------------------\n"
            "  Quantity                         Model        Paper\n"
            "  ------------------------------------------------------------------\n"
            f"  Oxide voltage share (tox=3um)    {frac*100:4.0f}%        > 50%  (Fig.4)\n"
            f"  Peak BV, tox=1um                 {peaks[1][0]:4.0f} V      ~380 V (Fig.6)\n"
            f"  Peak BV, tox=2um                 {peaks[2][0]:4.0f} V      ~480 V (Fig.6)\n"
            f"  Peak BV, tox=3um                 {peaks[3][0]:4.0f} V      ~580 V (Fig.6)\n"
            "  Optimum dose (RESURF)            ~1e12        ~1e12  cm^-2\n"
            "  Structure A needs thicker oxide  yes(2-D)     5um vs 3um (Fig.7)\n"
            "  ------------------------------------------------------------------\n\n"
            "VERDICT\n"
            "  Figs. 4, 5, 6 (structure B, the optimized device) are reproduced\n"
            "  quantitatively.  Fig. 7 (structure A) is a 2-D effect the vertical\n"
            "  model represents only in trend -- see the Fig.7 page.\n"
        )
        _text_page(pdf, TITLE, cover,
                   footer="Reproduction report - generated from code/ in this folder.")

        # ---- figures ----
        _image_page(
            pdf, "Fig. 4  -  Potential distribution, 3-D view (structure B, 500 V)",
            f"{FIGDIR}/fig4_potential_3d.png",
            "Potential as height over the cross-section, drawn to match the paper's\n"
            "Fig.4 wireframe: a 500 V plateau at the n-layer/drain sloping to the p+\n"
            "side, with a steep face over the oxide film at the front-left -- the\n"
            "signature of more than half the voltage dropping across the buried oxide.")

        _image_page(
            pdf, "Fig. 4  -  Potential distribution, 2-D map (structure B, 500 V)",
            f"{FIGDIR}/fig4_potential.png",
            "Same 2-D Poisson solution as a contour map. Equipotentials crowd into the\n"
            "buried oxide; the vertical cut under the drain (lower panel) shows 54% of\n"
            "the 500 V drops across the oxide, reproducing the paper's central claim.")

        _image_page(
            pdf, "Fig. 5  -  Field along the symmetry axis AA'",
            f"{FIGDIR}/fig5_field.png",
            "Field vs depth for four shallow-layer doses (structure B, tox=3um, 500 V).\n"
            "In the silicon the field is low and ordered D>C>B>A: more donor dose gives\n"
            "stronger shielding (lower bulk field). At the Si/oxide interface the field\n"
            "jumps up ~3x (eps ratio) and is constant across the charge-free oxide.")

        _image_page(
            pdf, "Fig. 6  -  Breakdown voltage vs. dose (structure B)",
            f"{FIGDIR}/fig6_bv.png",
            "Model curves (lines) vs digitized paper points (circles) for tox=1,2,3um.\n"
            "The RESURF optimum -- BV rising with dose, peaking near 1e12 cm^-2, then\n"
            "falling -- is reproduced, with peak BV ~ 400/480/565 V matching the paper's\n"
            "~380/480/580 V. This is the paper's headline quantitative result.")

        _image_page(
            pdf, "Fig. 7  -  Breakdown voltage vs. dose (structure A, tox=5um)",
            f"{FIGDIR}/fig7_bv.png",
            "Structures A and B share the same vertical stack, so a vertical model\n"
            "cannot separate them -- their difference is 2-D. The model (B-type\n"
            "sharing) lies above structure A's data; the gap is the signature of A's\n"
            "less complete field-sharing, consistent with the paper's finding that A\n"
            "needs a 5um oxide where B needs 3um to block 500 V.")

        # ---- closing notes ----
        notes = (
            "MODEL ASSUMPTIONS\n"
            "  * Full depletion of the 20 um silicon layer at blocking voltage.\n"
            "  * Background n- doping N_D0 = 1e14 cm^-3 (not stated in the paper;\n"
            "    chosen to centre the Fig.6 RESURF peak; physically standard for a\n"
            "    500 V-class 20 um drift layer).\n"
            "  * Shallow N layer 2 um thick ('a few microns', paper).\n"
            "  * Chynoweth Si ionization: alpha = 7.03e5 * exp(-1.231e6 / E).\n\n"
            "WHAT MATCHES\n"
            "  * Oxide carries > half the voltage (Fig.4): model 54-60%.\n"
            "  * Field shielding by the donor sheet (Fig.5): correct ordering and\n"
            "    the ~3x interface jump.\n"
            "  * BV(dose) RESURF curves and magnitudes (Fig.6): within ~5%.\n\n"
            "KNOWN LIMITATION\n"
            "  * Structure A (Fig.7) is a 2-D effect; the vertical model over-predicts\n"
            "    its absolute BV. Reproducing it quantitatively needs the full 2-D\n"
            "    device with the actual asymmetric drain/sinker doping.\n\n"
            "FILES\n"
            "  code/model.py        1-D electrostatics + avalanche\n"
            "  code/solver2d.py     2-D Poisson (Fig.4)\n"
            "  code/fig4..fig7 *.py figure scripts\n"
            "  derivation.md        full analytic derivation\n"
        )
        _text_page(pdf, "Notes, assumptions and limitations", notes)

    print("wrote", OUT)


if __name__ == "__main__":
    main()
