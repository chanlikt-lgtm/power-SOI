"""matplotlib-PdfPages version of the ENDIF report (no LaTeX needed).

Mirrors zhang2009_report.tex: cover+table, theory (mathtext eqns), 5 side-by-side
figure pages (paper crop | reproduction), notes + erratum.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import image as mpimg

import zhang2009_full as z

FIG = "../figs"
PAP = "../figs_paper"
OUT = "../report/zhang2009_report_pdfpages.pdf"
A4 = (8.27, 11.69)

TITLE = "Reproduction of Zhang, Li, Hu & Luo (IEEE T-ED 2009)"
SUB = ('"Field Enhancement for Dielectric Layer of HV Devices on SOI"\n'
       'IEEE Trans. Electron Devices 56(10), pp. 2327-2334 -- the ENDIF principle')


def text_page(pdf, title, blocks, footer=None):
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.95, title, fontsize=15, weight="bold", va="top")
    y = 0.89
    for kind, content in blocks:
        if kind == "h":
            y -= 0.008
            fig.text(0.08, y, content, fontsize=11.5, weight="bold", va="top")
            y -= 0.028
        elif kind == "m":
            fig.text(0.11, y, content, fontsize=13, va="top")
            y -= 0.055
        else:
            n = content.count("\n") + 1
            fig.text(0.08, y, content, fontsize=10.2, va="top", family="monospace",
                     linespacing=1.45)
            y -= 0.020 * n + 0.010
    if footer:
        fig.text(0.08, 0.035, footer, fontsize=8, color="gray")
    pdf.savefig(fig); plt.close(fig)


def compare_page(pdf, title, paper_img, model_img, caption):
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.955, title, fontsize=13.5, weight="bold", va="top")
    axp = fig.add_axes([0.04, 0.44, 0.46, 0.46])
    axp.imshow(mpimg.imread(paper_img)); axp.axis("off")
    axp.set_title("Paper (original)", fontsize=10)
    axm = fig.add_axes([0.52, 0.44, 0.46, 0.46])
    axm.imshow(mpimg.imread(model_img)); axm.axis("off")
    axm.set_title("This reproduction", fontsize=10)
    fig.text(0.08, 0.34, caption, fontsize=9.8, va="top", family="monospace",
             linespacing=1.4)
    pdf.savefig(fig); plt.close(fig)


def main():
    with PdfPages(OUT) as pdf:
        # cover
        cover = (
            f"{SUB}\n\n"
            "IDEA (ENDIF)\n"
            "  Vertical breakdown is limited by the SILICON critical field E_S,C,\n"
            "  but the buried DIELECTRIC holds a much larger field E_I set by\n"
            "  displacement continuity: eps_I E_I = eps_S E_S,C + q sigma_in.\n"
            "  Three knobs raise E_I (and V_B,V): thin Si (higher E_S,C), low-k\n"
            "  buried dielectric, and interface charge sigma_in.\n\n"
            "KEY RESULTS  (model vs paper -- all match to rounding)\n"
            "  ----------------------------------------------------------------\n"
            "  Quantity                              Model      Paper\n"
            "  ----------------------------------------------------------------\n"
            f"  E_S,C at t_S=0.1um                     {z.E_SC(0.1)/1e4:5.1f}      141   V/um\n"
            f"  E_I  (t_S=1, eps=3.9, sig=0)           {z.E_I(1.0)/1e4:5.0f}      161   V/um\n"
            f"  E_I  (t_S=1, eps=2.0, low-k)           {z.E_I(1.0,2.0)/1e4:5.0f}      319   V/um\n"
            f"  E_I  at t_S=0.1um                      {z.E_I(0.1)/1e4:5.0f}      423   V/um\n"
            f"  V_B,V at t_S=1um                       {z.V_BV(1.0):5.0f}      348   V\n"
            f"  V_B,V at t_S=0.1um                     {z.V_BV(0.1):5.0f}      852   V\n"
            f"  V_B,V (t_S=1, sig=4e12)                {z.V_BV(1.0,z.EPS_OX,4e12):5.0f}      719   V\n"
            "  ----------------------------------------------------------------\n\n"
            "SCOPE\n"
            "  The paper's analytical core (Figs 2-5). Figs 6-8 are specific 2-D\n"
            "  device structures (step-thickness / VK-LK / trench SOI) vs MEDICI +\n"
            "  experiment -- out of scope for this analytical reproduction."
        )
        text_page(pdf, TITLE, [("t", cover)],
                  footer="matplotlib-PdfPages version. Generated from code/.")

        # theory
        text_page(pdf, "Theory -- the ENDIF model", [
            ("h", "1.  Displacement continuity (why E_I is enhanced)"),
            ("t", "Normal displacement is continuous up to an interface sheet charge:"),
            ("m", r"$\varepsilon_I E_I=\varepsilon_S E_{S,C}+q\sigma_{in}"
                  r"\;\Rightarrow\;E_I=\frac{\varepsilon_S}{\varepsilon_I}E_{S,C}"
                  r"+\frac{q\sigma_{in}}{\varepsilon_I}$"),
            ("t", "Since eps_S/eps_I ~ 11.7/3.9 ~ 3, the dielectric already carries ~3x\n"
                  "the silicon field; low-k (smaller eps_I) or interface charge raise it."),
            ("h", "2.  Silicon critical field vs thickness (eq.10, t_S in um)"),
            ("m", r"$E_{S,C}(t_S)=9.783\,\frac{21.765-\ln t_S}{\ln t_S+3.975}"
                  r"\times10^{4}\ \mathrm{V/cm}$"),
            ("t", "rises sharply for thin silicon -- 141 V/um at t_S=0.1um, far above\n"
                  "the conventional constant 20 V/um. alpha_S,C = A exp(-B/E_S,C) is the\n"
                  "Chynoweth rate at E_S,C (A=7.03e5/cm, B=1.231e6 V/cm)."),
            ("h", "3.  Vertical breakdown voltage (eq.13)"),
            ("m", r"$V_{B,V}=\left\{\left[9.783\frac{21.765-\ln t_S}{\ln t_S+3.975}\right]"
                  r"(\varepsilon_S t_I+\frac{1}{2}\varepsilon_I t_S)"
                  r"+q\sigma_{in}t_I\times10^{-4}\right\}/\varepsilon_I$"),
            ("t", "U-shaped in t_S: high for thin silicon (large E_S,C) and again for\n"
                  "thick silicon; larger sigma_in and smaller eps_I lift the whole family."),
            ("h", "4.  Enhancement factors (eq.12, E_I,Con=60 V/um)"),
            ("m", r"$\Delta\gamma_1=\frac{\varepsilon_S\Delta E_{S,C}}"
                  r"{\varepsilon_{SiO_2}E_{I,Con}},\;\;"
                  r"\Delta\gamma_2=\frac{\Delta\varepsilon_I}"
                  r"{\varepsilon_{SiO_2}-\Delta\varepsilon_I},\;\;"
                  r"\Delta\gamma_3=\frac{q\sigma_{in}}{\varepsilon_{SiO_2}E_{I,Con}}$"),
            ("t", "each reaches 15.67 (E_I=1000 V/um) alone at t_S<50nm, eps_I=0.25eps0,\n"
                  "or sigma_in=2e13 -- so any ONE knob suffices; interface charge is the\n"
                  "most practical."),
        ])

        # figures
        compare_page(pdf, "Fig. 2 -- E_S,C and alpha_S,C vs t_S",
            f"{PAP}/fig2_paper.png", f"{FIG}/fig2_esc_alpha.png",
            "Critical field rises steeply for thin silicon (141 V/um at 0.1um); the\n"
            "critical ionisation rate alpha_S,C (right axis, Chynoweth at E_S,C) grows\n"
            "as t_S shrinks.")
        compare_page(pdf, "Fig. 3(a) -- dielectric field E_I vs t_S, low-k and sigma_in",
            f"{PAP}/fig3a_paper.png", f"{FIG}/fig3a_EI.png",
            "Upper family: sigma_in=1..20e12 (eps_I=3.9eps0). Lower black curves:\n"
            "eps_I=3.9/2.0/1.5 eps0 (sigma_in=0). Both knobs push E_I toward the\n"
            "1000 V/um dashed line.")
        compare_page(pdf, "Fig. 3(b) -- enhancement factors",
            f"{PAP}/fig3b_paper.png", f"{FIG}/fig3b_enhancement.png",
            "Each approach reaches delta-gamma_i=15.67 (E_I=1000 V/um) on its own:\n"
            "thin t_S, low eps_I, or interface charge (crossing at ~2e13, matching\n"
            "the paper).")
        compare_page(pdf, "Fig. 4 -- vertical breakdown voltage V_B,V vs t_S",
            f"{PAP}/fig4_paper.png", f"{FIG}/fig4_vbv.png",
            "U-shaped family (t_I=2um); higher sigma_in / lower eps_I lift V_B,V. The\n"
            "conventional E_S,C=20 V/um curve (bottom) is far lower for thin silicon.\n"
            "Model 348/852 V at t_S=1/0.1um match the paper exactly.")
        compare_page(pdf, "Fig. 5 -- vertical field under the drain vs depth",
            f"{PAP}/fig5_paper.png", f"{FIG}/fig5_vertical_field.png",
            "Each profile: silicon field limited to E_S,C, then jumps to the plateau\n"
            "E_I across the buried dielectric. The two envelopes (E_I and E_S,C vs\n"
            "t_S) rise sharply for thin silicon.")

        # notes + erratum
        text_page(pdf, "Notes, erratum and files", [
            ("h", "What matches"),
            ("t", "  The closed-form ENDIF quantities (E_S,C, alpha_S,C, E_I, V_B,V,\n"
                  "  enhancement factors) reproduce Figs 2-5 and every stated number\n"
                  "  to within rounding."),
            ("h", "Erratum found in the paper's text (p.2329)"),
            ("t", "  \"t_S=1um, eps_I=3.9eps0: E_I = 161 and 681 V/um for sigma_in of\n"
                  "  0 and 4e12\" is inconsistent:\n"
                  "    * eq.(11) gives E_I = 161 + 186 = 347 V/um\n"
                  "    * the paper's own V_B,V=719 V implies the same 347\n"
                  "    * the paper's own Fig.3(a) sigma=4 curve reads ~350 V/um\n"
                  "  Origin: using the LOW-K eps_I=2.0eps0 gives 319+362 = 681 exactly\n"
                  "  -- the LK value was carried into the interface-charge example.\n"
                  "  Text typo only; equations and figures are mutually consistent.\n"
                  "  (See the slide deck zhang2009_erratum.pptx and verify_681.py.)"),
            ("h", "Files"),
            ("t", "  code/zhang2009_full.py        model -> Figs 2-5 + self-check\n"
                  "  code/verify_681.py            three-way erratum verification\n"
                  "  code/build_pdfpages_report.py this report\n"
                  "  report/zhang2009_report.tex   LaTeX version (same content)\n"
                  "  report/zhang2009_erratum.pptx annotated erratum slides"),
        ])

    print("wrote", OUT)


if __name__ == "__main__":
    main()
