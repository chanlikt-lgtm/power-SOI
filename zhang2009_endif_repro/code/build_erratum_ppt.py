"""Build the erratum slide deck: Zhang 2009 ENDIF text typo (681 V/um).

Proves the error and shows the correct logic + proof plots.
  1 title
  2 the claim (sentence, 681 circled) + the governing equation (correct logic)
  3 PROOF plot: E_I vs sigma_in -- 681 lies on the wrong (low-k) line
  4 correct logic: decomposition 161+186=347, and 3 self-consistency checks
  5 root cause: 681 = low-k 2.0eps0 arithmetic (+ paper's own Fig.3a)
  6 verdict
"""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

FP = r"E:\MATERIAL\SOI\zhang2009_endif_repro\figs_paper"
OUT = r"E:\MATERIAL\SOI\zhang2009_endif_repro\report\zhang2009_erratum.pptx"

NAVY = RGBColor(0x1F, 0x2D, 0x5A)
RED  = RGBColor(0xC0, 0x00, 0x00)
GREEN= RGBColor(0x0A, 0x7D, 0x0A)
GRAY = RGBColor(0x55, 0x55, 0x55)
WHITE= RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = prs.slide_width, prs.slide_height


def tbox(slide, l, t, w, h):
    tb = slide.shapes.add_textbox(l, t, w, h)
    tb.text_frame.word_wrap = True
    return tb.text_frame


def para(tf, text, size=18, bold=False, color=None, align=PP_ALIGN.LEFT,
         first=False, space=6, font="Calibri"):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    p.space_after = Pt(space)
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold; r.font.name = font
    if color is not None:
        r.font.color.rgb = color
    return p


def bar(slide):
    s = slide.shapes.add_shape(1, 0, 0, SW, Inches(0.12))
    s.fill.solid(); s.fill.fore_color.rgb = NAVY; s.line.fill.background()


def box(slide, l, t, w, h, color, line=None):
    s = slide.shapes.add_shape(1, l, t, w, h)
    s.fill.solid(); s.fill.fore_color.rgb = color
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line; s.line.width = Pt(1.5)
    s.shadow.inherit = False
    return s


# ---------------- Slide 1: title ----------------
s = prs.slides.add_slide(BLANK)
box(s, 0, Inches(2.3), SW, Inches(2.9), NAVY)
tf = tbox(s, Inches(0.8), Inches(2.55), Inches(11.7), Inches(2.5))
para(tf, "An Erratum in Zhang et al. (IEEE T-ED 2009)", 34, True, WHITE,
     PP_ALIGN.CENTER, first=True)
para(tf, "The “Eᵢ = 681 V/µm” claim is a text typo — proven, "
         "with the correct value derived and plotted", 20, False,
     RGBColor(0xDD, 0xE3, 0xF0), PP_ALIGN.CENTER)
para(tf, "“Field Enhancement for Dielectric Layer of HV Devices on SOI”, "
         "IEEE T-ED 56(10), p.2329", 14, False, RGBColor(0xBB, 0xC6, 0xE0),
     PP_ALIGN.CENTER)
tf2 = tbox(s, Inches(0.8), Inches(5.6), Inches(11.7), Inches(1.0))
para(tf2, "ENDIF reproduction — verification note", 16, False, GRAY,
     PP_ALIGN.CENTER, first=True)

# ---------------- Slide 2: the claim + governing equation ----------------
s = prs.slides.add_slide(BLANK); bar(s)
tf = tbox(s, Inches(0.5), Inches(0.28), Inches(12.3), Inches(0.7))
para(tf, "The claim, and the equation that governs it", 26, True, NAVY, first=True)
s.shapes.add_picture(rf"{FP}\erratum_sentence_circled.png",
                     Inches(0.55), Inches(1.25), width=Inches(7.0))
tf = tbox(s, Inches(0.6), Inches(4.5), Inches(7.0), Inches(1.0))
para(tf, "Paper: “tₛ=1µm, εᵢ=3.9ε₀  →  "
         "Eᵢ = 161 and 681 V/µm for σᵢₙ = 0 and 4×10¹² cm⁻²”",
     15, False, GRAY, first=True)
# governing eq box
box(s, Inches(8.0), Inches(1.5), Inches(4.9), Inches(2.5), RGBColor(0xEE, 0xF1, 0xF8))
tf = tbox(s, Inches(8.2), Inches(1.65), Inches(4.6), Inches(2.3))
para(tf, "Governing equation (displacement continuity, eq.11):", 15, True, NAVY, first=True)
para(tf, "Eᵢ = εₛ·Eₛ,ᴄ/εᵢ + q·σᵢₙ/εᵢ",
     19, True, RGBColor(0x22,0x22,0x22), space=8)
para(tf, "Only two terms: a field term and a charge term. "
         "Both divide by the SAME εᵢ.", 14, False, GRAY)
tf = tbox(s, Inches(8.0), Inches(4.3), Inches(4.9), Inches(2.6))
para(tf, "With εᵢ = 3.9ε₀:", 16, True, NAVY, first=True)
para(tf, "field term  = εₛ·Eₛ,ᴄ/εᵢ      = 161 V/µm", 15, False, GRAY)
para(tf, "charge term = q·(4×10¹²)/(3.9ε₀) = 186 V/µm", 15, False, GRAY)
para(tf, "Eᵢ = 161 + 186 = 347 V/µm", 18, True, GREEN, space=4)
para(tf, "→ 681 is 2× too large.", 17, True, RED)

# ---------------- Slide 3: PROOF PLOT ----------------
s = prs.slides.add_slide(BLANK); bar(s)
tf = tbox(s, Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.7))
para(tf, "Proof: plot Eᵢ vs σᵢₙ — 681 sits on the WRONG line", 25, True, NAVY, first=True)
s.shapes.add_picture(rf"{FP}\proof_EI_vs_sigma.png",
                     Inches(0.4), Inches(1.15), height=Inches(5.8))
tf = tbox(s, Inches(8.5), Inches(1.4), Inches(4.5), Inches(5.4))
para(tf, "Read the plot:", 18, True, NAVY, first=True, space=8)
para(tf, "●  Green line = correct εᵢ=3.9ε₀.\n"
         "   At σᵢₙ=4×10¹² it gives 347 V/µm.",
     15, False, GREEN, space=10)
para(tf, "✖  Red X = the paper's 681.\n"
         "   It lands exactly on the dashed\n"
         "   εᵢ=2.0ε₀ (low-k) line.", 15, False, RED, space=10)
para(tf, "So 681 was computed with the low-k\npermittivity, not the stated 3.9ε₀.",
     15, True, GRAY, space=10)
para(tf, "(To reach 681 at 3.9ε₀ you would\nneed σᵢₙ ≈ 11×10¹², not 4.)",
     13, False, GRAY)

# ---------------- Slide 4: correct logic + self-consistency ----------------
s = prs.slides.add_slide(BLANK); bar(s)
tf = tbox(s, Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.7))
para(tf, "Correct value 347 V/µm — confirmed three independent ways", 23, True, NAVY, first=True)
s.shapes.add_picture(rf"{FP}\proof_decomposition.png",
                     Inches(0.4), Inches(1.15), height=Inches(5.7))
tf = tbox(s, Inches(7.3), Inches(1.3), Inches(5.7), Inches(5.6))
para(tf, "1.  Direct from eq.(11)", 17, True, NAVY, first=True)
para(tf, "161 + 186 = 347 V/µm  (left bar).", 14, False, GRAY, space=10)
para(tf, "2.  The paper's own Vʙ,ᴠ agrees", 17, True, NAVY)
para(tf, "Same paragraph: Vʙ,ᴠ = 348 → 719 V (eq.13,\nreproduced exactly). "
         "ΔV = 371 V over tᵢ=2µm\n→ ΔEᵢ = 186 V/µm → Eᵢ = 347 (not 681).",
     14, False, GRAY, space=10)
para(tf, "3.  The paper's own Fig.3(a) agrees", 17, True, NAVY)
para(tf, "σᵢₙ=4 curve at tₛ=1 reads ~35×10⁵ V/cm\n(~350 V/µm). 681 would be on the σᵢₙ=10 curve.",
     14, False, GRAY, space=10)
para(tf, "All three give ~347 — the paper's own numbers\ncontradict its own 681.",
     15, True, GREEN)

# ---------------- Slide 5: root cause + paper Fig 3a ----------------
s = prs.slides.add_slide(BLANK); bar(s)
tf = tbox(s, Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.7))
para(tf, "Root cause: the low-k εᵢ was carried over by mistake", 24, True, NAVY, first=True)
tf = tbox(s, Inches(0.7), Inches(1.15), Inches(6.4), Inches(5.6))
para(tf, "The sentence just before discusses the LOW-K case (εᵢ = 2.0ε₀). "
         "Using that εᵢ in BOTH terms:", 17, False, RGBColor(0x33,0x33,0x33),
     first=True, space=12)
para(tf, "Eᵢ = εₛ·Eₛ,ᴄ/2.0ε₀ + q·(4×10¹²)/(2.0ε₀)", 16, True, NAVY, space=4)
para(tf, "     = 319 + 362  =  681 V/µm", 20, True, RED, space=12)
para(tf, "319 is the paper's own rounded low-k Eᵢ (σᵢₙ=0). "
         "So both terms used 2.0ε₀ while the sentence states 3.9ε₀.",
     15, False, GRAY, space=12)
para(tf, "→ A copy/paste of the low-k permittivity into the interface-charge "
         "example. Physics is fine; only the printed number is wrong.",
     15, True, NAVY)
s.shapes.add_picture(rf"{FP}\erratum_fig3a_annotated.png",
                     Inches(7.4), Inches(1.35), width=Inches(5.6))
tf = tbox(s, Inches(7.4), Inches(6.55), Inches(5.6), Inches(0.6))
para(tf, "Paper's Fig.3(a):  ○ 347 (correct) vs ✖ 681 (on the wrong curve).",
     12, False, GRAY, PP_ALIGN.CENTER, first=True)

# ---------------- Slide 6: verdict ----------------
s = prs.slides.add_slide(BLANK); bar(s)
tf = tbox(s, Inches(0.5), Inches(0.3), Inches(12.3), Inches(0.8))
para(tf, "Verdict", 28, True, NAVY, first=True)
tf = tbox(s, Inches(0.8), Inches(1.4), Inches(11.7), Inches(5.4))
para(tf, "✓  The paper's TEXT has a typo: "
         "Eᵢ(σᵢₙ=4×10¹², εᵢ=3.9ε₀) = 347 V/µm, not 681.",
     19, True, GREEN, first=True, space=12)
para(tf, "✓  The paper's EQUATIONS and FIGURES are correct and mutually consistent "
         "— eq.(11), eq.(13), Fig.3(a), Fig.4 and Vʙ,ᴠ=719 V all agree on 347.",
     19, True, NAVY, space=12)
para(tf, "✓  Origin identified and plotted: the low-k εᵢ=2.0ε₀ was used instead "
         "of 3.9ε₀  (681 = 319 + 362).", 19, True, NAVY, space=12)
para(tf, "→  A copy-editing slip, not a physics error. The reproduction matches every "
         "other number in the paper exactly.", 18, False, GRAY, space=12)
para(tf, "Proof code: code/verify_681.py  ·  code/proof_plot.py", 15, False, GRAY)

prs.save(OUT)
print("wrote", OUT)
