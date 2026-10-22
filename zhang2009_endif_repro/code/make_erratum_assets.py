"""Create annotated crops for the erratum PPT:
  1. paper sentence with '681' circled in red
  2. paper Fig.3(a) zoom with markers where 347 (correct) vs 681 (claimed) sit
  3. paper V_B,V paragraph (719 V) crop -- the self-contradicting evidence
"""
import fitz
from PIL import Image, ImageDraw

SRC = r"E:\MATERIAL\SOI\Field Enhancement for Dielectric Layer of  HVIC  for SOI.pdf"
OUT = r"E:\MATERIAL\SOI\zhang2009_endif_repro\figs_paper"
DPI = 300
Z = DPI / 72.0

d = fitz.open(SRC)
p = d[2]                                     # page 2329 (index 2)

# ---- 1. the erroneous sentence, "681" circled --------------------------------
r681 = p.search_for("681")[0]
# crop a text block around the sentence (left column)
clip = fitz.Rect(40, r681.y0 - 58, 300, r681.y1 + 46)
pix = p.get_pixmap(dpi=DPI, clip=clip)
pix.save(rf"{OUT}\erratum_sentence.png")

img = Image.open(rf"{OUT}\erratum_sentence.png").convert("RGB")
dr = ImageDraw.Draw(img)
# 681 rect in crop-pixel coords
x0 = (r681.x0 - clip.x0) * Z; y0 = (r681.y0 - clip.y0) * Z
x1 = (r681.x1 - clip.x0) * Z; y1 = (r681.y1 - clip.y1 + clip.height) * Z
pad = 14
dr.ellipse([x0 - pad - 6, y0 - pad, x1 + pad + 6, y1 + pad], outline=(220, 0, 0), width=6)
img.save(rf"{OUT}\erratum_sentence_circled.png")
print("1. sentence crop:", img.size)

# ---- 2. V_B,V=719 sentence (the self-contradiction, on page 2330) -------------
p4 = d[3]
r719 = p4.search_for("719")[0]
clip2 = fitz.Rect(40, r719.y0 - 34, 300, r719.y1 + 34)
pix = p4.get_pixmap(dpi=DPI, clip=clip2)
pix.save(rf"{OUT}\erratum_719.png")
img = Image.open(rf"{OUT}\erratum_719.png").convert("RGB")
dr = ImageDraw.Draw(img)
x0 = (r719.x0 - clip2.x0) * Z; y0 = (r719.y0 - clip2.y0) * Z
x1 = (r719.x1 - clip2.x0) * Z; y1 = y0 + r719.height * Z
dr.ellipse([x0 - 18, y0 - 12, x1 + 18, y1 + 12], outline=(0, 120, 0), width=6)
img.save(rf"{OUT}\erratum_719_circled.png")
print("2. 719 crop:", img.size)

# ---- 3. Fig 3(a) zoom with correct-vs-claimed markers -------------------------
# figure occupies roughly x in [0.515,0.97]*W, y in [0.075,0.32]*H of the page
import numpy as np
r = p.rect
clip3 = fitz.Rect(r.x0 + 0.515 * r.width, r.y0 + 0.075 * r.height,
                  r.x0 + 0.97 * r.width, r.y0 + 0.32 * r.height)
pix = p.get_pixmap(dpi=DPI, clip=clip3)
pix.save(rf"{OUT}\erratum_fig3a.png")
img = Image.open(rf"{OUT}\erratum_fig3a.png").convert("RGB")
W, H = img.size
dr = ImageDraw.Draw(img)

# plot-axes pixel calibration (fractions measured from the 400-dpi zoom render):
# x-axis: 0.01 -> 100 um log scale ; y-axis: 0 -> 180 (x1e5 V/cm) linear
AXL, AXR = 0.115 * W, 0.874 * W          # x pixels of ts=0.01 and ts=100
AYB, AYT = 0.870 * H, 0.065 * H          # y pixels of E=0 and E=180
def px(ts):  return AXL + (np.log10(ts) - np.log10(0.01)) / 4.0 * (AXR - AXL)
def py(e5):  return AYB + (e5 - 0.0) / 180.0 * (AYT - AYB)

x1um = px(1.0)
# correct value: 34.6e5 (347 V/um) -- green circle on the sigma=4 curve
cx, cy = x1um, py(34.6)
dr.ellipse([cx - 26, cy - 26, cx + 26, cy + 26], outline=(0, 150, 0), width=7)
# claimed value: 68.1e5 (681 V/um) -- red X where the text says it should be
ex, ey = x1um, py(68.1)
dr.ellipse([ex - 26, ey - 26, ex + 26, ey + 26], outline=(220, 0, 0), width=7)
dr.line([ex - 20, ey - 20, ex + 20, ey + 20], fill=(220, 0, 0), width=7)
dr.line([ex - 20, ey + 20, ex + 20, ey - 20], fill=(220, 0, 0), width=7)
img.save(rf"{OUT}\erratum_fig3a_annotated.png")
print("3. fig3a annotated:", img.size)
print("done")
