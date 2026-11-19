"""Draw the Paramecium evolution diagram from the manuscript.

Run with python3 paramecium_diagram.py. Install Matplotlib once if needed:
    python3 -m pip install matplotlib
"""

from math import cos, sin, tau
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as Curve
from matplotlib.transforms import Affine2D


# 1. Figure settings ----------------------------------------------------------
# The figure is saved beside this script. Change .pdf to .svg for editable SVG.
OUTPUT_FILE = Path(__file__).with_name("Paramecium_Diagram.pdf")
PAGE_WIDTH, PAGE_HEIGHT = 1752, 897
FONT = "DejaVu Serif"
TITLE_SIZES = (25, 21)
HEADING_SIZE = 20
BODY_SIZE = 16
INK, NAVY = "#181717", "#153E78"
TEAL, RUST, GRAY = "#087D88", "#A94324", "#868889"
RED_TEXT = "#79261F"

# Each box has (left, top, width, height) on a 1752 × 897 page. Y increases down.
BOXES = {
    "treatments": (14, 282, 443, 274),
    "traits":     (635, 120, 523, 226),
    "home":       (673, 496, 455, 166),
    "away":       (1332, 301, 396, 241),
    "legend":     (14, 768, 393, 116),
    "takeaway":   (438, 768, 1264, 111),
}
#                       Fill       Border
BOX_COLORS = {
    "treatments": ("#E0EBD8", "#315A3F"),
    "traits":     ("#F9DFAB", "#B47A1E"),
    "home":       ("#D8E5F1", "#193E57"),
    "away":       ("#F6CCB6", "#9E4528"),
    "legend":     ("#FAFAFA", "#B5B5B5"),
    "takeaway":   ("#E7DFF0", "#A39AB5"),
}


# 2. Figure text --------------------------------------------------------------
# Use \n for a new line and $...$ for a mathematical expression.
TITLE = (
    "Evolutionary Selection, Evolved Phenotypes, and",
    r"Temperature-Dependent Demographic Responses in $\mathit{Paramecium\ caudatum}$",
)
HEADINGS = {
    "treatments": ("Evolutionary Treatments", "(Long-term selection)"),
    "traits": ("Evolved Phenotypic Traits", "(PC1 & PC2)"),
    "home": ("Realized Life History", "(Home Environment)"),
    "away": ("Home-to-Away\nDemographic Change", "(Alternative Temperature)"),
}
TREATMENTS = [("Temperature", "(cool vs. warm)"), ("Predation", "(copepod ±)")]
TRAITS = ["Cell volume", "Aspect ratio", "Swimming speed", "Turning variability"]
HOME_MEASURES = [r"$r$   Intrinsic population growth", r"$K$   Carrying capacity"]
AWAY_EQUATIONS = [
    r"$\Delta r = r_{\mathrm{away}} - r_{\mathrm{home}}$",
    r"$\Delta K = K_{\mathrm{away}} - K_{\mathrm{home}}$",
]
LEGEND_TITLE = "Relationship (model support)"
# Each entry gives the label, arrow color, and whether the arrow is dashed.
LEGEND = [
    ("Significant effect or association\n" + r"($p < 0.05$)", INK, False),
    ("Limited or no association\n" + r"($p \geq 0.05$)", GRAY, True),
]
TAKEAWAY = (
    "Take-home message: Warming and predation altered phenotypic traits and demographic "
    "characteristics, while evolved traits were more informative about demographic "
    "responses when populations experienced alternative temperatures."
)
TAKEAWAY_WRAP = 86  # Approximate number of characters per line


# 3. Connections and their labels --------------------------------------------
# For each connection, give the box name and a point within that box.
# (0, 0) is its top left, (1, 1) its bottom right. Label offsets are in page units.
def draw_connections(ax):
    def connect(source, start, target, end, color, *, dashed=False, curved=False):
        a, b = point(source, *start), point(target, *end)
        arrow(ax, a, b, color, dashed=dashed, curved=curved)
        return (a[0] + b[0]) / 2, (a[1] + b[1]) / 2

    mid = connect("treatments", (1, .13), "traits", (0, .44), TEAL)
    note(ax, mid, (-120, -106), "Altered traits", "(strong effect)", color=NAVY)

    mid = connect("treatments", (1, .81), "home", (0, .40), TEAL)
    note(ax, mid, (-22, 29), "Realized effects", "(at home)",
         "Warming ↑ $r$\nWarming & predation ↑ $K$\n(additive models)",
         color=NAVY, detail_size=13)

    mid = connect("traits", (.42, 1), "home", (.40, 0), GRAY, dashed=True)
    note(ax, mid, (152, -63), "Limited association", "(at home)",
         "PC2 ↑ $K$ only\n(others n.s.)")

    mid = connect("traits", (1, .40), "away", (0, .09), RUST)
    note(ax, mid, (175, -123), "Associated with change", "(stronger at away)",
         "PC1 → $Δr$\nPC2 → $ΔK$", color=RED_TEXT)

    mid = connect("home", (1, .55), "away", (0, .78), GRAY, dashed=True)
    note(ax, mid, (-12, 60), qualifier="Not a primary\npredictor", color="#505054")

    mid = connect("treatments", (.60, 1), "away", (.40, 1), RUST, curved=True)
    note(ax, mid, (0, 124), "Treatment history also affected change",
         details=r"Predation ↑ $\Delta r$; warming & predation ↓ $\Delta K$",
         color=RED_TEXT)


# 4. Place text and draw the boxes --------------------------------------------
def point(box, x, y):
    left, top, width, height = BOXES[box]
    return left + x * width, top + y * height


def label(ax, position, words, size=BODY_SIZE, *, color=INK, bold=False,
          italic=False, align="center", valign="center", spacing=1.18):
    ax.text(*position, words, fontsize=size, color=color, ha=align, va=valign,
            fontweight="bold" if bold else "normal",
            fontstyle="italic" if italic else "normal", linespacing=spacing, zorder=8)


def note(ax, midpoint, offset, heading="", qualifier="", details="", *,
         color=INK, detail_size=15):
    x, y = midpoint[0] + offset[0], midpoint[1] + offset[1]
    if heading:
        label(ax, (x, y), heading, BODY_SIZE, bold=True, color=color, valign="top")
        y += 29
    if qualifier:
        label(ax, (x, y), qualifier, BODY_SIZE, italic=True, color=color, valign="top")
        y += 33
    if details:
        label(ax, (x, y), details, detail_size, valign="top", spacing=1.32)


def arrow(ax, start, end, color, *, dashed=False, curved=False, legend=False):
    geometry = dict(posA=start, posB=end, shrinkA=1, shrinkB=3)
    if curved:  # Route the treatment-to-away arrow just above the footer.
        dx = end[0] - start[0]
        bend_y = BOXES["takeaway"][1] + 56
        controls = [(start[0] + .01 * dx, bend_y), (end[0] - .06 * dx, bend_y)]
        geometry = dict(path=Curve([start, *controls, end], [Curve.MOVETO] + [Curve.CURVE4] * 3))
    ax.add_patch(FancyArrowPatch(**geometry, arrowstyle="-|>", color=color,
                                 linewidth=3 if legend else 4.6,
                                 mutation_scale=15 if legend else 26,
                                 linestyle="--" if dashed else "-",
                                 zorder=4 if legend else 1))


def draw_boxes(ax):
    for name, (left, top, width, height) in BOXES.items():
        face, edge = BOX_COLORS[name]
        radius = min(42, height / 4)
        ax.add_patch(FancyBboxPatch((left, top), width, height, facecolor=face,
                     edgecolor=edge, linewidth=2, zorder=2,
                     boxstyle=f"round,pad=0,rounding_size={radius}"))
        if name in HEADINGS:
            ax.add_patch(FancyBboxPatch((left + 7, top + 7), width - 14, height - 14,
                         facecolor="none", edgecolor=edge, alpha=.5, linewidth=.6,
                         zorder=3, boxstyle=f"round,pad=0,rounding_size={radius - 7}"))


def draw_panel_text(ax):
    for name, (heading, subtitle) in HEADINGS.items():
        x, top = point(name, .5, 0)
        label(ax, (x, top + 17), heading, HEADING_SIZE, bold=True, valign="top")
        subtitle_y = top + 63 + 34 * heading.count("\n")
        label(ax, (x, subtitle_y), subtitle, BODY_SIZE)

    for y, (factor, levels) in zip((.42, .80), TREATMENTS):
        label(ax, point("treatments", .67, y), factor, BODY_SIZE + 1, bold=True)
        label(ax, point("treatments", .67, y + .11), levels)
    label(ax, point("treatments", .67, .63), "×", 21, bold=True)

    for index, trait in enumerate(TRAITS):
        y = .41 + .54 * (index + .5) / len(TRAITS)
        label(ax, point("traits", .56, y), "•  " + trait, 14, align="left")

    for name, lines in [("home", HOME_MEASURES), ("away", AWAY_EQUATIONS)]:
        for index, line in enumerate(lines):
            y = .58 + .37 * (index + .5) / len(lines)
            label(ax, point(name, .5, y), line)

    label(ax, point("legend", .045, .17), LEGEND_TITLE, 14, bold=True, align="left")
    for index, (words, color, dashed) in enumerate(LEGEND):
        y = .29 + .68 * (index + .5) / len(LEGEND)
        arrow(ax, point("legend", .05, y), point("legend", .25, y),
              color, dashed=dashed, legend=True)
        label(ax, point("legend", .29, y), words, 11, align="left")

    label(ax, point("takeaway", .5, .5), fill(TAKEAWAY, width=TAKEAWAY_WRAP),
          15, bold=True, spacing=1.42)


# 5. Organism drawings and thermometer ----------------------------------------
# Each illustration uses its own local coordinates so it can be moved as a unit.
def paramecium(ax, center, scale=1):
    transform = Affine2D().scale(scale).rotate_deg(-24).translate(*center) + ax.transData
    edge = "#59746B"
    for width, height, face in [(212, 92, "#D8E6D3"), (198, 78, "none")]:
        ax.add_patch(Ellipse((0, 0), width, height, transform=transform,
                             facecolor=face, edgecolor=edge, linewidth=.9, zorder=5))
    for i in range(98):  # Short cilia around the membrane
        angle = tau * i / 98
        x, y = cos(angle), sin(angle)
        ax.plot([106 * x, 117 * x - 2 * y], [46 * y, 56 * y + 2 * x],
                transform=transform, color=edge, linewidth=.45, zorder=4)
    for x, y, width, height in [(-14, 0, 40, 21), (25, -12, 12, 9)]:
        ax.add_patch(Ellipse((x, y), width, height, transform=transform,
                             facecolor="#AEC4B6", edgecolor=edge, linewidth=.9, zorder=6))
    for x, y, radius in [(-67, -9, 10), (-43, 16, 7), (10, 20, 5), (49, 12, 9), (66, -12, 6)]:
        ax.add_patch(Circle((x, y), radius, transform=transform, facecolor="#EFF0D9",
                            edgecolor=edge, linewidth=.6, zorder=6))
        ax.add_patch(Circle((x - 1, y - 1), radius / 4, transform=transform,
                            color="#B5C7B1", zorder=7))
    # Cytoplasmic granules
    ax.scatter([-78, -63, -48, -33, -22, -7, 8, 17, 32, 46, 62, 75],
               [6, 18, -20, -26, 24, -21, 10, -27, 17, -9, 5, 13],
               s=1.8, color="#839B88", linewidths=0, transform=transform, zorder=6)


def copepod(ax, center, scale):
    transform = Affine2D().scale(scale).rotate_deg(-8).translate(*center) + ax.transData
    edge = "#45584C"

    def stroke(xs, ys, width=.8):
        ax.plot(xs, ys, color=edge, linewidth=width, transform=transform, zorder=5)

    def segment(x, y, width, height, face):
        ax.add_patch(Ellipse((x, y), width, height, facecolor=face, edgecolor=edge,
                             linewidth=.85, transform=transform, zorder=6))

    for side in (-1, 1):
        stroke([side * x for x in (7, 28, 53, 72)], [-54, -82, -102, -113])
        for y in (-27, -9, 9, 26):
            stroke([side * x for x in (10, 31, 52)], [y, y + 15, y + 28])
            stroke([side * 31, side * 54], [y + 15, y + 20], .55)
        stroke([3, 3 + side * 8, 3 + side * 23], [74, 91, 109])
        stroke([3 + side * 8, 3 + side * 11], [91, 113], .55)
    segment(0, -17, 44, 80, "#D6DECC")
    segment(-2, -50, 34, 25, "#E8EADB")
    for i in range(5):
        segment(3, 16 + 13 * i, 32 - 5 * i, 17, "#C8D1BE")
    ax.scatter([-8, 5, -3, 9], [-28, -34, -13, -5], s=12 * scale**2,
               color="#899C8A", linewidths=0, transform=transform, zorder=7)
    ax.add_patch(Circle((-3, -52), 2.8, transform=transform, color=edge, zorder=7))


def thermometer(ax, center):
    x, y = center
    ax.add_patch(FancyBboxPatch((x - 6, y - 32), 12, 39, boxstyle="round,pad=0,rounding_size=5",
                                facecolor="white", edgecolor=INK, linewidth=1.4, zorder=5))
    ax.add_patch(Circle((x, y + 9), 10, facecolor="white", edgecolor=INK, linewidth=1.4, zorder=5))
    ax.plot([x, x], [y - 22, y + 9], color="#B93232", linewidth=3, zorder=6)
    ax.add_patch(Circle((x, y + 9), 6, color="#B93232", zorder=6))


# 6. Build and save -----------------------------------------------------------
def make_figure():
    fig = plt.figure(figsize=(PAGE_WIDTH / 100, PAGE_HEIGHT / 100), facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, PAGE_WIDTH), ylim=(PAGE_HEIGHT, 0), aspect="equal")
    ax.set_axis_off()
    draw_boxes(ax)
    draw_connections(ax)
    draw_panel_text(ax)
    for words, y, size in zip(TITLE, (36, 81), TITLE_SIZES):
        label(ax, (PAGE_WIDTH / 2, y), words, size, bold=True)
    paramecium(ax, point("traits", .29, .66))
    copepod(ax, point("treatments", .20, .66), .60)
    copepod(ax, point("treatments", .36, .74), .21)
    thermometer(ax, point("treatments", .38, .52))
    return fig


if __name__ == "__main__":
    with plt.rc_context({"font.family": FONT, "mathtext.fontset": "dejavuserif",
                         "pdf.fonttype": 42, "svg.fonttype": "none"}):
        figure = make_figure()
        figure.savefig(OUTPUT_FILE, facecolor="white")
        plt.close(figure)
    print(f"Saved {OUTPUT_FILE.resolve()}")
