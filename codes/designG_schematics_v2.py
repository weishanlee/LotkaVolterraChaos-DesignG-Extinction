#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproducible Design-G schematics (manuscript Figures 1, 3 and 9), designed at 6.5 inches wide.

architecture_methods_v5.png: shared multiplicative coupling template and the
forward-invariant state square, logistic prey edge, arcsine density, transverse
multiplier/exponent and computed coexistence point. Sources: eq:designG,
eq:template, eq:bounds, eq:logistic, eq:srb-density-mu4, eq:linear-transverse,
eq:lambda-mu4, eq:coexistence-fp and eq:designG-jury-band.

mechanism_methods_v5.png: exact logarithmic growth curves above five aligned,
independent property strips. Sources: eq:arcsine-moments, eq:jensen-gap,
eq:lambda-mu4, thm:ergodic-opt, thm:global (including eq:equality-case),
thm:basin-positive, prop:property-table and tab:properties. The worst exponent
is zero at 4/3: strict invariant-measure contraction excludes this endpoint,
whereas the nonlinear all-orbit extinction theorem includes it.

workflow_methods_v5.png: branching analytical dependencies, the constructive
positive-area/rare-predator branch, the ergodic-optimization/all-orbit branch,
and subordinate causality reductions and numerical checks. Sources:
sec:model, thm:transverse-criterion, eq:arcsine-moments,
eq:log-sine, eq:jensen-gap, eq:lambda-mu4,
thm:basin-positive, cor:rare-predator, thm:ergodic-opt, thm:global,
prop:property-table, tab:properties, thm:mod1, thm:skew and Section S12 of
the Supplementary Material. Arrows describe analytical operations, not ecological
forward causation. The dashed numerical node has no outgoing arrow.

Run in the bundle's codes directory: python designG_schematics_v2.py
(run_designG_figs.sh does this). Only NumPy and Matplotlib are third-party
dependencies. The manuscript sets the three figures at EMBED_LINEWIDTH of its
column (0.85, 0.731 and 0.714 of the 6.27-inch A4 text width, i.e. 5.33, 4.58
and 4.48 in); the mechanism and workflow figures therefore draw their type at
TYPE_SCALE = 1.15. The os/hashlib standard
library use preserves the v1 output and manifest contract. No data files,
external fonts, TeX executable, random sampling or network access are used.
Every render checks all visible nonempty Text extents for canvas containment,
pairwise non-overlap, a drawn size of at least 8 pt, and a printed size of at
least MIN_PRINTED_PT at the manuscript's width (the smallest printed type is
about 7.0, 6.9 and 6.7 pt in Figures 1, 9 and 3); every arrow label is also
checked against the node frames and against the arrows.
Analytic curve samples, density moments, fixed-point equations, Jacobian
eigenvalues, threshold alignment and interpolation error are checked before
output is saved. Finite viewing windows are disclosed in the artwork: the arcsine
density has its own calibrated axis and divergent endpoints; all logarithmic
curves diverge to minus infinity as nu1 approaches zero from above.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.path import Path
from matplotlib.text import Text

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT  = os.path.dirname(_HERE)
FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))


def _manifest_append(path):
    """Append a producer-authored ``| artifact | version | sha256 |`` row to DESIGNG_MANIFEST (if set)."""
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    import hashlib
    with open(path, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {os.path.basename(__file__)} | {digest} |\n")

os.makedirs(FIGDIR, exist_ok=True)

DPI = 220
WIDTH = 6.5
INK = "#20252a"
GRAY = "#61676b"
GRID = "#d6dade"
BLUE = "#0072b2"
ORANGE = "#b65b00"
GREEN = "#007f63"
PURPLE = "#85577c"
JURY = float(8 * np.sqrt(3) - 12)
# Type scale applied by text(): the mechanism and workflow figures, which the
# manuscript sets narrower than the architecture figure, draw their type larger.
TYPE_SCALE = 1.0
# Widths at which the manuscript includes the figures, as fractions of its text
# width COLUMN_IN (A4 with 1-inch margins), and the smallest printed type allowed.
COLUMN_IN = 6.27
EMBED_LINEWIDTH = {"architecture_methods_v5": 0.85, "mechanism_methods_v5": 0.731,
                   "workflow_methods_v5": 0.714}
MIN_PRINTED_PT = 6.5

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "mathtext.fontset": "dejavusans", "text.usetex": False,
    "text.color": INK, "axes.labelcolor": INK,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt",
    "path.simplify": False,
})


def density_samples():
    """Exact density evaluations; finite samples never include singular ends."""
    x = np.unique(np.r_[np.linspace(1e-6, 1-1e-6, 5000), .25, .5, .75])
    return x, 1 / (np.pi*np.sqrt(x*(1-x)))


def growth_samples():
    """Analytic functions, with all labelled thresholds sampled exactly."""
    nu = np.unique(np.r_[np.geomspace(1e-8, 4, 2400), 4/3, JURY, 2, 4])
    return nu, np.array([np.log(3*nu/4), np.log(nu/2), np.log(nu/4)])


def _interpolation_error(x, y, function, pixels_per_unit, visible_limits):
    """Sample nine interior points per segment; report visible error in pixels."""
    fractions = np.arange(1, 10)[:, None] / 10
    xx = x[:-1] + fractions * np.diff(x)
    exact = function(xx)
    interpolated = y[:-1] + fractions * np.diff(y)
    visible = (exact >= visible_limits[0]) & (exact <= visible_limits[1])
    return float(np.max(np.abs(interpolated[visible]-exact[visible])) * pixels_per_unit)


def audit_analytic_content():
    """Independent identities and substitutions, not simulation-based proof."""
    nu, growth = growth_samples()
    assert np.all(np.diff(growth, axis=1) > 0)
    assert np.all(growth[0] > growth[1]) and np.all(growth[1] > growth[2])
    gap_error = float(np.max(np.abs(growth[1]-growth[2]-np.log(2))))
    assert gap_error < 1e-13
    for row, crossing in enumerate([4/3, 2, 4]):
        assert abs(growth[row, np.flatnonzero(nu == crossing)[0]]) < 1e-14
    log_errors = [
        _interpolation_error(nu, values, function, DPI*2.12/4.45, (-3.2, 1.25))
        for values, function in zip(growth, [lambda v: np.log(3*v/4),
                                             lambda v: np.log(v/2),
                                             lambda v: np.log(v/4)])
    ]
    assert max(log_errors) < .1

    xx, rho = density_samples()
    assert np.all(rho > 0)
    assert abs(rho[np.flatnonzero(xx == .5)[0]] - 2/np.pi) < 1e-14
    density_error = _interpolation_error(
        xx, rho, lambda x: 1/(np.pi*np.sqrt(x*(1-x))), .125*DPI, (0, 5))
    assert density_error < .1
    # x=sin(theta)^2 transforms the arcsine measure to uniform theta.
    # Check the Jacobian-weighted normalisation separately, then the moments.
    moments = []
    for n in [16384, 32768]:
        theta = (np.arange(n)+.5)*np.pi/(2*n)
        x = np.sin(theta)**2
        weight = np.sin(2*theta)/(np.pi*np.sqrt(x*(1-x)))
        normalisation = float(np.mean(weight)*np.pi/2)
        assert abs(normalisation-1) < 1e-8
        moments.append((float(np.mean(x)), float(np.mean(np.log(x)))))
    log_moment = 2*moments[1][1]-moments[0][1]
    assert abs(moments[1][0]-.5) < 1e-13
    assert abs(log_moment+2*np.log(2)) < 1e-10

    # Differentiate the two update rules and compare against the manuscript's
    # independently expressed trace/determinant; also check the strict band.
    values = np.unique(np.r_[np.linspace(4/3+1e-5, 4-1e-5, 1001),
                             JURY-1e-6, JURY+1e-6, 2.1])
    x = 4/(4+values)
    y = (3*values-4)/(4*values)
    a = 4*(1-2*x)*(1-y)
    b = -4*x*(1-x)
    c = values*y*(1-y)
    d = values*x*(1-2*y)
    trace, det = a+d, a*d-b*c
    assert np.allclose(trace, -(values-4)**2/(values*(values+4)), rtol=0, atol=1e-13)
    assert np.allclose(det, (values**2+12*values-32)/(values*(values+4)), rtol=0, atol=1e-13)
    disc = np.sqrt((trace**2-4*det).astype(complex))
    radius = np.maximum(np.abs((trace+disc)/2), np.abs((trace-disc)/2))
    assert np.array_equal(radius < 1, values > JURY)
    residual = max(float(np.max(np.abs(4*x*(1-x)*(1-y)-x))),
                   float(np.max(np.abs(values*x*y*(1-y)-y))))
    assert residual < 1e-13
    for value, expected in [(JURY, [-1, 1/np.sqrt(3)]), (4, [-1j, 1j])]:
        xp, yp = 4/(4+value), (3*value-4)/(4*value)
        jacobian = np.array([[4*(1-2*xp)*(1-yp), -4*xp*(1-xp)],
                             [value*yp*(1-yp), value*xp*(1-2*yp)]])
        actual = np.linalg.eigvals(jacobian)
        assert all(np.min(np.abs(actual-target)) < 1e-12 for target in expected)
    print(f"Math QA: Jensen-gap error {gap_error:.3g}; "
          f"fixed-point residual {residual:.3g}; Jury endpoint and band checks passed")
    print(f"Math QA: density normalisation {normalisation:.12f}; "
          f"E[x]={moments[1][0]:.12f}; E[ln x]={log_moment:.12f} "
          "(transformed quadrature with endpoint-error extrapolation)")
    print(f"Math QA: sampled interpolation errors: density {density_error:.4f} px; "
          f"log curves {max(log_errors):.4f} px maximum")


def canvas(height):
    fig = plt.figure(figsize=(WIDTH, height), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, WIDTH), ylim=(0, height))
    ax.axis("off")
    # The transparent extent anchor makes bbox_inches='tight' retain the exact
    # designed canvas, so 8.5 pt in the drawing is 8.5 pt at width=6.5 inches.
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fc="none", ec="none", lw=0))
    fig._node_checks = []
    fig._arrow_checks = []
    fig._alignment_checks = []
    return fig, ax


def text(ax, x, y, value, size=9, color=INK, ha="left", weight="normal",
         va="center", rotation=0, background=False):
    return ax.text(x, y, value, fontsize=size*TYPE_SCALE, color=color, ha=ha, va=va,
                   weight=weight, rotation=rotation, linespacing=1.35,
                   bbox=({"facecolor": "white", "edgecolor": "none",
                          "pad": 1.5} if background else None), zorder=8)


def line(ax, xs, ys, color=GRAY, lw=1, style="-", zorder=2):
    return ax.plot(xs, ys, color=color, lw=lw, ls=style,
                   clip_on=False, zorder=zorder)[0]


def arrow(ax, points, label, at, color=INK, dashed=False, size=8.5,
          ha="center", rotation=0, double=False):
    """Every directed graphical edge has an explicit, separately placed label."""
    path = Path(points, [Path.MOVETO] + [Path.LINETO] * (len(points) - 1))
    patch = FancyArrowPatch(path=path, arrowstyle="<|-|>" if double else "-|>",
                            mutation_scale=10, lw=1.1, color=color,
                            linestyle="--" if dashed else "-", zorder=3)
    ax.add_patch(patch)
    artist = text(ax, *at, label, size=size, color=color, ha=ha,
                  rotation=rotation, background=True)
    ax.figure._arrow_checks.append((ax, [tuple(p) for p in points], artist))
    return artist


def node(ax, x, y, w, h, lines, color=INK, dashed=False):
    """Rectangular nodes on an inch-based grid; line tuples are text/size/weight."""
    patch = Rectangle((x-w/2, y-h/2), w, h, facecolor="white",
                      edgecolor=color, linewidth=1.15,
                      linestyle=(0, (4, 3)) if dashed else "-", zorder=4)
    ax.add_patch(patch)
    step = min(.28, (h - .20) / max(len(lines)-1, 1))
    ys = [y + (len(lines)-1)*step/2 - i*step for i in range(len(lines))]
    artists = [text(ax, x, yy, s, size=fs, weight=wt, color=color,
                    ha="center") for yy, (s, fs, wt) in zip(ys, lines)]
    ax.figure._node_checks.append((patch, artists))
    return patch


def _segments_hit_box(points, box, pad):
    """True if the polyline through `points` (display coordinates) enters `box` grown by `pad`."""
    x0, y0, x1, y1 = box.x0 - pad, box.y0 - pad, box.x1 + pad, box.y1 + pad
    for (ax_, ay_), (bx_, by_) in zip(points[:-1], points[1:]):
        for t in np.linspace(0.0, 1.0, 201):
            px, py = ax_ + t * (bx_ - ax_), ay_ + t * (by_ - ay_)
            if x0 <= px <= x1 and y0 <= py <= y1:
                return True
    return False


def _check_and_save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for data_axes, data_points, overlay_axes, overlay_points in fig._alignment_checks:
        actual = data_axes.transData.transform(data_points)
        intended = overlay_axes.transData.transform(overlay_points)
        assert np.max(np.abs(actual-intended)) < 1e-7, "Misaligned quantitative axes"
    artists = [t for t in fig.findobj(match=Text)
               if t.get_visible() and t.get_text().strip()]
    extents = [(t, t.get_window_extent(renderer)) for t in artists]
    problems = []
    for t, b in extents:
        if b.x0 < 0 or b.y0 < 0 or b.x1 > fig.bbox.x1 or b.y1 > fig.bbox.y1:
            problems.append("outside canvas: " + repr(t.get_text()))
        if t.get_fontsize() < 8:
            problems.append("under 8 pt: " + repr(t.get_text()))
    for i, (a, ba) in enumerate(extents):
        for b, bb in extents[i+1:]:
            if ba.overlaps(bb):
                problems.append("text overlap: " + repr(a.get_text()) +
                                " / " + repr(b.get_text()))
    for patch, labels in fig._node_checks:
        outer = patch.get_window_extent(renderer)
        for t in labels:
            b = t.get_window_extent(renderer)
            if (b.x0 < outer.x0 or b.x1 > outer.x1 or
                    b.y0 < outer.y0 or b.y1 > outer.y1):
                problems.append("outside node: " + repr(t.get_text()))
    frames = [patch.get_window_extent(renderer) for patch, _ in fig._node_checks]
    for _, _, label in fig._arrow_checks:
        b = label.get_window_extent(renderer)
        for frame in frames:   # an arrow label lies between nodes: it may not touch a node frame
            if b.overlaps(frame):
                problems.append("arrow label over a node frame: " + repr(label.get_text()))
        for arrow_axes, points, _ in fig._arrow_checks:   # nor cover an arrow, its own included
            pts = arrow_axes.transData.transform(np.asarray(points, dtype=float))
            if _segments_hit_box(pts, b, 1.0):
                problems.append("arrow label over an arrow: " + repr(label.get_text()))
    if problems:
        raise RuntimeError(stem + " layout check failed:\n" + "\n".join(sorted(set(problems))))
    out = os.path.join(FIGDIR, stem + ".png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0)
    _manifest_append(out)
    # PNG IHDR dimensions; no image-processing dependency is imported.
    with open(out, "rb") as fh:
        header = fh.read(24)
    pixels = (int.from_bytes(header[16:20], "big"),
              int.from_bytes(header[20:24], "big"))
    printed_in = EMBED_LINEWIDTH[stem] * COLUMN_IN
    effective_min = min(t.get_fontsize() for t in artists) * printed_in * DPI / pixels[0]
    assert effective_min >= MIN_PRINTED_PT, "Printed type is below MIN_PRINTED_PT at the manuscript's width"
    print(f"QA: {stem}: {len(artists)} text artists; 0 overlaps; "
          f"0 out of bounds; min printed type {effective_min:.2f} pt at "
          f"{EMBED_LINEWIDTH[stem]:.3f} x {COLUMN_IN} in; {pixels[0]} x {pixels[1]} px")
    plt.close(fig)
    print(f"Saved: {out}")


def make_architecture():
    fig, ax = canvas(6.15)
    text(ax, .22, 5.93, "(a)  Shared template, opposite partner factors", 10.5,
         weight="bold")
    text(ax, 3.25, 5.61,
         r"rate $\times\ z(1-z)\ \times$ partner factor", 10, ha="center")

    node(ax, 1.36, 4.86, 2.20, 1.02,
         [(r"Prey  $x$", 10.5, "bold"),
          (r"$x'=\mu_0\,x(1-x)(1-y)$", 10, "normal"),
          (r"self-limitation: $x(1-x)$", 8.5, "normal")], BLUE)
    node(ax, 5.14, 4.86, 2.20, 1.02,
         [(r"Predator  $y$", 10.5, "bold"),
          (r"$y'=\nu_1\,x\,y(1-y)$", 10, "normal"),
          (r"self-limitation: $y(1-y)$", 8.5, "normal")], ORANGE)
    arrow(ax, [(2.50, 5.02), (4.00, 5.02)],
          "feeds predators\n" + r"partner $x$: increasing in $x$", (3.25, 5.23), BLUE)
    arrow(ax, [(4.00, 4.49), (2.50, 4.49)],
          "harms prey\n" + r"partner $(1-y)$:" + "\ndecreasing in $y$", (3.25, 4.75), ORANGE)
    line(ax, [.22, 6.28], [4.13, 4.13], GRID, .8)
    text(ax, .22, 3.91, "(b)  State space and invasion geometry", 10.5,
         weight="bold")

    # A metric square: both physical dimensions are 2.30 inches.
    sx, sy, side = .48, 1.27, 2.30
    ax.add_patch(Rectangle((sx, sy), side, side, fill=False, ec=GRAY, lw=1.1))
    line(ax, [sx, sx+side], [sy, sy], BLUE, 2.7)
    text(ax, sx+side/2, 3.33, "Forward invariant", 9, ha="center")
    text(ax, sx+side/2, 3.06, r"$[0,1]^2,\quad 0<\mu_0,\nu_1\leq4$", 9,
         ha="center")
    text(ax, .29, sy+side, r"$1$", 8.5, ha="center")
    text(ax, sx, sy-.17, r"$0$", 8.5, ha="center")
    text(ax, sx+side, sy-.17, r"$1$", 8.5, ha="center")
    text(ax, .30, sy+side/2, r"$y$", 10, ha="center")
    text(ax, 2.94, sy, r"$x$", 10, ha="center")
    text(ax, sx+side/2, sy-.17, r"$\mathcal{M}_0=\{y=0\}$", 9,
         color=BLUE, ha="center")
    arrow(ax, [(.86, sy+.02), (.86, 2.58)], "transverse",
          (.72, 2.02), color=ORANGE, rotation=90)

    nu_plot = 2.1
    xp = 4 / (4 + nu_plot)
    yp = (3*nu_plot - 4) / (4*nu_plot)
    # Check the displayed point is an actual fixed point of the plotted map.
    assert np.allclose([4*xp*(1-xp)*(1-yp), nu_plot*xp*yp*(1-yp)],
                       [xp, yp], rtol=0, atol=1e-13)
    ax.plot(sx+side*xp, sy+side*yp, "o", ms=5, color=PURPLE, zorder=6)
    text(ax, sx+side*xp+.14, sy+side*yp+.16, r"$(x^*,y^*)$", 9, color=PURPLE)
    text(ax, 1.83, 2.53, rf"$\mu_0=4,\ \nu_1={nu_plot:.1f}$", 8.5,
         ha="center", color=PURPLE)

    # Quantitatively calibrated auxiliary density axis, separate from state y.
    # Plot the analytic curve through the frame; no artificial finite endcaps.
    xx, rho = density_samples()
    base, density_scale = .35, .125
    density = fig.add_axes([sx/WIDTH, base/6.15, side/WIDTH,
                            (5*density_scale)/6.15])
    density.set_zorder(-1)
    density.set(xlim=(0, 1), ylim=(0, 5))
    density.axis("off")
    density.plot(xx, rho, color=BLUE, lw=1.65)
    calibration = np.array([[0, 0], [.5, 2], [1, 4]])
    fig._alignment_checks.append((density, calibration, ax,
                                  np.column_stack((sx+side*calibration[:, 0],
                                                   base+density_scale*calibration[:, 1]))))
    line(ax, [sx, sx+side], [base, base], GRID, .7)
    line(ax, [sx, sx], [base, base+5*density_scale], GRAY, .7)
    for value in [0, 2, 4]:
        yy = base+density_scale*value
        line(ax, [sx-.035, sx], [yy, yy], GRAY, .7)
        text(ax, sx-.075, yy, str(value), 8.5, ha="right")
    text(ax, sx+side/2, .18, "Arcsine density: separate vertical scale", 8.5, ha="center")
    text(ax, .19, .69, r"$\rho_{\mathrm{SRB}}$", 9, ha="center", color=BLUE)

    text(ax, 3.15, 3.53, "Coexistence at " + r"$\mu_0=4$", 9.5, weight="bold")
    text(ax, 3.15, 3.17,
         r"$(x^*,y^*)=\left(\frac{4}{4+\nu_1},\frac{3\nu_1-4}{4\nu_1}\right)$", 10)
    text(ax, 3.15, 2.85, r"Interior iff $\nu_1>4/3$", 9)
    text(ax, 3.15, 2.57, r"Jury-stable iff $8\sqrt{3}-12<\nu_1<4$", 8.5)
    text(ax, 3.15, 2.20,
         r"$\left.\frac{\partial G}{\partial y}\right|_{y=0}=\nu_1x$", 10,
         color=ORANGE)
    text(ax, 3.15, 1.78,
         r"$\lambda_\perp=\int\ln(\nu_1x)\,d\mu_4=\ln(\nu_1/4)$", 9,
         color=ORANGE)
    text(ax, 3.15, 1.38, "Exact prey-edge dynamics", 9.5, weight="bold")
    text(ax, 3.15, 1.08, r"$f_{\mu_0}(x)=\mu_0x(1-x)$", 10, color=BLUE)
    text(ax, 3.15, .77, r"At $\mu_0=4$: arcsine SRB measure $\mu_4$", 8.5)
    text(ax, 3.15, .38, r"$\rho_{\mathrm{SRB}}(x)=\frac{1}{\pi\sqrt{x(1-x)}}$,  $0<x<1$", 10,
         color=BLUE)
    text(ax, 3.15, .13, r"$\rho_{\mathrm{SRB}}(x)\to\infty$ as $x\to0^+$ or $1^-$", 8.5, color=BLUE)
    print(f"Computed coexistence example: nu1={nu_plot}, x*={xp:.10f}, y*={yp:.10f}")
    _check_and_save(fig, "architecture_methods_v5")


def make_mechanism():
    fig, ax = canvas(7.80)
    text(ax, .20, 7.58, "(a)  Analytic growth curves", 10.5, weight="bold")
    text(ax, 6.24, 7.58, r"$\mu_0=4$", 10, ha="right")
    legends = [
        (7.24, INK, "-", r"Worst invariant measure $\delta_{3/4}$:  exponent $\ln(3\nu_1/4)$"),
        (6.97, ORANGE, "--", r"Arithmetic mean $\nu_1\mathbb{E}_{\mu_4}[x]=\nu_1/2$:  $\ln(\nu_1/2)$"),
        (6.70, BLUE, "-.", r"Geometric mean under $\mu_4$ $=\nu_1/4$:  $\lambda_\perp=\ln(\nu_1/4)$"),
    ]
    for yy, col, ls, label in legends:
        line(ax, [.29, .72], [yy, yy], col, 1.6, ls)
        text(ax, .84, yy, label, 9)

    x0, x1 = 1.78, 5.99
    plot_y, plot_h = 4.29, 2.12
    g = fig.add_axes([x0/WIDTH, plot_y/7.80, (x1-x0)/WIDTH, plot_h/7.80])
    # Curve axes sit behind the labelled inch-based overlay, including the
    # Jensen-gap arrow and the explanation inside the highlighted interval.
    g.set_zorder(-1)
    g.set(xlim=(0, 4), ylim=(-3.2, 1.25))
    g.axis("off")
    g.axvspan(0, 4/3, fc="#f0f2f3", ec="none", zorder=0)
    g.axvspan(2, 4, fc="#fbf0df", ec="none", zorder=0)
    for threshold in [4/3, JURY, 2, 4]:
        g.axvline(threshold, color=GRID, lw=.7, ls=":" if threshold == JURY else "--")
    g.axhline(0, color=GRAY, lw=.9)
    nu, growth = growth_samples()
    for yy, col, ls in zip(growth, [INK, ORANGE, BLUE], ["-", "--", "-."]):
        g.plot(nu, yy, color=col, lw=1.7, ls=ls)
    for xv, col, marker in [(4/3, INK, "^"), (2, ORANGE, "s"), (4, BLUE, "o")]:
        g.plot(xv, 0, marker=marker, ms=5, mfc="white", mec=col,
               mew=1.2, clip_on=False, zorder=6)
    def X(v):
        return x0 + (x1-x0)*v/4
    def Y(v):
        return plot_y + plot_h*(v+3.2)/4.45
    calibration = np.array([[4/3, 0], [JURY, 0], [2, 0], [4, 0]])
    fig._alignment_checks.append((g, calibration, ax,
                                  np.column_stack((X(calibration[:, 0]),
                                                   Y(calibration[:, 1])))))
    line(ax, [x0, x0], [plot_y, plot_y+plot_h], GRAY, .9)
    line(ax, [x0, x1], [plot_y, plot_y], GRAY, .9)
    for val in [-3, -2, -1, 0, 1]:
        line(ax, [x0-.04, x0], [Y(val), Y(val)], GRAY, .8)
        text(ax, x0-.12, Y(val), rf"${val}$", 8.5, ha="right")
    text(ax, .24, 6.18, "Log growth\nper step", 8.5, weight="bold")
    text(ax, .24, 5.50, "Mean-driven:\nevery invariant\nmeasure\ncontracts",
         8.5)
    text(ax, .24, 4.90, r"$0<\nu_1<4/3$", 9)
    text(ax, .24, 4.49,
         r"As $\nu_1\to0^+$," + "\n" + r"all curves $\to-\infty$" + "\n(tails clipped)",
         8.5)
    text(ax, X(3), Y(-2.08), "Fluctuation penalty decides", 8.5,
         ha="center", color=ORANGE)
    text(ax, X(3), Y(-2.65), r"$2<\nu_1<4$", 9, ha="center", color=ORANGE)
    n_gap = 3.15
    arrow(ax, [(X(n_gap), Y(np.log(n_gap/4))),
               (X(n_gap), Y(np.log(n_gap/2)))],
          r"$\ln2$", (X(n_gap)+.18, Y(np.log(n_gap/4)+np.log(2)/2)),
          color=ORANGE, double=True, ha="left")
    for val, label in [(0, "0"), (4/3, r"4/3"), (2, "2"), (4, "4")]:
        text(ax, X(val), 4.11, "$"+label+"$", 8.5, ha="center")
    text(ax, 6.20, 4.11, r"$\nu_1$", 10, ha="center")
    text(ax, .24, 3.82,
         r"At $4/3$: worst exponent $=0$; self-limitation still proves extinction.", 8.5)
    text(ax, .24, 3.56,
         r"Jury endpoints: at $\nu_1^J$, multiplier $-1$; at $4$, $\det J=1$, multipliers $\pm i$.", 8.5)

    text(ax, .20, 3.23, "(b)  Five separate properties", 10.5, weight="bold")
    text(ax, 6.24, 3.23, rf"$\nu_1^J=8\sqrt{{3}}-12\approx {JURY:.3f}$", 8.5, ha="right")
    rows = [2.80, 2.25, 1.70, 1.15, .60]
    for val in [0, 4/3, JURY, 2, 4]:
        line(ax, [X(val), X(val)], [.47, 2.98], GRID, .65,
             ":" if val == JURY else "--")

    def endpoint(v, yy, color, closed):
        ax.plot(X(v), yy, "o", ms=4.4, mew=1.05, mec=color,
                mfc=color if closed else "white", zorder=6)
    def band(lo, hi, yy, color, label, closed=(False, False), style="-", label_y=.13, ends=(True, True)):
        line(ax, [X(lo), X(hi)], [yy, yy], color, 2.25, style, 4)
        if ends[0]:
            endpoint(lo, yy, color, closed[0])
        if ends[1]:
            endpoint(hi, yy, color, closed[1])
        text(ax, X((lo+hi)/2), yy+label_y, label, 8.5, color=color, ha="center")
    labels = ["Interior fixed point\nexists", "Jury-stable\nfixed point",
              r"Sign of $\lambda_\perp$", "All-orbit\nextinction",
              "Extinction basin\n" + r"$\operatorname{Leb}_2(W^s)>0$"]
    for yy, label in zip(rows, labels):
        text(ax, .20, yy, label, 8.5)
    yy = rows[0]
    band(0, 4/3, yy, GRAY, "no", style="--")
    band(4/3, 4, yy, PURPLE, "yes", closed=(False, True))
    text(ax, X(4/3), yy-.17, "on the edge", 8.5, ha="center")
    yy = rows[1]
    band(0, 4/3, yy, GRAY, "—", closed=(False, True), style=":")
    band(4/3, JURY, yy, GRAY, "no", style="--", ends=(False, True))   # 4/3 belongs to the band on its left
    band(JURY, 4, yy, PURPLE, "yes")
    text(ax, X(JURY), yy-.17, "neutral", 8.5, ha="center")
    text(ax, X(4), yy-.17, "neutral", 8.5, ha="center")
    yy = rows[2]
    band(0, 4, yy, BLUE, r"$<0$")
    text(ax, X(4), yy+.14, r"$=0$", 8.5, ha="center")
    yy = rows[3]
    band(0, 4/3, yy, GREEN, "PROVED", closed=(False, True))
    line(ax, [X(4/3)+.04, X(4)], [yy, yy], GRAY, 1.2, "--")
    endpoint(4, yy, GRAY, True)
    text(ax, X(8/3), yy+.14, "impossible: interior fixed orbit", 8.5, ha="center")
    yy = rows[4]
    band(0, 4/3, yy, GREEN, "full area", closed=(False, True))
    line(ax, [X(4/3)+.04, X(4)], [yy, yy], BLUE, 2.25)
    endpoint(4, yy, BLUE, False)
    text(ax, X(2.62), yy+.14, "PROVED positive area", 8.5, ha="center", color=BLUE)
    text(ax, X(4), yy+.14, "unresolved", 8.5, ha="center")
    # A shared metric x scale; the nearby Jury label is staggered below 2.
    line(ax, [x0, x1], [.32, .32], GRAY, .9)
    for val, label in [(0, "0"), (4/3, "4/3"), (2, "2"), (4, "4")]:
        line(ax, [X(val), X(val)], [.29, .35], GRAY, .8)
        text(ax, X(val), .17, "$"+label+"$", 8.5, ha="center")
    text(ax, X(JURY)-.10, .15, r"$\nu_1^J$", 8.5, ha="right")
    line(ax, [X(JURY), X(JURY)-.08], [.31, .20], GRAY, .7)
    text(ax, 6.22, .17, r"$\nu_1$", 10, ha="center")
    print(f"Computed Jury threshold: {JURY:.12f}")
    _check_and_save(fig, "mechanism_methods_v5")


def make_workflow():
    fig, ax = canvas(8.00)
    text(ax, .22, 7.79, "Analytical workflow", 11, weight="bold")
    text(ax, 6.24, 7.79, "Arrows denote analytical steps", 8.5, ha="right", color=GRAY)
    node(ax, 3.25, 7.33, 3.40, .50,
         [(r"The map at $\mu_0=4$", 10, "bold")], INK)
    node(ax, 1.79, 6.32, 2.96, .84,
         [(r"Prey chaos on $\mathcal{M}_0$", 9.5, "bold"),
          (r"$T=f_4$; arcsine SRB $\mu_4$", 9, "normal"),
          (r"$\mathbb{E}_{\mu_4}[\ln x]=-2\ln2$", 9, "normal")], BLUE)
    node(ax, 4.83, 6.32, 2.72, .84,
         [("Log-sine integral", 8.5, "bold"),
          (r"$2\int_0^1\ln|\sin\pi u|\,du=-\ln4$", 9, "normal"),
          (r"$x=\sin^2(\pi\theta/2)$ gives the moment", 8.5, "normal")], BLUE)
    arrow(ax, [(2.20, 7.06), (1.79, 6.77)],
          r"restrict to $y=0$", (1.30, 6.94), BLUE)
    arrow(ax, [(4.29, 7.06), (4.83, 6.77)],
          "evaluate the moment", (4.74, 6.97), BLUE, ha="left")

    node(ax, 1.79, 5.24, 2.96, .60,
         [(r"Transverse multiplier $\nu_1x$", 9.5, "bold"),
          ("Linearised criterion", 9, "normal")], ORANGE)
    node(ax, 4.83, 5.24, 2.72, .66,
         [("Jensen comparison", 9, "bold"),
          (r"$\lambda_\perp-\ln(\nu_1/2)=-\ln2$", 9, "normal")], ORANGE)
    arrow(ax, [(1.79, 5.87), (1.79, 5.57)], "linearise", (1.28, 5.72), ORANGE)
    arrow(ax, [(4.83, 5.87), (4.83, 5.60)], "use the moments", (4.95, 5.74), ORANGE, ha="left")

    node(ax, 3.25, 4.24, 5.30, .62,
         [("Closed-form transverse exponent", 9.5, "bold"),
          (r"$\lambda_\perp=\ln\nu_1+\mathbb{E}_{\mu_4}[\ln x]=\ln(\nu_1/4)$", 10, "normal")], ORANGE)
    arrow(ax, [(1.79, 4.91), (1.79, 4.58)],
          r"average under $\mu_4$", (1.10, 4.75), ORANGE)
    arrow(ax, [(4.83, 4.58), (4.83, 4.88)],
          "quantify the penalty", (4.95, 4.73), ORANGE, ha="left")

    node(ax, 1.79, 3.18, 2.96, .70,
         [("Positive-area extinction basin", 9, "bold"),
          (r"$\operatorname{Leb}_2(W^s)>0$,  $0<\nu_1<4$", 9, "normal")], GREEN)
    node(ax, 4.83, 3.18, 2.72, .70,
         [("Ergodic optimization", 9.5, "bold"),
          (r"$\sup_m\int\ln x\,dm=\ln(3/4)$", 9, "normal"),
          (r"attained at $\delta_{3/4}$", 8.5, "normal")], PURPLE)
    arrow(ax, [(1.79, 3.90), (1.79, 3.56)],
          "construct the basin", (1.04, 3.73), GREEN)
    arrow(ax, [(4.83, 3.90), (4.83, 3.56)],
          "vary invariant law", (4.95, 3.73), PURPLE, ha="left")

    node(ax, 1.79, 2.18, 2.96, .64,
         [(r"Rare-predator consequence on $Z_+$", 8.5, "bold"),
          ("Prey chaotic itinerary retained", 8.5, "normal")], GREEN)
    node(ax, 4.83, 2.18, 2.72, .64,
         [("All-orbit extinction", 9.5, "bold"),
          (r"iff $0<\nu_1\leq4/3$", 10, "normal")], PURPLE)
    arrow(ax, [(1.79, 2.80), (1.79, 2.53)],
          "track the prey orbit", (1.06, 2.67), GREEN)
    arrow(ax, [(4.83, 2.80), (4.83, 2.53)],
          "sub-action", (4.71, 2.665), PURPLE, ha="right")
    text(ax, 4.95, 2.665, "+ self-limitation", 8.5, color=PURPLE, ha="left")

    node(ax, 3.25, 1.22, 3.45, .57,
         [(r"Properties by $\nu_1$", 10, "bold"),
          ("Separately established properties", 8.5, "normal")], INK)
    arrow(ax, [(1.79, 1.83), (2.36, 1.54)], "tabulate", (1.47, 1.63), GREEN)
    arrow(ax, [(4.83, 1.83), (4.16, 1.54)], "tabulate", (5.19, 1.63), PURPLE)

    # Subordinate, explicitly separate nodes: no transfer of the Design-G
    # basin result to reductions and no numerical arrow into an analytic node.
    text(ax, 2.17, .81, r"Causality reductions at $\mu_0=4$", 8.5,
         color=GRAY, ha="center")
    node(ax, 1.15, .40, 1.88, .60,
         [("Specialist escapement", 8.5, "bold"),
          (r"$\lambda_\perp=\ln(\nu_1/4)$", 8.5, "normal")], GRAY, True)
    node(ax, 3.25, .40, 1.94, .60,
         [("One-way skew product", 8.5, "bold"),
          (r"$\lambda_\perp=\ln(\nu_1/4)$", 8.5, "normal")], GRAY, True)
    node(ax, 5.37, .40, 1.88, .60,
         [("Numerical checks", 8.5, "bold"),
          ("of the theorems (S12)", 8.5, "normal")], GRAY, True)
    _check_and_save(fig, "workflow_methods_v5")


if __name__ == "__main__":
    print(f"Writing figures to: {FIGDIR}")
    audit_analytic_content()
    TYPE_SCALE = 1.0; make_architecture()
    TYPE_SCALE = 1.15; make_mechanism()
    TYPE_SCALE = 1.15; make_workflow()
    print("All done.")
