#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G three-panel extinction bifurcation diagram (prey x, predator y, and the orbit's two Lyapunov exponents) over a
mu0 sweep at fixed nu1 = 2.1, seed (x0, y0) = (0.2, 1e-3): 75,000-iterate burn-in, 400 retained
iterates per mu0. CURRENT PROTOCOL (2026-10-17): every mu0 grid is rounded to the nearest double
of its decimal value (_grid) so that this script and designG_onset_ensemble_v1.py evaluate
identical parameters; detect_onsets uses the SAME "sustained" rule as the ensemble (chaos onset =
least mu0 beyond which lambda_max > 0 at EVERY later swept point; extinction onset = least mu0
beyond which the retained predator stays below Y_EXT at every later swept point); the largest
Lyapunov exponent is the larger of two finite-time QR (Gram-Schmidt) estimates propagated from the
orthonormal frame ((1, 0), (0, 1)), as in the ensemble (at finite length the larger accumulated rate
is not in general the growth rate of the largest singular value of the whole product). Shipped outcome
(designG_extinction_bif_values.json): chaos onset = extinction onset = 3.96 on the 0.01 grid.

WHAT THE FIGURE SHOWS AND DOES NOT SHOW. A rare predator invades once the coexistence stock is
positive (mu0 > 21/11) and regulates the prey to the Jury-stable coexistence equilibrium over most
of the sampled range; before the sustained onset the classification on the 0.01 detector grid
alternates (this seed meets the extinction criterion at 3.85, inside the period-three window of the
prey, and coexists again at 3.86-3.95; on the 4e-4 display grid it is extinct already at 3.912 and
at about half of the points of [3.94, 3.96)), and from 3.96 on it is extinct at every later sampled point. The ensemble
(designG_onset_ensemble_v1.py, 24 seeds x 3 resolutions) finds the sustained onset seed-dependent
(median 3.96, range [3.94, 4.00] at 0.01), resolution-dependent within one coarse grid step
(3.96 / 3.955 / 3.958 for this seed), and sensitive to the last bits of mu0 (3.950000000000003 vs
3.9499999999999797 changed the fate of the retained iterates); 2 of 24 seeds have no sustained
onset (extinct at earlier sampled mu0, recaptured by the coexistence equilibrium at mu0 = 4);
in all 22 seeds with a sustained onset the two onsets coincide exactly at the three resolutions of
the ensemble (on the 4e-4 display grid of Figure 10 the extinction onset of this seed, 3.9596, comes
before its sustained-chaos onset, 3.9612), so the coincidence is no independent evidence that chaos
drives the extinction -- that evidence is the
on-axis exponent lambda_perp = ln(nu1/4) < 0, the controlled comparisons and the positive-area
theorem of the manuscript. No period-doubling cascade is visible while the predator persists.
Replaces the old single-panel designG_bifurcation_mu0.png (which seeded y0=0 and so could only
show the prey logistic cascade on M0={y=0}).

FOLDED (2026-07-02): the nu1-sweep diagnostics from the now-retired designG_bifurcation_v1.py
have been folded into this script, so this SINGLE script now produces ALL Design-G bifurcation
figures -- the mu0 extinction pair (designG_bifurcation_mu0.png + designG_lyapunov_mu0.png), the
nu1-sweep pair across the coexistence equilibrium's unit-circle locus nu1 = 2 mu0/(mu0-2)
(designG_bifurcation_nu1.png + designG_lyapunov_nu1.png; purely observational -- no normal-form
coefficient or criticality classification is computed or claimed, see the manuscript's
Proposition prop:ns-locus and its stated limitation), and the M0-logistic diagnostic pair
(designG_bifurcation_mu0_M0logistic.png + designG_lyapunov_mu0_M0logistic.png), plus both JSONs
(designG_extinction_bif_values.json + designG_bifurcation_values.json).

It also regenerates the companion transverse-exponent figure designG_lyapunov_mu0.png at the SAME
nu1, so the pair is coherent: lambda_perp=<ln(nu1 x)> on M0 is negative on the fixed-point branch
below mu0 = 21/11 ~ 1.909, positive from there through the periodic cascade and much of the chaotic
range, and of both signs at high mu0: negative inside the period-three window (e.g. -0.12 at
mu0 = 3.83), positive where the threshold nu1*(mu0) lies below nu1 (for example at mu0 = 3.9), and
negative again towards mu0 = 4, where it equals the closed form ln(nu1/4) ~ -0.64, the
exponent that forbids a rare predator from invading the fully chaotic prey.

Map (Design-G):  x' = mu0 x (1-x)(1-y),   y' = nu1 x y (1-y).
Stacked-panel format follows v5/mappings/mapping_V503.py::plotFiguresSameX, the same layout as the
old-paper reference Extinction_LogisticMap.png.

WHY nu1=2.1 (STRICTLY inside the OPEN chaos-driven band 2<nu1<4 -- nu1=2 is the EXCLUDED
endpoint): at nu1=2.1, mu0=4 the arithmetic-mean multiplier nu1<x> = 2.1 * 1/2 = 1.05 > 1, so
a mean-field reading would let the predator grow; the chaotic spread drives the geometric mean
nu1/4 = 0.525 < 1 (lambda_perp = ln(2.1/4) = -0.644), so the extinction is unambiguously the work
of the Jensen / chaotic-spread penalty, not a mean deficit -- genuinely chaos-driven. This seed's
sustained onset is 3.96 / 3.955 / 3.958 at the three ensemble resolutions (see the head of this
docstring; the ensemble median over its 22 seeds with an onset is 3.96); the chaotic band up to 4.0 is visible. nu1=2.1 is the
smallest clean choice: for nu1 >~ 2.3-2.5 the rare-invasion seed is captured by the (overlapping)
Jury-stable coexistence basin, so the sweep no longer shows clean extinction. Re-scan with
DESIGNG_NU1 if a different representative inside (2,4) is wanted.

numba JIT. Run in the bundle's codes directory: python designG_bifurcation_v2.py (run_designG_figs.sh does this).
2026-11-08 (v30): Figure 10 gains panel (c), the two Lyapunov exponents of the plotted orbit (lyap_spectrum_2d, 75,000-iterate burn-in, 20,000 averaged iterations); it is written only when its layout check finds no violation, otherwise the existing figure is left unchanged and the script exits 1.
"""
from __future__ import annotations
import json, os, sys
from datetime import datetime, timezone
import numpy as np
from numba import njit

_HERE = os.path.dirname(os.path.abspath(__file__)); _OUT = os.path.dirname(_HERE)
FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))


def _grid(lo, hi, step):
    """mu0 (or nu1) sweep grid rounded to the nearest double of each DECIMAL grid value, so that
    every script sweeping the same parameter evaluates identical float64 values: a difference
    of ~2e-14 in mu0 (np.arange from 1.0 versus from 3.0) changes the fate of the sweep seed
    near the collapse."""
    return np.round(np.arange(lo, hi + 1e-9, step), 10)


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

DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))

NU1 = float(os.environ.get("DESIGNG_NU1", "2.1"))   # strictly inside chaos-driven band (2,4)
X0 = float(os.environ.get("DESIGNG_X0", "0.2"))
# Predator seeded RARE (invasion from low density): it grows where lambda_perp>0 and stays extinct
# where lambda_perp<0; on M0 lambda_perp is negative below mu0 = 21/11, positive through the cascade
# and much of the chaotic range, and negative again at high mu0, periodic windows included (the companion
# M0 Lyapunov figure of the manuscript).
# A large y0 instead lands in the coexistence basin (a non-chaotic stable fixed point).
Y0 = float(os.environ.get("DESIGNG_Y0", "1e-3"))
N_BURN = 75_000    # paper-consistent burn-in
N_KEEP = 400
Y_EXT = 1e-3       # predator deemed extinct when the whole kept window stays below this
# Manuscript Figure 10 is laid out ZOOM-as-main + FULL-RANGE-as-inset: the collapse (chaos burst +
# extinction) lives only at mu0~3.94-4.0 (~2% of [1,4]), so the main panels
# zoom to [ZOOM_LO,4] and the full mu0 in [1,4] view is demoted to a locator inset.
ZOOM_LO = float(os.environ.get("DESIGNG_ZOOM_LO", "3.80"))


@njit(cache=True)
def sweep_mu0_xy(mu0_grid, nu1, x0, y0, n_burn, n_keep):
    """Post-burn-in (x,y) samples of the full 2-D map across a mu0 grid (fixed rare seed)."""
    N = mu0_grid.size
    P = np.empty(N * n_keep); X = np.empty(N * n_keep); Y = np.empty(N * n_keep)
    cnt = 0
    for i in range(N):
        mu0 = mu0_grid[i]
        x = x0; y = y0; escaped = False
        for _ in range(n_burn):
            xn = mu0 * x * (1.0 - x) * (1.0 - y)
            yn = nu1 * x * y * (1.0 - y)
            x = xn; y = yn
            if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
                escaped = True; break
        if escaped:
            continue
        for _ in range(n_keep):
            xn = mu0 * x * (1.0 - x) * (1.0 - y)
            yn = nu1 * x * y * (1.0 - y)
            x = xn; y = yn
            if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
                break
            P[cnt] = mu0; X[cnt] = x; Y[cnt] = y; cnt += 1
    return P[:cnt], X[:cnt], Y[:cnt]


@njit(cache=True)
def lyap_spectrum_2d(mu0, nu1, x0, y0, n_burn, n_iter):
    """Both Lyapunov exponents of the full 2-D orbit (Benettin two-vector Gram-Schmidt).

    Returns (lam1, lam2). lam1 comes from the first tangent vector started at (1, 0) with the
    arithmetic of the previous single-vector helper, operation for operation, so it is
    bit-identical to it on every orbit on which no norm is exactly zero; lam2 comes from a
    second vector started at (0, 1) and orthogonalised against the renormalised first vector
    after each step. On the invariant axis y = 0 the Jacobian is upper triangular: the first
    vector stays tangential and the second recovers the transverse rate, so the largest
    exponent there is max(lam1, lam2), not lam1 alone (since
    2026-10-21; oracle mu0 = 1.8, nu1 = 2.1, x0 = 4/9, y0 = 0: multipliers 0.2 and 14/15).
    Exact-zero convention, as in m0_lyapunov: a first norm exactly 0 gives lam1 = -inf and
    the run continues with the second vector alone; a second norm exactly 0 gives lam2 = -inf.
    Returns (nan, nan) if the orbit escapes the square during burn-in or iteration."""
    x = x0; y = y0
    for _ in range(n_burn):
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x = xn; y = yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            return np.nan, np.nan
    vx = 1.0; vy = 0.0; s1 = 0.0; alive1 = True   # the previous single-vector arithmetic
    ux = 0.0; uy = 1.0; s2 = 0.0; alive2 = True   # the second vector (transverse on the axis)
    for _ in range(n_iter):
        a = mu0 * (1.0 - 2.0 * x) * (1.0 - y); b = -mu0 * x * (1.0 - x)
        c = nu1 * y * (1.0 - y);               d = nu1 * x * (1.0 - 2.0 * y)
        if alive1:
            nvx = a * vx + b * vy; nvy = c * vx + d * vy
            nrm = np.sqrt(nvx * nvx + nvy * nvy)
            if nrm == 0.0:
                alive1 = False
            else:
                s1 += np.log(nrm); vx = nvx / nrm; vy = nvy / nrm
        if alive2:
            nux = a * ux + b * uy; nuy = c * ux + d * uy
            if alive1:
                dot = nux * vx + nuy * vy
                nux = nux - dot * vx; nuy = nuy - dot * vy
            nrm2 = np.sqrt(nux * nux + nuy * nuy)
            if nrm2 == 0.0:
                alive2 = False
            else:
                s2 += np.log(nrm2); ux = nux / nrm2; uy = nuy / nrm2
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x = xn; y = yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            return np.nan, np.nan
    lam1 = s1 / n_iter if alive1 else -np.inf
    lam2 = s2 / n_iter if alive2 else -np.inf
    return lam1, lam2


@njit(cache=True)
def lyap_max(mu0, nu1, x0, y0, n_burn, n_iter):
    """Largest Lyapunov exponent of the full 2-D orbit: max of the two exponents of
    lyap_spectrum_2d, which is the largest exponent also for orbits on the invariant axis y = 0
    (where a single tangent started at (1, 0) would report only the tangential rate). NaN if the
    orbit escapes the square."""
    l1, l2 = lyap_spectrum_2d(mu0, nu1, x0, y0, n_burn, n_iter)
    if np.isnan(l1) and np.isnan(l2):
        return np.nan
    return max(l1, l2)


@njit(cache=True)
def m0_lyapunov(mu0, nu1, x0, n_burn, n_iter):
    """On the extinction manifold M0={y=0} the prey is the logistic map x'=mu0 x(1-x).
    Returns (lambda_par, lambda_perp): tangential prey exponent <ln|mu0(1-2x)|> and the
    transverse predator exponent lambda_perp=<ln(nu1 x)> (the y-multiplier at y=0), both
    averaged over the full retained window of n_iter steps.

    Boundary: an orbit that reaches x = 0 or 1 during burn-in returns (nan, nan) -- there is no
    retained window to average. One that reaches it inside the retained window is followed
    exactly (1 -> 0, and 0 is fixed), so the remaining steps contribute ln(nu1 * 1) and then
    ln(nu1 * 0) = -inf to the transverse sum and ln|mu0| to the tangential sum: a vanishing
    multiplier makes lambda_perp -inf (the earlier version averaged the shortened
    prefix, e.g. +ln 1.05 for 1/2 -> 1 -> 0 at mu0 = 4, nu1 = 2.1). x outside [0, 1] returns nan."""
    x = x0
    for _ in range(n_burn):
        x = mu0 * x * (1.0 - x)
        if x <= 0.0 or x >= 1.0:
            return np.nan, np.nan
    if n_iter <= 0:
        return np.nan, np.nan
    spar = 0.0; sperp = 0.0
    hit_critical = False   # an exact zero derivative makes the tangential exponent -inf
    hit_zero = False       # a vanishing transverse multiplier (x == 0) makes lambda_perp -inf
    for _ in range(n_iter):
        dpar = abs(mu0 * (1.0 - 2.0 * x))
        if dpar > 0.0:
            spar += np.log(dpar)
        else:
            hit_critical = True   # the orbit sits exactly on the critical point x = 1/2
        if nu1 * x > 0.0:
            sperp += np.log(nu1 * x)
        else:
            hit_zero = True       # the orbit has reached the absorbing boundary point x = 0
        x = mu0 * x * (1.0 - x)
        if x < 0.0 or x > 1.0:
            return np.nan, np.nan
    # An exact critical hit (e.g. the superstable fixed point x = 1/2 at mu0 = 2, where
    # f'(1/2) = 0) is propagated as -inf, not skipped: the earlier version dropped the term
    # and reported a finite (zero) tangential exponent for that orbit.
    # Display code clips -inf for plotting; the value itself is not clipped.
    lam_par = -np.inf if hit_critical else spar / n_iter
    lam_perp = -np.inf if hit_zero else sperp / n_iter
    return lam_par, lam_perp


def detect_onsets(nu1, x0, y0):
    """chaos onset = the least mu0 beyond which lambda_max > 0 at EVERY swept point (the
    "sustained" rule of designG_onset_ensemble_v1.py, so the figure and the ensemble use
    one detector); extinction onset = the least mu0 beyond which the kept predator window
    stays below Y_EXT at every later swept point."""
    gl = _grid(1.0, 4.0, 0.01)
    lam = np.array([lyap_max(m, nu1, x0, y0, 20_000, 20_000) for m in gl])
    ymax = np.empty(gl.size)
    for i, m in enumerate(gl):
        _, _, Y = sweep_mu0_xy(np.array([m]), nu1, x0, y0, N_BURN, N_KEEP)
        ymax[i] = Y.max() if Y.size else 0.0
    chaos_onset = None
    for i in range(gl.size):
        if np.all(np.isfinite(lam[i:])) and np.all(lam[i:] > 0.0):
            chaos_onset = float(gl[i]); break
    ext_onset = None
    below = ymax < Y_EXT
    for i in range(gl.size):
        if below[i] and below[i:].all():
            ext_onset = float(gl[i]); break
    return chaos_onset, ext_onset


# --- nu1-sweep diagnostics folded in from designG_bifurcation_v1 (retired 2026-07-02); observational only ---
def designG_step(x, y, mu0, nu1):
    return mu0*x*(1.0-x)*(1.0-y), nu1*x*y*(1.0-y)
def jacobian(x, y, mu0, nu1):
    return np.array([[mu0*(1-2*x)*(1-y), -mu0*x*(1-x)],
                     [nu1*y*(1-y),        nu1*x*(1-2*y)]])
def coexistence_fixed_point(mu0, nu1):
    return mu0/(mu0+nu1), 1.0 - 1.0/mu0 - 1.0/nu1
def fp_seed(mu0, nu1):
    x, y = coexistence_fixed_point(mu0, nu1); return x+1e-3, y+1e-3

def in_box(x, y): return (0.0 <= x <= 1.0) and (0.0 <= y <= 1.0)

def lyapunov_spectrum(x0, y0, mu0, nu1, n_burn=N_BURN, n_iter=20000):
    """Benettin/QR Lyapunov spectrum [lambda1, lambda2] of the 2-D map."""
    x, y = x0, y0
    for _ in range(n_burn):
        x, y = designG_step(x, y, mu0, nu1)
        if not in_box(x, y): return np.array([np.nan, np.nan])
    Q = np.eye(2); s = np.zeros(2)
    for _ in range(n_iter):
        Z = jacobian(x, y, mu0, nu1) @ Q
        Q, R = np.linalg.qr(Z)
        d = np.sign(np.diag(R)); d[d == 0] = 1.0
        Q = Q*d; s += np.log(np.abs(np.diag(R)*d))
        x, y = designG_step(x, y, mu0, nu1)
        if not in_box(x, y): return np.array([np.nan, np.nan])
    return s/n_iter

def bifurcation_points(param_name, params, fixed, x0, y0, n_burn=N_BURN, n_keep=400):
    """Post-burn-in x samples. param_name in {'mu0','nu1'} (the other is `fixed`)."""
    P, X = [], []
    for p in params:
        mu0, nu1 = (p, fixed) if param_name == "mu0" else (fixed, p)
        x, y = x0, y0; ok = True
        for _ in range(n_burn):
            x, y = designG_step(x, y, mu0, nu1)
            if not in_box(x, y): ok = False; break
        if not ok: continue
        for _ in range(n_keep):
            x, y = designG_step(x, y, mu0, nu1)
            if not in_box(x, y): break
            P.append(p); X.append(x)
    return np.array(P), np.array(X)

def nu1_sweep(mu0, nu1_grid, n_keep=400):
    """Forward nu1-sweep at fixed mu0, each nu1 seeded from fp_seed (post-burn-in x samples)."""
    P, X = [], []
    for p in nu1_grid:
        x, y = fp_seed(mu0, p); ok = True
        for _ in range(N_BURN):
            x, y = designG_step(x, y, mu0, p)
            if not in_box(x, y): ok = False; break
        if not ok: continue
        for _ in range(n_keep):
            x, y = designG_step(x, y, mu0, p)
            if not in_box(x, y): break
            P.append(p); X.append(x)
    return np.array(P), np.array(X)

def _save_bif(ax, P, X, xlabel):
    ax.plot(P, X, "k.", ms=0.2, alpha=0.5)
    ax.set_xlabel(xlabel); ax.set_ylabel("$x$"); ax.grid(True, alpha=0.3)


# ---- Manuscript Figure 10 (2026-11-08): three panels sharing the mu0 axis -- prey x, predator y and the two
# Lyapunov exponents of the same orbit -- each with a full-range locator inset placed in room free of data, text
# and lines. Sizes are fixed in inches, and every letter is drawn large enough to print at >= MIN_PT when the
# 9-in figure is set at the manuscript's text width PRINT_W; fig10_violations measures both after drawing.
N_AVG = 20_000                     # iterations averaged for the exponents of panel (c), after the N_BURN burn-in
W, H, DPI = 9.0, 9.2, 200          # size (in) and resolution of designG_bifurcation_mu0.png
PRINT_W = 159.25 / 25.4            # the text width of the manuscript (main.pdf), 159.25 mm, in inches
R = PRINT_W / W                    # printed pt = drawn pt x R
MIN_PT = 7.0                       # the smallest printed letter size the layout check accepts
FS = dict(tick=11.0, label=13.0, title=12.0, letter=12.0, legend=10.5, lam=13.0, ins=10.5)
M_IN = 0.10                        # clearance between an inset's outer edge and data, text or frame, inches
L_IN, R_IN, T_IN, B_IN, GAP_IN = 1.10, 0.30, 0.80, 0.80, 0.14   # margins and the gap between panels, inches
RATIOS = (1.0, 0.8, 1.2)           # (b) has the most empty room; (c) needs height for an inset above lambda = 0


def orbit_spectrum(grid, nu1=NU1, x0=X0, y0=Y0, n_burn=N_BURN, n_avg=N_AVG):
    """Both Lyapunov exponents of the plotted orbit at every mu0 of grid, as (largest, smallest):
    lyap_spectrum_2d from the figure's seed after n_burn iterations, averaged over n_avg iterations.
    An orbit that escapes the square gives nan in both, kept as nan."""
    L = np.array([lyap_spectrum_2d(m, nu1, x0, y0, n_burn, n_avg) for m in grid])
    return np.max(L, axis=1), np.min(L, axis=1)


def fig10_data(P, X, Y, chaos_onset, ext_onset):
    """The arrays Figure 10 draws: the full-range sweep (P, X, Y) of main() on gf for the insets of (a) and (b),
    a finer sweep (Pz, Xz, Yz) on gz = [ZOOM_LO, 4] for the main panels, and the orbit's two exponents on both
    grids (l1f, l2f for the inset of (c); l1z, l2z for panel (c)), with the two onsets of detect_onsets."""
    gf = _grid(1.0, 4.0, 0.003)
    gz = _grid(ZOOM_LO, 4.0, 0.0004)
    Pz, Xz, Yz = sweep_mu0_xy(gz, NU1, X0, Y0, N_BURN, N_KEEP)
    l1z, l2z = orbit_spectrum(gz)
    l1f, l2f = orbit_spectrum(gf)
    return dict(P=P, X=X, Y=Y, Pz=Pz, Xz=Xz, Yz=Yz, gf=gf, gz=gz, l1f=l1f, l2f=l2f, l1z=l1z, l2z=l2z,
                chaos_onset=chaos_onset, ext_onset=ext_onset)


class Fig10NoRoom(Exception):
    """The panels of Figure 10 leave no room for the insets (raised by draw_fig10); main() reports it as a layout
    violation: Figure 10 is not written, the rest is."""


def draw_fig10(d, inset_bounds=None):
    """Draw Figure 10 from the dict of fig10_data without writing it; returns (fig, [ax_a, ax_b, ax_c],
    {"a": inset, "b": inset, "c": inset}). Each inset sits in the widest run of mu0 columns left of the onset line
    where the orbit rests at the coexistence fixed point, with a column of data on each side (the only data-free room
    in panel (a)): under the branch in (a) and (b), above the exponents and the zero line in (c). Raises Fig10NoRoom,
    leaving no figure open, when there is no such run or it is too small to hold the insets (for instance with
    DESIGNG_ZOOM_LO at or above the onset). inset_bounds = {panel: [x0, y0, w, h]} (that panel's axes fractions)
    replaces the computed place of that panel's inset. FS, R and MIN_PT are read here, when it is called."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import AutoMinorLocator
    from matplotlib.transforms import Bbox

    P, X, Y, Pz, Xz, Yz = (d[k] for k in ("P", "X", "Y", "Pz", "Xz", "Yz"))
    gf, gz, l1f, l2f, l1z, l2z = (d[k] for k in ("gf", "gz", "l1f", "l2f", "l1z", "l2z"))
    chaos_onset, ext_onset = d["chaos_onset"], d["ext_onset"]
    inset_bounds = inset_bounds or {}
    n = len(gz)
    assert len(Pz) == n * N_KEEP and np.all(Pz.reshape(n, N_KEEP) == gz[:, None]), "sweep layout is not one row per mu0"
    Xm, Ym = Xz.reshape(n, N_KEEP), Yz.reshape(n, N_KEEP)
    ext = np.all(Ym < Y_EXT, axis=1)
    calm = (~ext) & (np.ptp(Xm, axis=1) < 1e-9) & (np.ptp(Ym, axis=1) < 1e-9)   # predator alive, orbit at one point

    def clip(v):                   # display only: an exponent of -inf is drawn below the axis floor
        return np.where(np.isfinite(v), v, -9.0)

    # widest run of consecutive calm columns left of the dashed line with a column of data on each side: the only span
    # free in panel (a) below the line (a run at an end of the zoom window has no column there to bound it)
    edge = chaos_onset if chaos_onset is not None else np.inf
    best, i = (-1, 0, -1), 0
    while i < n:
        if calm[i] and gz[i] < edge - 1e-9:
            j = i
            while j + 1 < n and calm[j + 1] and gz[j + 1] < edge - 1e-9:
                j += 1
            if i > 0 and j < n - 1:
                best = max(best, (j - i, i, j))
            i = j + 1
        else:
            i += 1
    _, I0, I1 = best
    if I1 < 0:
        raise Fig10NoRoom(f"no room for the insets: no run of mu0 columns at the coexistence fixed point with a column "
                          f"of data on each side left of the onset line in the zoom window [{gz[0]:.4f}, {gz[-1]:.4f}]")
    print(f"free span: calm columns mu0 {gz[I0]:.4f}..{gz[I1]:.4f} ({I1 - I0 + 1} columns); "
          f"neighbours {gz[I0 - 1]:.4f} and {gz[I1 + 1]:.4f}")

    # fixed layout in inches
    AX_W = W - L_IN - R_IN
    AX_HS = [(H - T_IN - B_IN - 2 * GAP_IN) * q / sum(RATIOS) for q in RATIOS]
    fig = plt.figure(figsize=(W, H), dpi=DPI)
    ax = []
    for k in range(3):
        y0_in = B_IN + sum(AX_HS[j] + GAP_IN for j in range(k + 1, 3))
        rect = (L_IN / W, y0_in / H, AX_W / W, AX_HS[k] / H)
        ax.append(fig.add_axes(rect, sharex=ax[0] if ax else None))

    ax[0].scatter(Pz, Xz, color="k", s=0.2, marker=".", rasterized=True)
    ax[0].set_ylabel(r"$x$  (prey)", size=FS["label"])
    ax[0].set_ylim(0, 1.02)
    ax[1].scatter(Pz, Yz, color="k", s=0.2, marker=".", rasterized=True)
    ax[1].set_ylabel(r"$y$  (predator)", size=FS["label"])
    ax[2].plot(gz, clip(l1z), color="tab:blue", lw=0.8)
    ax[2].plot(gz, clip(l2z), color="tab:red", lw=0.8)
    ax[2].axhline(0.0, color="k", lw=0.6)
    ax[2].text(3.975, 0.80, r"$\lambda_1$", color="tab:blue", size=FS["lam"])
    ax[2].text(3.975, -0.66, r"$\lambda_2$", color="tab:red", size=FS["lam"])
    ax[2].set_ylabel(r"$\lambda$  (this orbit)", size=FS["label"])
    ax[2].set_xlabel(r"$\mu_0$", size=FS["label"])
    ax[2].set_ylim(-1.25, 1.2)
    ax[2].set_yticks([-1.0, -0.5, 0.0, 0.5])
    ax[2].set_xlim(ZOOM_LO, 4.0)
    # mark the sustained prey-chaos onset and the predator-extinction onset of THIS seed.
    # Both are read off the same event (the seed's exit from the coexistence basin as mu0
    # grows), so when they coincide one line is drawn; the label states the seed dependence
    # and asserts no theorem (the onset ensemble of Section 6 and Section S12 gives the spread over seeds).
    for a in ax:
        if chaos_onset is not None:
            a.axvline(chaos_onset, color="r", ls="--", lw=1.1, alpha=0.8)
        a.grid(True, alpha=0.3)
        a.yaxis.set_minor_locator(AutoMinorLocator(5))
        a.xaxis.set_minor_locator(AutoMinorLocator(5))
        a.tick_params(labelsize=FS["tick"])
    for a in ax[:2]:
        a.tick_params(labelbottom=False)
    for a, t in zip(ax, "abc"):
        a.text(-0.085, 1.0, f"({t})", transform=a.transAxes, size=FS["letter"], fontweight="bold", va="top")
    if chaos_onset is not None:
        same = ext_onset is not None and abs(ext_onset - chaos_onset) < 1e-9
        if same:
            lbl = r"sustained prey-chaos onset $=$ predator-extinction onset, this seed ($\mu_0\!\approx\!%.2f$)" % chaos_onset
        else:
            lbl = r"sustained prey-chaos onset, this seed ($\mu_0\!\approx\!%.2f$)" % chaos_onset
        handles = [Line2D([0], [0], color="r", ls="--", lw=1.1, label=lbl)]
        if ext_onset is not None and not same:
            for a in ax:
                a.axvline(ext_onset, color="b", ls=":", lw=1.1, alpha=0.8)
            handles.append(Line2D([0], [0], color="b", ls=":", lw=1.1,
                                  label=r"predator-extinction onset, this seed ($\mu_0\!\approx\!%.2f$)" % ext_onset))
        ax[0].legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.02), borderaxespad=0.0,
                     frameon=False, fontsize=FS["legend"])
    fig.text(0.5, 1 - 0.22 / H, rf"Predator extinction from prey chaos ($\nu_1={NU1:g}$)", ha="center", va="center",
             size=FS["title"])

    def draw_inset(parent, kind, bounds_axes):     # full-range locator inset (mu0 in [1,4]); zoom window shaded red
        ins = parent.inset_axes(bounds_axes)
        if kind == "a":
            ins.scatter(P, X, color="k", s=0.06, marker=".", rasterized=True)
            ins.set_ylim(0, 1.02)
            ins.set_yticks([0, 0.5, 1.0])
        elif kind == "b":
            ins.scatter(P, Y, color="k", s=0.06, marker=".", rasterized=True)
            ins.set_ylim(0, 0.30)
            ins.set_yticks([0, 0.1, 0.2, 0.3])
        else:
            ins.plot(gf, clip(l1f), color="tab:blue", lw=0.6)
            ins.plot(gf, clip(l2f), color="tab:red", lw=0.6)
            ins.axhline(0.0, color="k", lw=0.4)
            ins.set_ylim(-2.5, 1.0)
            ins.set_yticks([-2, -1, 0, 1])
        ins.axvspan(ZOOM_LO, 4.0, color="tab:red", alpha=0.15)
        if chaos_onset is not None:
            ins.axvline(chaos_onset, color="r", ls="--", lw=0.7, alpha=0.7)
        ins.set_xlim(1.0, 4.0)
        ins.set_xticks([1, 2, 3, 4])
        ins.tick_params(labelsize=FS["ins"], length=2.5, pad=1.5, width=0.6)
        ins.set_title(r"full range $\mu_0\in[1,4]$", size=FS["ins"], pad=2)
        for s in ins.spines.values():
            s.set_linewidth(0.6)
        return ins

    def to_axes(parent, box):      # display-pixel Bbox -> parent axes-fraction bounds
        (x0, y0), (x1, y1) = parent.transAxes.inverted().transform([(box.x0, box.y0), (box.x1, box.y1)])
        return [x0, y0, x1 - x0, y1 - y0]

    px = DPI                       # pixels per inch
    m = M_IN * px
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    xl = ax[0].transData.transform((gz[I0 - 1], 0))[0] + m          # free span, display x
    xr = ax[0].transData.transform((gz[I1 + 1], 0))[0] - m

    def free_rect(k):
        a = ax[k]
        if k == 0:
            top = a.transData.transform((0, Xm[I0:I1 + 1].min()))[1] - m
            return Bbox([[xl, a.bbox.y0 + m], [xr, top]])
        if k == 1:
            top = a.transData.transform((0, Ym[I0:I1 + 1].min()))[1] - m
            return Bbox([[xl, a.bbox.y0 + m], [xr, top]])
        lo = a.transData.transform((0, max(0.0, np.max(clip(l1z[I0:I1 + 1])))))[1] + m   # above the zero line and the span's exponents
        return Bbox([[xl, lo], [xr, a.bbox.y1 - m]])

    def no_room(size):             # the free span cannot hold the insets: close the figure, report it
        plt.close(fig)
        return Fig10NoRoom(f"no room for the insets: the calm columns mu0 {gz[I0]:.4f}..{gz[I1]:.4f} leave {size}")

    # pass 1: a provisional box per panel, to measure how far tick labels and title reach beyond the box
    pads = []
    for k, kind in enumerate("abc"):
        R0 = free_rect(k)
        box = Bbox([[R0.x0 + 0.45 * px, R0.y0 + 0.30 * px], [R0.x1, R0.y1 - 0.30 * px]])
        if box.width <= 0 or box.height <= 0:
            raise no_room(f"a {box.width / px:.2f} in x {box.height / px:.2f} in box in panel ({kind})")
        ins = draw_inset(ax[k], kind, to_axes(ax[k], box))
        fig.canvas.draw()
        T = ins.get_tightbbox(r)
        pads.append((box.x0 - T.x0, box.y0 - T.y0, T.x1 - box.x1, T.y1 - box.y1))
        ins.remove()
    pl = max(p[0] for p in pads) + 1
    pb = max(p[1] for p in pads) + 1
    pr = max(p[2] for p in pads) + 1
    pt = max(p[3] for p in pads) + 1
    box_w = (xr - xl) - pl - pr
    h_a = free_rect(0).height - pb - pt
    h_b = free_rect(1).height - pb - pt
    box_h = min(h_a, h_b)
    if box_w <= 0 or box_h <= 0:
        raise no_room(f"an inset box of {box_w / px:.2f} in x {box_h / px:.2f} in once its labels are placed")
    # raise the top of panel (c) until its free rectangle holds a box of the same height
    need = box_h + pb + pt + m
    for _ in range(30):
        Rc = free_rect(2)
        if Rc.height >= need - 0.5:
            break
        lo_data, hi_data = ax[2].get_ylim()
        per_px = (hi_data - lo_data) / ax[2].bbox.height
        ax[2].set_ylim(lo_data, hi_data + (need - Rc.height) * per_px * 1.05)
        fig.canvas.draw()
    print(f"panel (c) y range {ax[2].get_ylim()[0]:.2f}..{ax[2].get_ylim()[1]:.2f}; inset box "
          f"{box_w / px:.2f} x {box_h / px:.2f} in drawn = {box_w / px * R:.2f} x {box_h / px * R:.2f} in printed")

    insets = {}
    for k, kind in enumerate("abc"):
        Rk = free_rect(k)
        if k < 2:   # just under the line or branch the inset must not touch
            box = Bbox([[xl + pl, Rk.y1 - pt - box_h], [xl + pl + box_w, Rk.y1 - pt]])
        else:       # just above the zero line
            box = Bbox([[xl + pl, Rk.y0 + pb], [xl + pl + box_w, Rk.y0 + pb + box_h]])
        insets[kind] = draw_inset(ax[k], kind, inset_bounds[kind] if kind in inset_bounds else to_axes(ax[k], box))
    fig.canvas.draw()
    return fig, ax, insets


def fig10_violations(fig, axes, insets):
    """Measure the drawn Figure 10 and return one string per violation (an empty list: the layout passes). Checked:
    each inset (box, tick labels, title) inside its panel and clear of data, lines and text; no overlap among the
    insets, the panel texts, tick labels and axis labels, the legend and the title, nor among an inset's own title
    and tick labels; nothing outside the canvas; the data-placed labels of (c) inside the panel and off the curves;
    the legend above panel (a) and off the data; every letter printing at >= MIN_PT at the width PRINT_W. Prints
    one line per check, and one per violation found."""
    from matplotlib.text import Text
    r = fig.canvas.get_renderer()
    ax = list(axes)
    ins_all = [insets[k] for k in "abc"]

    def seg_hits(xy, B):
        """True when a segment of the display polyline xy meets the Bbox B (Liang-Barsky)."""
        for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:]):
            dx, dy, t0, t1, ok = x1 - x0, y1 - y0, 0.0, 1.0, True
            for p, q in ((-dx, x0 - B.x0), (dx, B.x1 - x0), (-dy, y0 - B.y0), (dy, B.y1 - y0)):
                if p == 0:
                    if q < 0:
                        ok = False
                        break
                else:
                    t = q / p
                    if p < 0:
                        t0 = max(t0, t)
                    else:
                        t1 = min(t1, t)
                    if t0 > t1:
                        ok = False
                        break
            if ok:
                return True
        return False

    def data_hits(a, B):
        hits = []
        for coll in a.collections:
            pts = coll.get_offset_transform().transform(coll.get_offsets())
            k = int(np.sum((pts[:, 0] >= B.x0) & (pts[:, 0] <= B.x1) & (pts[:, 1] >= B.y0) & (pts[:, 1] <= B.y1)))
            if k:
                hits.append(f"{k} data points")
        for ln in a.lines:
            xy = ln.get_transform().transform(ln.get_path().vertices)
            if seg_hits(xy, B):
                hits.append(f"line {ln.get_label() or ln.get_color()}")
        return hits

    def visible_ticklabels(a):
        out = []
        lo, hi = sorted(a.get_xlim())
        out += [t for t in a.get_xticklabels() if t.get_visible() and t.get_text() and lo <= t.get_position()[0] <= hi]
        lo, hi = sorted(a.get_ylim())
        out += [t for t in a.get_yticklabels() if t.get_visible() and t.get_text() and lo <= t.get_position()[1] <= hi]
        return out

    bad = []
    items = []    # (name, Bbox) of every text-like thing outside the insets, and of each inset as one block
    for k, (a, ins) in enumerate(zip(ax, ins_all)):
        p = "abc"[k]
        T = ins.get_tightbbox(r).padded(1)
        inner = a.bbox.padded(-2)
        inside = T.x0 >= inner.x0 and T.x1 <= inner.x1 and T.y0 >= inner.y0 and T.y1 <= inner.y1
        hits = data_hits(a, T) + [f"text {t.get_text()!r}" for t in a.texts if t.get_window_extent(r).overlaps(T)]
        print(f"inset ({p}): inside the panel frame {inside}; covers {hits or 'nothing'}")
        bad += ([] if inside else [f"inset ({p}): outside the panel frame"]) + [f"inset ({p}): covers {h}" for h in hits]
        items.append((f"inset ({p})", T))
        items += [(f"({p}) text {t.get_text()!r}", t.get_window_extent(r)) for t in a.texts]
        items += [(f"({p}) tick {t.get_text()!r}", t.get_window_extent(r)) for t in visible_ticklabels(a)]
        items += [(f"({p}) axis label", lab.get_window_extent(r)) for lab in (a.xaxis.label, a.yaxis.label)
                  if lab.get_text()]
    leg = ax[0].get_legend()
    if leg is not None:
        items.append(("legend", leg.get_window_extent(r)))
    items += [(f"figure text {t.get_text()[:20]!r}", t.get_window_extent(r)) for t in fig.texts]
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i][1].overlaps(items[j][1]):
                msg = f"OVERLAP {items[i][0]} <-> {items[j][0]}"
                print(msg)
                bad.append(msg)
    for k, ins in enumerate(ins_all):          # inside each inset: its title against its own tick labels
        own = [("title", ins.title.get_window_extent(r))] + \
              [(f"tick {t.get_text()!r}", t.get_window_extent(r)) for t in visible_ticklabels(ins)]
        for i in range(len(own)):
            for j in range(i + 1, len(own)):
                if own[i][1].overlaps(own[j][1]):
                    msg = f"OVERLAP inside inset ({'abc'[k]}): {own[i][0]} <-> {own[j][0]}"
                    print(msg)
                    bad.append(msg)
    for name, B in items:                      # nothing may leave the canvas (a cut-off tick label is a defect)
        if not (B.x0 >= fig.bbox.x0 and B.x1 <= fig.bbox.x1 and B.y0 >= fig.bbox.y0 and B.y1 <= fig.bbox.y1):
            msg = f"OUTSIDE THE CANVAS {name}"
            print(msg)
            bad.append(msg)
    for k, a in enumerate(ax):                 # data-placed labels stay inside their panel
        for t in a.texts:
            if t.get_transform() is a.transData and not a.bbox.contains(*t.get_window_extent(r).p0) \
                    or t.get_transform() is a.transData and not a.bbox.contains(*t.get_window_extent(r).p1):
                msg = f"LABEL OUTSIDE PANEL ({'abc'[k]}) {t.get_text()!r}"
                print(msg)
                bad.append(msg)
    if leg is not None:
        legB = leg.get_window_extent(r)
        leg_in_axes = legB.overlaps(ax[0].bbox.padded(-1))
        lh = data_hits(ax[0], legB)
        print(f"legend: inside panel (a) {leg_in_axes}; covers {lh or 'nothing'}")
        bad += (["legend: inside panel (a)"] if leg_in_axes else []) + [f"legend: covers {h}" for h in lh]
    for t in ax[2].texts:                      # the direct labels must not sit on the curves
        h = data_hits(ax[2], t.get_window_extent(r).padded(1))
        print(f"label {t.get_text()!r}: covers {h or 'nothing'}")
        bad += [f"label {t.get_text()!r}: covers {x}" for x in h]
    small = sorted({(round(t.get_fontsize() * R, 2), t.get_text()[:28]) for t in fig.findobj(Text)
                    if t.get_visible() and t.get_text().strip() and t.get_fontsize() * R < MIN_PT - 1e-9})
    sizes = sorted({round(t.get_fontsize() * R, 2) for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()})
    print(f"printed sizes present (pt): {sizes}; under {MIN_PT}: {small or 'none'}")
    bad += [f"printed size {pt} pt under {MIN_PT}: {txt!r}" for pt, txt in small]
    return bad


def save_fig10(fig, path, violations):
    """Write Figure 10 to path only when violations is empty, and say which happened; returns True when written.
    With violations, the file at path (the figure already in the manuscript, when path is in figures/) is left
    unchanged."""
    if violations:
        print(f"Figure 10 layout check: {len(violations)} violation(s)")
        for v in violations:
            print(f"  - {v}")
        print(f"Figure 10 NOT written: {path} left unchanged")
        return False
    print("Figure 10 layout check: 0 violations")
    fig.savefig(path, dpi=DPI)
    _manifest_append(path)
    print(f"written {path}")
    return True


def main():
    os.makedirs(FIGDIR, exist_ok=True); os.makedirs(DATADIR, exist_ok=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    mu0_grid = _grid(1.0, 4.0, 0.003)
    P, X, Y = sweep_mu0_xy(mu0_grid, NU1, X0, Y0, N_BURN, N_KEEP)
    chaos_onset, ext_onset = detect_onsets(NU1, X0, Y0)
    print(f"chaos onset mu0 ~ {chaos_onset}; predator extinction onset mu0 ~ {ext_onset}")

    # ---- Manuscript Figure 10: three panels (prey, predator, this orbit's exponents), ZOOM-as-main + FULL-RANGE-as-inset.
    #      Main panels zoom to mu0 in [ZOOM_LO,4] (where the chaos burst + this seed's
    #      predator extinction live -- only ~2% of [1,4]); a locator inset in each panel
    #      shows the full mu0 in [1,4] range with the zoom window shaded. Presentation
    #      only: same nu1, same data as the full-range sweep above (P,X,Y feeds the insets of (a)
    #      and (b); a finer zoom sweep feeds the main panels; panel (c) and its inset are the two
    #      Lyapunov exponents of the same orbits). The figure is written only when its layout
    #      check finds no violation (no room for the insets is one); otherwise the file already in
    #      FIGDIR is left unchanged, the rest is written as usual and main() returns 1. ----
    d = fig10_data(P, X, Y, chaos_onset, ext_onset)
    path10 = os.path.join(FIGDIR, "designG_bifurcation_mu0.png")
    try:
        fig, axes, insets = draw_fig10(d)
    except Fig10NoRoom as e:
        fig10_written = save_fig10(None, path10, [str(e)])
    else:
        fig10_written = save_fig10(fig, path10, fig10_violations(fig, axes, insets))
        plt.close(fig)

    # ---- Supplementary Figure S3: companion M0 Lyapunov exponents (lambda_par, lambda_perp) at the SAME nu1 ----
    gl = _grid(1.0, 4.0, 0.01)
    lams = np.array([m0_lyapunov(m, NU1, 0.2, 20_000, 40_000) for m in gl])
    # Display only: a -inf tangential exponent (exact critical hit) is drawn at the
    # axis floor; the JSON/returned values keep -inf. With x0 = 0.2 on this grid no exact
    # hit occurs, so the plotted curve is unchanged by the clip.
    lam_par_disp = np.where(np.isfinite(lams[:, 0]), lams[:, 0], -1.5)
    figL, axL = plt.subplots(figsize=(9, 3.6))
    axL.plot(gl, lam_par_disp, "b-", lw=0.9, label=r"$\lambda_\parallel$ (prey on $\mathcal{M}_0$)")
    axL.plot(gl, lams[:, 1], "r-", lw=0.9, label=r"$\lambda_\perp=\langle\ln(\nu_1 x)\rangle$ (transverse / predator)")
    axL.axhline(0.0, color="k", lw=0.5)
    axL.axhline(np.log(NU1 / 4.0), color="r", ls=":", lw=0.8,
                label=r"$\ln(\nu_1/4)=%.3f$ at $\mu_0{=}4$" % np.log(NU1 / 4.0))
    axL.set_xlabel(r"$\mu_0$", size=12); axL.set_ylabel(r"$\lambda$", size=12)
    axL.set_ylim(-1.5, 1.0); axL.grid(True, alpha=0.3); axL.legend(fontsize=9, loc="upper left")   # no curve passes there
    axL.set_title(rf"Exponents on $\mathcal{{M}}_0$ at $\nu_1={NU1:g}$: prey $\lambda_\parallel$ and the on-axis rare-invader rate $\lambda_\perp$ "
                  rf"($\to\ln(\nu_1/4)$ at $\mu_0{{=}}4$)" + "" + r"",
                  size=10, wrap=True)
    figL.tight_layout()
    _p = os.path.join(FIGDIR, "designG_lyapunov_mu0.png"); figL.savefig(_p, dpi=200, bbox_inches="tight", pad_inches=0.02); _manifest_append(_p)  # 2026-10-17: trim the white margin
    plt.close(figL)

    # ---- Slice 1: mu0-sweep at nu1=2.5 (bounded box) ----
    mu0_grid = _grid(1.0, 4.0, 0.005)
    P1, X1 = bifurcation_points("mu0", mu0_grid, 2.5, 0.2, 0.0)
    mu0_lyap = _grid(1.0, 4.0, 0.025)   # coarser grid for the costlier QR Lyapunov
    lam1 = np.array([lyapunov_spectrum(0.2, 0.0, m, 2.5, n_iter=8000) for m in mu0_lyap])
    fig, ax = plt.subplots(figsize=(9, 4)); _save_bif(ax, P1, X1, r"$\mu_0$")
    ax.set_title(r"Prey on the extinction manifold $M_0$ ($y{=}0$): logistic period-doubling route to chaos"); fig.tight_layout()
    # Distinct name: the canonical designG_bifurcation_mu0.png is the three-panel Figure 10;
    # this v1 slice is only the M0-logistic diagnostic (kept, not overwriting v2).
    _p = os.path.join(FIGDIR, "designG_bifurcation_mu0_M0logistic.png"); fig.savefig(_p, dpi=200); _manifest_append(_p); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(mu0_lyap, lam1[:, 0], "b-", lw=0.8, label=r"$\lambda_\parallel$ (prey on $M_0$)")
    ax.plot(mu0_lyap, lam1[:, 1], "r-", lw=0.8, label=r"$\lambda_\perp$ (transverse / predator)")
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel(r"$\mu_0$"); ax.set_ylabel(r"$\lambda$")
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_title(r"On $M_0$ ($\nu_1=2.5$): prey chaos ($\lambda_\parallel>0$); transverse $\lambda_\perp=\langle\ln(\nu_1 x)\rangle$; closed form $\ln(\nu_1/4)<0$ at $\mu_0{=}4$")
    ax.set_ylim(-1.5, 1.0)
    fig.tight_layout(); _p = os.path.join(FIGDIR, "designG_lyapunov_mu0_M0logistic.png"); fig.savefig(_p, dpi=200); _manifest_append(_p); plt.close(fig)
    # ---- Slice 2: nu1-sweeps across the unit-circle locus nu1 = 2 mu0/(mu0-2) at mu0 = 3 (locus 6)
    # and mu0 = 3.5 (locus 14/3). The labels record what the sweeps SHOW (a hard onset with 1:5
    # locking at mu0 = 3; a soft onset into a quasiperiodic circle at mu0 = 3.5); they are not a
    # criticality classification -- no normal-form coefficient is computed here.
    cases = [(3.0, 6.0, "observed: hard onset, 1:5 locking", (4.5, 6.6)),
             (3.5, 14.0/3.0, "observed: soft onset, quasiperiodic circle", (4.0, 5.8))]
    figB, axB = plt.subplots(1, 2, figsize=(13, 4.2))
    figL, axL = plt.subplots(1, 2, figsize=(13, 4.2))
    win = {}
    for k, (mu0, nu1ns, label, (lo, hi)) in enumerate(cases):
        P, X = nu1_sweep(mu0, _grid(lo, hi, 0.01))
        win[str(mu0)] = [float(P.min()) if P.size else None, float(P.max()) if P.size else None]
        _save_bif(axB[k], P, X, r"$\nu_1$")
        axB[k].axvline(nu1ns, color="b", ls="--", lw=1, alpha=0.7)
        axB[k].set_title(f"$\\mu_0={mu0}$ ({label})")
        gl = _grid(lo, hi, 0.025)
        lam = np.array([lyapunov_spectrum(*fp_seed(mu0, p), mu0, p, n_iter=8000) for p in gl])
        axL[k].plot(gl, lam[:, 0], "b-", lw=0.8, label=r"$\lambda_1$")
        axL[k].plot(gl, lam[:, 1], "r-", lw=0.8, label=r"$\lambda_2$")
        axL[k].axhline(0, color="k", lw=0.5); axL[k].axvline(nu1ns, color="b", ls="--", lw=1, alpha=0.7)
        axL[k].set_xlabel(r"$\nu_1$"); axL[k].set_ylabel(r"$\lambda$"); axL[k].legend(); axL[k].grid(True, alpha=0.3)
        axL[k].set_title(f"$\\mu_0={mu0}$ ({label})")
    figB.suptitle(r"Design-G $\nu_1$-sweeps across the unit-circle locus $\nu_1=2\mu_0/(\mu_0-2)$ at $\mu_0=3$ and $3.5$ (observed behaviour)")
    figB.tight_layout(); _p = os.path.join(FIGDIR, "designG_bifurcation_nu1.png"); figB.savefig(_p, dpi=200); _manifest_append(_p); plt.close(figB)
    figL.suptitle(r"Lyapunov spectrum: hard/locked ($\mu_0=3$) vs $\lambda_{\max}\approx0$ quasiperiodic ($\mu_0=3.5$)")
    figL.tight_layout(); _p = os.path.join(FIGDIR, "designG_lyapunov_nu1.png"); figL.savefig(_p, dpi=200); _manifest_append(_p); plt.close(figL)
    json.dump({"utc": datetime.now(timezone.utc).isoformat(),
               "nu1_window_mu0_3": win.get("3.0"), "nu1_window_mu0_3.5": win.get("3.5")},
              open(os.path.join(DATADIR, "designG_bifurcation_values.json"), "w"), indent=2)

    # The transverse exponent is negative inside periodic windows
    # of the high-mu0 region too, not "only in the chaotic band"; the values the manuscript
    # quotes at mu0 = 3.83 and 3.84 (the period-three window) are recorded from the figure's
    # own axis protocol (x0 = 0.2, 20,000 burn-in, 40,000 retained) so they are value-checked.
    window_check = {}
    for mu0_w in (3.83, 3.84):
        lp, lt = m0_lyapunov(round(mu0_w, 10), NU1, 0.2, 20_000, 40_000)
        window_check[f"{mu0_w:.2f}"] = {"lambda_par": float(lp), "lambda_perp": float(lt),
                                        "prey_periodic": bool(lp < 0.0), "transverse_negative": bool(lt < 0.0)}
    json.dump({"utc": datetime.now(timezone.utc).isoformat(), "nu1": NU1,
               "seed": [X0, Y0], "n_burn": N_BURN, "n_keep": N_KEEP,
               "chaos_onset_mu0": chaos_onset, "extinction_onset_mu0": ext_onset,
               "lambda_perp_mu0_4_closed_form": float(np.log(NU1 / 4.0)),
               "periodic_window_check": window_check},
              open(os.path.join(DATADIR, "designG_extinction_bif_values.json"), "w"), indent=2)
    if fig10_written:
        print("OK three-panel extinction bifurcation (Figure 10) + companion M0 Lyapunov (Figure S3) written")
    else:
        print("Figure 10 NOT written (layout violations above); companion M0 Lyapunov (Figure S3) written")
    return 0 if fig10_written else 1


if __name__ == "__main__":
    sys.exit(main())
