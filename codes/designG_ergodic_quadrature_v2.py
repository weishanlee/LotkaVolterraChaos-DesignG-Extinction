#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G ergodic / quadrature / enumeration validation suite (mu0=4 layer).

Numerical CONSISTENCY CHECKS (never proofs) for the closed-form / ergodic-
optimization content of ``output/methods_v6.tex`` attached to the symmetric
multiplicative ("Design G") map.  On the prey-only manifold M0 the dynamics is
the logistic map ``x_{n+1} = mu0 x(1-x)``; at ``mu0 = 4`` its physical measure
is the arcsine (Chebyshev) measure ``dmu4 = dx / (pi sqrt(x(1-x)))``.

This is one of two CPU-only scripts.  It covers spec section 2B:

  (a) arcsine moments via Gauss--Chebyshev quadrature  -- eq:arcsine-moments, eq:jensen-gap
  (b) equilibrium-measure identity int ln|x-r| dmu4    -- thm:equilibrium-measure, eq:equilibrium-potential
  (c) ergodic-optimization sub-action H(x) <= ln(3/4)  -- thm:ergodic-opt, eq:subaction
      + exact-angle enumeration of EVERY interior periodic orbit through period 16
        (fractions.Fraction) with float64 ln x orbit averages -- rmk:ergodic-numeric
  (d) I(mu0) family + threshold nu1*(mu0)=exp(-I)      -- cor:mu0family, rmk:Icurve

The companion script ``designG_dynamics_v4.py`` covers the 2-D map
simulation layer (spec 2A).

Key closed-form targets (mu0 = 4 specific):
    E_mu4[ln x] = -2 ln2 = -1.386294 ;  E_mu4[x] = 1/2
    Jensen gap  = E[ln x] - ln E[x] = -ln2 = -0.693147
    int ln|x-r| dmu4 = -ln4 = -1.386294  for r in [0,1]  (~ -0.519 @ r=1.2, ~ +0.376 @ r=2.0)
    sup over invariant measures of int ln x = ln(3/4) = -0.287682, attained at x*=3/4
    H(x) = ln x + (4/9)(3x - 4x^2) <= ln(3/4), equality at 3/4, H''(x) < 0 on (0,1]
    I(4) = -2 ln2 -> nu1*(4) = 4.000 ; all-orbit uniform bound mu0/(mu0-1), min 4/3 at mu0=4.

Honesty notes:
  * Every number is a finite-grid / finite-precision consistency check, NOT a
    proof.  The sub-action telescoping argument (thm:ergodic-opt) is the proof;
    the period-<=16 enumeration only corroborates it.
  * a.e./SRB statements are NOT all-orbit statements.
  * The closed forms (arcsine density, I(4) = -2 ln2, nu1*(4) = 4) are SPECIFIC
    to mu0 = 4; for mu0 < 4 the SRB measure is not arcsine.  I(mu0) is still
    known in closed form where the prey's attractor is the fixed point
    (1 < mu0 <= 3: I = ln(1 - 1/mu0)) or the two-cycle (3 < mu0 <= 1 + sqrt 6:
    I = ln((mu0 + 1)/mu0^2)/2); nu1_star_closed below uses both.  At the sampled
    parameters beyond that range (3.60 to 3.99: chaotic, or in a periodic
    window) no closed form is used, so nu1*(mu0) = exp(-I) is a per-parameter
    orbit average there, not monotone in mu0 (the sampled values fall from 3.60
    to 3.70, both chaotic; periodic windows interrupt the chaotic set), not a
    smooth global curve.

Engineering follows the canonical siblings ``codes/forward_causal_*.py``:
portable ``_HERE/_OUT`` paths, ``matplotlib.use("Agg")``, Numba ``@njit`` hot
loops, provenance + values JSON.  CPU only; run in the bundle's codes directory
(run_designG_figs.sh does this; README.md gives the commands)::

    python designG_ergodic_quadrature_v2.py

The PNGs (manuscript Figures 7, 8 and S1) are written unless ``--no-figs`` is supplied.
"""
# === Revision v2 (2026-06-04) =================================================
# Changes:
#   the r=0.5 singular-quadrature comment, provenance fields,
#   the run command, a --max-period domain guard, a capacity-sign typo,
#   an off-support legend entry.
# No algorithm / theorem / equation / numerical result changed.
#
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
from fractions import Fraction
from datetime import datetime, timezone
from typing import Dict, List, Tuple

import numpy as np
try:
    import numba
    from numba import njit
except ModuleNotFoundError:  # numba is optional: keep the module import-safe in a numpy-only env
    import types as _types
    def njit(*args, **kwargs):
        # usable as both @njit and @njit(...): no-op (pure-Python) fallback when numba is absent
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda f: f
    numba = _types.SimpleNamespace(__version__="not installed (pure-Python fallback)")

# -----------------------------------------------------------------------------
# Portable paths (mirror the canonical forward_causal_*.py idiom).
# -----------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))   # .../output/codes

def _rel(p):
    """Record paths RELATIVE to the bundle/authoring root (the parent of codes/), never an
    absolute author-local path; non-string or empty values pass through."""
    if not isinstance(p, str) or not p:
        return p
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        return os.path.relpath(os.path.abspath(p), root)
    except ValueError:
        return os.path.basename(p)

_OUT = os.path.dirname(_HERE)                         # .../output
DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))
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


LN2 = math.log(2.0)
LN_THREE_QUARTER = math.log(0.75)   # = -0.287682072...
MU_INF = 3.5699456718               # Feigenbaum accumulation point (matches sibling scripts)


# =============================================================================
# Closed forms of the parameter box 0 < mu0, nu1 <= 4 (manuscript Figure 8)
# =============================================================================
#
# The map is (x, y) -> (mu0 x(1-x)(1-y), nu1 x y(1-y)). Nothing here is computed by iteration.
#   nu1_all_orbit(mu0)            mu0/(mu0-1): every interior orbit loses its predator if and only if
#                                 nu1 <= this value, i.e. 1/mu0 + 1/nu1 >= 1 (Theorem 5.2);
#   coexistence_point(mu0, nu1)   (x*, y*) = (mu0/(mu0+nu1), 1 - 1/mu0 - 1/nu1), in the open square iff y* > 0;
#   coexistence_jacobian(...)     the Jacobian of the map at that point;
#   nu1_flip(mu0)                 the curve on which that point has a multiplier -1: the positive root of
#                                 (5-mu0) nu1^2 + mu0(mu0+2) nu1 - 3 mu0^2 = 0; the point is a sink above it;
#   nu1_star_closed(mu0)          exp(-<ln x>) on the prey's attracting cycle where it is known in closed form:
#                                 mu0/(mu0-1) on the fixed point (1 < mu0 <= 3), mu0/sqrt(mu0+1) on the
#                                 two-cycle (3 < mu0 <= 1 + sqrt 6), 4 at mu0 = 4 (arcsine law); nan elsewhere.
TWO_CYCLE_END = 1.0 + math.sqrt(6.0)           # the prey's two-cycle attracts for 3 < mu0 <= 1 + sqrt 6
JURY_EDGE_AT_4 = 8.0 * math.sqrt(3.0) - 12.0   # nu1_flip(4), the lower Jury edge of the main text
# The theorem numbers printed inside Figure 8 are those of the manuscript; they must follow its numbering.
FIG8_LABEL_ALL_ORBIT_THEOREM = "Theorem 5.2"
FIG8_LABEL_BASIN_THEOREM = "Theorem 3.2"


def nu1_all_orbit(mu0):
    """mu0/(mu0-1) for mu0 > 1; inf for mu0 <= 1 (there every interior orbit loses its predator for every nu1)."""
    mu0 = np.asarray(mu0, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(mu0 > 1.0, mu0 / (mu0 - 1.0), np.inf)


def nu1_flip(mu0):
    """The positive root of (5-mu0) nu1^2 + mu0 (mu0+2) nu1 - 3 mu0^2 = 0, for 0 < mu0 <= 4."""
    mu0 = np.asarray(mu0, dtype=float)
    return mu0 * (np.sqrt(mu0 * mu0 - 8.0 * mu0 + 64.0) - mu0 - 2.0) / (2.0 * (5.0 - mu0))


def nu1_star_closed(mu0):
    """Invasion threshold of the prey's attracting cycle where it has a closed form (see the block comment)."""
    mu0 = np.asarray(mu0, dtype=float)
    out = np.full(mu0.shape, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where((mu0 > 1.0) & (mu0 <= 3.0), mu0 / (mu0 - 1.0), out)
        out = np.where((mu0 > 3.0) & (mu0 <= TWO_CYCLE_END), mu0 / np.sqrt(mu0 + 1.0), out)
    return np.where(mu0 == 4.0, 4.0, out)


def coexistence_point(mu0, nu1):
    """(x*, y*); it lies in the open square exactly when y* > 0, i.e. 1/mu0 + 1/nu1 < 1."""
    return mu0 / (mu0 + nu1), 1.0 - 1.0 / mu0 - 1.0 / nu1


def coexistence_jacobian(mu0, nu1):
    """Jacobian of (x, y) -> (mu0 x(1-x)(1-y), nu1 x y(1-y)) at the coexistence point."""
    xs, ys = coexistence_point(mu0, nu1)
    return np.array([[mu0 * (1.0 - 2.0 * xs) * (1.0 - ys), -mu0 * xs * (1.0 - xs)],
                     [nu1 * ys * (1.0 - ys), nu1 * xs * (1.0 - 2.0 * ys)]])


def coexistence_is_sink(mu0, nu1):
    """True when the coexistence point exists (y* > 0) and both its multipliers lie in the open unit disc."""
    _, ys = coexistence_point(mu0, nu1)
    return bool(ys > 0.0 and np.all(np.abs(np.linalg.eigvals(coexistence_jacobian(mu0, nu1))) < 1.0))


# =============================================================================
# (a) Arcsine moments via Gauss--Chebyshev quadrature
# =============================================================================
#
# For the arcsine weight w(x) = 1 / (pi sqrt(x(1-x))) on (0,1), substitute
# x = (1 - cos theta)/2 = sin^2(theta/2), theta in (0, pi).  Then
#   int_0^1 g(x) dmu4(x) = (1/pi) int_0^pi g(sin^2(theta/2)) dtheta.
# Gauss--Chebyshev (first kind) with nodes theta_k = pi (k - 1/2)/M, k=1..M,
# integrates (1/pi) int_0^pi f(theta) dtheta as a uniform average f(theta_k)/M,
# which is the midpoint rule and is spectrally accurate for smooth integrands.
# The only singularity of ln x is the integrable log at x=0 (theta=0), a null
# point of the measure; the midpoint nodes avoid the endpoints exactly.


def arcsine_quadrature_moment(func, n_nodes: int) -> float:
    """int_0^1 func(x) dmu4(x) via the Gauss--Chebyshev / midpoint rule.

    Nodes x_k = sin^2(theta_k/2), theta_k = pi (k - 0.5)/n_nodes, weights 1/n_nodes.
    """
    k = np.arange(1, n_nodes + 1, dtype=np.float64)
    theta = math.pi * (k - 0.5) / n_nodes
    x = np.sin(theta / 2.0) ** 2
    vals = func(x)
    return float(np.mean(vals))


def check_arcsine_moments(n_nodes: int) -> Tuple[List["Check"], Dict[str, object]]:
    """(a) Arcsine moments + Jensen gap (eq:arcsine-moments, eq:jensen-gap)."""
    rows: List[Check] = []
    e_lnx = arcsine_quadrature_moment(lambda x: np.log(x), n_nodes)
    e_x = arcsine_quadrature_moment(lambda x: x, n_nodes)
    jensen_gap = e_lnx - math.log(e_x)
    rec = {
        "n_nodes": n_nodes,
        "E_lnx": e_lnx,
        "E_lnx_closed_form": -2.0 * LN2,
        "E_x": e_x,
        "E_x_closed_form": 0.5,
        "jensen_gap": jensen_gap,
        "jensen_gap_closed_form": -LN2,
    }
    rows.append(
        Check("E_mu4[ln x]", -2.0 * LN2, e_lnx, 1e-3, "= -2 ln2 = -1.386294 (eq:arcsine-moments)")
    )
    rows.append(Check("E_mu4[x]", 0.5, e_x, 1e-3, "= 1/2 (eq:arcsine-moments)"))
    rows.append(
        Check("Jensen gap E[ln x]-ln E[x]", -LN2, jensen_gap, 1e-3, "= -ln2 = -0.693147 (eq:jensen-gap)")
    )
    return rows, rec


# =============================================================================
# (b) Equilibrium-measure identity  int ln|x - r| dmu4(x) = -ln4  for r in [0,1]
# =============================================================================
#
# The arcsine measure is the logarithmic equilibrium measure of [0,1] (capacity
# cap = 1/4), so the energy integral int ln|x-r| dmu4 equals ln(cap) = ln(1/4)
# = -ln4 on the support (thm:equilibrium-measure, eq:equilibrium-potential).
# (The standard logarithmic potential U(r) = int ln(1/|x-r|) dmu4 = -ln(cap) =
# +ln4; the sign differs because this code integrates ln|x-r|, not its negative.)
# For r INSIDE [0,1] the
# integrand has an integrable log-singularity at x=r; for r OUTSIDE [0,1] there
# is no singularity and the value departs from -ln4.
#
# Off-support reference values: there is a closed form
#   int ln|x - r| dmu4(x) = ln( (|2r - 1| + 2 sqrt(r^2 - r)) / 4 )   for r not in [0,1],
# obtained from Jensen's formula for the quadratic with the arcsine substitution
# (printed in the caption of Figure S1; Table S10 lists both values).  We use it to corroborate r=1.2 (~ -0.519)
# and r=2.0 (~ +0.376).


def _potential_offsupport_closed_form(r: float) -> float:
    """Closed-form logarithmic potential of mu4 at r OUTSIDE [0,1]."""
    disc = r * r - r  # r^2 - r > 0 for r<0 or r>1
    return math.log((abs(2.0 * r - 1.0) + 2.0 * math.sqrt(disc)) / 4.0)


def _potential_quadrature(r: float, n_nodes: int) -> float:
    """int ln|x - r| dmu4(x) via arcsine quadrature.

    For r in (0,1) the midpoint nodes x_k = sin^2(theta_k/2) generally miss x=r,
    but this is NOT guaranteed: for r=0.5 and ODD n the central node lands on
    x=0.5 in exact arithmetic (saved here only by float64 round-off, |x_k-0.5| ~
    1e-16, plus the 1e-300 floor below).  Because the log-singularity is
    INTEGRABLE the midpoint rule still converges; we further damp the
    slowly-converging singular tail by averaging two staggered node counts
    (n and n+1), one of which always misses the singularity.  Empirically this
    matches -ln4 to ~5.6e-5 at the default node count (cross-checked against the
    off-support closed form and a node-doubling test), i.e. convergent-by-
    integrability, NOT a knife-edge.  A principled singular quadrature would
    remove the round-off dependence but is unnecessary at the accuracy this
    consistency check requires.
    """
    def potential_with(n: int) -> float:
        k = np.arange(1, n + 1, dtype=np.float64)
        theta = math.pi * (k - 0.5) / n
        x = np.sin(theta / 2.0) ** 2
        d = np.abs(x - r)
        d = np.where(d < 1e-300, 1e-300, d)
        return float(np.mean(np.log(d)))

    return 0.5 * (potential_with(n_nodes) + potential_with(n_nodes + 1))


def check_equilibrium_measure(n_nodes: int) -> Tuple[List["Check"], Dict[str, object]]:
    """(b) Equilibrium-measure identity (thm:equilibrium-measure)."""
    rows: List[Check] = []
    rec: Dict[str, object] = {"n_nodes": n_nodes, "on_support": {}, "off_support": {}}
    # On-support r values: identity predicts constant -ln4.
    on_support_rs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for r in on_support_rs:
        if r <= 0.0 or r >= 1.0:
            # endpoints: no interior singularity; use a high-node midpoint estimate.
            val = _potential_quadrature(r, max(n_nodes, 200000))
        else:
            val = _potential_quadrature(r, n_nodes)
        rec["on_support"][f"{r:.2f}"] = val
        rows.append(
            Check(
                f"int ln|x-r| dmu4, r={r:.2f}",
                predicted=-2.0 * LN2,  # -ln4
                computed=val,
                tol=1e-3 if 0.0 < r < 1.0 else 5e-3,
                note="= -ln4 = -1.386294 on [0,1] (eq:equilibrium-potential)",
            )
        )
    # Off-support reference values (no singularity): r = 1.2 (~ -0.519), r = 2.0 (~ +0.376).
    for r, target in [(1.2, -0.519), (2.0, 0.376)]:
        val = _potential_quadrature(r, n_nodes)
        cf = _potential_offsupport_closed_form(r)
        rec["off_support"][f"{r:.2f}"] = {"quadrature": val, "closed_form": cf, "spec_target": target}
        rows.append(
            Check(
                f"int ln|x-r| dmu4, r={r:.2f} (off support)",
                predicted=cf,
                computed=val,
                tol=1e-3,
                note=f"closed form {cf:.4f}; spec target ~{target:+.3f} (departs from -ln4)",
            )
        )
    return rows, rec


# =============================================================================
# (c) Ergodic optimization: sub-action H(x) and period-<=16 orbit enumeration
# =============================================================================
#
# Sub-action (eq:subaction) at mu0 = 4, g(x) = 4x/9 (a constant added to g changes nothing below):
#   H(x) = ln x + g(T(x)) - g(x) = ln x + (4/9)(3x - 4x^2),  T(x) = 4x(1-x).
#   H'(x) = 1/x + (4/9)(3 - 8x),  H''(x) = -1/x^2 - 32/9 < 0 on (0,1].
#   H(x) <= ln(3/4) with equality at x = 3/4.


def check_subaction(n_grid: int) -> Tuple[List["Check"], Dict[str, object]]:
    """(c.1) Sub-action H(x) <= ln(3/4), equality at 3/4, strict concavity (eq:subaction)."""
    rows: List[Check] = []
    xs = np.linspace(1e-6, 1.0, n_grid)
    H = np.log(xs) + (4.0 / 9.0) * (3.0 * xs - 4.0 * xs * xs)
    Hpp = -1.0 / (xs * xs) - 32.0 / 9.0
    max_H = float(np.max(H))
    argmax_H = float(xs[int(np.argmax(H))])
    H_at_three_quarter = math.log(0.75) + (4.0 / 9.0) * (3.0 * 0.75 - 4.0 * 0.75 * 0.75)
    rec = {
        "n_grid": n_grid,
        "max_H": max_H,
        "argmax_H": argmax_H,
        "H_at_3/4": H_at_three_quarter,
        "ln_3/4": LN_THREE_QUARTER,
        "Hpp_max": float(np.max(Hpp)),  # must be < 0
        "Hpp_all_negative": bool(np.all(Hpp < 0.0)),
    }
    rows.append(
        Check("sup_x H(x)", LN_THREE_QUARTER, max_H, 1e-3, "<= ln(3/4) = -0.287682 (eq:subaction)")
    )
    rows.append(
        Check("H(3/4) value", LN_THREE_QUARTER, H_at_three_quarter, 1e-12, "equality at x=3/4")
    )
    rows.append(
        Check("argmax H(x)", 0.75, argmax_H, 2e-3, "interior maximiser x*=3/4")
    )
    # Strict concavity. TWO rows, because one was doing the work of two and lying about it.
    #
    # DEFECT FIXED 2026-07-17. This row previously read
    #     Check("max H''(x) (concavity)", 0.0, max(0.0, float(np.max(Hpp))), ...)
    # and so REPORTED 0 under the label "max H''(x)". But 0 is not max H''; it is
    # max(0, max H'') -- a CLIPPED VIOLATION RESIDUAL, whose 0 means "no violation
    # found". The actual maximum is attained at the right endpoint,
    #     H''(x) = -1/x^2 - 32/9,   max_{(0,1]} H'' = H''(1) = -1 - 32/9 = -41/9 ~ -4.5556,
    # which Section 5 and Section S7 already quote correctly as <= -4.556. The table row was the
    # only place that said 0, and it said it because it was named after the quantity we
    # wanted rather than the one it computed. The honest value was even already in the
    # record beside it, as `Hpp_max`. Split into: (1) the violation residual, under its
    # real name; (2) the physical maximum, against its closed form.
    rows.append(
        Check("nonconcavity violation max(H'',0)", 0.0, max(0.0, float(np.max(Hpp))), 1e-12,
              "clipped residual: 0 == no nonconcave grid point. NOT max H'' -- see next row")
    )
    rows.append(
        Check("max H''(x) on (0,1]", -41.0 / 9.0, float(np.max(Hpp)), 1e-6,
              "closed form H''(1) = -1 - 32/9 = -41/9 ~ -4.5556 (strictly negative => concave)")
    )
    return rows, rec


def _mobius(n: int) -> int:
    """Moebius function mu(n) (exact integer arithmetic)."""
    if n == 1:
        return 1
    result = 1
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0  # squared prime factor
            result = -result
        p += 1
    if m > 1:
        result = -result
    return result


def _divisors(n: int) -> List[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def _least_period_points(n: int) -> int:
    """Number of points of LEAST period n for T(x)=4x(1-x) (== count for tent map).

    Fix(T^n) = 2^n, so by Moebius inversion the number of least-period-n points is
        sum_{d|n} mu(n/d) 2^d .
    The number of orbits of least period n is this divided by n (Lyndon count).
    These are EXACT integer counts -- no floating-point iteration is involved, so
    they do not suffer the exponential precision loss of forward orbit iteration.
    """
    return sum(_mobius(n // d) * (2 ** d) for d in _divisors(n))


def _periodic_angles(max_period: int) -> set:
    """Exact rational angles ``r`` (0 <= r <= 1/2) of EVERY periodic point of T(x)=4x(1-x)
    of least period <= ``max_period``.

    Under the conjugacy x = sin^2(pi r), with r identified with 1 - r so that
    0 <= r <= 1/2, the mu0 = 4 logistic map is the angle-doubling map
        D(r) = min(2r mod 1, 1 - (2r mod 1)),
    and D^n(r) = r iff 2^n r = +/- r (mod 1), i.e. r = j/(2^n - 1) or r = j/(2^n + 1).
    Every member is a fractions.Fraction: no floating-point number enters the set, so
    membership, orbit closure and least periods below are decided exactly.
    """
    angles: set = set()
    for n in range(1, max_period + 1):
        for den in (2 ** n - 1, 2 ** n + 1):
            angles.update(Fraction(j, den) for j in range(den // 2 + 1))
    return angles


def _doubling(r: Fraction) -> Fraction:
    """Exact angle-doubling D(r) = min(2r mod 1, 1 - (2r mod 1)) on fractions.Fraction."""
    t = (2 * r) % 1
    return min(t, 1 - t)


def enumerate_periodic_orbits(max_period: int) -> Tuple[Dict[str, object], List["Check"]]:
    """Exact enumeration of EVERY interior periodic orbit of T(x)=4x(1-x) of least
    period <= ``max_period`` and the float64 orbit average of ln x on each one
    (rmk:ergodic-numeric).

    Two layers, kept separate on purpose:

    1. EXACT LAYER (fractions.Fraction; no floating point).  The periodic points are the
       rational angles of ``_periodic_angles``; each orbit is closed by exact doubling and
       its least period is the exact cycle length.  The per-period orbit and point counts
       so MEASURED are cross-checked against the independent Moebius/Lyndon formula
       (#least-period-n points = sum_{d|n} mu(n/d) 2^d, orbits = that / n).  The single
       boundary periodic point r = 0 (x = 0, ln x = -inf) is excluded; the angle 1/2
       (x = 1) is a pre-image of 0, not periodic, and never occurs (the denominators are
       odd).
    2. FLOAT LAYER.  For every interior orbit the average of ln x over its points
       x_k = sin^2(pi r_k) is evaluated in double precision (math.fsum).  These averages are
       consistency values, not interval-certified bounds.  The analytic subaction of
       thm:ergodic-opt is the PROOF that the supremum over invariant measures is ln(3/4),
       attained at the fixed point x* = 3/4 (angle 1/3); this routine only checks that no
       periodic orbit through period ``max_period`` beats it, with coverage reported per
       least period and asserted equal to the exact count independently of the maximiser.

    History: through 2026-10-17 this routine rounded the sin^2 values to nine decimals and
    re-identified each candidate by unstable float64 iteration with a return tolerance of
    1e-9; under the pinned environment that search accepted only twenty starting points
    through period 12 while the manuscript described it as exhaustive through period 16.
    It was replaced, not patched.
    """
    # ---- (1) exact combinatorial counts (Moebius/Lyndon), independent of the enumeration
    points_by_period = {n: _least_period_points(n) for n in range(1, max_period + 1)}
    orbits_by_period = {n: points_by_period[n] // n for n in range(1, max_period + 1)}
    total_points = sum(points_by_period.values())
    total_orbits = sum(orbits_by_period.values())
    interior_points = total_points - 1      # minus the boundary fixed point x = 0
    interior_orbits = total_orbits - 1

    # ---- (2) exact enumeration of the orbits themselves ---------------------------
    angles = _periodic_angles(max_period)
    unseen = set(angles)
    enum_orbits = {n: 0 for n in range(1, max_period + 1)}
    enum_points = {n: 0 for n in range(1, max_period + 1)}
    eval_orbits = {n: 0 for n in range(1, max_period + 1)}
    max_by_period: Dict[int, float] = {}
    best = (-math.inf, 0, Fraction(0))        # (mean ln x, least period, smallest angle)
    runner_up = (-math.inf, 0, Fraction(0))
    boundary_orbits = 0
    while unseen:
        start = unseen.pop()
        cycle = [start]
        r = _doubling(start)
        while r != start:
            if r not in unseen:   # closed, disjoint cycles: a repeat would be a logic error
                raise AssertionError(f"exact doubling left the periodic set: {start} -> {r}")
            unseen.remove(r)
            cycle.append(r)
            r = _doubling(r)
        n = len(cycle)
        if n > max_period:
            raise AssertionError(f"cycle of length {n} > max_period {max_period} from angle {start}")
        if n == 1 and cycle[0] == 0:
            boundary_orbits += 1          # x = 0: the sole boundary periodic point, excluded
            continue
        enum_orbits[n] += 1
        enum_points[n] += n
        # ---- (3) float64 evaluation of the orbit average of ln x --------------------
        mean = math.fsum(math.log(math.sin(math.pi * float(q)) ** 2) for q in cycle) / n
        eval_orbits[n] += 1
        if mean > max_by_period.get(n, -math.inf):
            max_by_period[n] = mean
        cand = (mean, n, min(cycle))
        if mean > best[0]:
            best, runner_up = cand, best
        elif mean > runner_up[0]:
            runner_up = cand

    # The formula counts include the boundary fixed point x = 0 (least period 1); the
    # enumeration counts are interior, so the n = 1 reference is reduced by that one point.
    mismatches = sum(
        1 for n in range(1, max_period + 1)
        if enum_orbits[n] != orbits_by_period[n] - (1 if n == 1 else 0)
        or enum_points[n] != points_by_period[n] - (1 if n == 1 else 0)
    )
    # Largest average over least periods 2..max_period: every one must stay strictly below
    # ln(3/4). From period 3 on they increase with the period (long orbits shadowing the
    # repelling fixed point 3/4 spend most of their time near it; period 2 lies above period 3),
    # so the closest competitor is a long orbit, not the period-2 orbit; the
    # gap to ln(3/4) is recorded, never assumed.
    n_periods_at_or_above = sum(1 for n, v in max_by_period.items() if n >= 2 and v >= LN_THREE_QUARTER)
    closest_gap = min((LN_THREE_QUARTER - v for n, v in max_by_period.items() if n >= 2), default=math.inf)
    n_enum_orbits = sum(enum_orbits.values())
    n_enum_points = sum(enum_points.values())
    n_eval_orbits = sum(eval_orbits.values())
    coverage_complete = bool(
        mismatches == 0 and boundary_orbits == 1
        and n_enum_orbits == interior_orbits == n_eval_orbits
        and n_enum_points == interior_points
    )
    best_x0 = math.sin(math.pi * float(best[2])) ** 2
    runner_x0 = math.sin(math.pi * float(runner_up[2])) ** 2

    rec = {
        "max_period": max_period,
        "enumeration_method": (
            "exact rational angle doubling (fractions.Fraction) on the conjugacy x = sin^2(pi r); "
            "least periods and counts exact; the ln x orbit averages are float64 (math.fsum)"
        ),
        "points_by_least_period": points_by_period,            # Moebius/Lyndon formula (exact)
        "orbits_by_least_period": orbits_by_period,
        "enumerated_orbits_by_least_period": enum_orbits,       # measured by the exact enumeration
        "enumerated_points_by_least_period": enum_points,
        "evaluated_orbits_by_least_period": eval_orbits,        # orbits whose ln x average was evaluated
        "max_orbit_avg_lnx_by_least_period": {n: max_by_period[n] for n in sorted(max_by_period)},
        "total_periodic_points": int(total_points),
        "total_orbits": int(total_orbits),
        "interior_periodic_points_exact": int(interior_points),
        "interior_orbits_exact": int(interior_orbits),
        "interior_orbits_enumerated": int(n_enum_orbits),
        "interior_points_enumerated": int(n_enum_points),
        "interior_orbits_evaluated": int(n_eval_orbits),
        "boundary_orbits_excluded": int(boundary_orbits),
        "least_period_count_mismatches": int(mismatches),
        "coverage_complete": coverage_complete,
        "max_orbit_avg_lnx": float(best[0]),
        "argmax_orbit_period": int(best[1]),
        "argmax_orbit_angle": str(best[2]),
        "argmax_orbit_x0": float(best_x0),
        "runner_up_avg_lnx": float(runner_up[0]),
        "runner_up_period": int(runner_up[1]),
        "runner_up_angle": str(runner_up[2]),
        "runner_up_x0": float(runner_x0),
        "period_2_orbit_avg_lnx_closed_form": 0.5 * math.log(5.0 / 16.0),
        "n_least_periods_2_to_max_with_avg_at_or_above_ln34": int(n_periods_at_or_above),
        "closest_gap_below_ln34_over_periods_2_to_max": float(closest_gap),
        "target_max_avg_ln_3/4": LN_THREE_QUARTER,
        "float_evaluation_note": (
            "counts and periods are exact; the sin/log orbit averages are double precision, "
            "not interval-certified; the analytic subaction (thm:ergodic-opt) is the proof"
        ),
    }
    rows: List[Check] = []
    # Load-bearing consistency check: the largest orbit average is ln(3/4) at x* = 3/4.
    rows.append(Check(
        f"max orbit-avg ln x (all interior orbits, period<= {max_period})",
        predicted=LN_THREE_QUARTER, computed=float(best[0]), tol=1e-9,
        note="= ln(3/4) = -0.287682 at x*=3/4, angle 1/3 (rmk:ergodic-numeric); float64 evaluation of every orbit",
    ))
    rows.append(Check(
        "argmax orbit point x0", predicted=0.75, computed=float(best_x0), tol=1e-9,
        note="maximiser is the interior fixed point x*=3/4 (least period 1)",
    ))
    rows.append(Check(
        "argmax least period", predicted=1.0, computed=float(best[1]), tol=0.0,
        note="the maximising orbit is the fixed point (period 1)",
    ))
    # Independent oracle for the float layer: the unique period-2 orbit {(5 -+ sqrt5)/8}
    # has x1 x2 = 20/64, so its ln x average is ln(5/16)/2 = -0.581575 exactly. It applies
    # only when period 2 was requested (max_period=1 used to fail it).
    if max_period >= 2:
        rows.append(Check(
            "period-2 orbit-avg ln x (closed form ln(5/16)/2)",
            predicted=0.5 * math.log(5.0 / 16.0), computed=float(max_by_period.get(2, math.nan)), tol=1e-9,
            note="the single least-period-2 orbit, evaluated by the same float64 path as every other orbit",
        ))
    rows.append(Check(
        f"least periods 2..{max_period} whose largest orbit average reaches ln(3/4)",
        predicted=0.0, computed=float(n_periods_at_or_above), tol=0.0,
        note=f"all strictly below; closest gap {closest_gap:.4f} (period {runner_up[1]}, an orbit shadowing x*=3/4)",
    ))
    # Coverage checks, independent of the maximising value.
    if max_period == 16:
        pred_orbits, pred_points = 8799.0, 130485.0
        orbit_note = "exact Lyndon/Moebius reference for the period-16 protocol"
        point_note = "exact point count minus the boundary fixed point x=0"
    else:
        pred_orbits, pred_points = float(interior_orbits), float(interior_points)
        orbit_note = f"exact Lyndon/Moebius count for period<= {max_period} (no calibrated reference)"
        point_note = f"exact point count for period<= {max_period} (no calibrated reference)"
    rows.append(Check(
        f"interior-orbit count (exact enumeration, period<= {max_period})",
        predicted=pred_orbits, computed=float(n_enum_orbits), tol=0.5, note=orbit_note,
    ))
    rows.append(Check(
        f"interior periodic-point count (exact enumeration, period<= {max_period})",
        predicted=pred_points, computed=float(n_enum_points), tol=0.5, note=point_note,
    ))
    rows.append(Check(
        "least-period count mismatches (enumeration vs Moebius/Lyndon)",
        predicted=0.0, computed=float(mismatches), tol=0.0,
        note="every least period 1..max_period: enumerated orbits/points == formula",
    ))
    rows.append(Check(
        f"orbits with an evaluated ln x average (coverage, period<= {max_period})",
        predicted=float(interior_orbits), computed=float(n_eval_orbits), tol=0.5,
        note="coverage measured, not inferred from the maximising value",
    ))
    rows.append(Check(
        "boundary periodic orbits excluded", predicted=1.0, computed=float(boundary_orbits), tol=0.0,
        note="exactly one: the fixed point x=0 (angle 0)",
    ))
    return rec, rows


# =============================================================================
# (d) I(mu0) family and threshold nu1*(mu0) = exp(-I(mu0))
# =============================================================================


@njit(cache=True, fastmath=False)
def _logistic_I_average(mu0: float, x0: float, n_burn: int, n_keep: int) -> float:
    """I(mu0) = (1/N) sum ln x_n over the logistic SRB orbit at mu0.

    Used for mu0 < 4 (single long orbit).  At mu0 = 4 the orbit average is
    unreliable (rare deep x->0 excursions corrupt the floating-point sum), so
    the caller substitutes the EXACT arcsine value -2 ln2 instead.
    """
    x = x0
    for _ in range(n_burn):
        x = mu0 * x * (1.0 - x)
    s = 0.0
    for _ in range(n_keep):
        x = mu0 * x * (1.0 - x)
        xc = x if x > 1e-300 else 1e-300
        s += math.log(xc)
    return s / n_keep


@njit(cache=True, fastmath=False)
def _logistic_I_multi_ic(
    mu0: float, n_ic: int, n_burn: int, n_keep: int
) -> Tuple[float, float]:
    """Many-IC SRB average of ln x at mu0 (rmk:Icurve), plus prey Lyapunov exponent.

    Averages over ``n_ic`` interior initial conditions to stabilise the estimate
    against single-orbit pathologies.  Returns (mean I(mu0), mean prey-lyap).
    The prey Lyapunov exponent lambda_prey = mean ln|mu0 (1-2x)| flags whether
    mu0 sits in a chaotic band (>0) or a periodic window (<=0).
    """
    s_I = 0.0
    s_ly = 0.0
    total = 0
    for j in range(n_ic):
        x0 = (j + 0.5) / n_ic  # spread ICs across (0,1)
        x = x0
        for _ in range(n_burn):
            x = mu0 * x * (1.0 - x)
        for _ in range(n_keep):
            d = abs(mu0 * (1.0 - 2.0 * x))
            s_ly += math.log(d) if d > 1e-300 else -50.0
            x = mu0 * x * (1.0 - x)
            xc = x if x > 1e-300 else 1e-300
            s_I += math.log(xc)
            total += 1
    return s_I / total, s_ly / total


def check_I_family(n_burn: int, n_keep: int, n_ic: int) -> Tuple[List["Check"], Dict[str, object]]:
    """(d) I(mu0) family + threshold nu1*(mu0) = exp(-I) (cor:mu0family, rmk:Icurve)."""
    rows: List[Check] = []
    rec: Dict[str, object] = {"n_burn": n_burn, "n_keep": n_keep, "n_ic": n_ic, "family": {}}

    # mu0 = 4: EXACT arcsine quadrature, NOT orbit average (rmk:Icurve).
    I4_quad = arcsine_quadrature_moment(lambda x: np.log(x), 200000)
    nu1star_4 = math.exp(-I4_quad)
    rec["I_4_quadrature"] = I4_quad
    rec["I_4_closed_form"] = -2.0 * LN2
    rec["nu1star_4"] = nu1star_4
    rows.append(
        Check("I(4) (arcsine quadrature)", -2.0 * LN2, I4_quad, 1e-3, "= -2 ln2 = -1.386294 (rmk:Icurve)")
    )
    rows.append(
        Check("nu1*(4) = exp(-I(4))", 4.0, nu1star_4, 1e-2,
              "closed form exp(2 ln 2) = 4.000 at the band edge mu0=4 (rmk:Icurve); no maximality over the family is asserted")
    )

    # mu0 < 4 family: many-IC SRB averages. Targets (design specification, section 5):
    #   nu1* ~ 1.653 @ 3.60, ~1.592 @ 3.70, ~1.765 @ 3.74, ~2.038 @ 3.90.
    # These are not monotone in mu0 (they fall from 3.60 to 3.70). Tolerance is wide (consistency).
    family_targets = {
        3.60: 1.653,
        3.70: 1.592,
        3.74: 1.765,
        3.90: 2.038,
        3.95: None,
        3.99: None,
    }
    for mu0 in sorted(family_targets.keys()):
        I_mu, lyap = _logistic_I_multi_ic(mu0, n_ic, n_burn, n_keep)
        nu1star = math.exp(-I_mu)
        all_orbit_bound = mu0 / (mu0 - 1.0)
        rec["family"][f"{mu0:.2f}"] = {
            "I": I_mu,
            "nu1star_ae": nu1star,
            "prey_lyap": lyap,
            "chaotic": bool(lyap > 5e-3),
            "all_orbit_bound_mu0_over_mu0m1": all_orbit_bound,
            "spec_target_nu1star": family_targets[mu0],
        }
        tgt = family_targets[mu0]
        if tgt is not None:
            rows.append(
                Check(
                    f"nu1*({mu0:.2f}) = exp(-I)",
                    predicted=tgt,
                    computed=nu1star,
                    tol=0.25,  # wide: orbit-average + periodic-window sensitivity
                    note="non-monotone family threshold (rmk:Icurve); consistency, not certified",
                )
            )

    # All-orbit uniform bound mu0/(mu0-1), minimised at mu0=4 as 4/3.
    bound_at_4 = 4.0 / (4.0 - 1.0)
    rec["all_orbit_bound_min_at_4"] = bound_at_4
    rows.append(
        Check("all-orbit bound mu0/(mu0-1) at mu0=4", 4.0 / 3.0, bound_at_4, 1e-12, "min = 4/3 = 1.3333 (cor:mu0family)")
    )
    return rows, rec


# =============================================================================
# Validation-table machinery
# =============================================================================


class Check:
    """A single PASS/FAIL row of the validation table."""

    __slots__ = ("name", "predicted", "computed", "tol", "abs_err", "passed", "note")

    def __init__(self, name: str, predicted: float, computed: float, tol: float, note: str = "") -> None:
        self.name = name
        self.predicted = float(predicted)
        self.computed = float(computed)
        self.tol = float(tol)
        self.abs_err = abs(self.predicted - self.computed)
        self.passed = self.abs_err <= self.tol
        self.note = note

    def as_dict(self) -> Dict[str, object]:
        return {
            "quantity": self.name,
            "predicted": self.predicted,
            "computed": self.computed,
            "abs_err": self.abs_err,
            "tol": self.tol,
            "pass": bool(self.passed),
            "note": self.note,
        }


# =============================================================================
# Optional diagnostic figures (P2; guarded by --no-figs)
# =============================================================================


def _draw_parameter_box(family: Dict[str, Dict[str, object]], out_png: str, dpi: int = 300) -> str:
    """Manuscript Figure 8, two panels 6.27 in wide (the manuscript's line width; the figure is included at that
    width, so the sizes below are the printed sizes): (a) the box 1 <= mu0 <= 4, 0 < nu1 <= 4; (b) the window
    3 <= mu0 <= 4, 1.25 <= nu1 <= 2.30. Every curve is a closed form (the functions above); the markers are the
    sampled thresholds nu1star_ae of the I_family block (filled: chaotic parameter; open: periodic window), and
    the star is the closed-form value nu1*(4) = 4, as in the earlier form of this figure (the arcsine quadrature's
    3.99999 is reported in the values JSON and in Table S10, not plotted as if it were exact)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    sampled = sorted((float(k), float(v["nu1star_ae"]), bool(v["chaotic"])) for k, v in family.items())
    fs = 7.4            # smallest lettering, in printed points
    rc = {"font.size": fs, "axes.labelsize": 8.0, "xtick.labelsize": fs, "ytick.labelsize": fs,
          "mathtext.fontset": "dejavusans", "axes.linewidth": 0.7, "hatch.linewidth": 0.45}
    grey, dark, edge = "#d6d6d6", "#000000", "#3c3c3c"
    win = (3.0, 4.0, 1.25, 2.30)
    ls_two = (0, (4, 2))
    ls_flip = (0, (6, 1.6, 1.2, 1.6))
    with plt.rc_context(rc):
        fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.27, 3.10), dpi=dpi)
        fig.subplots_adjust(left=0.066, right=0.962, bottom=0.262, top=0.935, wspace=0.25)

        def markers(axis, pts, ms):
            for mu0, v, chaotic in pts:
                axis.plot([mu0], [v], marker="o", ms=ms + 1.6, mfc="white", mec="white", mew=0, ls="none", zorder=6, clip_on=False)
                axis.plot([mu0], [v], marker="o", ms=ms, mfc=dark if chaotic else "white", mec=dark, mew=0.8, ls="none",
                          zorder=7, clip_on=False)

        # ---- panel (a): the box
        m = np.linspace(1.0, 4.0, 1201)
        top = np.minimum(4.0, nu1_all_orbit(np.maximum(m, 1.0 + 1e-9)))
        ax.fill_between(m, 0.0, top, color=grey, lw=0, zorder=1)
        mc = np.linspace(4.0 / 3.0, 4.0, 801)
        ax.plot(mc, nu1_all_orbit(mc), color=dark, lw=1.5, zorder=4)
        m2 = np.linspace(3.0, TWO_CYCLE_END, 200)
        ax.plot(m2, nu1_star_closed(m2), color=dark, lw=1.1, ls=ls_two, zorder=4)
        ax.plot([4.0, 4.0], [0.0, 4.0], color=edge, lw=2.4, zorder=5, clip_on=False, solid_capstyle="butt")
        markers(ax, sampled, 3.3)
        ax.plot([4.0], [4.0], marker="*", ms=9.5, mfc=dark, mec="white", mew=0.5, ls="none", zorder=8, clip_on=False)
        ax.add_patch(Rectangle((win[0], win[2]), win[1] - win[0], win[3] - win[2], fill=False, ec=dark, lw=0.7,
                               ls=(0, (1.5, 1.5)), zorder=5))
        ax.set_xlim(1.0, 4.0)
        ax.set_ylim(0.0, 4.0)
        ax.set_xticks([1.0, 4.0 / 3.0, 2.0, 3.0, 4.0])
        ax.set_xticklabels(["1", "4/3", "2", "3", "4"])
        ax.set_yticks([0.0, 1.0, 4.0 / 3.0, 2.0, 3.0, 4.0])
        ax.set_yticklabels(["0", "1", "4/3", "2", "3", "4"])
        ax.set_xlabel(r"prey growth rate $\mu_0$", labelpad=1.5)
        ax.set_ylabel(r"conversion rate $\nu_1$", labelpad=1.5)
        ax.tick_params(length=2.5, width=0.6, pad=2)
        thm = (" (%s)" % FIG8_LABEL_ALL_ORBIT_THEOREM) if FIG8_LABEL_ALL_ORBIT_THEOREM else ""
        ax.text(2.45, 0.64, "every interior orbit loses its predator\n" + r"$1/\mu_0+1/\nu_1\geq 1$" + thm, ha="center",
                va="center", fontsize=fs, linespacing=1.3, zorder=9)
        ax.text(2.80, 3.30, "the coexistence fixed point exists\n" + r"$1/\mu_0+1/\nu_1<1$", ha="center", va="center",
                fontsize=fs, linespacing=1.3, zorder=9)
        ax.text(3.03, 2.345, "window of (b)", ha="left", va="bottom", fontsize=fs, zorder=9)
        ax.text(0.0, 1.03, "(a)", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.4, fontweight="bold")

        # ---- panel (b): the window
        mb = np.linspace(win[0], win[1], 801)
        allb = nu1_all_orbit(mb)
        flb = nu1_flip(mb)
        bx.fill_between(mb, win[2], allb, color=grey, lw=0, zorder=1)
        bx.fill_between(mb, flb, win[3], facecolor="none", edgecolor="#a6a6a6", hatch="////", lw=0, zorder=1)
        bx.plot(mb, allb, color=dark, lw=1.5, zorder=4)
        bx.plot(mb, flb, color=dark, lw=1.1, ls=ls_flip, zorder=4)
        bx.plot(m2, nu1_star_closed(m2), color=dark, lw=1.1, ls=ls_two, zorder=4)
        inside = [pt for pt in sampled if win[0] <= pt[0] <= win[1] and win[2] <= pt[1] <= win[3]]
        for mu0, v, chaotic in inside:      # the band at each sampled parameter: from the all-orbit curve to nu1*
            bx.plot([mu0, mu0], [float(nu1_all_orbit(mu0)), v], color="#6e6e6e", lw=0.6, ls=(0, (1, 1.6)), zorder=3)
        bx.plot([4.0, 4.0], [win[2], win[3]], color=edge, lw=2.4, zorder=5, clip_on=False, solid_capstyle="butt")
        markers(bx, inside, 3.9)
        bx.set_xlim(win[0], win[1])
        bx.set_ylim(win[2], win[3])
        bx.set_xticks([3.0, 3.2, TWO_CYCLE_END, 3.6, 3.8, 4.0])
        bx.set_xticklabels(["3", "3.2", r"$1{+}\sqrt{6}$", "3.6", "3.8", "4"])
        bx.set_yticks([4.0 / 3.0, 1.5, JURY_EDGE_AT_4, 2.0, 2.2])
        bx.set_yticklabels(["4/3", "3/2", r"$8\sqrt{3}{-}12$", "2", "2.2"])
        bx.set_xlabel(r"prey growth rate $\mu_0$", labelpad=1.5)
        bx.tick_params(length=2.5, width=0.6, pad=2)
        box = dict(boxstyle="square,pad=0.2", fc="white", ec="none")
        bx.text(3.36, 2.09, "the coexistence point is a sink", ha="center", va="center", fontsize=fs, bbox=box, zorder=9)
        xa = 3.33                           # the band, between the all-orbit curve and nu1* (two-cycle range)
        bx.annotate("", xy=(xa, float(nu1_all_orbit(xa)) + 0.004), xytext=(xa, float(nu1_star_closed(xa)) - 0.004),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color=dark, shrinkA=0, shrinkB=0, mutation_scale=5), zorder=6)
        bx.text(xa + 0.022, 0.5 * (float(nu1_all_orbit(xa)) + float(nu1_star_closed(xa))) - 0.012, "band", ha="left",
                va="center", fontsize=fs, zorder=9)
        bx.text(3.50, 1.283, "every interior orbit loses its predator", ha="center", va="center", fontsize=fs, zorder=9)
        basin = (": %s" % FIG8_LABEL_BASIN_THEOREM) if FIG8_LABEL_BASIN_THEOREM else ""
        bx.text(1.028, 0.5, r"$\mu_0=4$" + basin, transform=bx.transAxes, rotation=90, ha="left", va="center", fontsize=fs,
                zorder=9, clip_on=False)
        bx.text(0.0, 1.03, "(b)", transform=bx.transAxes, ha="left", va="bottom", fontsize=8.4, fontweight="bold")

        handles = [
            Line2D([], [], color=dark, lw=1.5, label=r"all-orbit curve $1/\mu_0+1/\nu_1=1$"),
            Line2D([], [], color=dark, lw=1.1, ls=ls_two, label=r"$\nu_1^{\ast}=\mu_0/\sqrt{\mu_0+1}$ (attracting two-cycle)"),
            Line2D([], [], color=dark, lw=1.1, ls=ls_flip, label=r"$\nu_1^{\mathrm{flip}}(\mu_0)$"),
            Line2D([], [], marker="o", ms=3.9, mfc=dark, mec=dark, mew=0.8, ls="none", label=r"$\nu_1^{\ast}$ sampled, chaotic prey"),
            Line2D([], [], marker="o", ms=3.9, mfc="white", mec=dark, mew=0.8, ls="none", label=r"$\nu_1^{\ast}$ sampled, periodic window"),
            Line2D([], [], marker="*", ms=8.5, mfc=dark, mec=dark, ls="none", label=r"$\nu_1^{\ast}(4)=4$ (closed form)"),
        ]
        fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=fs, frameon=False, handlelength=2.5,
                   handletextpad=0.6, columnspacing=1.6, labelspacing=0.45, borderaxespad=0.15, bbox_to_anchor=(0.5, 0.0))
        fig.savefig(out_png, dpi=dpi)
        plt.close(fig)
    return out_png


def make_figures(values: Dict[str, object]) -> List[str]:
    """Write optional diagnostic PNGs to FIGDIR. Returns list of written paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIGDIR, exist_ok=True)
    written: List[str] = []

    # Manuscript Figure 7: the function H of the sub-action inequality, H(x) <= ln(3/4).
    xs = np.linspace(1e-3, 1.0, 1000)
    H = np.log(xs) + (4.0 / 9.0) * (3.0 * xs - 4.0 * xs * xs)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(xs, H, color="tab:blue", lw=1.5, label=r"$H(x)=\ln x+\frac{4}{9}(3x-4x^2)$")
    ax.axhline(LN_THREE_QUARTER, color="tab:red", ls="--", lw=1.0, label=r"$\ln(3/4)=-0.2877$")
    ax.plot([0.75], [LN_THREE_QUARTER], "ko", ms=6)
    ax.annotate(r"$x=3/4$", (0.75, LN_THREE_QUARTER), xytext=(6, -14), textcoords="offset points", fontsize=9)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$H(x)$")
    ax.set_ylim(-3.0, 0.2)
    ax.set_title("Sub-action inequality: $H(x)\\leq\\ln(3/4)$")
    ax.legend(fontsize=8, loc="lower right")   # clear of the maximum at x = 3/4
    p = os.path.join(FIGDIR, "designG_subaction.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    _manifest_append(p)
    plt.close(fig)
    written.append(p)

    # Supplementary Figure S1: equilibrium-measure potential vs r.
    eq = values["equilibrium_measure"]
    rs_on = sorted(float(k) for k in eq["on_support"].keys())
    vals_on = [eq["on_support"][f"{r:.2f}"] for r in rs_on]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(rs_on, vals_on, "o-", ms=5, color="tab:green", label=r"$\int\ln|x-r|\,d\mu_4$ (on support)")
    ax.axhline(-2.0 * LN2, color="tab:red", ls="--", lw=1.0, label=r"$-\ln4=-1.3863$")
    _off_label = r"$\int\ln|x-r|\,d\mu_4$ (off support)"
    for r in ("1.20", "2.00"):
        if r in eq["off_support"]:
            ax.plot([float(r)], [eq["off_support"][r]["quadrature"]], "ks", ms=6, label=_off_label)
            _off_label = None  # label only the first off-support marker (single legend entry)
    ax.set_xlabel("$r$")
    ax.set_ylabel(r"$\int\ln|x-r|\,d\mu_4(x)$")
    ax.set_title("Equilibrium-measure identity: constant $-\\ln4$ on $[0,1]$")
    ax.legend(fontsize=8)
    p = os.path.join(FIGDIR, "designG_equilibrium_potential.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    _manifest_append(p)
    plt.close(fig)
    written.append(p)

    # Manuscript Figure 8: the map of the parameter box (closed forms, and the sampled thresholds of the I_family block).
    p = os.path.join(FIGDIR, "designG_I_family.png")
    _draw_parameter_box(values["I_family"]["family"], p)
    _manifest_append(p)
    written.append(p)

    return written


# =============================================================================
# Driver
# =============================================================================


def _provenance(args: argparse.Namespace) -> Dict[str, object]:
    return {
        "generated_by": "output/codes/" + os.path.basename(__file__),

        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model": "logistic prey axis of Design G at mu0=4 (arcsine measure)",
        "honesty": (
            "Numerics are CONSISTENCY CHECKS, not proofs (the sub-action telescoping argument is the proof). "
            "a.e./SRB != all-orbit. The arcsine density, I(4)=-2ln2 and nu1*(4)=4 hold at mu0=4 only; nu1*(mu0) "
            "has a closed form on an attracting cycle of the prey (mu0/(mu0-1) for 1<mu0<=3, mu0/sqrt(mu0+1) for "
            "3<mu0<=1+sqrt6); the six values at mu0=3.60 to 3.99 are orbit averages, non-monotone (periodic windows)."
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "numba": numba.__version__,
        "output_data_dir": _rel(DATADIR),
        "output_fig_dir": _rel(FIGDIR),
        "no_figs": bool(args.no_figs),
        "protocol": {
            "quadrature_nodes": args.n_nodes,
            "subaction_grid": args.sub_grid,
            "max_period": args.max_period,
            "enumeration": "exact rational angle doubling (fractions.Fraction); every interior orbit evaluated",
            "I_burn_in": args.n_burn,
            "I_retain": args.n_keep,
            "I_n_ic": args.n_ic,
        },
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Design-G ergodic / quadrature / enumeration validation suite (spec 2B)."
    )
    parser.add_argument("--no-figs", action="store_true", help="skip optional diagnostic PNGs (P2)")
    parser.add_argument("--n-nodes", type=int, default=200000, help="Gauss--Chebyshev quadrature nodes")
    parser.add_argument("--sub-grid", type=int, default=200001, help="sub-action grid points")
    parser.add_argument("--max-period", type=int, default=16,
                        help="max least period of the exact-angle periodic-orbit enumeration (every orbit evaluated)")
    parser.add_argument("--n-burn", type=int, default=75000, help="I(mu0) burn-in")
    parser.add_argument("--n-keep", type=int, default=50000, help="I(mu0) retained iterations")
    parser.add_argument("--n-ic", type=int, default=32, help="I(mu0) many-IC count for mu0<4")
    args = parser.parse_args(argv)

    # Guard the CLI domain. max_period < 1 makes the combinatorial count loop empty, so
    # the interior counts (total minus the boundary fixed point) go negative yet still
    # self-assert and "pass". The default period-16 path is unaffected; a non-default
    # 1 <= max_period != 16 is honestly an exact-count-vs-itself consistency check (no
    # independent reference), acceptable but it must be sane.
    if args.max_period < 1:
        parser.error(f"--max-period ({args.max_period}) must be >= 1.")

    os.makedirs(DATADIR, exist_ok=True)

    print("=" * 96)
    print("Design-G ergodic / quadrature / enumeration suite (spec 2B) -- CONSISTENCY CHECKS, NOT PROOFS")
    print("  logistic prey axis at mu0=4 (arcsine measure dmu4 = dx / (pi sqrt(x(1-x))))")
    print("=" * 96)

    all_rows: List[Check] = []
    values: Dict[str, object] = {"provenance": _provenance(args)}

    print("\n[a] Arcsine moments + Jensen gap (eq:arcsine-moments, eq:jensen-gap) ...")
    rows, rec = check_arcsine_moments(args.n_nodes)
    all_rows += rows
    values["arcsine_moments"] = rec

    print("\n[b] Equilibrium-measure identity int ln|x-r| dmu4 = -ln4 (thm:equilibrium-measure) ...")
    rows, rec = check_equilibrium_measure(args.n_nodes)
    all_rows += rows
    values["equilibrium_measure"] = rec

    print("\n[c.1] Sub-action H(x) <= ln(3/4), strict concavity (eq:subaction) ...")
    rows, rec = check_subaction(args.sub_grid)
    all_rows += rows
    values["subaction"] = rec

    print(f"\n[c.2] Exact-angle enumeration of every interior periodic orbit, period<= {args.max_period},"
          f" max avg ln x (rmk:ergodic-numeric) ...")
    rec, rows = enumerate_periodic_orbits(args.max_period)
    all_rows += rows
    values["periodic_enumeration"] = rec
    print(f"    interior orbits enumerated = {rec['interior_orbits_enumerated']} (exact Lyndon/Moebius"
          f" {rec['interior_orbits_exact']}) / interior points = {rec['interior_points_enumerated']}"
          f" (exact {rec['interior_periodic_points_exact']}); evaluated = {rec['interior_orbits_evaluated']};"
          f" coverage_complete = {rec['coverage_complete']}")
    print("    evaluated orbits by least period: "
          + ", ".join(f"{n}:{rec['evaluated_orbits_by_least_period'][n]}" for n in sorted(rec['evaluated_orbits_by_least_period'])))
    print(f"    max avg ln x = {rec['max_orbit_avg_lnx']:.6f} at x0={rec['argmax_orbit_x0']:.6f}"
          f" (angle {rec['argmax_orbit_angle']}, period {rec['argmax_orbit_period']}); target ln(3/4) = {LN_THREE_QUARTER:.6f};"
          f" runner-up {rec['runner_up_avg_lnx']:.6f} (period {rec['runner_up_period']})")

    print("\n[d] I(mu0) family + threshold nu1*(mu0)=exp(-I) (cor:mu0family, rmk:Icurve) ...")
    rows, rec = check_I_family(args.n_burn, args.n_keep, args.n_ic)
    all_rows += rows
    values["I_family"] = rec

    # ---- PASS/FAIL validation table --------------------------------------
    print("\n" + "=" * 96)
    print("VALIDATION TABLE  (Design-G ergodic/quadrature; consistency checks, not proofs)")
    print("-" * 96)
    header = f"{'quantity':<46} {'predicted':>13} {'computed':>13} {'abs-err':>11}  result"
    print(header)
    print("-" * 96)
    n_pass = 0
    for c in all_rows:
        status = "PASS" if c.passed else "FAIL"
        if c.passed:
            n_pass += 1
        print(f"{c.name:<46} {c.predicted:>13.6f} {c.computed:>13.6f} {c.abs_err:>11.3e}  {status}")
    print("-" * 96)
    n_total = len(all_rows)
    print(f"SUMMARY: {n_pass}/{n_total} checks PASS.")
    print("NOTE: the sub-action telescoping argument is the PROOF; the exact-angle enumeration only corroborates it.")
    print("      At mu0 = 4 the invariant density is the arcsine law, which gives I(4) = -2 ln 2 and nu1*(4) = 4 in\n"
          "      closed form; the arcsine law holds at mu0 = 4 only. nu1*(mu0) also has a closed form on an attracting\n"
          "      cycle of the prey: mu0/(mu0 - 1) for 1 < mu0 <= 3 (fixed point; this is the all-orbit curve of\n"
          "      Figure 8) and mu0/sqrt(mu0 + 1) for 3 < mu0 <= 1 + sqrt 6 (two-cycle; the dashed curve of Figure 8).\n"
          "      The six values at mu0 = 3.60 to 3.99 come from orbit averages and are not monotone in mu0.")
    print("=" * 96)

    values["validation_table"] = [c.as_dict() for c in all_rows]
    values["summary"] = {"n_pass": n_pass, "n_total": n_total, "all_pass": bool(n_pass == n_total)}

    # ---- optional figures -------------------------------------------------
    if not args.no_figs:
        try:
            written = make_figures(values)
            values["figures"] = [_rel(w) for w in written]   # bundle-relative
            for p in written:
                print("  wrote figure", p)
        except Exception as exc:
            # The values JSON is still written below, but the run FAILS (see the return): these PNGs are
            # manuscript figures, and a runner that got exit 0 would keep the previous PNG and report success.
            print(f"  [ERROR] figure generation FAILED: {exc!r}")
            values["figures_error"] = repr(exc)
    else:
        print("  (--no-figs: diagnostic PNGs skipped)")

    # ---- JSON deposit -----------------------------------------------------
    out_json = os.path.join(DATADIR, os.path.splitext(os.path.basename(__file__))[0] + "_values.json")
    with open(out_json, "w") as f:
        json.dump(values, f, indent=2)
    print("wrote", out_json)

    return 0 if (n_pass == n_total and "figures_error" not in values) else 1


if __name__ == "__main__":
    raise SystemExit(main())
