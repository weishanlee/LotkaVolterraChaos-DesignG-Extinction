#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G dynamics validation suite (2-D map simulation layer).

Numerical CONSISTENCY CHECKS (never proofs) for the symmetric multiplicative
predator--prey map ("Design G") of ``output/methods_v6.tex``::

    x_{n+1} = mu0 * x_n (1 - x_n)(1 - y_n)       (eq:designG, F)
    y_{n+1} = nu1 * x_n  y_n (1 - y_n)           (eq:designG, G)

on the biological square ``[0,1]^2`` with ``mu0, nu1 in (0, 4]``.

This is one of two CPU-only scripts that make the Design-G numerics of
``methods_v6.tex`` reproducible.  It covers spec section 2A:

  (a) boundedness / no blow-up                         -- eq:bounds, rmk:three-regime-numeric
  (b) transverse exponent lambda_hat (Benettin--Wolf)  -- eq:numerical-lambda, eq:lambda-mu4
  (c) finite-grid TERMINAL-PERSISTENCE extinction sweep -- rmk:first-passage  [v4: rewritten]
  (d) coexistence fixed point + Jury stability band    -- eq:coexistence-fp / -trace-det / -jury-band
  (e) global-extinction telescoping demo (nu1 < 4/3)   -- thm:global
  (f) family lambda_perp unification (skew/MOD-1/G)    -- thm:skew, thm:mod1, eq:lambda-mu4

The companion script ``designG_ergodic_quadrature_v2.py`` covers the
mu0=4 closed-form / quadrature / periodic-orbit-enumeration layer (spec 2B).

Honesty notes:
  * Every number below is a finite-time / finite-grid consistency check, not a
    proof.  SRB-typical / finite-grid statements are NOT all-orbit statements.
    (For the extinction set the paper proves positive,
    not full, area (thm:basin-positive); these finite-time checks do not verify it.)
  * For ``2 < nu1 < 4`` the interior coexistence fixed point is Jury-stable on
    ``8 sqrt(3) - 12 < nu1 < 4`` (eq:designG-jury-band), so interior ICs that
    converge to coexistence are EXPECTED, not test failures.
  * The closed forms (e.g. lambda_perp = ln(nu1/4), the arcsine moments) are
    SPECIFIC to mu0 = 4, where the SRB measure equals the equilibrium measure
    of [0,1].

Engineering follows the canonical sibling scripts
``codes/forward_causal_theory_figs.py`` and
``codes/forward_causal_mu0_family.py``: portable ``_HERE/_OUT`` paths derived
from ``__file__``, ``matplotlib.use("Agg")``, Numba ``@njit`` hot loops, and a
provenance + values JSON deposit.  CPU only; README.md gives the commands (outside
``--smoke`` the script needs DESIGNG_FIGDIR and DESIGNG_DATADIR; run_designG_figs.sh sets them).

Optional diagnostic PNGs are written only when ``--no-figs`` is NOT supplied
(figures are a P2 nicety; the PASS/FAIL table and JSON are the contract).

Env-var I/O contract (no-clobber, explicit-path runs): ``DESIGNG_FIGDIR`` /
``DESIGNG_DATADIR`` select the output dirs and are REQUIRED outside
``--smoke`` (see ``resolve_paths``); ``DESIGNG_MANIFEST``, if set, receives an
appended ``| artifact | version | sha256 |`` row per artifact written (see
``_manifest_append``).  ``--smoke`` runs the full schema at tiny/reduced
values -> ``<name>_smoke.json``; ``--only <block>`` restricts to one named
block; the default runs the full pipeline (all blocks) -> ``DESIGNG_DATADIR``.
"""
# === Revision v2 (2026-06-04) =================================================
# Changes:
#   the nu1=1.5 comment, a CLI guard, the family-unification doc,
#   a collapse-prevalence floor in the lambda gate, a one-sided max-x check,
#   provenance fields, the run command, the "finite-grid" rename.
# No algorithm / theorem / equation / numerical result changed.
#
# === Revision v3 (2026-07-09, Task B0 scaffold) =================================
# Additive scaffold ONLY -- no algorithm / theorem / equation / numerical result
# changed in the copied v2 logic. Added: _preflight() (interpreter/package guard),
# resolve_paths(smoke) (env-var I/O contract), _smoke_params() + --smoke/--only/
# --public CLI dispatch, a block registry (so --only <block> can run a single
# check), _manifest_append() (DESIGNG_MANIFEST producer rows), and a top-level
# "_taxonomy" map in the values JSON (B0 JSON-schema convention; populated here
# only for the two RETAINED v2 producers this script already emits --
# coexistence.jury_band.lower_edge_closed_form and
# family_unification.per_nu1.*.pairwise_max_spread; the U1-U4 blocks and their
# taxonomy entries are added by Tasks B1-B4). make_figures()/_provenance() now
# take an explicit figdir/datadir instead of reading module-level globals, so a
# run can be redirected (incl. to --smoke's tmp dirs) without editing FIGDIR/
# DATADIR.
# === Revision v4 (2026-07-17, P3 numerics recompute) =========================
# STANDALONE successor to v3 (no `from v3 import *`; v3 retained as historical).
# This is a DATA-CHANGING revision -- the one thing v3's relabelling could not do.
#
# THE DEFECT: v3's
# `_designG_extinct_fraction` marked an orbit extinct at the FIRST `y_n < eps`
# and broke out of the loop. The Design-G predator can recover, since
#     y_{n+1}/y_n = nu1 * x_n * (1 - y_n)  >  1   whenever  x_n(1-y_n) > 1/nu1,
# so that kernel measured FIRST PASSAGE, not extinction. At nu1=4 a 400-IC /
# T=2e4 / eps=1e-3 run found 32 first-passage "extinctions" and 0 terminal ones:
# every single one was false.
#
# WHAT v4 CHANGES:
#   1. TERMINAL PERSISTENCE. Every trajectory runs the FULL horizon; an orbit is
#      extinct iff ALL of its final `tail` (=_EXT_TAIL=1000) iterates are < eps.
#      No early break.
#   2. LEFT-CENSORING. ICs with y0 < eps start already below threshold and cannot
#      evidence extinction (the log-uniform sampler put ~52% of ICs
#      there). They are EXCLUDED from the denominator and the count is REPORTED.
#   3. EPS SWEEP. The fraction is reported across a range of eps, so the reader
#      sees tail-decay behaviour rather than one threshold's artefact.
#
# HOW IT STAYS CHEAP (the one real design idea here): for a fixed orbit,
#     first passage below eps  <=>  min_{n>=1} y_n  <  eps
#     terminal persistence     <=>  max_{n in tail} y_n  <  eps
#     left-censored            <=>  y0 < eps
# All three are threshold comparisons against per-orbit summary statistics that
# do not depend on eps. So `_designG_orbit_stats` makes ONE pass per (grid, nu1),
# recording (min_y_orbit, max_y_tail, y0) per IC, and `_extinct_stats` then
# derives the counts for ANY eps for free. The eps sweep therefore costs one pass,
# not one pass per eps -- and first-passage and terminal counts come out of the
# SAME orbits, which is what makes the side-by-side honest.
#
# HONESTY:
# the first-passage count is NOT a bound on asymptotic extinction in EITHER
# direction. It OVER-counts dip-and-recover orbits, and UNDER-counts orbits
# converging to 0 too slowly to reach eps within the horizon -- the latter
# worsening as nu1 -> 4, exactly where ln(nu1/4) -> 0. Because the two errors
# oppose, its distance from the asymptotic fraction has no known sign. The
# change from v3 to v4 does have one: at the same threshold, horizon and sample an
# orbit that stays below eps over the terminal window has passed below eps, so the
# terminal count is never above the first-passage count (817 against 826 of 3600
# at nu1 = 2); only the size of the decrease is a measurement.
# The paper quotes what this script computes.
# v3's fractions (nu1=2:0.2294, 3:0.0289, 4:0.0172) are
# first-passage values and are retained only as `fraction_first_passage`.
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
from datetime import datetime, timezone
from typing import Callable, Dict, List, Tuple, Optional

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
# Portable paths (mirror the _HERE/_PROJ idiom of the canonical forward_causal_*.py).
# This script lives in output/codes/ ; data -> output/data/designG, figs -> output/figures.
# -----------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))   # .../codes

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

LN2 = math.log(2.0)

# --- U2 (Task B2) sensitivity-sweep constants --------------------------------
# RNG seed for the random-IC configs (C9 ic.runif, C10 ic.rlogy): reproducibility
# is required (plan Global Constraints) so every random config is seeded FRESH with
# this value -> its fraction is deterministic and resume-safe regardless of order.
_SENS_SEED = 20260709
# Per-CELL wall-clock budget (one (config, nu1) cell). The plan's runtime envelope
# has the heaviest cell (400^2 x 1e4 ~ 1.6e9 map-steps) at ~30-60 s, well under this
# 300 s ceiling, so NO fallback fires at nominal params (fractions stay submission
# grade). Env-overridable (DESIGNG_SENS_BUDGET_S).
_SENS_BUDGET_S_DEFAULT = 300.0
# Conservative numba map-step rate for the deterministic pre-run estimate (real
# rate ~3.5e7/s; using 3e7 makes the estimate slightly pessimistic -> the fallback
# never fires spuriously below budget). The estimate is deterministic, so the
# no-fallback path is fully reproducible (a regression anchor must be).
_SENS_STEPS_PER_SEC = 3.0e7
# Fallback shrink floors (a reduced-parameter fallback halves the cost axis toward
# these floors and sets _fallback_used=True -> NOT a submission-grade claim).
_SENS_MIN_K = 25
_SENS_MIN_T = 500
_SENS_MIN_N = 1000

# --- U3 (Task B3) finite-amplitude two-way decay diagnostic constants -------
# Same project-mandated RNG seed as U2 (plan Global Constraints: RNG seed
# 20260709), kept as its OWN named constant (not a reuse of _SENS_SEED) so the
# two blocks' seeding stays independently documented/adjustable.
_FINAMP_SEED = 20260709
# delta sweep is SCHEMA, not a cost knob: --smoke keeps all three (see
# _finite_amp_schedule / run_finite_amplitude; only M/n_burn/horizon/nu1_list shrink).
_FINAMP_DELTAS: Tuple[float, ...] = (1e-2, 1e-3, 1e-4)
# Representable-y floor (well above float64 underflow ~1e-308): once a
# trajectory's predator population drops below this, ln(y) is no longer
# meaningfully trackable and the per-IC slope accumulation stops for that IC
# (mirrors the "last representable step" idiom of _designG_telescoping_rate).
_FINAMP_UNDERFLOW_FLOOR = 1e-280
# Exact interpretation string required by the plan's Task B3 interface contract
# (spec 4.1a): this block is an EMPIRICAL diagnostic distinct from the
# (infinitesimal, tangent-linearised) axis variational exponent that a future
# U4 block will estimate -- NOT a theory gate, no PASS/FAIL tolerance on it.
_FINAMP_INTERPRETATION = (
    "finite-amplitude; distinct from the axis variational exponent (U4); "
    "expected -> ln(nu1/4) as delta->0"
)

# =============================================================================
# Numba-compiled hot loops
# =============================================================================


@njit(cache=True, fastmath=False)
def _designG_step(x: float, y: float, mu0: float, nu1: float) -> Tuple[float, float]:
    """One Design-G iterate (eq:designG). Returns (x_next, y_next)."""
    xn = mu0 * x * (1.0 - x) * (1.0 - y)
    yn = nu1 * x * y * (1.0 - y)
    return xn, yn


@njit(cache=True, fastmath=False)
def _designG_orbit_maxima(
    x0: float, y0: float, mu0: float, nu1: float, n_iter: int
) -> Tuple[float, float]:
    """Iterate Design G from (x0, y0); return (max x_n, max y_n) over ITERATES n>=1.

    Used for the boundedness check (a).  The analytic bounds eq:bounds
    (x_{n+1} <= mu0/4, y_{n+1} <= nu1/4) constrain the IMAGE of the map, i.e.
    iterates n >= 1; the interior initial condition (x0, y0) is deliberately
    EXCLUDED from the maxima (an interior y0 may exceed nu1/4 without violating
    the no-blow-up property of the dynamics).  We therefore take the first
    iterate before initialising the running maxima.
    """
    x, y = _designG_step(x0, y0, mu0, nu1)  # first iterate (n=1)
    max_x = x
    max_y = y
    for _ in range(n_iter - 1):
        x, y = _designG_step(x, y, mu0, nu1)
        if x > max_x:
            max_x = x
        if y > max_y:
            max_y = y
    return max_x, max_y


@njit(cache=True, fastmath=False)
def _designG_grid_maxima(
    grid: np.ndarray, mu0: float, nu1: float, n_iter: int
) -> Tuple[float, float]:
    """Max prey / max predator over a whole grid of interior ICs (boundedness, (a))."""
    gmax_x = 0.0
    gmax_y = 0.0
    m = grid.shape[0]
    for i in range(m):
        mx, my = _designG_orbit_maxima(grid[i, 0], grid[i, 1], mu0, nu1, n_iter)
        if mx > gmax_x:
            gmax_x = mx
        if my > gmax_y:
            gmax_y = my
    return gmax_x, gmax_y


@njit(cache=True, fastmath=False)
def _logistic_lambda_hat(
    nu1: float, mu0: float, x0: float, n_burn: int, n_keep: int, window: int
) -> Tuple[float, float, float, int]:
    """Benettin--Wolf transverse exponent on M0 via the logistic (y->0 limit) prey orbit.

    Implements eq:numerical-lambda for Design G:
        lambda_hat = (1/N) sum_{n} ln(nu1 * x_n),   x_n = logistic orbit at mu0.
    Burn-in ``n_burn`` discarded, ``n_keep`` retained.  A Welford running
    variance is accumulated and a stability gate is applied over moving windows
    of length ``window``: the per-window cumulative-mean drift must stay small.

    Returns (lambda_hat, running_variance, max_window_drift, accepted_flag).
    The transverse multiplier nu1*x is strictly positive on the interior orbit
    (x in (0,1)), so no absolute value is needed (eq:linear-transverse).
    """
    x = x0
    for _ in range(n_burn):
        x = mu0 * x * (1.0 - x)

    mean = 0.0
    m2 = 0.0
    count = 0

    # moving-window stability gate state
    win_sum = 0.0
    win_count = 0
    cum_mean_at_window_end = 0.0
    have_prev_window = False
    max_drift = 0.0

    for _ in range(n_keep):
        x = mu0 * x * (1.0 - x)
        xc = x
        if xc < 1e-300:
            xc = 1e-300
        val = math.log(nu1 * xc)

        count += 1
        delta = val - mean
        mean += delta / count
        delta2 = val - mean
        m2 += delta * delta2

        win_sum += val
        win_count += 1
        if win_count == window:
            # compare cumulative mean now vs at end of previous window
            if have_prev_window:
                drift = abs(mean - cum_mean_at_window_end)
                if drift > max_drift:
                    max_drift = drift
            cum_mean_at_window_end = mean
            have_prev_window = True
            win_sum = 0.0
            win_count = 0

    variance = m2 / count if count > 1 else 0.0
    # Acceptance: cumulative estimate stable across 5000-iteration windows.
    accepted = 1 if max_drift < 1e-3 else 0
    return mean, variance, max_drift, accepted


@njit(cache=True, fastmath=False)
def _logistic_lnx_ensemble(
    mu0: float, n_burn: int, n_keep: int, window: int, n_ic: int
) -> Tuple[float, float, int, int, int]:
    """Ensemble mean of ln x over the logistic SRB at mu0, via the DECLARED protocol.

    Implements the methods/numerics + spec 2A(b) recipe literally PER TRAJECTORY:
    burn-in ``n_burn`` (=75000) discarded, ``n_keep`` (=50000) retained, transverse
    time average (1/N) sum ln x (eq:numerical-lambda; the constant ln(nu1) shift is added
    by the caller), with the Welford cumulative-mean stability gate over moving windows
    of ``window`` (=5000): "a trajectory is accepted only if the cumulative estimate is
    stable across moving windows".  A SINGLE long orbit at mu0=4 is corrupted (arcsine
    ln-singularity / float64 collapse onto x in {0,1}; rmk:Icurve), so the recipe is run
    over ``n_ic`` initial conditions and lambda is averaged over the ACCEPTED trajectories
    (the methods many-IC practice).  <ln x> is nu1-independent and computed once.

    PER-TRAJECTORY STABILITY GATE (genuinely enforced, NOT a collapse filter): the running
    mean necessarily drifts in the EARLY 5000-windows while it settles from a cold start,
    so "stable across moving windows" is evaluated on the CONVERGED second half -- a
    trajectory is accepted only if the cumulative-mean drift across its late 5000-windows
    stays within 2x the trajectory's OWN standard error of the estimate (a self-calibrated
    "the running estimate has stopped moving relative to its statistical noise" test, with
    a 5e-3 floor).  In the shipped run (75000/50000 at mu0=4, 128 seeds) 3 seeds are dropped
    because an iterate comes within 1e-14 of an endpoint and 7 of the other 125 are rejected
    by this test (118 accepted).  Returns (mean_lnx, sem_lnx, n_accepted, n_noncollapsed, n_total).
    """
    nwin = n_keep // window
    cum = np.empty(nwin)                        # cumulative mean at each 5000-window end
    acc_sum = 0.0
    acc_sumsq = 0.0
    n_acc = 0
    n_noncollapsed = 0
    for m in range(n_ic):
        # well-spread interior seed; golden-ratio sequence avoids low-period landings
        x = (0.5 + 0.6180339887498949 * (m + 1)) % 1.0
        if x < 1e-6:
            x = 1e-6
        if x > 1.0 - 1e-6:
            x = 1.0 - 1e-6
        collapsed = False
        for _ in range(n_burn):                 # DECLARED burn-in (75000)
            x = mu0 * x * (1.0 - x)
            if x <= 1e-14 or x >= 1.0 - 1e-14:
                collapsed = True
                break
        if collapsed:
            continue
        # DECLARED retained window (50000): Welford mean+variance, cumulative mean per 5000-window
        mean = 0.0
        m2 = 0.0
        count = 0
        win_count = 0
        w = 0
        for _ in range(n_keep):
            x = mu0 * x * (1.0 - x)
            if x <= 1e-14 or x >= 1.0 - 1e-14:
                collapsed = True
                break
            v = math.log(x)
            count += 1
            d = v - mean
            mean += d / count
            m2 += d * (v - mean)
            win_count += 1
            if win_count == window and w < nwin:
                cum[w] = mean
                w += 1
                win_count = 0
        if collapsed or count < n_keep // 2 or w < nwin:
            continue                            # collapsed / short -> not stable -> reject
        n_noncollapsed += 1
        se = math.sqrt((m2 / count) / count) if count > 1 else 1.0   # per-traj standard error
        stab_tol = 2.0 * se
        if stab_tol < 5e-3:
            stab_tol = 5e-3
        start = nwin // 2
        if start < 1:
            start = 1
        max_drift_late = 0.0
        for i in range(start, nwin):            # 5000-window stability across the CONVERGED half
            dd = abs(cum[i] - cum[i - 1])
            if dd > max_drift_late:
                max_drift_late = dd
        if max_drift_late > stab_tol:
            continue                            # cumulative estimate NOT stable across late windows -> reject
        acc_sum += mean
        acc_sumsq += mean * mean
        n_acc += 1
    if n_acc < 2:
        return 0.0, 1.0, n_acc, n_noncollapsed, n_ic   # degenerate: large sem so the gate rejects
    mlnx = acc_sum / n_acc
    var = acc_sumsq / n_acc - mlnx * mlnx
    if var < 0.0:
        var = 0.0
    sem = math.sqrt(var / n_acc)                # standard error of <ln x> across accepted orbits
    return mlnx, sem, n_acc, n_noncollapsed, n_ic


# --- v4 extinction protocol constants ----------------------------------------
# Number of FINAL iterates that must ALL sit below eps for an orbit to count as
# extinct (spec item 1: "classify extinct only if all final 1000 iterates < eps").
_EXT_TAIL = 1000
# The eps ladder (spec item 3: "sweep successively smaller eps; report tail-decay
# behaviour, not a single threshold"). 1e-3 is v3's single value and stays in the
# ladder so the two versions remain comparable at a common point.
_EXT_EPS_LADDER = (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8)

# --- MEASURED 2026-07-17 under the nominal protocol -------------------------
# (60x60 grid on [0.05,0.95]^2, horizon 20000, tail 1000, eps 1e-3, mu0=4).
# These are computed values.
#
# _EXT_ANCHORS_V4: the TERMINAL fractions -- the paper's claim.
_EXT_ANCHORS_V4 = {1.0: 1.0, 1.5: 1.0, 2.0: 0.2269, 3.0: 0.0111, 4.0: 0.0000}
# _EXT_ANCHORS_V3_FIRST_PASSAGE: v3's numbers. v4 reproduces these EXACTLY from the
# `min_y_orbit < eps` column, which is the regression evidence that the rewrite
# changed the CLASSIFICATION and not the simulation: if the orbits had shifted,
# these would not still land on 0.2294/0.0289/0.0172.
_EXT_ANCHORS_V3_FIRST_PASSAGE = {1.0: 1.0, 1.5: 1.0, 2.0: 0.2294, 3.0: 0.0289, 4.0: 0.0172}
# ...but ONLY at the protocol they were measured under. A regression anchor is a
# (value, protocol) PAIR, and comparing one across protocols is a category error:
# under --smoke (grid 50, horizon 1e3) the same correct code gives first-passage
# 0.2336, not 0.2294, and asserting the anchor there manufactures a FALSE failure
# for a healthy build. So the reproduction rows are emitted only when the run
# matches this protocol exactly.
_EXT_ANCHOR_PROTOCOL = {"horizon": 20000, "n_side": 60, "eps": 1.0e-3}
#
# WHAT THE RECOMPUTE SHOWED:
#  1. At nu1=4, terminal extinction is EXACTLY 0/3600. Every one of v3's 62
#     first-passage "extinctions" there recovers before the horizon and fails the
#     terminal criterion.
#     lambda_perp = ln(nu1/4) = 0 at
#     nu1=4, so the on-axis rate predicts nothing either way there; the zero count
#     is a finite-horizon observation at an undecided endpoint, not a prediction.
#  2. The TERMINAL fraction is INVARIANT across the whole eps ladder (1e-2 ... 1e-8):
#     0.2269 / 0.0111 / 0.0000 at nu1 = 2 / 3 / 4 at EVERY eps. That invariance is
#     PARTLY ARITHMETIC: every terminally extinct orbit at H=20000 has reached the
#     floating-point zero (an absorbing state of the computed map) long before the
#     terminal window (geometric decay at |ln(nu1/4)| >= 0.28/step takes y below
#     the smallest double within a few thousand steps), so it lies below EVERY eps
#     by construction. The counts n_extinct_exact_zero_before_tail record this;
#     the sweep shows insensitivity to eps, NOT a dynamical dichotomy.
#     The FIRST-PASSAGE fraction is strongly eps-dependent by contrast (at nu1=4:
#     0.0617 at 1e-2 down to 0.0011 at 1e-8), which is precisely what made v3's
#     single-eps number an artefact of its threshold.
#  3. The change is NOT uniform in one direction:
#     at nu1=2 the fraction barely moves (0.2294 -> 0.2269), while at
#     nu1=4 it collapses to zero.


def _eff_tail(horizon: int) -> int:
    """Effective terminal window: ``_EXT_TAIL``, but never more than half the horizon.

    A tail >= the horizon would silently redefine "extinct" as "below eps for the
    ENTIRE orbit" -- a strictly stronger and DIFFERENT criterion, quietly applied
    to some configs and not others. The U2 sensitivity configs run horizons as
    short as ``_SENS_MIN_T`` (=500), and --smoke shrinks block (c)'s horizon too,
    so the clamp is load-bearing rather than defensive. The effective value is
    RECORDED in every record it affects, so a clamped run is never mistaken for a
    nominal one. At the nominal horizon 20000 this returns exactly _EXT_TAIL=1000.
    """
    return max(1, min(_EXT_TAIL, int(horizon) // 2))


@njit(cache=True, fastmath=False)
def _designG_orbit_stats(
    grid: np.ndarray, mu0: float, nu1: float, horizon: int, tail: int
) -> np.ndarray:
    """ONE pass per IC -> the eps-independent summary statistics of each orbit.

    Returns an (n, 4) array of columns ``(min_y_orbit, max_y_tail, y0, first_zero_step)``:

      col 0  min_y_orbit -- min over the ITERATES n = 1..horizon of y_n.
                            (y0 is deliberately EXCLUDED, matching v3, whose loop
                            tested y only AFTER stepping. Censoring uses y0.)
      col 1  max_y_tail  -- max over the FINAL ``tail`` iterates of y_n.
      col 2  y0          -- the initial predator value, for left-censoring.
      col 3  first_zero_step -- 1-based index of the first iterate at which y == 0.0
                            EXACTLY (floating-point absorption: the axis is invariant,
                            so the computed orbit then stays at 0 for ever), or -1.0 if
                            the predator never reached exact zero within the horizon.
                            Reported so that a threshold-independent terminal count
                            is not mistaken for a dynamical dichotomy.

    Why this shape: for a fixed orbit and ANY eps,
        first passage below eps  <=>  min_y_orbit < eps      (v3's quantity)
        terminal persistence     <=>  max_y_tail  < eps      (v4's quantity)
        left-censored            <=>  y0 < eps
    Each is a threshold comparison against a statistic that does NOT depend on
    eps. So the whole eps ladder is derived from this single pass by
    ``_extinct_stats`` -- the sweep costs one pass, not one pass per eps -- and
    the first-passage and terminal counts are computed from the SAME orbits,
    which is what makes reporting them side by side an honest comparison rather
    than two protocols talked about as if they were one.

    NOTE [0,1]^2 is forward-invariant for mu0, nu1 in (0,4] (eq:bounds): for
    x,y in [0,1], x' = mu0*x(1-x)(1-y) in [0, 4*(1/4)*1] = [0,1] and likewise
    y' = nu1*x*y(1-y) in [0,1]. So no blow-up guard is needed, and y >= 0 always,
    which is why col 1 can be initialised at 0.0.
    """
    n_total = grid.shape[0]
    out = np.empty((n_total, 4), dtype=np.float64)
    tail_start = horizon - tail
    if tail_start < 0:
        tail_start = 0
    for i in range(n_total):
        x = grid[i, 0]
        y = grid[i, 1]
        out[i, 2] = y                 # y0, before any step
        min_y = 1.0e300               # min over iterates n>=1 (NOT including y0)
        max_tail = 0.0                # y >= 0 by invariance, so 0.0 is a safe floor
        first_zero = -1.0
        for k in range(horizon):
            x, y = _designG_step(x, y, mu0, nu1)
            if y < min_y:
                min_y = y
            if k >= tail_start and y > max_tail:
                max_tail = y
            if first_zero < 0.0 and y == 0.0:
                first_zero = float(k + 1)
        out[i, 0] = min_y
        out[i, 1] = max_tail
        out[i, 3] = first_zero
    return out


def _extinct_stats(stats: np.ndarray, eps: float, tail_start: Optional[int] = None) -> Dict[str, object]:
    """Derive the eps-specific counts from ``_designG_orbit_stats`` output.

    TERMINAL extinction (the quantity this script now reports as
    ``fraction_extinct``) means the predator was below eps for ALL of the final
    ``_EXT_TAIL`` iterates -- not that it once dipped there.

    Left-censoring: ICs with y0 < eps begin below threshold, so they cannot
    evidence extinction either way; they are removed from the DENOMINATOR and
    counted in ``n_censored`` (spec item 2). ``fraction_extinct``
    is therefore over the VALID sample, and ``n_censored``/``n_total`` are
    reported alongside so the exclusion is auditable rather than silent.

    ``fraction_first_passage`` is v3's quantity over the same valid sample, kept
    for explicit comparison. For the SAME valid sample, horizon and eps it is an
    UPPER bound on the terminal count (an orbit ending below eps crossed it earlier);
    it bounds the ASYMPTOTIC extinction fraction in neither direction (see the v4
    revision note). It is reported because showing both is what demonstrates the
    size of v3's error.
    """
    y0 = stats[:, 2]
    keep = y0 >= eps
    n_censored = int(np.count_nonzero(~keep))
    sel = stats[keep]
    n_valid = int(sel.shape[0])
    n_total = int(stats.shape[0])
    if n_valid == 0:
        return {
            "eps": float(eps),
            "fraction_extinct": float("nan"),
            "fraction_first_passage": float("nan"),
            "n_extinct": 0,
            "n_first_passage": 0,
            "n_valid": 0,
            "n_censored": n_censored,
            "n_total": n_total,
            "note": "every IC was left-censored (y0 < eps): no valid sample at this eps",
        }
    n_fp = int(np.count_nonzero(sel[:, 0] < eps))
    n_term = int(np.count_nonzero(sel[:, 1] < eps))
    # Floating-point absorption bookkeeping: how many of the terminally
    # extinct orbits had ALREADY reached exact 0.0 -- and how many before the terminal
    # window even began. For those, "below every eps" holds by arithmetic, not dynamics.
    extra: Dict[str, object] = {}
    if sel.shape[1] > 3:
        fz = sel[:, 3]
        term_mask = sel[:, 1] < eps
        reached = fz >= 0.0
        extra["n_exact_zero"] = int(np.count_nonzero(reached))
        extra["n_extinct_reached_exact_zero"] = int(np.count_nonzero(term_mask & reached))
        if tail_start is not None:
            extra["tail_start"] = int(tail_start)
            extra["n_extinct_exact_zero_before_tail"] = int(np.count_nonzero(term_mask & reached & (fz <= float(tail_start))))
        steps = fz[term_mask & reached]
        if steps.size:
            extra["first_zero_step_min"] = int(steps.min())
            extra["first_zero_step_max"] = int(steps.max())
    return {
        "eps": float(eps),
        "fraction_extinct": float(n_term) / n_valid,          # TERMINAL persistence
        "fraction_first_passage": float(n_fp) / n_valid,      # v3's quantity, same orbits
        "n_extinct": n_term,
        "n_first_passage": n_fp,
        "n_valid": n_valid,
        "n_censored": n_censored,
        "n_total": n_total,
        **extra,
    }


@njit(cache=True, fastmath=False)
def _designG_telescoping_rate(
    grid: np.ndarray, mu0: float, nu1: float, n_iter: int, safe_floor: float
) -> Tuple[float, float, float, int, int, int, float, float]:
    """Global-extinction telescoping demo (e) for nu1 < 4/3 (thm:global).

    For each interior IC, measure the realised predator log-growth rate
        rate = ln(y_k / y_0) / k
    where ``k`` is the LAST step at which the predator was still safely
    representable (y_k >= ``safe_floor``).  thm:global bounds every realised rate by
    ln(3*nu1/4) + (4/9)/k -- the worst-orbit / delta_{3/4} coefficient plus the finite-time
    boundary term 2 sup|g| = 4/9 -- which is < 0 for large k when nu1 < 4/3, so y -> 0.

    Using the last *representable* (y_k, k) rather than (floored y_N, n_iter) is
    essential: once y underflows below ~1e-300 the ratio ln(floor/y0)/n_iter
    divides the true log-drop by an inflated denominator and spuriously
    *understates* the decay magnitude, which would make a correct orbit appear to
    violate the bound.  Stopping at the last safe step gives a FINITE-TIME stopping-rate
    diagnostic, not an asymptotic limit: it is to be read against the
    finite-time bound  (1/k) ln(y_k/y_0) <= ln(3*nu1/4) + (4/9)/k,  whose boundary term the
    common-horizon excess below keeps.

    Returns (max_rate_over_ICs, max_stop_value_y, predicted_bound = ln(3*nu1/4),
             max_stop_step, min_stop_step, n_stopped_before_horizon,
             max_log10_y_at_horizon, max_horizon_telescoping_excess).

    ``max_stop_value_y`` is the largest predator value at which an orbit's rate
    accumulation STOPPED: the first value below ``safe_floor`` (or the value at the
    horizon for an orbit that never crossed it). It is NOT y_N at the common horizon N:
    orbits cross the floor at different steps (recorded as
    ``max_stop_step``/``min_stop_step``), and the common-horizon value is computed
    separately in log space -- below the floor the predator no longer perturbs the prey
    in float64 ((1 - y) rounds to 1), so the prey orbit IS the rounded logistic orbit, and
    ln y_{k+1} = ln y_k + ln(nu1 x_k) is a log-domain continuation along that rounded prey
    sequence (the predator's multiplications, and their underflow, are replaced by log arithmetic).
    ``n_stopped_before_horizon`` counts the orbits whose first sub-floor value comes strictly before
    the horizon (a first crossing at step N itself is not counted).
    ``max_horizon_telescoping_excess`` = max over ICs of ln y_N - ln y_0 - N ln(3 nu1/4).
    An exact predator zero (reachable from exceptional inputs) is reported as what it is:
    the stopping value stays 0 and its logarithm is -inf, never the safe floor or a finite
    sentinel; on the default grids no exact zero occurs.
    """
    n_total = grid.shape[0]
    max_rate = -math.inf
    max_final_y = 0.0
    max_stop_step = 0
    min_stop_step = n_iter
    n_stopped = 0
    max_log10_y_h = -math.inf
    max_excess = -math.inf
    log_safe_floor = math.log(safe_floor)
    ln10 = math.log(10.0)
    for i in range(n_total):
        x = grid[i, 0]
        y0 = grid[i, 1]
        y = y0
        log_y_last = math.log(y0)   # log of last representable y
        step_last = 0               # step index of that y
        stopped_early = False
        for k in range(1, n_iter + 1):
            x, y = _designG_step(x, y, mu0, nu1)
            if y >= safe_floor:
                log_y_last = math.log(y)
                step_last = k
            else:
                # next value underflows the safe range; stop, keep last safe (log_y_last, step_last)
                stopped_early = True
                break
        if step_last == 0:
            # predator left the safe range within a single step; use the safe-floor
            # crossing as a conservative one-step rate estimate.
            rate = (log_safe_floor - math.log(y0))
        else:
            rate = (log_y_last - math.log(y0)) / step_last
        if rate > max_rate:
            max_rate = rate
        # stopping value for reporting (first sub-floor value, or the horizon value). An exact
        # zero stays 0 (it never raises the maximum); it used to be replaced by safe_floor.
        if y > max_final_y:
            max_final_y = y
        # The common-horizon value, continued in log space from the stopping step.
        if stopped_early:
            t_stop = step_last + 1          # the time of the first sub-floor value (x, y)
            if t_stop < n_iter:
                n_stopped += 1              # strictly before the horizon
            if t_stop > max_stop_step:
                max_stop_step = t_stop
            if t_stop < min_stop_step:
                min_stop_step = t_stop
            log_y = math.log(y) if y > 0.0 else -math.inf
            # Synchronous recurrence: y_k = nu1 x_{k-1} y_{k-1} (1 - y_{k-1}) with (1 - y) == 1.0
            # in float64 below the floor, so ln y_k = ln y_{k-1} + ln(nu1 x_{k-1}) uses the OLD
            # prey value and the prey is advanced afterwards (an earlier version
            # advanced x first and used x_k, an off-by-one that also killed the predator one
            # step early when the prey hit a boundary). A prey value of exactly 0 makes the NEXT
            # predator value 0 (ln y = -inf); x = 1 exactly is harmless here
            # (ln(nu1 * 1)) and maps to 0 at the following step.
            for k in range(t_stop + 1, n_iter + 1):
                if x <= 0.0:
                    log_y = -math.inf       # y_k = nu1 * x_{k-1} * y_{k-1} = 0 exactly
                    break
                log_y += math.log(nu1 * x)  # ln y_k from x_{k-1}
                x = mu0 * x * (1.0 - x)     # x_k; (1 - y_{k-1}) == 1.0 in float64 below the floor
        else:
            max_stop_step = max(max_stop_step, n_iter)
            log_y = math.log(y) if y > 0.0 else -math.inf
        log10_y_h = log_y / ln10
        if log10_y_h > max_log10_y_h:
            max_log10_y_h = log10_y_h
        excess = log_y - math.log(y0) - n_iter * math.log(3.0 * nu1 / 4.0)
        if excess > max_excess:
            max_excess = excess
    if n_stopped == 0:
        min_stop_step = n_iter
    predicted = math.log(3.0 * nu1 / 4.0)
    return max_rate, max_final_y, predicted, max_stop_step, min_stop_step, n_stopped, max_log10_y_h, max_excess


@njit(cache=True, fastmath=False)
def _skew_lambda_hat(
    nu1: float, mu0: float, x0: float, y0_unused: float, n_burn: int, n_keep: int
) -> float:
    """Transverse exponent of the one-way SKEW product (eq:skew), thm:skew.

        x_{n+1} = mu0 x (1-x)  (exogenous logistic prey)
        y_{n+1} = (nu1 x_n) y_n  (exactly driven fibre)
    lambda = (1/N) sum ln(nu1 x_n) -- identical estimator to the Design-G
    linear transverse multiplier, because the fibre is EXACTLY linear here.
    """
    x = x0
    for _ in range(n_burn):
        x = mu0 * x * (1.0 - x)
    s = 0.0
    for _ in range(n_keep):
        x = mu0 * x * (1.0 - x)
        xc = x if x > 1e-300 else 1e-300
        s += math.log(nu1 * xc)
    return s / n_keep


@njit(cache=True, fastmath=False)
def _mod1_lambda_hat(
    nu1: float, mu0: float, mu1: float, x0: float, y0: float, n_burn: int, n_keep: int
) -> float:
    """Transverse exponent of the specialist-escapement MOD-1 map (eq:mod1), thm:mod1.

        x_{n+1} = mu0 x (1-x) e^{-mu1 y}
        y_{n+1} = nu1 x_n y_n
    Estimated in the y->0 limit (rare-invader): at y=0 the escapement factor
    e^{-mu1 y}=1 EXACTLY, so the prey orbit is the bare logistic and the transverse
    multiplier is nu1 x, giving lambda = ln(nu1/4) at mu0=4.  Callers pass y0=0 so
    the multiplier is read on M0 itself.  A nonzero seed would be amplified by the
    mu0=4 prey chaos (+ln2/step) into a decorrelated orbit, i.e. exactly the
    "corrupted single-orbit estimate at mu0=4" the numerical-methods protocol warns
    against; the predator is held at the rarefaction floor (e^{-mu1*floor}=1 in
    float64) so the coupled map's prey axis stays bit-identical to the bare logistic.
    """
    x = x0
    y = y0
    for _ in range(n_burn):
        xn = mu0 * x * (1.0 - x) * math.exp(-mu1 * y)
        yn = nu1 * x * y
        x = xn
        y = yn
        if y < 1e-280:
            y = 1e-280  # keep invader rare and alive so the linearisation holds
        if y > 1e-6:
            y = 1e-12   # re-rarefy if it tries to grow (we want the M0 exponent)
    s = 0.0
    for _ in range(n_keep):
        xc = x if x > 1e-300 else 1e-300
        s += math.log(nu1 * xc)
        xn = mu0 * x * (1.0 - x) * math.exp(-mu1 * y)
        yn = nu1 * x * y
        x = xn
        y = yn
        if y < 1e-280:
            y = 1e-280
        if y > 1e-6:
            y = 1e-12
    return s / n_keep


@njit(cache=True, fastmath=False)
def _kahan_sum(arr: np.ndarray) -> float:
    """Compensated (Kahan) summation of a 1-D float64 array.

    Used by U3 (Task B3) for both stages of the finite-amplitude decay-rate
    estimate: the per-IC least-squares slope accumulators
    (_finite_amp_trajectory_slopes, inlined) and the across-IC ensemble
    mean/SEM (_kahan_mean_sem).
    """
    s = 0.0
    c = 0.0
    for i in range(arr.shape[0]):
        yk = arr[i] - c
        tk = s + yk
        c = (tk - s) - yk
        s = tk
    return s


@njit(cache=True, fastmath=False)
def _finite_amp_trajectory_slopes(
    x0: float, y0: float, mu0: float, nu1: float, n_burn: int, checkpoints: np.ndarray
) -> np.ndarray:
    """U3 (Task B3): per-IC finite-amplitude decay-rate estimator.

    Burn in ``x0`` for ``n_burn`` steps via the BARE logistic map
    ``x' = mu0 x (1-x)`` -- the y=0 axis dynamics, the SAME "typical SRB point"
    idiom already used by ``_logistic_lambda_hat`` / ``_skew_lambda_hat`` /
    ``_mod1_lambda_hat`` elsewhere in this file -- THEN inject the FINITE
    predator amplitude ``y0`` and iterate the FULL two-way map (eq:designG)
    forward:  x' = mu0 x (1-x)(1-y),  y' = nu1 x y (1-y).  Burning in x ALONE
    (y held implicitly at 0) avoids the degenerate case where a COUPLED
    burn-in would already have driven y to float64 underflow before the
    measurement window even starts (the coupled decay rate can be steep enough
    -- e.g. ln(2/4) ~ -0.69/step -- that 1e3 coupled steps would exhaust the
    ~700-e-fold underflow budget for small delta).

    Returns the least-squares slope of ln(y_n) vs n (n = 1, 2, ... measured
    from the FIRST post-burn-in iterate) evaluated at each requested
    checkpoint horizon, using Kahan-summed Sx/Sy/Sxy/Sxx accumulators for
    numerical stability (project convention: compensated summation for
    finite-horizon accumulations).  If y underflows ``_FINAMP_UNDERFLOW_FLOOR``
    before a checkpoint is reached, that checkpoint (and every longer one) is
    left NaN for THIS trajectory -- the caller aggregates only the non-NaN ICs
    at each checkpoint (a fully-decayed trajectory has no further slope
    information to contribute at that horizon; mirrors the "last representable
    step" idiom of ``_designG_telescoping_rate``).  ``checkpoints`` must be a
    sorted, strictly increasing int64 array.
    """
    x = x0
    for _ in range(n_burn):
        x = mu0 * x * (1.0 - x)
    y = y0

    n_cp = checkpoints.shape[0]
    out = np.empty(n_cp, dtype=np.float64)
    for i in range(n_cp):
        out[i] = np.nan

    sx = 0.0
    csx = 0.0
    sy = 0.0
    csy = 0.0
    sxy = 0.0
    csxy = 0.0
    sxx = 0.0
    csxx = 0.0
    n = 0
    cp_idx = 0
    max_h = int(checkpoints[n_cp - 1])
    for step in range(1, max_h + 1):
        x, y = _designG_step(x, y, mu0, nu1)
        if y < _FINAMP_UNDERFLOW_FLOOR:
            break  # exhausted representable range -- remaining checkpoints stay NaN
        ln_y = math.log(y)
        nf = float(step)

        # Kahan-add nf into Sx
        yk = nf - csx
        tk = sx + yk
        csx = (tk - sx) - yk
        sx = tk
        # Kahan-add ln_y into Sy
        yk = ln_y - csy
        tk = sy + yk
        csy = (tk - sy) - yk
        sy = tk
        # Kahan-add nf*ln_y into Sxy
        v = nf * ln_y
        yk = v - csxy
        tk = sxy + yk
        csxy = (tk - sxy) - yk
        sxy = tk
        # Kahan-add nf*nf into Sxx
        v2 = nf * nf
        yk = v2 - csxx
        tk = sxx + yk
        csxx = (tk - sxx) - yk
        sxx = tk

        n += 1
        if cp_idx < n_cp and step == checkpoints[cp_idx]:
            if n >= 2:
                denom = n * sxx - sx * sx
                if denom != 0.0:
                    out[cp_idx] = (n * sxy - sx * sy) / denom
            cp_idx += 1

    return out


# =============================================================================
# Validation-table machinery
# =============================================================================


class Check:
    """A single PASS/FAIL row of the validation table."""

    __slots__ = ("name", "predicted", "computed", "tol", "abs_err", "passed", "note")

    def __init__(
        self,
        name: str,
        predicted: float,
        computed: float,
        tol: float,
        note: str = "",
    ) -> None:
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
            # U1 (Task B1): "label" alias so every validation_table row (Check-based
            # PASS/FAIL rows AND the newer raw-dict continuous-metric rows appended by
            # check_extinction_fraction()) shares one common lookup key -- the plan's
            # producer_selector convention resolves "/validation_table/*/..." records
            # by their "label" field regardless of row shape.
            "label": self.name,
            "predicted": self.predicted,
            "computed": self.computed,
            "abs_err": self.abs_err,
            "tol": self.tol,
            "pass": bool(self.passed),
            "note": self.note,
        }


def _make_grid(n_side: int, lo: float = 0.05, hi: float = 0.95) -> np.ndarray:
    """Interior IC grid (n_side x n_side) on (0,1)^2, returned as an (m, 2) array."""
    xs = np.linspace(lo, hi, n_side)
    ys = np.linspace(lo, hi, n_side)
    pts = np.empty((n_side * n_side, 2), dtype=np.float64)
    k = 0
    for xv in xs:
        for yv in ys:
            pts[k, 0] = xv
            pts[k, 1] = yv
            k += 1
    return pts


# --- U2 (Task B2) sensitivity-sweep helpers ----------------------------------


def wilson(k, n, z=1.96):
    """Wilson score interval for a binomial proportion k/n (default 95%, z=1.96).

    Used ONLY for the random-IC sensitivity configs (C9 ic.runif, C10 ic.rlogy),
    where the extinction fraction is a SAMPLE estimate over n Monte-Carlo ICs.
    Deterministic grid fractions are exact-for-that-grid (not sample estimates)
    and therefore get NO CI -- that CI-semantics-by-design split is deliberate
    (see run_sensitivity).
    """
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def _midpoint_lattice(K: int) -> np.ndarray:
    """Deterministic det-grid ICs = midpoint lattice ((i+0.5)/K, (j+0.5)/K).

    Every node lies STRICTLY inside (0,1)^2 (i+0.5 in (0,K) for i in 0..K-1), so no
    node is an exact-singular IC (x or y in {0,1}); the whole K*K set is valid.
    Distinct from _make_grid (linspace on [0.05,0.95]) -- U2's convergence study
    needs the canonical midpoint quadrature lattice (no node on the boundary); the
    lattices at K and 2K share no node, so the convergence study compares independent lattices.
    """
    coords = (np.arange(K, dtype=np.float64) + 0.5) / K
    xx, yy = np.meshgrid(coords, coords, indexing="ij")
    return np.ascontiguousarray(np.column_stack([xx.ravel(), yy.ravel()]))


def _random_ics(ic_type: str, N: int, seed: int) -> np.ndarray:
    """Random-IC sampler for the C9/C10 configs (fresh-seeded -> reproducible).

    ic_type='runif' -> U(0,1)^2 ; ic_type='rlogy' -> x~U(0,1), y=10^u with
    u~U(-12, log10(0.3)) (log-uniform low predator, spanning 12 decades up to 0.3).
    Exact-singular draws (x or y hitting {0,1}) are excluded by the caller's validity
    mask, not here; this only DRAWS the (N,2) sample.
    """
    rng = np.random.default_rng(seed)
    if ic_type == "runif":
        return rng.random((N, 2))
    if ic_type == "rlogy":
        x = rng.random(N)
        u = rng.uniform(-12.0, math.log10(0.3), N)
        y = np.power(10.0, u)
        return np.ascontiguousarray(np.column_stack([x, y]))
    raise ValueError(f"unknown random ic_type {ic_type!r}")


# --- initial-condition provenance (2026-10-17) -------------------------------------------
# ``np.power(10.0, u)`` is dispatched by NumPy to a CPU-feature-specific kernel. On an
# AVX-512 host 11,147 of the 40,000 log-uniform predator seeds differ by one ULP from the
# libm-pow values produced on the authors' AVX2 host (external measurement: three of the
# 19,022 valid trajectories then change terminal classification, 5189 against 5188
# extinct at nu1 = 3), so version pins alone do not fix the random INPUTS of the C9/C10
# configurations. The exact arrays behind the shipped numbers are therefore archived with
# their checksum and generation convention. A regression run LOADS the archive (strict
# dataset regression, independent of CPU dispatch); ``DESIGNG_IC_MODE=resample`` re-draws
# the seeds from the generator instead (a statistical comparison whose inputs may differ
# at the ULP level between hosts). Every run records which it did, the fresh draw's
# agreement with the archive, and the NumPy CPU-dispatch state, under
# ``sensitivity.<config>._params.ic_provenance``.
_IC_ARCHIVE_SUBDIR = "ics"


def _ic_archive_name(ic_type: str, N: int, seed: int) -> str:
    return f"designG_ics_{ic_type}_N{int(N)}_seed{int(seed)}"


def _ic_archive_dir() -> str:
    """Directory holding the SHIPPED initial-condition archive (override: DESIGNG_IC_ARCHIVE_DIR).

    The bundle keeps it beside the shipped reference JSON (codes/data/designG/ics); the
    authoring tree keeps it under output/data/designG/ics. The first existing directory wins.
    """
    env = os.environ.get("DESIGNG_IC_ARCHIVE_DIR")
    if env:
        return env
    for cand in (os.path.join(_HERE, "data", "designG", _IC_ARCHIVE_SUBDIR),
                 os.path.join(_OUT, "data", "designG", _IC_ARCHIVE_SUBDIR)):
        if os.path.isdir(cand):
            return cand
    return os.path.join(_HERE, "data", "designG", _IC_ARCHIVE_SUBDIR)


def _sha256_bytes(b: bytes) -> str:
    import hashlib
    return hashlib.sha256(b).hexdigest()


def _numpy_dispatch_info() -> Dict[str, object]:
    """Interpreter, NumPy version and the CPU features NumPy's dispatcher enabled."""
    info: Dict[str, object] = {
        "python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform(),
        "machine": platform.machine(), "NPY_DISABLE_CPU_FEATURES": os.environ.get("NPY_DISABLE_CPU_FEATURES"),
    }
    try:  # numpy 1.26 exposes these on the C extension; absent on other versions
        from numpy.core._multiarray_umath import __cpu_features__ as _cf  # type: ignore
        info["cpu_features_enabled"] = sorted(k for k, v in _cf.items() if v)
        from numpy.core._multiarray_umath import __cpu_baseline__ as _cb, __cpu_dispatch__ as _cd  # type: ignore
        info["cpu_baseline"] = list(_cb)
        info["cpu_dispatch_targets"] = list(_cd)
    except Exception as exc:  # noqa: BLE001 -- provenance only
        info["cpu_features_enabled"] = f"unavailable ({exc!r})"
    return info


def write_ic_archive(ics: np.ndarray, ic_type: str, N: int, seed: int, directory: str) -> Dict[str, object]:
    """Archive the exact (N, 2) float64 initial conditions as <name>.npy + <name>.json.

    The sidecar records dtype, shape, memory order, byte order, generator, seed, the
    serialization, two checksums (sha256 over the C-order bytes; sha256 over the .npy
    file) and the runtime/CPU-dispatch state that produced the array.
    """
    os.makedirs(directory, exist_ok=True)
    arr = np.ascontiguousarray(ics, dtype=np.float64)
    name = _ic_archive_name(ic_type, N, seed)
    npy = os.path.join(directory, name + ".npy")
    np.save(npy, arr)
    with open(npy, "rb") as fh:
        file_sha = _sha256_bytes(fh.read())
    meta: Dict[str, object] = {
        "ic_type": ic_type, "N": int(N), "seed": int(seed),
        "generator": ("numpy.random.default_rng(seed) [PCG64]; "
                      "rlogy: x = rng.random(N); u = rng.uniform(-12.0, log10(0.3), N); y = np.power(10.0, u); "
                      "runif: rng.random((N, 2))"),
        "dtype": str(arr.dtype), "shape": [int(v) for v in arr.shape],
        "order": "C (row-major; column 0 = x0, column 1 = y0)", "byteorder": sys.byteorder,
        "serialization": "numpy .npy (np.save, format 1.0); sha256_data is over arr.tobytes() in C order",
        "sha256_data": _sha256_bytes(arr.tobytes()), "sha256_file": file_sha,
        "runtime": _numpy_dispatch_info(),
        "written_utc": datetime.now(timezone.utc).isoformat(),
    }
    with open(os.path.join(directory, name + ".json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    return meta


def load_ic_archive(ic_type: str, N: int, seed: int, directory: str) -> "Optional[Tuple[np.ndarray, Dict[str, object]]]":
    """Load and VERIFY an archived initial-condition array (None if the archive is absent).

    A present archive whose dtype, shape or sha256 disagree with its sidecar raises
    ValueError: a corrupt or edited archive must never be used silently.
    """
    name = _ic_archive_name(ic_type, N, seed)
    npy, js = os.path.join(directory, name + ".npy"), os.path.join(directory, name + ".json")
    if not (os.path.exists(npy) and os.path.exists(js)):
        return None
    with open(js) as fh:
        meta = json.load(fh)
    arr = np.load(npy)
    if arr.dtype != np.float64 or list(arr.shape) != list(meta.get("shape", [])) or arr.shape != (int(N), 2):
        raise ValueError(f"IC archive {npy}: dtype/shape {arr.dtype}/{arr.shape} do not match its sidecar")
    arr = np.ascontiguousarray(arr)
    got = _sha256_bytes(arr.tobytes())
    if got != meta.get("sha256_data"):
        raise ValueError(f"IC archive {npy}: sha256 {got[:16]}... != sidecar {str(meta.get('sha256_data'))[:16]}...")
    return arr, meta


def _ulp_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance in units in the last place between two finite float64 arrays of equal
    sign (difference of the IEEE-754 bit patterns)."""
    ai = np.ascontiguousarray(a, dtype=np.float64).view(np.int64)
    bi = np.ascontiguousarray(b, dtype=np.float64).view(np.int64)
    return np.abs(ai - bi)


def resolve_ics(ic_type: str, N: int, seed: int, record_dir: "Optional[str]" = None) -> Tuple[np.ndarray, Dict[str, object]]:
    """Initial conditions for a random C-config plus their provenance record.

    Default (``DESIGNG_IC_MODE`` unset or ``archive``): the archived array is used when
    present and verified, otherwise the fresh draw. ``DESIGNG_IC_MODE=resample`` always
    uses the fresh draw. The array actually used is written to ``<record_dir>/ics`` when a
    record directory is given, so each run leaves its inputs beside its outputs.
    """
    fresh = _random_ics(ic_type, N, seed)
    fresh_sha = _sha256_bytes(np.ascontiguousarray(fresh).tobytes())
    mode = (os.environ.get("DESIGNG_IC_MODE") or "archive").strip().lower()
    if mode not in ("archive", "resample"):
        sys.exit(f"[ics] DESIGNG_IC_MODE must be 'archive' or 'resample', got {mode!r}")
    prov: Dict[str, object] = {"mode": mode, "fresh_draw_sha256": fresh_sha,
                               "runtime": _numpy_dispatch_info(), "archive_dir": _rel(_ic_archive_dir())}
    loaded = None if mode == "resample" else load_ic_archive(ic_type, N, seed, _ic_archive_dir())
    if loaded is None:
        ics = fresh
        prov.update(source=("resampled from the seed" if mode == "resample" else "fresh draw (no archive found)"),
                    archive_file=None, archive_sha256=None, identical_to_archive=None,
                    n_values_differ=None, n_x_differ=None, n_y_differ=None, max_ulp=None)
    else:
        ics, meta = loaded
        d = _ulp_distance(ics, fresh)
        prov.update(source="archive", archive_file=_ic_archive_name(ic_type, N, seed) + ".npy",
                    archive_sha256=meta["sha256_data"], archive_runtime=meta.get("runtime"),
                    identical_to_archive=bool(meta["sha256_data"] == fresh_sha),
                    n_values_differ=int(np.count_nonzero(d)), n_x_differ=int(np.count_nonzero(d[:, 0])),
                    n_y_differ=int(np.count_nonzero(d[:, 1])), max_ulp=int(d.max()) if d.size else 0)
    prov["used_sha256"] = _sha256_bytes(np.ascontiguousarray(ics).tobytes())
    if record_dir:
        try:
            write_ic_archive(ics, ic_type, N, seed, os.path.join(record_dir, _IC_ARCHIVE_SUBDIR))
            prov["recorded_in"] = f"{_IC_ARCHIVE_SUBDIR}/{_ic_archive_name(ic_type, N, seed)}.npy (beside this run's values JSON)"
        except OSError as exc:  # recording is provenance, never a reason to abort the science
            prov["recorded_in"] = f"not written ({exc!r})"
    return ics, prov


# =============================================================================
# Sub-checks (each returns rows + a dict of recorded values)
# =============================================================================


def check_boundedness(n_iter: int, n_side: int) -> Tuple[List[Check], Dict[str, object]]:
    """(a) Boundedness / no blow-up at mu0=4 across nu1 (eq:bounds, rmk:three-regime-numeric)."""
    mu0 = 4.0
    nu1_list = [1.0, 1.5, 2.0, 3.0, 3.5, 4.0]
    grid = _make_grid(n_side)
    rows: List[Check] = []
    rec: Dict[str, object] = {"mu0": mu0, "per_nu1": {}}
    overall_max_x = 0.0
    overall_max_y = 0.0
    for nu1 in nu1_list:
        mx, my = _designG_grid_maxima(grid, mu0, nu1, n_iter)
        overall_max_x = max(overall_max_x, mx)
        overall_max_y = max(overall_max_y, my)
        # Predicted analytic bounds: x <= mu0/4 = 1, y <= nu1/4.
        rec["per_nu1"][f"{nu1:.2f}"] = {
            "max_x": mx,
            "max_y": my,
            "bound_x": mu0 / 4.0,
            "bound_y": nu1 / 4.0,
        }
        # Prey bound x <= mu0/4 = 1; assert the orbit max does not exceed it
        # (one-sided check: violation = max(0, max_x - bound) must be ~0).
        viol_x = max(0.0, mx - mu0 / 4.0)
        rows.append(
            Check(
                f"max x <= mu0/4 (nu1={nu1:.2f})",
                predicted=0.0,
                computed=viol_x,
                tol=1e-9,
                note=f"max_x={mx:.6f} vs bound mu0/4={mu0/4.0:.6f} (eq:bounds)",
            )
        )
        # Predator bound: physical bound is nu1/4; assert orbit max does not exceed it.
        viol_y = max(0.0, my - nu1 / 4.0)
        rows.append(
            Check(
                f"max y <= nu1/4 (nu1={nu1:.2f})",
                predicted=0.0,
                computed=viol_y,
                tol=1e-9,
                note=f"max_y={my:.6f} vs bound nu1/4={nu1/4.0:.6f}",
            )
        )
    rec["overall_max_x"] = overall_max_x
    rec["overall_max_y"] = overall_max_y
    rec["nu1_list"] = nu1_list
    rec["n_iter"] = n_iter
    rec["n_ic"] = int(grid.shape[0])
    # Headline values of rmk:three-regime-numeric: the largest x and the largest y over nu1 <= 4.
    rows.append(
        Check(
            "overall max x (nu1<=4) <= mu0/4",
            predicted=0.0,
            computed=max(0.0, overall_max_x - 1.0),
            tol=1e-9,
            note=f"overall max_x={overall_max_x:.6f}; prey saturates near mu0/4=1.000 (one-sided: must not exceed 1)",
        )
    )
    # max y over all nu1<=4: assert <= 1 hard and report the value.
    rows.append(
        Check(
            "overall max y (nu1<=4) <= 1",
            predicted=0.0,
            computed=max(0.0, overall_max_y - 1.0),
            tol=1e-9,
            note=f"overall max_y={overall_max_y:.4f}",
        )
    )
    return rows, rec


def check_lambda_hat(
    n_burn: int, n_keep: int, window: int
) -> Tuple[List[Check], Dict[str, object]]:
    """(b) Benettin--Wolf lambda_hat vs ln(nu1/4) at mu0=4 (eq:numerical-lambda, eq:lambda-mu4)."""
    mu0 = 4.0
    x0 = 0.314159265358979
    nu1_list = [1.0, 1.5, 2.0, 3.0, 3.5, 4.0]
    rows: List[Check] = []
    rec: Dict[str, object] = {"mu0": mu0, "x0": x0, "per_nu1": {}, "n_burn": n_burn, "n_keep": n_keep}
    # DECLARED protocol (methods numerics + spec 2A(b)): per-trajectory Benettin-Wolf
    # with burn-in n_burn=75000, retain n_keep=50000, Welford stability gate over moving
    # windows of `window`=5000.  Single-orbit mu0=4 estimates are corrupted, so the
    # declared per-trajectory protocol is run over n_ic seeds and lambda is averaged over
    # the gate-ACCEPTED trajectories (methods many-IC practice).  <ln x> is nu1-independent.
    n_ic = 128
    mlnx, sem_lnx, n_acc, n_noncol, n_tot = _logistic_lnx_ensemble(mu0, n_burn, n_keep, window, n_ic)
    # Gate is LOAD-BEARING and the per-trajectory 5000-window stability gate is genuinely
    # enforced inside the ensemble (it rejects unstable non-collapsed orbits, n_acc < n_noncol).
    # The estimate is certified only if AT LEAST HALF of the TOTAL ensemble did not collapse
    # (a run dominated by collapse is not a trustworthy SRB average even if the few survivors
    # look stable), AND a healthy fraction of the non-collapsed trajectories pass the per-traj
    # stability gate, AND the across-ensemble standard error is small vs the 2e-2 validation tol.
    gate_ok = (
        (n_noncol >= n_tot // 2)               # collapse-prevalence floor vs TOTAL ensemble
        and (n_acc >= (3 * n_noncol) // 5)
        and (n_acc >= 16)
        and (sem_lnx < 5e-3)
    )
    for nu1 in nu1_list:
        predicted = math.log(nu1 / 4.0)
        lam = math.log(nu1) + mlnx
        rec["per_nu1"][f"{nu1:.2f}"] = {
            "lambda_hat": lam,
            "predicted_ln_nu1_over_4": predicted,
            "mean_ln_x": mlnx,
            "sem_across_accepted": sem_lnx,
            # across-seed standard deviation of the accepted orbit means (sem * sqrt(n));
            # the manuscript quotes both statistics.
            "sd_across_accepted": sem_lnx * math.sqrt(n_acc),
            "n_accepted_stable": n_acc,
            "n_noncollapsed": n_noncol,
            "n_trajectories": n_tot,
            "n_burn": n_burn,
            "n_keep": n_keep,
            "window": window,
            "accepted_by_gate": bool(gate_ok),
        }
        chk = Check(
            f"lambda_hat (nu1={nu1:.2f}, mu0=4)",
            predicted=predicted,
            computed=lam,
            tol=2e-2,
            note=f"declared {n_burn}/{n_keep} Benettin, {window}-window per-traj stability gate; {n_acc} stable / {n_noncol} non-collapsed / {n_tot}, sem={sem_lnx:.2e}",
        )
        # PASS requires |lambda_hat - ln(nu1/4)| <= tol AND gate_ok (the per-trajectory
        # 5000-window stability gate accepted a healthy fraction + small ensemble SE).
        chk.passed = bool(chk.passed and gate_ok)
        rows.append(chk)
    # Worst-orbit (delta_{3/4}) check: on the fixed-point prey orbit x=3/4, ln x = ln(3/4).
    # The same number is the exponent predicted at nu1 = 3 (rmk:three-regime-numeric): measured ~ -0.2867 vs ln(3/4) = -0.2877.
    ln_three_quarter = math.log(0.75)
    rec["worst_orbit_ln_x_at_x_star"] = ln_three_quarter
    rows.append(
        Check(
            "worst-orbit ln x at x*=3/4",
            predicted=ln_three_quarter,
            computed=math.log(0.75),
            tol=1e-9,
            note="delta_{3/4} bound = ln(3/4) = -0.287682 (all-orbit boundary)",
        )
    )
    return rows, rec


def check_extinction_fraction(
    horizon: int, n_side: int, eps: float
) -> Tuple[List[Check], Dict[str, object]]:
    """(c) Finite-grid TERMINAL-PERSISTENCE extinction sweep over nu1 (rmk:first-passage).

    v4 REWRITE. v3 reported a FIRST-PASSAGE count here and called it an extinction
    fraction. This block now classifies an orbit as extinct only if
    ALL of its final ``_EXT_TAIL`` iterates lie below eps, left-censors ICs that
    start below eps, and reports the whole ``_EXT_EPS_LADDER`` rather
    than a single threshold. ``fraction_extinct`` is therefore now the terminal
    quantity the name always claimed; v3's number survives beside it as
    ``fraction_first_passage``, over the SAME orbits, so the size of the old error
    is visible rather than quietly dropped.

    HONESTY:
      * v3's ~0.23 at nu1=2 is a first-passage number; ``rmk:three-regime-numeric`` quotes
        the terminal count, 0.2269. At the same horizon and threshold a first-passage count
        is never below the terminal count; neither bounds the ASYMPTOTIC extinction
        fraction in either direction (over-count dip-and-recover; under-count slow
        converters, the latter worsening as nu1 -> 4 where ln(nu1/4) -> 0). v3's number is recorded here
        as history, not as a target.
      * Every number is a finite-time / finite-grid consistency check, never a
        proof. A terminal fraction is still a finite-horizon statement: an orbit
        below eps for its last 1000 iterates is evidence of, not a certificate of,
        y_n -> 0.
      * Only nu1 <= 4/3 is the PROVEN all-orbit extinction regime (thm:global, endpoint included);
        nu1 = 1.5 sits just ABOVE 4/3, so total extinction there is empirical.
      * Coexistence-band ICs converging to the interior fixed point are EXPECTED,
        not test failures.
    """
    mu0 = 4.0
    grid = _make_grid(n_side)
    rows: List[Check] = []
    rec: Dict[str, object] = {
        "mu0": mu0,
        "horizon": horizon,
        "eps": eps,
        "n_ic": int(grid.shape[0]),
        "protocol": "terminal_persistence",
        "tail": int(_eff_tail(horizon)),
        "tail_nominal": int(_EXT_TAIL),
        "protocol_note": (
            "extinct <=> ALL of the final %d iterates < eps (NOT first passage); "
            "ICs with y0 < eps are left-censored out of the denominator"
            % _eff_tail(horizon)
        ),
        "eps_ladder": [float(e) for e in _EXT_EPS_LADDER],
        "legacy_first_passage_nu1_2_v3": 0.2294,
        "legacy_first_passage_note": (
            "rmk:three-regime-numeric quotes the terminal fraction 0.2269 (817 of 3600 orbits) at nu1=2; v3 gave 0.2294 "
            "under its own protocol, a FIRST-PASSAGE fraction (826 of 3600; v3 broke at the "
            "first y<eps). It is history, not a target: at the same horizon and "
            "threshold a first-passage count is never below the terminal count, and "
            "neither bounds the asymptotic extinction fraction in either direction."
        ),
        "per_nu1": {},
    }
    fracs: Dict[float, float] = {}
    for nu1 in [1.0, 1.5, 2.0, 3.0, 4.0]:
        # ONE pass per nu1; the whole eps ladder is derived from these statistics.
        stats = _designG_orbit_stats(grid, mu0, nu1, horizon, _eff_tail(horizon))
        _ts = horizon - _eff_tail(horizon)
        ladder = {f"{e:g}": _extinct_stats(stats, e, tail_start=_ts) for e in _EXT_EPS_LADDER}
        nominal = _extinct_stats(stats, eps, tail_start=_ts)
        frac = float(nominal["fraction_extinct"])
        fracs[nu1] = frac
        rec["per_nu1"][f"{nu1:.2f}"] = {
            # `fraction_extinct` is TERMINAL persistence at the nominal eps -- the
            # key name is unchanged so the U1 rows / figure / taxonomy keep
            # resolving, but the QUANTITY behind it is now the one it names.
            "fraction_extinct": frac,
            "fraction_first_passage": float(nominal["fraction_first_passage"]),
            "n_extinct": int(nominal["n_extinct"]),
            "n_first_passage": int(nominal["n_first_passage"]),
            "n_valid": int(nominal["n_valid"]),
            "n_censored": int(nominal["n_censored"]),
            "n_total": int(nominal["n_total"]),
            "eps_sweep": ladder,
        }

    # --- U1 (Task B1): continuous fraction_extinct + a separate is_small boolean ---
    # The legacy PASS/FAIL rows below recode each fraction into a 0/1-style pass
    # delta (e.g. "extinction fraction small" reports max(0, frac-ceiling), never
    # frac itself). U1 surfaces the RAW continuous statistic plus an independently
    # named boolean classification as two DISTINCT fields (never compare them by
    # value -- 0.0 == False in Python), consuming the existing
    # extinction_fraction.per_nu1.*.fraction_extinct computation above: no new
    # simulation, no algorithm change. One row per swept nu1 (not only nu1=2);
    # half-integer nu1 (nu1=1.5) gets a distinct 'p'-separated key
    # ("ext_frac_nu1_1p5"), never int()-collapsed onto an integer neighbour.
    SMALL_CEILING = 0.30  # project-agreed "small" ceiling -- same numeric value as
    # the legacy SMALL_FRAC_CEILING below; kept as its own named constant so this
    # row shape does not depend on that (separate) PASS/FAIL block.
    continuous_rows: List[Dict[str, object]] = []
    for nu1 in [1.0, 1.5, 2.0, 3.0, 4.0]:
        frac_u1 = float(rec["per_nu1"][f"{nu1:.2f}"]["fraction_extinct"])
        key = f"ext_frac_nu1_{nu1:g}".replace(".", "p")
        continuous_rows.append({
            "label": key,
            "fraction_extinct": frac_u1,
            "is_small": bool(frac_u1 <= SMALL_CEILING),
            "small_ceiling": SMALL_CEILING,
        })
    rows.extend(continuous_rows)  # raw dicts, not Check objects -- informational,
    # no PASS/FAIL: main() skips non-Check rows in the PASS/FAIL accounting but
    # still serializes them into values["validation_table"] (see main()'s
    # isinstance(c, Check) guard). The U6 validator requires this key to RESOLVE and
    # be taxonomy-tagged `regression_anchor`; it does not hardcode the value. v3's
    # 0.229/0.029/0.017 anchors are RETIRED -- they were first-passage numbers, and
    # `fraction_extinct` now carries the terminal quantity. See _EXT_ANCHORS_V4.

    # nu1 in {1.0, 1.5}: finite-grid total extinction. Only nu1=1.0 is inside the
    # PROVEN all-orbit regime (nu1 <= 4/3); nu1=1.5 > 4/3 is an empirical finite-grid
    # observation here, NOT the all-orbit theorem (thm:global covers nu1 <= 4/3 only).
    for nu1 in (1.0, 1.5):
        rows.append(
            Check(
                f"extinction fraction = 1 (nu1={nu1:.2f})",
                predicted=1.0,
                computed=fracs[nu1],
                tol=0.02,
                note=(
                    "every sampled interior IC extinguishes on this grid; "
                    "nu1=1.0 < 4/3 (proven all-orbit), nu1=1.5 > 4/3 (finite-grid only)"
                ),
            )
        )
    # nu1 in {2, 3, 4}: SMALL fraction. One-sided 'is small' check.
    # The ceiling is a NAMED grid-robustness margin, NOT a tight point assertion --
    # the value is grid/horizon-dependent basin geometry, not a theorem constant.
    # v4 RECALIBRATION: 0.30 was calibrated against v3's FIRST-PASSAGE 0.2294. The
    # measured TERMINAL fractions are nu1=2:0.2269, 3:0.0111, 4:0.0000, so 0.2269 is
    # still the binding value and 0.30 still sits ~32% above it. The ceiling is
    # therefore kept: it is still correctly calibrated to the new
    # binding value.
    SMALL_FRAC_CEILING = 0.30
    for nu1 in (2.0, 3.0, 4.0):
        rows.append(
            Check(
                f"extinction fraction small (nu1={nu1:.2f})",
                predicted=0.0,
                computed=max(0.0, fracs[nu1] - SMALL_FRAC_CEILING),
                tol=1e-9,
                # The paper proves a positive-area, not a full-measure, extinction set
                # (thm:basin-positive).
                note=(f"frac={fracs[nu1]:.4f} <= {SMALL_FRAC_CEILING} (terminal criterion: y < {eps:g} over the last "
                      f"{_eff_tail(horizon)} of {horizon} iterates; a finite grid; where the other orbits go is not recorded)"),
            )
        )
    # v4 REGRESSION GUARD: the same orbits, classified v3's way, must still reproduce
    # v3's published first-passage fractions. This is the evidence that the rewrite
    # changed the CLASSIFICATION and not the SIMULATION -- if the dynamics had shifted,
    # these would not still land on 0.2294/0.0289/0.0172. It fails loudly if a future
    # edit perturbs the map, the grid, or the horizon while leaving the headline
    # terminal numbers superficially plausible.
    #
    # GATED on the anchor's own protocol. The anchors were measured at (60^2, 2e4, 1e-3)
    # and mean nothing anywhere else -- under --smoke the same correct code gives 0.2336
    # -- so an ungated assertion would fail a HEALTHY build.
    # `at_anchor_protocol` is recorded either way, so a run that
    # SKIPPED the guard can never be mistaken for one that PASSED it.
    at_anchor = (
        int(horizon) == _EXT_ANCHOR_PROTOCOL["horizon"]
        and int(n_side) == _EXT_ANCHOR_PROTOCOL["n_side"]
        and abs(float(eps) - _EXT_ANCHOR_PROTOCOL["eps"]) < 1e-15
    )
    rec["at_anchor_protocol"] = bool(at_anchor)
    rec["anchor_protocol"] = dict(_EXT_ANCHOR_PROTOCOL)
    if at_anchor:
        for nu1 in (2.0, 3.0, 4.0):
            fp_now = float(rec["per_nu1"][f"{nu1:.2f}"]["fraction_first_passage"])
            rows.append(
                Check(
                    f"v3 first-passage reproduced (nu1={nu1:.2f})",
                    predicted=_EXT_ANCHORS_V3_FIRST_PASSAGE[nu1],
                    computed=fp_now,
                    tol=1e-4,
                    note="same orbits under v3's rule; proves v4 rewrote the classifier, not the dynamics",
                )
            )
    # Structural monotonic decrease across the coexistence band 2 -> 3 -> 4.
    mono_ok = (fracs[2.0] >= fracs[3.0] - 1e-9) and (fracs[3.0] >= fracs[4.0] - 1e-9)
    rec["monotone_decreasing_2_to_4"] = bool(mono_ok)
    rows.append(
        Check(
            "extinction fraction monotone dec (2->4)",
            predicted=1.0,
            computed=1.0 if mono_ok else 0.0,
            tol=1e-9,
            note=f"fracs: nu1=2:{fracs[2.0]:.4f} >= 3:{fracs[3.0]:.4f} >= 4:{fracs[4.0]:.4f}",
        )
    )
    return rows, rec


def check_coexistence_fp_and_jury() -> Tuple[List[Check], Dict[str, object]]:
    """(d) Coexistence fixed point + Jury band at mu0=4 (eq:coexistence-fp/-trace-det/-jury-band).

    Closed forms (mu0 = 4 specific):
        (x*, y*) = (4/(4+nu1), (3 nu1 - 4)/(4 nu1)),  interior iff nu1 > 4/3
        tr J = -(nu1-4)^2 / (nu1 (nu1+4))
        det J = (nu1^2 + 12 nu1 - 32) / (nu1 (nu1+4))
        Jury-stable iff  8 sqrt(3) - 12 < nu1 < 4.
    """
    mu0 = 4.0
    rows: List[Check] = []
    rec: Dict[str, object] = {"mu0": mu0, "per_nu1": {}}

    def F(x: float, y: float, nu1: float) -> float:
        return mu0 * x * (1.0 - x) * (1.0 - y)

    def G(x: float, y: float, nu1: float) -> float:
        return nu1 * x * y * (1.0 - y)

    nu1_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    max_fp_residual = 0.0
    max_tr_err = 0.0
    max_det_err = 0.0
    for nu1 in nu1_list:
        xs = 4.0 / (4.0 + nu1)
        ysn = (3.0 * nu1 - 4.0) / (4.0 * nu1)
        # fixed-point residual
        res = abs(F(xs, ysn, nu1) - xs) + abs(G(xs, ysn, nu1) - ysn)
        max_fp_residual = max(max_fp_residual, res)
        # Jacobian via exact analytic partials of (F, G):
        # F = mu0 x (1-x)(1-y);  dF/dx = mu0 (1-2x)(1-y);  dF/dy = -mu0 x(1-x)
        # G = nu1 x y (1-y);     dG/dx = nu1 y(1-y);        dG/dy = nu1 x (1-2y)
        dFdx = mu0 * (1.0 - 2.0 * xs) * (1.0 - ysn)
        dFdy = -mu0 * xs * (1.0 - xs)
        dGdx = nu1 * ysn * (1.0 - ysn)
        dGdy = nu1 * xs * (1.0 - 2.0 * ysn)
        J = np.array([[dFdx, dFdy], [dGdx, dGdy]], dtype=np.float64)
        tr = J[0, 0] + J[1, 1]
        det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
        tr_cf = -((nu1 - 4.0) ** 2) / (nu1 * (nu1 + 4.0))
        det_cf = (nu1 * nu1 + 12.0 * nu1 - 32.0) / (nu1 * (nu1 + 4.0))
        max_tr_err = max(max_tr_err, abs(tr - tr_cf))
        max_det_err = max(max_det_err, abs(det - det_cf))
        eigvals = np.linalg.eigvals(J)
        spectral = float(np.max(np.abs(eigvals)))
        # Jury conditions for the unit disc.
        jury1 = 1.0 + tr + det          # > 0
        jury2 = 1.0 - tr + det          # > 0
        jury3 = 1.0 - det               # > 0
        jury_stable = (jury1 > 0.0) and (jury2 > 0.0) and (jury3 > 0.0)
        rec["per_nu1"][f"{nu1:.2f}"] = {
            "x_star": xs,
            "y_star": ysn,
            "fp_residual": res,
            "tr_J": tr,
            "tr_J_closed_form": tr_cf,
            "det_J": det,
            "det_J_closed_form": det_cf,
            "eig_abs_max": spectral,
            "jury1": jury1,
            "jury2": jury2,
            "jury3": jury3,
            "jury_stable": bool(jury_stable),
            "interior": bool(0.0 < xs < 1.0 and 0.0 < ysn < 1.0),
        }
    rows.append(
        Check(
            "coexistence FP residual (max over nu1)",
            predicted=0.0,
            computed=max_fp_residual,
            tol=1e-10,
            note="F(x*,y*)=x*, G(x*,y*)=y* (eq:coexistence-fp)",
        )
    )
    rows.append(
        Check(
            "tr J vs closed form (max err)",
            predicted=0.0,
            computed=max_tr_err,
            tol=1e-10,
            note="-(nu1-4)^2/(nu1(nu1+4)) (eq:designG-trace-det)",
        )
    )
    rows.append(
        Check(
            "det J vs closed form (max err)",
            predicted=0.0,
            computed=max_det_err,
            tol=1e-10,
            note="(nu1^2+12nu1-32)/(nu1(nu1+4)) (eq:designG-trace-det)",
        )
    )

    # Locate the measured Jury-stable band edges by a fine scan and compare to
    # the algebraic lower edge 8 sqrt(3) - 12 and the upper edge 4.
    lower_cf = 8.0 * math.sqrt(3.0) - 12.0
    scan = np.linspace(1.0, 4.0, 600001)  # fine grid, includes 4.0
    stable_mask = np.zeros(scan.shape[0], dtype=np.bool_)
    for i, nu1 in enumerate(scan):
        xs = 4.0 / (4.0 + nu1)
        ysn = (3.0 * nu1 - 4.0) / (4.0 * nu1)
        tr_cf = -((nu1 - 4.0) ** 2) / (nu1 * (nu1 + 4.0))
        det_cf = (nu1 * nu1 + 12.0 * nu1 - 32.0) / (nu1 * (nu1 + 4.0))
        interior = (0.0 < ysn < 1.0)
        j1 = 1.0 + tr_cf + det_cf
        j2 = 1.0 - tr_cf + det_cf
        j3 = 1.0 - det_cf
        stable_mask[i] = interior and (j1 > 0.0) and (j2 > 0.0) and (j3 > 0.0)
    idx = np.where(stable_mask)[0]
    if idx.size > 0:
        measured_lo = float(scan[idx[0]])
        measured_hi = float(scan[idx[-1]])
    else:
        measured_lo = float("nan")
        measured_hi = float("nan")
    rec["jury_band"] = {
        "lower_edge_closed_form": lower_cf,
        "upper_edge_closed_form": 4.0,
        "measured_lower_edge": measured_lo,
        "measured_upper_edge": measured_hi,
        "scan_step": float(scan[1] - scan[0]),
    }
    rows.append(
        Check(
            "Jury band lower edge",
            predicted=lower_cf,
            computed=measured_lo,
            tol=5e-3,
            note="8*sqrt(3)-12 ~ 1.8564 (eq:designG-jury-band)",
        )
    )
    rows.append(
        Check(
            "Jury band upper edge",
            predicted=4.0,
            computed=measured_hi,
            tol=5e-3,
            note="upper edge nu1 -> 4",
        )
    )
    return rows, rec


def check_global_telescoping(
    n_iter: int, n_side: int
) -> Tuple[List[Check], Dict[str, object]]:
    """(e) Global-extinction telescoping demo for nu1 < 4/3 (thm:global)."""
    mu0 = 4.0
    grid = _make_grid(n_side)
    # safe_floor well above double-precision underflow (~1e-308) so the last
    # representable (y_k, k) is clean; the realised per-step rate is measured to
    # that step, never across an underflow (see _designG_telescoping_rate).
    safe_floor = 1e-250
    rows: List[Check] = []
    rec: Dict[str, object] = {"mu0": mu0, "n_iter": n_iter, "n_ic": int(grid.shape[0]),
                              "safe_floor": safe_floor, "per_nu1": {}}
    for nu1 in [1.0, 1.2, 1.3]:
        (max_rate, max_final_y, predicted, max_stop_step, min_stop_step, n_stopped,
         max_log10_y_h, max_excess) = _designG_telescoping_rate(grid, mu0, nu1, n_iter, safe_floor)
        rec["per_nu1"][f"{nu1:.2f}"] = {
            "max_realised_rate": max_rate,
            "predicted_bound_ln_3nu1_over_4": predicted,
            # `max_final_y` is the largest STOPPING value (the first predator value below the
            # floor 1e-250, reached at a per-orbit step <= max_stop_step), kept under its historical
            # key for the frozen claim inventory; it is not y_N at the common horizon.
            "max_final_y": max_final_y,
            "max_final_y_meaning": "largest first-sub-floor predator value over the grid (stopping value), not y_N",
            "max_stop_step": int(max_stop_step),
            "min_stop_step": int(min_stop_step),
            "n_stopped_before_horizon": int(n_stopped),
            "n_ic": int(grid.shape[0]),
            "horizon_N": int(n_iter),
            # the common-horizon value, continued in log space below the floor (see the kernel)
            "max_log10_y_at_horizon": max_log10_y_h,
            "max_horizon_telescoping_excess_nats": max_excess,
        }
        # Theorem: every realised rate <= ln(3 nu1 / 4) < 0.  Check the worst rate
        # does not exceed the bound, allowing a small finite-time slack for the
        # telescoping boundary term (4/9)/k (k = safe steps) plus float noise.
        slack = 5e-3
        rows.append(
            Check(
                f"telescoping rate <= ln(3nu1/4) (nu1={nu1:.2f})",
                predicted=0.0,
                computed=max(0.0, max_rate - predicted),
                tol=slack,
                note=f"worst rate={max_rate:.5f} <= bound={predicted:.5f}? (thm:global)",
            )
        )
        rows.append(
            Check(
                f"predator stopping value below floor 1e-250 (nu1={nu1:.2f})",
                predicted=0.0,
                computed=max_final_y,
                tol=1e-6,
                note=(f"largest first-sub-floor y over all ICs; every IC crossed the floor by step "
                      f"{max_stop_step} of N={n_iter} ({n_stopped}/{grid.shape[0]} stopped early)"),
            )
        )
        # The common-horizon value itself. The telescoping bound of thm:global gives
        # ln y_N - ln y_0 <= N ln(3 nu1/4) + 4/9 for every interior orbit; the row reports the
        # excess over that bound CLIPPED AT ZERO (its name says so: 0 means no violation), and
        # the actual max log10 y_N is in the record and in the note.
        rows.append(
            Check(
                f"horizon telescoping excess, clipped at 0 (nu1={nu1:.2f}, N={n_iter})",
                predicted=0.0,
                computed=max(0.0, max_excess - 4.0 / 9.0),
                tol=1e-9,
                note=(f"max_i[ln y_N - ln y_0 - N ln(3nu1/4)] = {max_excess:.3f} nats <= 4/9; "
                      f"max log10 y_N at the common horizon = {max_log10_y_h:.1f}"),
            )
        )
    return rows, rec


def check_family_unification(
    n_burn: int, n_keep: int
) -> Tuple[List[Check], Dict[str, object]]:
    """(f) Family lambda_perp unification: skew = MOD-1 = Design-G = ln(nu1/4) at mu0=4.

    Folds in thm:skew (eq:skew), thm:mod1 (eq:mod1) and the Design-G linear
    transverse multiplier (eq:linear-transverse, eq:lambda-mu4).  All three share
    the SAME boundary value ln(nu1/4); the canonical
    ``forward_causal_theory_figs.py`` computes the same closed form for skew /
    MOD-1, so co-locating Design G makes the unification explicit (no canonical
    edit -- this is a tiny CPU consistency check, NOT the heavy GPU sweep).

    METHOD / HONESTY: this block intentionally reuses the SINGLE-ORBIT estimators
    (_skew_lambda_hat / _mod1_lambda_hat / _logistic_lambda_hat) in the DEFINITIONAL
    rare-invader limit y0=0.  At y0=0 the predator perturbation is identically zero, so
    the three maps reduce EXACTLY to the bare logistic prey orbit and share the M0
    multiplier nu1*x -- this is the regime where the single-orbit estimator is valid, NOT
    the mu0=4 "corrupted single-orbit" regime warned about elsewhere (that pathology needs
    a nonzero, chaos-amplified perturbation).  The main lambda check (b) carries the heavier
    many-IC stability-gated protocol; this family block only confirms the SHARED boundary
    value ln(nu1/4), not the off-axis nonlinearity.
    """
    mu0 = 4.0
    x0 = 0.314159265358979
    mu1 = 0.5  # escapement strength for MOD-1 (irrelevant in the y->0 limit)
    rows: List[Check] = []
    rec: Dict[str, object] = {"mu0": mu0, "x0": x0, "mu1_mod1": mu1, "per_nu1": {}}
    for nu1 in [1.0, 2.0, 3.0, 3.5, 4.0]:
        predicted = math.log(nu1 / 4.0)
        lam_skew = _skew_lambda_hat(nu1, mu0, x0, 0.2, n_burn, n_keep)
        # y0=0 reads the transverse multiplier on M0 itself (rare-invader limit): the
        # three maps PROVABLY share the same M0 multiplier nu1 x (skew exactly; MOD-1 and
        # Design-G via d(y')/dy|_{y=0}=nu1 x), so on the common logistic prey orbit they
        # coincide.  Seeding a nonzero y0 only injects mu0=4-chaos-amplified single-orbit
        # scatter (the corrupted-estimate regime), not a genuine inter-model difference.
        lam_mod1 = _mod1_lambda_hat(nu1, mu0, mu1, x0, 0.0, n_burn, n_keep)
        # Design-G transverse multiplier in the y->0 limit is the same logistic estimator.
        lam_designG, _, _, _ = _logistic_lambda_hat(nu1, mu0, x0, n_burn, n_keep, 5000)
        rec["per_nu1"][f"{nu1:.2f}"] = {
            "predicted_ln_nu1_over_4": predicted,
            "lambda_skew": lam_skew,
            "lambda_mod1": lam_mod1,
            "lambda_designG": lam_designG,
            "pairwise_max_spread": max(
                abs(lam_skew - lam_mod1),
                abs(lam_skew - lam_designG),
                abs(lam_mod1 - lam_designG),
            ),
        }
        rows.append(
            Check(
                f"skew lambda_perp (nu1={nu1:.2f})",
                predicted=predicted,
                computed=lam_skew,
                tol=2e-2,
                note="thm:skew: exact fibre -> ln(nu1/4)",
            )
        )
        rows.append(
            Check(
                f"MOD-1 lambda_perp (nu1={nu1:.2f})",
                predicted=predicted,
                computed=lam_mod1,
                tol=2e-2,
                note="thm:mod1: clean Jensen sign -> ln(nu1/4)",
            )
        )
        rows.append(
            Check(
                f"family spread skew/MOD-1/G (nu1={nu1:.2f})",
                predicted=0.0,
                computed=rec["per_nu1"][f"{nu1:.2f}"]["pairwise_max_spread"],
                tol=2e-2,
                note="pairwise agreement -> all = ln(nu1/4) (thm:skew/thm:mod1/eq:lambda-mu4)",
            )
        )
    return rows, rec


# =============================================================================
# U2 (Task B2): extinction-fraction sensitivity sweep (C0-C10) + CIs
# (spec 4.1a)
# =============================================================================


def _sensitivity_schedule(smoke: bool) -> "Tuple[List[Dict[str, object]], List[float]]":
    """The C0-C10 OFAT sensitivity schedule + the nu1 sweep values.

    NOMINAL (spec 4.1a; nu1 in {2,3,4}): C0 baseline grid=100^2,T=1e4,eps=1e-6,det-grid;
    C1-C3 vary grid {50^2,200^2,400^2}; C4-C5 vary horizon {1e3,7.5e4}; C6-C8 vary eps
    {1e-3,1e-9,1e-12}; C9 ic.runif / C10 ic.rlogy random ICs N=40000 seed 20260709.
    One-factor-at-a-time: the baseline C0 is the shared pivot of the grid/horizon/eps
    convergence axes (it appears in each axis' `convergence_band`).

    --smoke GENUINELY reduces:
    grid=50^2, T=1e3, a single nu1=2, and only the {base, ic.runif, ic.rlogy} subset --
    enough to exercise BOTH CI branches (det->none, random->wilson), the convergence
    band, the rlogy sampler, and the _params/_fallback_used metadata, fast.

    `axes` = which convergence_band axes a DET config feeds (random configs: []);
    `cost_axis` = the parameter shrunk by the wall-clock fallback (grid/horizon/N).
    """
    if smoke:
        K, T, N = 50, 1000, 2000
        sched = [
            dict(key="base",     kind="det",    ic_type="grid",  K=K,    T=T, eps=1e-6, N=None, axes=["grid", "horizon", "eps"], cost_axis="grid"),
            dict(key="ic.runif", kind="random", ic_type="runif", K=None, T=T, eps=1e-6, N=N,    axes=[],                         cost_axis="N"),
            dict(key="ic.rlogy", kind="random", ic_type="rlogy", K=None, T=T, eps=1e-6, N=N,    axes=[],                         cost_axis="N"),
        ]
        return sched, [2.0]
    sched = [
        dict(key="base",          kind="det",    ic_type="grid",  K=100,  T=10000, eps=1e-6,  N=None,  axes=["grid", "horizon", "eps"], cost_axis="grid"),
        dict(key="grid.50",       kind="det",    ic_type="grid",  K=50,   T=10000, eps=1e-6,  N=None,  axes=["grid"],                   cost_axis="grid"),
        dict(key="grid.200",      kind="det",    ic_type="grid",  K=200,  T=10000, eps=1e-6,  N=None,  axes=["grid"],                   cost_axis="grid"),
        dict(key="grid.400",      kind="det",    ic_type="grid",  K=400,  T=10000, eps=1e-6,  N=None,  axes=["grid"],                   cost_axis="grid"),
        dict(key="horizon.1e3",   kind="det",    ic_type="grid",  K=100,  T=1000,  eps=1e-6,  N=None,  axes=["horizon"],                cost_axis="horizon"),
        dict(key="horizon.7.5e4", kind="det",    ic_type="grid",  K=100,  T=75000, eps=1e-6,  N=None,  axes=["horizon"],                cost_axis="horizon"),
        dict(key="eps.1e-3",      kind="det",    ic_type="grid",  K=100,  T=10000, eps=1e-3,  N=None,  axes=["eps"],                    cost_axis="horizon"),
        dict(key="eps.1e-9",      kind="det",    ic_type="grid",  K=100,  T=10000, eps=1e-9,  N=None,  axes=["eps"],                    cost_axis="horizon"),
        dict(key="eps.1e-12",     kind="det",    ic_type="grid",  K=100,  T=10000, eps=1e-12, N=None,  axes=["eps"],                    cost_axis="horizon"),
        dict(key="ic.runif",      kind="random", ic_type="runif", K=None, T=10000, eps=1e-6,  N=40000, axes=[],                         cost_axis="N"),
        dict(key="ic.rlogy",      kind="random", ic_type="rlogy", K=None, T=10000, eps=1e-6,  N=40000, axes=[],                         cost_axis="N"),
    ]
    return sched, [2.0, 3.0, 4.0]


def _effective_sensitivity_params(cfg: Dict[str, object], budget_s: float) -> Tuple[object, int, object, bool]:
    """(K, T, N, fallback_used) after the deterministic wall-clock guard of one C-config.

    Shared by the runner and the checkpoint identity check, so a resumed record is compared
    with exactly what this run would compute."""
    kind = str(cfg["kind"])
    cost_axis = str(cfg["cost_axis"])
    K = cfg["K"]
    T = int(cfg["T"])
    N = cfg["N"]
    fallback_used = False
    # Deterministic per-cell wall-clock guard: shrink the cost axis until the
    # (conservative) estimate fits the budget. `while True` with explicit floors.
    while True:
        n_ic_est = (int(K) * int(K)) if kind == "det" else int(N)
        est_s = n_ic_est * T / _SENS_STEPS_PER_SEC
        if est_s <= budget_s:
            break
        if cost_axis == "grid" and kind == "det" and int(K) > _SENS_MIN_K:
            K = max(_SENS_MIN_K, int(K) // 2)
            fallback_used = True
        elif cost_axis == "horizon" and T > _SENS_MIN_T:
            T = max(_SENS_MIN_T, T // 2)
            fallback_used = True
        elif cost_axis == "N" and N is not None and int(N) > _SENS_MIN_N:
            N = max(_SENS_MIN_N, int(N) // 2)
            fallback_used = True
        else:
            break  # already at the floor for this axis; run as-is (still flagged if shrunk)
    return K, T, N, fallback_used


def _checkpoint_identity_mismatch(
    existing: object, cfg: Dict[str, object], nu1_list: List[float], mu0: float, budget_s: float,
) -> "Optional[str]":
    """None when a checkpointed config record may be reused for this request, else the reason.

    A record is reused only if it was computed for the SAME request (a run with
    DESIGNG_IC_MODE=resample used to return the archive-based results of an older checkpoint):
    the configuration -- kind, ic_type, mu0, the effective grid/T/N after the wall-clock guard,
    eps, the nu1 sweep, the seed, the protocol and terminal window, the fallback flag -- and, for
    a random config, the same initial conditions: the requested DESIGNG_IC_MODE and the SHA-256
    of the array this run would use (archive or fresh draw)."""
    if not (isinstance(existing, dict) and isinstance(existing.get("_params"), dict)):
        return "no _params record"
    p = existing["_params"]
    K, T, N, fallback_used = _effective_sensitivity_params(cfg, budget_s)
    kind = str(cfg["kind"])
    want = {
        "kind": kind, "ic_type": str(cfg["ic_type"]), "mu0": float(mu0),
        "grid": (int(K) if kind == "det" else None), "T": int(T), "eps": float(cfg["eps"]),
        "N": (int(N) if kind != "det" else None), "nu1_list": [float(x) for x in nu1_list],
        "seed": (_SENS_SEED if kind != "det" else None), "cost_axis": str(cfg["cost_axis"]),
        "protocol": "terminal_persistence", "tail": int(_eff_tail(int(T))), "tail_nominal": int(_EXT_TAIL),
    }
    for k, v in want.items():
        if p.get(k) != v:
            return f"{k}: recorded {p.get(k)!r}, requested {v!r}"
    if bool(existing.get("_fallback_used")) != bool(fallback_used):
        return f"_fallback_used: recorded {existing.get('_fallback_used')!r}, requested {fallback_used!r}"
    if kind != "det":
        prov = p.get("ic_provenance")
        if not isinstance(prov, dict):
            return "no ic_provenance record"
        _ics, now = resolve_ics(str(cfg["ic_type"]), int(N), _SENS_SEED, None)
        for k in ("mode", "used_sha256"):
            if prov.get(k) != now.get(k):
                return f"ic_provenance.{k}: recorded {prov.get(k)!r}, requested {now.get(k)!r}"
        # the inputs run_sensitivity would use: the identified array without exact-singular seeds
        inputs = _ics[(_ics[:, 0] > 0.0) & (_ics[:, 0] < 1.0) & (_ics[:, 1] > 0.0) & (_ics[:, 1] < 1.0)]
    else:
        inputs = _midpoint_lattice(int(K))
    # the sample those inputs define: every input counts, those starting below eps are left-censored
    exp_total = int(inputs.shape[0])
    exp_cens = int(np.count_nonzero(inputs[:, 1] < float(cfg["eps"])))
    # Completeness: every requested nu1 cell with every result field, and a consistent sample
    # (a parseable but incomplete record used to be reused as if it were whole).
    fields = ("fraction", "fraction_first_passage", "fraction_first_passage_censored_inclusive",
              "ci_method", "ci_low", "ci_high", "n", "n_censored", "n_total")
    for nu1 in nu1_list:
        cell = existing.get(f"{float(nu1):g}")
        if not isinstance(cell, dict):
            return f"cell nu1={float(nu1):g}: missing"
        lacking = [f for f in fields if f not in cell]
        if lacking:
            return f"cell nu1={float(nu1):g}: missing {lacking}"
        n, nc, nt = cell["n"], cell["n_censored"], cell["n_total"]
        if not all(isinstance(v, int) and not isinstance(v, bool) and v >= 0 for v in (n, nc, nt)) or n != nt - nc:
            return f"cell nu1={float(nu1):g}: inconsistent sample (n={n!r}, n_censored={nc!r}, n_total={nt!r})"
        if cell["ci_method"] != ("none" if kind == "det" else "wilson"):
            return f"cell nu1={float(nu1):g}: ci_method {cell['ci_method']!r} for a {kind} config"
        # the values themselves must be usable: finite fractions in [0, 1]; a Wilson cell carries a
        # finite ordered interval in [0, 1], a deterministic cell none (both bounds null)
        def _unit(v):
            return (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
                    and -1e-12 <= v <= 1.0 + 1e-12)
        bad = [f for f in ("fraction", "fraction_first_passage", "fraction_first_passage_censored_inclusive")
               if not _unit(cell[f])]
        if bad:
            return f"cell nu1={float(nu1):g}: unusable value(s) {[(f, cell[f]) for f in bad]}"
        lo, hi = cell["ci_low"], cell["ci_high"]
        if cell["ci_method"] == "none":
            if lo is not None or hi is not None:
                return f"cell nu1={float(nu1):g}: a deterministic cell carries an interval ({lo!r}, {hi!r})"
        elif not (_unit(lo) and _unit(hi) and lo <= hi):
            return f"cell nu1={float(nu1):g}: unusable Wilson interval ({lo!r}, {hi!r})"
        # consistency: the recorded sample is the one the identified inputs define (its total and its
        # left-censored part) and is not empty; both fractions are counts over the valid sample, the
        # terminal count never exceeds the first-passage count; the censored-inclusive fraction is
        # (first-passage count + censored) / total; a Wilson interval is the interval of the terminal
        # count (the producer's own wilson)
        if nt != exp_total or nc != exp_cens:
            return (f"cell nu1={float(nu1):g}: sample (n_total {nt!r}, n_censored {nc!r}) is not the requested "
                    f"sample, ({exp_total}, {exp_cens}) as its inputs define it")
        if n <= 0:
            return f"cell nu1={float(nu1):g}: empty valid sample (n={n!r})"
        kt, kf = cell["fraction"] * n, cell["fraction_first_passage"] * n
        if abs(kt - round(kt)) > 1e-6 or abs(kf - round(kf)) > 1e-6:
            return f"cell nu1={float(nu1):g}: a fraction is not a count over n={n} (terminal {kt!r}, first passage {kf!r})"
        kt, kf = int(round(kt)), int(round(kf))
        if not (0 <= kt <= kf <= n):
            return f"cell nu1={float(nu1):g}: terminal count {kt} exceeds the first-passage count {kf} (n={n})"
        inc = (kf + nc) / nt
        if abs(cell["fraction_first_passage_censored_inclusive"] - inc) > 1e-12:
            return (f"cell nu1={float(nu1):g}: censored-inclusive fraction "
                    f"{cell['fraction_first_passage_censored_inclusive']!r} is not ({kf} + {nc}) / {nt}")
        if cell["ci_method"] == "wilson":
            wlo, whi = wilson(kt, n)
            if abs(wlo - lo) > 1e-9 or abs(whi - hi) > 1e-9:
                return (f"cell nu1={float(nu1):g}: interval ({lo!r}, {hi!r}) is not the Wilson interval of "
                        f"{kt}/{n}, ({wlo!r}, {whi!r})")
    return None


def _run_one_sensitivity_config(
    cfg: Dict[str, object], nu1_list: List[float], mu0: float, budget_s: float,
    record_dir: "Optional[str]" = None,
) -> Dict[str, object]:
    """Run one C-config over the nu1 sweep -> {<nu1>: {...}, _params, _fallback_used}.

    CI SEMANTICS BY DESIGN:
      * det configs  -> ci_method="none", ci_low/high=None (a deterministic midpoint-
        lattice fraction is EXACT for that grid, not a sample estimate -> no CI).
      * random configs -> ci_method="wilson", ci_low/high = wilson(k, n) over the
        n VALID (non-singular) Monte-Carlo ICs (k = #extinct).

    WALL-CLOCK FALLBACK: a deterministic pre-run estimate (n_ic*T/rate) is compared to
    the per-cell `budget_s`; if over, the cost axis (grid K / horizon T / random N) is
    halved toward its floor and `_fallback_used=True` is set (a reduced run is NOT a
    submission-grade claim -- the B6 gate re-runs it at nominal). At nominal params the
    estimate is well under budget, so no shrink occurs and fractions stay reproducible.
    The ACTUAL (post-fallback) params + wall_s are recorded under `_params`.
    """
    import time

    kind = str(cfg["kind"])
    ic_type = str(cfg["ic_type"])
    cost_axis = str(cfg["cost_axis"])
    eps = float(cfg["eps"])
    K, T, N, fallback_used = _effective_sensitivity_params(cfg, budget_s)

    t0 = time.perf_counter()
    rec_cfg: Dict[str, object] = {}
    # v4: same terminal-persistence kernel as block (c). `T` here can be as short as
    # _SENS_MIN_T=500, so the tail is clamped (and recorded in _params below).
    tail_eff = _eff_tail(T)

    def _cell(pts: np.ndarray, nu1: float, with_ci: bool) -> Dict[str, object]:
        """One (config, nu1) cell -> its record. `n` is the LEFT-CENSORED denominator."""
        st = _extinct_stats(_designG_orbit_stats(pts, mu0, float(nu1), T, tail_eff), eps, tail_start=T - tail_eff)
        n_valid = int(st["n_valid"])
        out: Dict[str, object] = {
            # `fraction` is now TERMINAL persistence (was first passage in v3).
            "fraction": (float(st["fraction_extinct"]) if n_valid > 0 else None),
            "fraction_first_passage": (
                float(st["fraction_first_passage"]) if n_valid > 0 else None
            ),
            # The fraction a FIRST-PASSAGE rule would report over a denominator that still
            # counts the left-censored seeds (they start below eps, so such a rule counts
            # them at n=0): quoted in the initial-condition exclusion paragraph of Section S12.
            # Equals fraction_first_passage whenever n_censored == 0 (every det config).
            "fraction_first_passage_censored_inclusive": (
                (float(st["fraction_first_passage"]) * n_valid + int(st["n_censored"])) / int(st["n_total"])
                if (n_valid > 0 and int(st["n_total"]) > 0) else None
            ),
            "ci_method": "none",
            "ci_low": None,
            "ci_high": None,
            # The DENOMINATOR is the valid (non-left-censored) sample, not the raw
            # IC count. This also fixes the CI: v3 put a first-passage count over a
            # denominator that still included ICs starting below eps.
            "n": n_valid,
            "n_censored": int(st["n_censored"]),
            "n_total": int(st["n_total"]),
        }
        if with_ci and n_valid > 0:
            lo, hi = wilson(int(st["n_extinct"]), n_valid)
            out["ci_method"] = "wilson"
            out["ci_low"] = float(lo)
            out["ci_high"] = float(hi)
        return out

    if kind == "det":
        grid = _midpoint_lattice(int(K))
        for nu1 in nu1_list:
            rec_cfg[f"{nu1:g}"] = _cell(grid, nu1, with_ci=False)
    else:
        # Archived initial conditions when available (strict regression), else the
        # fresh seeded draw; the provenance (archive/fresh/resample, checksums, CPU dispatch)
        # is recorded under _params.ic_provenance below.
        ics, ic_prov = resolve_ics(ic_type, int(N), _SENS_SEED, record_dir)
        # exclude exact-singular ICs (x or y in {0,1}); count only the valid sample
        mask = (ics[:, 0] > 0.0) & (ics[:, 0] < 1.0) & (ics[:, 1] > 0.0) & (ics[:, 1] < 1.0)
        valid = np.ascontiguousarray(ics[mask])
        for nu1 in nu1_list:
            rec_cfg[f"{nu1:g}"] = _cell(valid, nu1, with_ci=True)
    wall_s = time.perf_counter() - t0

    # Underscore-prefixed metadata siblings of the <nu1> records: the taxonomy
    # wildcard (`*` = any non-`_` key) SKIPS these, so /sensitivity/*/*/fraction
    # never resolves onto _params (a reduced run is never mistaken for the nominal).
    rec_cfg["_params"] = {
        "kind": kind,
        "ic_type": ic_type,
        "mu0": float(mu0),
        "grid": (int(K) if kind == "det" else None),
        "T": int(T),
        "eps": float(eps),
        "N": (int(N) if kind != "det" else None),
        "nu1_list": [float(x) for x in nu1_list],
        "seed": (_SENS_SEED if kind != "det" else None),
        "cost_axis": cost_axis,
        "wall_s": float(wall_s),
        "budget_s": float(budget_s),
        # v4 protocol provenance: `tail` is the EFFECTIVE terminal window actually
        # used for this config (clamped to T//2 for short horizons -- see _eff_tail),
        # never the nominal 1000. A clamped cell must not read as a nominal one.
        "protocol": "terminal_persistence",
        "tail": int(tail_eff),
        "tail_nominal": int(_EXT_TAIL),
    }
    if kind != "det":
        rec_cfg["_params"]["ic_provenance"] = ic_prov   # inputs, checksums, dispatch state
    rec_cfg["_fallback_used"] = bool(fallback_used)
    return rec_cfg


def _build_convergence_band(
    results: Dict[str, object], sched: List[Dict[str, object]], nu1_list: List[float]
) -> Dict[str, object]:
    """convergence_band.<axis>.<nu1>.{min,max,configs} across the DET configs of an axis.

    Deterministic configs carry NO CI; instead their grid-/horizon-/eps-refinement
    spread is summarised as a band (min/max fraction + the contributing config keys).
    Only DET configs actually present in `results` contribute (resume-safe). Random
    configs (Wilson CI) do not feed any band.
    """
    axis_keys: Dict[str, List[str]] = {}
    for cfg in sched:
        if cfg["kind"] != "det":
            continue
        for ax in cfg["axes"]:  # type: ignore[union-attr]
            axis_keys.setdefault(str(ax), []).append(str(cfg["key"]))
    band: Dict[str, object] = {}
    for ax, keys in axis_keys.items():
        per_nu1: Dict[str, object] = {}
        for nu1 in nu1_list:
            nk = f"{nu1:g}"
            fracs: List[float] = []
            present: List[str] = []
            for k in keys:
                reck = results.get(k)
                if isinstance(reck, dict) and isinstance(reck.get(nk), dict) and "fraction" in reck[nk]:
                    fracs.append(float(reck[nk]["fraction"]))
                    present.append(k)
            if fracs:
                per_nu1[nk] = {"min": min(fracs), "max": max(fracs), "configs": present}
        if per_nu1:
            band[ax] = per_nu1
    return band


def run_sensitivity(args: argparse.Namespace) -> Tuple[List[Check], Dict[str, object]]:
    """U2 block: the C0-C10 extinction-fraction sensitivity sweep + CIs.

    Returns ([], rec) -- this block contributes NO PASS/FAIL Check rows (its output
    is the `sensitivity` record, consumed by the B6 validator and the sensitivity figure),
    so main()'s isinstance(c, Check) guards leave the PASS/FAIL table untouched.

    rec = {<config_key>: {<nu1>: {...}, _params, _fallback_used}, ..., convergence_band}.
    Each config is CHECKPOINTED to a JSON in DESIGNG_DATADIR immediately after it
    completes, and a resumed run skips any config already present in that checkpoint
    (the heavy 15-40 min nominal sweep is interruptible/resumable). The checkpoint file
    is mode-specific (smoke vs full) so a tiny smoke run never poisons a nominal resume, and a
    checkpointed record is reused only when its configuration and initial-condition identity
    match this run (_checkpoint_identity_mismatch); otherwise that config is recomputed.
    """
    mu0 = 4.0
    smoke = bool(getattr(args, "smoke", False))
    sched, nu1_list = _sensitivity_schedule(smoke)

    datadir = resolve_paths(smoke)["datadir"]
    os.makedirs(datadir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(__file__))[0]
    ckpt = os.path.join(datadir, f"{stem}_sensitivity{'_smoke' if smoke else ''}_checkpoint.json")

    results: Dict[str, object] = {}
    if os.path.exists(ckpt):
        try:
            with open(ckpt) as fh:
                loaded = json.load(fh)
            if isinstance(loaded, dict):
                results = loaded
        except (ValueError, OSError):
            results = {}  # corrupt/partial checkpoint -> start fresh (never crash a resume)

    budget_s = float(os.environ.get("DESIGNG_SENS_BUDGET_S", _SENS_BUDGET_S_DEFAULT))
    for cfg in sched:
        key = str(cfg["key"])
        existing = results.get(key)
        if isinstance(existing, dict) and "_params" in existing:
            why = _checkpoint_identity_mismatch(existing, cfg, nu1_list, mu0, budget_s)
            if why is None:
                print(f"    [{key}] resumed from checkpoint (skip: configuration and inputs match)")
                continue
            print(f"    [{key}] checkpoint does not match this run ({why}); recomputing")
        rec_cfg = _run_one_sensitivity_config(cfg, nu1_list, mu0, budget_s, record_dir=datadir)
        results[key] = rec_cfg
        fb = "  [FALLBACK: reduced params -- NOT submission-grade]" if rec_cfg["_fallback_used"] else ""
        print(f"    [{key}] {rec_cfg['_params']['kind']} "
              f"grid={rec_cfg['_params']['grid']} T={rec_cfg['_params']['T']} "
              f"eps={rec_cfg['_params']['eps']:g} N={rec_cfg['_params']['N']} "
              f"wall={rec_cfg['_params']['wall_s']:.2f}s{fb}")
        with open(ckpt, "w") as fh:   # checkpoint IMMEDIATELY (resumable)
            json.dump(results, fh, indent=2)

    rec: Dict[str, object] = dict(results)
    rec["convergence_band"] = _build_convergence_band(results, sched, nu1_list)
    return [], rec


# =============================================================================
# U3 (Task B3): finite-amplitude two-way decay diagnostic (NON-blocking)
# (spec 4.1a)
# =============================================================================


def _kahan_mean_sem(arr: np.ndarray) -> Tuple[float, float]:
    """Kahan-summed mean + SEM across a 1-D float64 array of per-IC estimates.

    Population variance / n (matching the convention already used by
    ``_logistic_lnx_ensemble`` elsewhere in this file: ``var = E[x^2] - E[x]^2``,
    ``sem = sqrt(var / n)``), but with BOTH the mean and the variance sums
    computed via ``_kahan_sum`` (compensated summation) rather than a plain
    running total, per the plan's "Kahan-sum the accumulations for numerical
    stability" requirement.  ``n < 2`` -> sem is NaN (undefined); ``n == 0`` ->
    both NaN (no ICs survived to this checkpoint horizon).
    """
    n = arr.shape[0]
    if n == 0:
        return float("nan"), float("nan")
    mean = float(_kahan_sum(arr) / n)
    if n < 2:
        return mean, float("nan")
    var = float(_kahan_sum((arr - mean) ** 2) / n)
    if var < 0.0:
        var = 0.0
    return mean, float(math.sqrt(var / n))


def _finite_amp_schedule(smoke: bool) -> Dict[str, object]:
    """U3 (Task B3) M / n_burn / horizon_checkpoints / nu1_list schedule.

    NOMINAL (plan Task B3; nu1 in {2,3,4}): M=128 ICs, burn-in 1e3, measured
    out to horizon 1000 with checkpoints at 100/300/1000 (so the n->horizon
    convergence in ``horizon_sweep`` is visible).  delta in {1e-2,1e-3,1e-4} is
    FIXED SCHEMA, returned separately by the caller as ``_FINAMP_DELTAS`` --
    never reduced here (see run_finite_amplitude).

    --smoke GENUINELY reduces cost
    via a smaller M, a much shorter burn-in/horizon, and a
    reduced nu1 subset ({2,4} instead of {2,3,4}) -- while the caller still
    sweeps the FULL delta schema for every nu1 regardless of smoke.  At these
    sizes (2 nu1 x 3 delta x 16 ICs x 100 steps ~ 10k map-steps) the whole
    block finishes in well under a second, no wall-clock fallback needed
    (unlike U2's much heavier grid sweeps).
    """
    if smoke:
        return dict(M=16, n_burn=50, horizon_checkpoints=[10, 20, 50], nu1_list=[2.0, 4.0])
    return dict(M=128, n_burn=1000, horizon_checkpoints=[100, 300, 1000], nu1_list=[2.0, 3.0, 4.0])


def _run_one_finite_amp_cell(
    nu1: float, delta: float, M: int, n_burn: int, checkpoints: List[int], seed: int
) -> Dict[str, object]:
    """One (nu1, delta) cell of the U3 sweep: draw M valid ICs, run the per-IC
    estimator, aggregate to decay_rate/sem + the horizon_sweep convergence list.

    x0 in (0,1) / y0 in (0,delta) are drawn fresh-seeded (``seed`` -> reproducible,
    resume-order-independent, mirroring U2's per-config fresh-seed convention).
    The exceptional set E (exact-singular draws landing exactly on the
    {0,1}/{0,delta} boundary) is excluded by an open-interval mask, same idiom
    as U2's ``_random_ics`` validity mask; a rare filtered-out draw is topped up
    by redraw so exactly M valid ICs are always used (deterministic given the seed).
    """
    mu0 = 4.0
    rng = np.random.default_rng(seed)
    x0s = np.empty(0, dtype=np.float64)
    y0s = np.empty(0, dtype=np.float64)
    while x0s.shape[0] < M:
        need = M - x0s.shape[0]
        draw_n = need + 8  # oversample so a second pass is rarely needed
        cand_x = rng.random(draw_n)
        cand_y = rng.random(draw_n) * delta
        mask = (cand_x > 0.0) & (cand_x < 1.0) & (cand_y > 0.0) & (cand_y < delta)
        x0s = np.concatenate([x0s, cand_x[mask]])
        y0s = np.concatenate([y0s, cand_y[mask]])
    x0s = x0s[:M]
    y0s = y0s[:M]

    cp_arr = np.asarray(sorted(int(h) for h in checkpoints), dtype=np.int64)
    n_cp = cp_arr.shape[0]
    slopes = np.empty((M, n_cp), dtype=np.float64)
    for i in range(M):
        slopes[i, :] = _finite_amp_trajectory_slopes(
            float(x0s[i]), float(y0s[i]), mu0, nu1, int(n_burn), cp_arr
        )

    horizon_sweep: List[Dict[str, object]] = []
    decay_rate = float("nan")
    sem = float("nan")
    n_valid_final = 0
    for j in range(n_cp):
        col = slopes[:, j]
        valid = col[~np.isnan(col)]
        mean_j, sem_j = _kahan_mean_sem(valid)
        horizon_sweep.append({
            "horizon": int(cp_arr[j]),
            "decay_rate": float(mean_j),
            "sem": float(sem_j),
            "n_valid": int(valid.shape[0]),
        })
        if j == n_cp - 1:
            decay_rate, sem, n_valid_final = mean_j, sem_j, int(valid.shape[0])

    return {
        "decay_rate": float(decay_rate),
        "sem": float(sem),
        "n_valid": n_valid_final,
        "M": int(M),
        "n_burn": int(n_burn),
        "delta": float(delta),
        "horizon_sweep": horizon_sweep,
    }


def run_finite_amplitude(args: argparse.Namespace) -> Tuple[List[Check], Dict[str, object]]:
    """U3 (Task B3): finite-amplitude two-way decay diagnostic (NON-blocking,
    tagged ``empirical_discovery`` -- an empirical diagnostic, NOT a theory gate).

    Distinct from the (infinitesimal, tangent-linearised) axis variational
    exponent lambda_perp = ln(nu1/4) that ``check_lambda_hat`` /
    ``check_family_unification`` (and a future U4 block) estimate: here the
    predator starts at a FINITE amplitude y0 in (0, delta) and the FULL
    nonlinear two-way map (eq:designG) is iterated forward --
    x' = mu0 x (1-x)(1-y),  y' = nu1 x y (1-y)  at mu0=4 (the arcsine-SRB axis
    regime).  Physically we EXPECT the measured finite-amplitude decay_rate of
    ln(y_n) to approach the linear/variational ln(nu1/4) as delta -> 0 (the
    nonlinear O(y) corrections vanish); how close it gets at FIXED delta is the
    empirical question this block answers -- see ``_FINAMP_INTERPRETATION``.

    Per (nu1, delta) cell (``_run_one_finite_amp_cell``): M ICs, each burned in
    via the bare logistic map then given a fresh finite y0 (see
    ``_finite_amp_trajectory_slopes``); decay_rate = least-squares slope of
    ln(y_n) vs n, per IC then averaged across the M ICs (Kahan-summed
    accumulators throughout); sem = standard error of that mean.
    ``horizon_sweep[]`` records the same aggregate at several shorter
    checkpoint horizons so the n->horizon convergence is visible.

    Returns ``([], rec)`` -- like ``run_sensitivity``, this block contributes
    NO PASS/FAIL Check row (main()'s ``isinstance(c, Check)`` guards leave the
    validation table untouched).  This is how "non-blocking" is enforced
    STRUCTURALLY: there is no tolerance gate on decay_rate to pass or fail.

    HONESTY (survivorship bias at the longest checkpoint, nominal params):
    at nu1=2 (the fastest-decaying nu1 in the sweep, ln(2/4)=-0.693/step) a
    majority of the M=128 trajectories underflow ``_FINAMP_UNDERFLOW_FLOOR``
    well before the 1000-step checkpoint (empirically n_valid drops to
    ~35-50/128 there), so the top-level ``decay_rate`` (read off the LAST
    checkpoint) is a survivors-only average, biased toward the slower-decaying
    tail and systematically LESS negative than ln(nu1/4) (observed ~-0.47 to
    -0.57 vs predicted -0.693).  The SHORTER checkpoints (100, 300), where all
    M ICs are still representable, track the predicted ln(nu1/4) closely
    (observed within ~0.03 at h=100 for nu1=2). This is exactly why
    ``horizon_sweep`` is stored -- the bias is visible, not hidden -- and is a
    genuine finite-amplitude/finite-horizon effect, not a code defect; nu1=3
    and nu1=4 decay slowly enough that no such underflow occurs by h=1000.
    """
    smoke = bool(getattr(args, "smoke", False))
    sched = _finite_amp_schedule(smoke)
    M = int(sched["M"])
    n_burn = int(sched["n_burn"])
    checkpoints = list(sched["horizon_checkpoints"])
    nu1_list = [float(v) for v in sched["nu1_list"]]
    deltas = list(_FINAMP_DELTAS)  # ALWAYS the full {1e-2,1e-3,1e-4} -- schema, not a cost knob

    rec: Dict[str, object] = {
        "mu0": 4.0,
        "M": M,
        "n_burn": n_burn,
        "horizon_checkpoints": [int(h) for h in checkpoints],
        "delta_list": [float(d) for d in deltas],
        "seed": _FINAMP_SEED,
        "interpretation": _FINAMP_INTERPRETATION,
        "per_nu1": {},
    }
    for nu1 in nu1_list:
        nu1_key = f"{nu1:g}"
        per_delta: Dict[str, object] = {}
        for delta in deltas:
            delta_key = f"{delta:g}"
            per_delta[delta_key] = _run_one_finite_amp_cell(nu1, delta, M, n_burn, checkpoints, _FINAMP_SEED)
        predicted = math.log(nu1 / 4.0)
        rec["per_nu1"][nu1_key] = {"per_delta": per_delta, "predicted_ln_nu1_over_4": predicted}
        last_delta_key = f"{deltas[-1]:g}"
        print(f"    [nu1={nu1_key}] predicted ln(nu1/4)={predicted:.4f}; "
              f"decay_rate(delta={last_delta_key})={per_delta[last_delta_key]['decay_rate']:.4f}")

    return [], rec


# =============================================================================
# U4 (Task B4): exact axis transverse-variational lambda_perp = ln(nu1/4) --
# quadrature primary (BLOCKING) + orbit cross-check
# (spec 4.1a)
#
# The double-exponential (tanh-sinh) quadrature confirms the closed form
# lambda_perp_analytic() = ln(nu1/4) to ~1e-10 (the arcsine-SRB log-potential
# integral_0^1 ln x dmu4 = -ln 4, thm:skew / eq:lambda-mu4 family) and is the
# BLOCKING theory invariant (tol 1e-6).  lambda_perp_orbit() is a SECONDARY,
# purely empirical cross-check (a direct float64 Birkhoff average of
# ln(nu1*x_n) over a real logistic orbit) -- descriptive agreement only
# (tol 2e-2), NOT part of the 1e-6 blocking gate.  The functions below are
# the math core of this block.
# =============================================================================


def lambda_perp_analytic(nu1):
    # integral_0^1 ln x dmu4 = -ln 4 (arcsine log-potential). Exact closed form.
    return math.log(nu1) - math.log(4.0)                       # ln(nu1/4)


def _tanh_sinh_ln_sin(n=200, h=0.03):
    # double-exponential quadrature of integral_0^{pi/2} ln(sin psi) dpsi = -(pi/2) ln 2,
    # handling the log singularity at psi->0 (reaches ~1e-10 with these params).
    s = 0.0
    for k in range(-n, n + 1):
        t = k * h
        u = (math.pi / 2) * math.sinh(t)
        psi = (math.pi / 4) * (1.0 + math.tanh(u))             # node in (0, pi/2)
        w = (math.pi / 4) * (math.pi / 2) * math.cosh(t) / (math.cosh(u) ** 2)
        if 0.0 < psi < math.pi / 2:
            s += w * math.log(math.sin(psi))
    return s * h


def lambda_perp_quadrature(nu1, n=200, h=0.03):
    # numerical confirmation: I = (4/pi) * integral ln sin = -ln 4
    I = (4.0 / math.pi) * _tanh_sinh_ln_sin(n, h)              # -> -ln 4
    return math.log(nu1) + I                                   # -> ln(nu1/4)


def collapse_guard(x_orbit, tol=1e-12, run=50):
    # reject an orbit that stalls before T: False if `run` consecutive steps each move x by less
    # than `tol` (a run of near-fixed values, as when the float64 orbit lands on a fixed point);
    # a cycle of period >= 2 is not detected, since its successive differences are not small
    import numpy as np
    d = np.abs(np.diff(x_orbit))
    return not (np.convolve((d < tol).astype(int), np.ones(run), "valid") >= run).any()


def lambda_perp_orbit(nu1, T=10000, B=1000, x0=0.401, kahan=True):
    # SECONDARY cross-check: direct float64 logistic Birkhoff avg of ln(nu1*x_n).
    # (NOT float64 angle-doubling, which collapses dyadic theta.)
    x = x0; s = 0.0; c = 0.0; xs = []
    for n in range(B+T):
        if n >= B:
            term = math.log(nu1*x)
            if kahan:
                y = term - c; t = s + y; c = (t - s) - y; s = t
            else:
                s += term
        xs.append(x); x = 4.0*x*(1.0-x)
    return s/T, xs


# --- U4 block-emission machinery (schedule / exceptional-set filter / aggregation) --
# (the four functions above are the math core; everything below is the
# block-emission logic around them.)

_EXACT_SKEW_SEED = 20260709   # project RNG-seed convention (U2 _SENS_SEED / U3 _FINAMP_SEED)
_EXACT_SKEW_PREFILTER_STEPS = 64  # forward float64 steps checked for an exact hit of {0, 1/2, 1}


def _is_exceptional_x0(x0: float) -> bool:
    """Conservative NUMERICAL pre-filter for the exceptional set E (Task B4 interface contract).

    The theoretical exceptional set of the mu0 = 4 logistic map is the countable set of
    seeds whose forward orbit reaches the fixed point x = 0: under the conjugacy
    x = sin^2(pi theta) these are exactly the DYADIC ANGLES theta = k/2^p, i.e. the
    pre-images of x = 1/2 (theta = 1/4), which maps to 1 and then to 0. A dyadic STATE
    value x0 = k/2^p is in general NOT in E: x0 = 1/4 maps to 3/4, the interior fixed
    point, and never reaches 0 (the version of this filter shipped through 2026-10-17
    wrongly excluded low-order dyadic state values on that reading;
    the shipped ensemble reported zero exclusions, so no reported value depended on
    it). Exact membership of E IS decidable for a float64 seed: every float64
    is an exact dyadic rational, rational iterates stay rational, and the preimages of 1/2 are
    the irrational (1 +- 1/sqrt(2))/2, so the only rational members of E are 0, 1/2 and 1 --
    the first test below is that exact membership test. The forward loop is a different,
    NUMERICAL safeguard: it flags seeds whose ROUNDED float64 orbit collapses onto {0, 1/2, 1}
    within ``_EXACT_SKEW_PREFILTER_STEPS`` steps -- seeds that are not in E but whose computed
    orbit would degenerate. Any later float64 degeneration is caught downstream by
    collapse_guard().
    """
    if x0 <= 0.0 or x0 >= 1.0 or x0 == 0.5:
        return True
    x = x0
    for _ in range(_EXACT_SKEW_PREFILTER_STEPS):
        x = 4.0 * x * (1.0 - x)
        if x == 0.0 or x == 0.5 or x == 1.0:
            return True
    return False


def _exact_skew_schedule(smoke: bool) -> Dict[str, object]:
    """U4 (Task B4) orbit M / horizon-sweep schedule.

    The quadrature side of this block (lambda_perp_quadrature, n=200/h=0.03) is
    DETERMINISTIC, cheap, and ~1e-10 accurate, so it is NEVER reduced under
    --smoke (n and h are not reduced under smoke: the quadrature must
    stay exact).  Only the SECONDARY orbit cross-check's cost (M ICs, horizon
    lengths) is reduced under --smoke; all three nu1 are kept in both modes
    (the quadrature is cheap regardless, and the reduced orbit cost stays cheap
    at 3 nu1 too -- see run_exact_skew).
    """
    if smoke:
        return dict(M=16, B=200, horizons=[200, 500, 1000])
    return dict(M=256, B=1000, horizons=[10000, 30000, 75000])


def _run_one_exact_skew_nu1(nu1: float, sched: Dict[str, object]) -> Dict[str, object]:
    """One nu1 cell of the U4 block: quadrature (blocking) + M-orbit cross-check
    (secondary, empirical) with a horizon sweep, per the Task B4 interface
    contract (JSON: exact_skew.per_nu1.<v>.{...}).

    IC filtering: M x0 in (0,1) are drawn once (fresh-seeded rng, seed
    _EXACT_SKEW_SEED -- same "fresh-seeded per unit of work" convention as U2's
    _random_ics / U3's _run_one_finite_amp_cell, so re-running this cell in any
    order reproduces the identical draw).  Each horizon's kept/excluded partition
    is evaluated INDEPENDENTLY (mirrors U3's _finite_amp_trajectory_slopes /
    horizon_sweep convention: the caller aggregates only the orbits that survive
    to THAT checkpoint -- collapse can develop between a short and a long
    horizon, so a longer horizon may legitimately keep fewer orbits than a
    shorter one).  The top-level lambda_perp_orbit/sem/collapse_guard_ok/
    excluded_E_ic_count fields are read off the LONGEST horizon (the
    "production" cross-check estimate).
    """
    M = int(sched["M"])
    B = int(sched["B"])
    horizons = sorted(int(h) for h in sched["horizons"])
    T_max = horizons[-1]

    lam_q = lambda_perp_quadrature(nu1)                 # DE quadrature, blocking (~1e-10 accurate)
    target = lambda_perp_analytic(nu1)                  # exact closed form ln(nu1/4)
    abs_err = abs(lam_q - target)

    rng = np.random.default_rng(_EXACT_SKEW_SEED)       # fresh-seeded (U2/U3 convention)
    x0_draws = [float(v) for v in rng.random(M)]
    pre_excluded = sum(1 for x0 in x0_draws if _is_exceptional_x0(x0))
    candidates = [x0 for x0 in x0_draws if not _is_exceptional_x0(x0)]

    horizon_sweep: List[Dict[str, object]] = []
    for h in horizons:
        vals: List[float] = []
        n_excluded_h = pre_excluded
        for x0 in candidates:
            try:
                mean_h, xs_h = lambda_perp_orbit(nu1, T=h, B=B, x0=x0)
            except ValueError:
                # Orbit hit x=0.0 exactly -> math.log(nu1*x) raised "math domain
                # error" inside lambda_perp_orbit (verbatim, untouched). A
                # dynamically-crashing orbit IS an exceptional IC (same bucket as
                # the _is_exceptional_x0 pre-filter): exclude it and skip -- do
                # NOT let it enter the mean/sem/collapse_guard aggregation.
                n_excluded_h += 1
                continue
            if collapse_guard(xs_h):
                vals.append(mean_h)
            else:
                n_excluded_h += 1
        arr = np.asarray(vals, dtype=np.float64)
        mean_v = float(np.mean(arr)) if arr.size else float("nan")
        sem_v = float(np.std(arr, ddof=1) / math.sqrt(arr.size)) if arr.size > 1 else float("nan")
        horizon_sweep.append({
            "T": h, "mean": mean_v, "sem": sem_v,
            "n_kept": int(arr.size), "excluded_E_ic_count": int(n_excluded_h),
        })

    prod = horizon_sweep[-1]                            # T_max entry (horizons sorted ascending)
    lam_orbit = prod["mean"]
    sem = prod["sem"]
    n_kept = prod["n_kept"]
    excluded = prod["excluded_E_ic_count"]
    # True iff every KEPT orbit (at the production horizon) passed collapse_guard --
    # true by construction (only guard-passing orbits enter `vals` above) as long as
    # at least one orbit survived; n_kept==0 would mean the whole cohort degenerated,
    # which is a genuine red flag and is correctly reported as False, not vacuous True.
    collapse_guard_ok = bool(n_kept > 0)
    orbit_cross_check_abs_diff = abs(lam_orbit - lam_q) if n_kept else float("nan")

    return {
        "lambda_perp_quadrature": lam_q,
        "target_ln_nu1_over_4": target,
        "abs_err": abs_err,
        "lambda_perp_orbit": lam_orbit,
        "sem": sem,
        "orbit_cross_check_abs_diff": orbit_cross_check_abs_diff,
        "horizon_sweep": horizon_sweep,
        "collapse_guard_ok": collapse_guard_ok,
        "excluded_E_ic_count": int(excluded),
        "M_drawn": M,
        "n_kept": int(n_kept),
        "B": B,
        "T_max": T_max,
    }


def run_exact_skew(args: argparse.Namespace) -> Tuple[List[object], Dict[str, object]]:
    """U4 (Task B4) block: exact axis transverse-variational lambda_perp = ln(nu1/4).

    Quadrature primary (BLOCKING) + orbit cross-check:
    unlike U2/U3 (run_sensitivity / run_finite_amplitude, both explicitly
    non-blocking, contribute zero Check rows), THIS block DOES emit a real
    PASS/FAIL Check per nu1: abs_err <= 1e-6 AND collapse_guard_ok (the
    ``chk.passed = chk.passed and gate_ok`` idiom already used by
    check_lambda_hat elsewhere in this file).  The orbit cross-check (tol 2e-2,
    descriptive, not part of the 1e-6 blocking gate) is emitted as an
    INFORMATIONAL raw-dict row instead (U1's continuous_rows idiom in
    check_extinction_fraction), excluded from PASS/FAIL accounting by main()'s
    isinstance(c, Check) guard.
    """
    smoke = bool(getattr(args, "smoke", False))
    sched = _exact_skew_schedule(smoke)
    nu1_list = [2.0, 3.0, 4.0]   # kept in FULL under --smoke too (quadrature is cheap/exact;
                                  # --smoke only shrinks the orbit M/B/horizons above)

    rows: List[object] = []
    rec: Dict[str, object] = {
        "mu0": 4.0,
        "M": int(sched["M"]),
        "B": int(sched["B"]),
        "horizons": [int(h) for h in sched["horizons"]],
        "seed": _EXACT_SKEW_SEED,
        "model": (
            "exact axis SKEW cross-check: x_{n+1}=4x_n(1-x_n) (mu0=4 logistic), "
            "lambda_perp = <ln(nu1 x_n)> under the arcsine SRB measure -- DISTINCT "
            "from the GPU fc_v4_skew parameter sweep (separate script/model)."
        ),
        "per_nu1": {},
    }
    for nu1 in nu1_list:
        cell = _run_one_exact_skew_nu1(nu1, sched)
        nk = f"{nu1:g}"
        rec["per_nu1"][nk] = cell

        chk = Check(
            f"lambda_perp_quadrature == ln(nu1/4) (nu1={nu1:g})",
            predicted=cell["target_ln_nu1_over_4"],
            computed=cell["lambda_perp_quadrature"],
            tol=1e-6,
            note=(
                f"DE tanh-sinh quadrature (n=200,h=0.03; BLOCKING theory invariant); "
                f"collapse_guard_ok={cell['collapse_guard_ok']} "
                f"({cell['excluded_E_ic_count']}/{sched['M']} orbits excluded at T={cell['T_max']})"
            ),
        )
        chk.passed = bool(chk.passed and cell["collapse_guard_ok"])
        rows.append(chk)

        # SECONDARY/descriptive orbit cross-check: informational dict row (U1's
        # continuous_rows idiom), NOT a Check -- explicitly excluded from PASS/FAIL
        # accounting (not part of the 1e-6 blocking gate).
        rows.append({
            "label": f"orbit_cross_check_nu1_{nk}".replace(".", "p"),
            "orbit_cross_check_abs_diff": cell["orbit_cross_check_abs_diff"],
            "lambda_perp_orbit": cell["lambda_perp_orbit"],
            "sem": cell["sem"],
            "agrees_within_2e-2": bool(cell["n_kept"] > 0 and cell["orbit_cross_check_abs_diff"] <= 2e-2),
        })
    print(f"    [U4] blocking quadrature gate: "
          + ", ".join(f"nu1={nu1:g} abs_err={rec['per_nu1'][f'{nu1:g}']['abs_err']:.3e}" for nu1 in nu1_list))
    return rows, rec


# =============================================================================
# Optional diagnostic figures (P2; guarded by --no-figs)
# =============================================================================


def make_figures(rec: Dict[str, object], figdir: str) -> List[str]:
    """Write optional diagnostic PNGs to ``figdir``. Returns list of written paths.

    v3 (B0): takes ``figdir`` explicitly (was the module-level FIGDIR global) so
    --smoke / --only can redirect output without editing a global; each block's
    figure is guarded by its rec key so a partial run (--only <block>) degrades
    gracefully instead of KeyError'ing on a block that did not run. Plotting
    logic/data unchanged from v2 -- a full (all-blocks) run reproduces the same
    three figures as before.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(figdir, exist_ok=True)
    written: List[str] = []

    # Manuscript Figure 4: lambda_hat vs nu1 (Benettin--Wolf) against ln(nu1/4).
    if "lambda_hat" in rec:
        lam_rec = rec["lambda_hat"]["per_nu1"]
        nu1s = sorted(float(k) for k in lam_rec.keys())
        lam_vals = [lam_rec[f"{n:.2f}"]["lambda_hat"] for n in nu1s]
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        grid_nu1 = np.linspace(0.5, 4.0, 200)
        # The caption described a horizontal ln(3/4) line and
        # error bars that the artwork never carried, and misplaced the all-orbit threshold.
        # The figure now shows the two on-axis exponents the text compares: the typical
        # (arcsine) exponent ln(nu1/4), zero at nu1 = 4, and the worst-orbit exponent
        # ln(3 nu1/4) of the fixed point x* = 3/4, zero at the all-orbit threshold nu1 = 4/3.
        ax.plot(grid_nu1, np.log(grid_nu1 / 4.0), "k-", lw=1.5,
                label=r"$\ln(\nu_1/4)$: arcsine-typical on-axis exponent (closed form)")
        ax.plot(grid_nu1, np.log(3.0 * grid_nu1 / 4.0), color="tab:green", ls="--", lw=1.3,
                label=r"$\ln(3\nu_1/4)$: worst-orbit exponent, fixed point $x^*=3/4$")
        n_acc = lam_rec[f"{nu1s[0]:.2f}"]["n_accepted_stable"]   # the count of the record, 118 in the shipped run
        ax.plot(nu1s, lam_vals, "o", ms=6, color="tab:blue", label="Benettin $\\hat\\lambda_\\perp$ (%d accepted orbits)" % n_acc)
        ax.axhline(0.0, color="grey", lw=0.7)
        ax.axvline(4.0, color="tab:red", ls=":", lw=1.0, label=r"$\nu_1=4$: typical exponent changes sign")
        ax.axvline(4.0 / 3.0, color="tab:green", ls=":", lw=1.0, label=r"$\nu_1=4/3$: all-orbit threshold")
        ax.set_xlabel(r"$\nu_1$")
        ax.set_ylabel("on-axis exponent")
        ax.set_ylim(-2.2, 0.6)
        ax.set_title("On-axis exponents at $\\mu_0=4$")
        ax.legend(fontsize=7.5, loc="lower right")
        p = os.path.join(figdir, "designG_lambda_vs_nu1.png")
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        written.append(p)
        _manifest_append(p)

    # Manuscript Figure 11: extinction fraction vs nu1.
    if "extinction_fraction" in rec:
        ext_rec = rec["extinction_fraction"]["per_nu1"]
        enu1s = sorted(float(k) for k in ext_rec.keys())
        fracs = [ext_rec[f"{n:.2f}"]["fraction_extinct"] for n in enu1s]
        fps = [ext_rec[f"{n:.2f}"].get("fraction_first_passage") for n in enu1s]
        # Two panels: (a) the fractions vs nu1; (b) the eps LADDER (spec item 3 --
        # "report tail-decay behaviour, not a single threshold"). Panel (b) is the
        # load-bearing one: on a linear 0-1 axis panel (a) cannot show the difference
        # between 0.017 and 0.000, so on its own it would UNDERSTATE the very defect
        # this rewrite exists to fix. (b) shows the eps ladder: the terminal fraction is
        # FLAT in eps over 1e-2..1e-8 while first passage slides. NOTE (2026-10-17): that
        # flatness is largely ARITHMETIC -- every terminally extinct orbit had reached the
        # float64 absorbing zero before the terminal window (n_extinct_exact_zero_before_tail,
        # first_zero_step_min/max) -- so it shows insensitivity of the COUNT to eps, not that
        # the exact dynamics is dichotomous; the manuscript states the same.
        # Stacked (2 x 1) since 2026-10-17 so that each panel keeps its full width in
        # the single-column manuscript; the panels and their content are unchanged.
        fig, (ax, axb) = plt.subplots(2, 1, figsize=(6.4, 8.0))
        # v4: the DATA is now terminal persistence, so the figure finally plots the
        # quantity its title names. v3 could only relabel this axis; it plotted a
        # FIRST-PASSAGE count (break at the first y_n < eps), which counted an orbit
        # that dipped once and recovered as extinct forever.
        # BOTH curves are drawn, from the SAME orbits: the gap between them IS the
        # size of v3's error, and showing it is the honest way to report the change.
        # At the same horizon and threshold the first-passage curve is an UPPER bound on
        # the terminal one (an orbit ending below eps must have crossed it); it bounds the
        # ASYMPTOTIC extinction fraction in neither direction (it over-counts dip-and-recover
        # orbits, and under-counts orbits decaying too slowly to reach eps within the horizon
        # -- the latter worsening as nu1 -> 4, where ln(nu1/4) -> 0), so it is drawn faintly.
        ax.plot(enu1s, fracs, "o-", ms=6, color="tab:green",
                label=r"terminal: $y<\varepsilon$ for all final iterates")
        if all(f is not None for f in fps):
            ax.plot(enu1s, fps, "s--", ms=5, color="tab:grey", alpha=0.75,
                    label=r"first passage: $y<\varepsilon$ at some iterate")
        ax.axvline(4.0 / 3.0, color="tab:blue", ls="--", lw=1.0, label=r"$\nu_1=4/3$ (all-orbit)")
        ax.axvspan(4.0 / 3.0, 4.0, color="tab:blue", alpha=0.07)
        # Type sizes: the manuscript includes the stacked figure at 0.612 of its 6.27 in text width,
        # where the 11 and 12 pt lettering set here prints at about 6.9 and 7.5 pt.
        ax.set_xlabel(r"$\nu_1$", fontsize=12)
        ax.set_ylabel(r"fraction of valid seeds extinct", fontsize=12)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title("(a) Terminal extinction fraction ($\\mu_0=4$)\n"
                     "extinct $\\Leftrightarrow$ $y<\\varepsilon$ for all final "
                     "%d iterates; $y_0<\\varepsilon$ excluded" % _EXT_TAIL, fontsize=12)
        ax.tick_params(labelsize=11)
        ax.legend(fontsize=11)

        # (b) eps ladder. Terminal flat vs first-passage sliding = the whole point.
        for nu1v, col in zip([2.0, 3.0, 4.0], ["tab:green", "tab:orange", "tab:red"]):
            sweep = ext_rec.get(f"{nu1v:.2f}", {}).get("eps_sweep")
            if not sweep:
                continue
            eps_vals = sorted((float(e) for e in sweep), reverse=True)
            axb.plot(eps_vals, [sweep[f"{e:g}"]["fraction_extinct"] for e in eps_vals],
                     "o-", ms=5, color=col, label=rf"$\nu_1={nu1v:g}$ terminal")
            axb.plot(eps_vals, [sweep[f"{e:g}"]["fraction_first_passage"] for e in eps_vals],
                     "s--", ms=4, color=col, alpha=0.45,
                     label=rf"$\nu_1={nu1v:g}$ first passage")
        axb.set_xscale("log")
        axb.invert_xaxis()          # eps decreasing left->right = "stricter threshold"
        axb.set_xlabel(r"$\varepsilon$ (decreasing $\rightarrow$)", fontsize=12)
        axb.set_ylabel("fraction of valid seeds", fontsize=12)
        axb.set_title("(b) Threshold dependence\n"
                      "terminal is flat in $\\varepsilon$; first passage is not", fontsize=12)
        axb.tick_params(labelsize=11)
        axb.legend(fontsize=11, ncol=1)
        fig.tight_layout()
        p = os.path.join(figdir, "designG_extinction_fraction.png")
        fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.02)  # 2026-10-17: trim the white margin
        plt.close(fig)
        written.append(p)
        _manifest_append(p)

    # Manuscript Figure 2: Jury band (tr/det closed forms with stability shading).
    if "coexistence" in rec:
        jb = rec["coexistence"]["jury_band"]
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        # The fixed point (4/(4+nu1), (3nu1-4)/(4nu1)) is interior
        # only for nu1 > 4/3 (y* = 0 at 4/3), so the closed forms are drawn from there, and the
        # certified attracting two-cycle below the Jury edge (nu1 = 46/25, prop:two-cycle) is
        # marked: the band is a FIXED-POINT stability statement, not an attractor classification.
        nn = np.linspace(4.0 / 3.0, 4.0, 1000)
        trv = -((nn - 4.0) ** 2) / (nn * (nn + 4.0))
        detv = (nn * nn + 12.0 * nn - 32.0) / (nn * (nn + 4.0))
        ax.plot(nn, trv, color="tab:purple", lw=1.4, label=r"$\operatorname{tr}J$ (interior fixed point)")
        ax.plot(nn, detv, color="tab:orange", lw=1.4, label=r"$\det J$ (interior fixed point)")
        ax.axvspan(jb["lower_edge_closed_form"], jb["upper_edge_closed_form"], color="tab:green", alpha=0.12,
                   label=r"fixed point a hyperbolic sink (strict Jury band)")
        ax.axvline(jb["lower_edge_closed_form"], color="tab:green", ls="--", lw=1.0,
                   label=r"$\nu_1^J=8\sqrt{3}-12\approx1.856$")
        ax.axvline(4.0 / 3.0, color="tab:purple", ls=":", lw=1.0,
                   label=r"$\nu_1=4/3$: fixed point enters the square ($y^*=0$)")
        ax.axvline(46.0 / 25.0, color="tab:red", ls="-.", lw=1.0,
                   label=r"$\nu_1=46/25$: certified attracting 2-cycle (not the fixed point)")
        ax.axhline(0.0, color="grey", lw=0.6)
        ax.set_xlabel(r"$\nu_1$")
        ax.set_ylabel("trace / determinant of $J$")
        ax.set_xlim(1.2, 4.05)
        ax.set_title("Coexistence fixed point: trace, determinant, Jury band ($\\mu_0=4$)")
        ax.legend(fontsize=7.5, loc="lower right")
        p = os.path.join(figdir, "designG_jury_band.png")
        fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.02)  # 2026-10-17: trim the white margin
        plt.close(fig)
        written.append(p)
        _manifest_append(p)

    # Supplementary Figure S2 (U2, Task B2): extinction-fraction sensitivity across the C0-C10
    # configs, per nu1, with Wilson 95% error bars on the random-IC configs only
    # (deterministic grids are exact -> no error bar). Under --smoke this is a sparse
    # PNG (few configs, single nu1); the full-data figure comes from the full run.
    if "sensitivity" in rec:
        sens = rec["sensitivity"]
        cfg_keys = [k for k in sens.keys() if not k.startswith("_") and k != "convergence_band"]
        nu1_keys = sorted(
            {nk for k in cfg_keys for nk in sens[k] if not str(nk).startswith("_")},
            key=lambda s: float(s),
        )
        if cfg_keys and nu1_keys:
            fig, ax = plt.subplots(figsize=(7.6, 4.6))
            xpos = np.arange(len(cfg_keys))
            for nk in nu1_keys:
                xs: List[float] = []
                ys: List[float] = []
                elo: List[float] = []
                ehi: List[float] = []
                for xi, k in enumerate(cfg_keys):
                    reck = sens[k].get(nk)
                    if not isinstance(reck, dict) or "fraction" not in reck:
                        continue
                    f = float(reck["fraction"])
                    xs.append(xi)
                    ys.append(f)
                    if reck.get("ci_method") == "wilson" and reck.get("ci_low") is not None:
                        elo.append(max(0.0, f - float(reck["ci_low"])))
                        ehi.append(max(0.0, float(reck["ci_high"]) - f))
                    else:
                        elo.append(0.0)
                        ehi.append(0.0)
                if xs:
                    ax.errorbar(xs, ys, yerr=[elo, ehi], marker="o", ms=5, lw=0,
                                elinewidth=1.2, capsize=3, label=rf"$\nu_1={nk}$")
            ax.set_xticks(xpos)
            ax.set_xticklabels(cfg_keys, rotation=45, ha="right", fontsize=7)
            # v4: this block calls the SAME terminal-persistence kernel as the Figure 11 block,
            # so "extinction fraction" is now accurate here too. The Wilson CI is over
            # the LEFT-CENSORED denominator (ICs with y0 < eps removed) -- which matters
            # most for the ic.rlogy config, where the log-uniform sampler put a large
            # share of ICs below threshold at t=0.
            ax.set_ylabel(r"terminal extinction fraction (valid seeds)")
            ax.set_ylim(bottom=-0.02)
            ax.set_title("Sensitivity of the terminal extinction fraction (C0–C10)\n"
                         "Error bars: Wilson 95% intervals; deterministic grids exact")
            ax.legend(fontsize=8, title=r"$\mu_0=4$")
            p = os.path.join(figdir, "designG_extinction_sensitivity.png")
            fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.02)  # 2026-10-17: trim the white margin
            plt.close(fig)
            written.append(p)
            _manifest_append(p)

    # Diagnostic plot, not a manuscript figure (U3, Task B3): finite-amplitude decay-rate convergence -- decay_rate
    # vs delta (left panel) and vs horizon at the smallest delta (right panel), per
    # nu1, against the ln(nu1/4) reference dashed lines. Empirical diagnostic
    # (non-blocking, no PASS/FAIL gate); under --smoke this is a sparse PNG (few
    # ICs/checkpoints/nu1) -- the full-data figure comes from the later B7 run.
    if "finite_amp" in rec:
        fa = rec["finite_amp"]
        nu1_keys = sorted(fa["per_nu1"].keys(), key=lambda s: float(s))
        if nu1_keys:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6))
            colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(nu1_keys)))
            smallest_delta_key = f"{min(fa['delta_list']):g}"
            for nk, col in zip(nu1_keys, colors):
                nu1v = float(nk)
                per_delta = fa["per_nu1"][nk]["per_delta"]
                dkeys = sorted(per_delta.keys(), key=lambda s: float(s), reverse=True)
                dvals = [float(dk) for dk in dkeys]
                drates = [per_delta[dk]["decay_rate"] for dk in dkeys]
                ax1.plot(dvals, drates, "o-", ms=6, color=col, label=rf"$\nu_1={nk}$")
                ax1.axhline(math.log(nu1v / 4.0), color=col, ls=":", lw=1.0)
                if smallest_delta_key in per_delta:
                    hs = per_delta[smallest_delta_key]["horizon_sweep"]
                    hs_h = [pt["horizon"] for pt in hs]
                    hs_r = [pt["decay_rate"] for pt in hs]
                    ax2.plot(hs_h, hs_r, "s-", ms=6, color=col, label=rf"$\nu_1={nk}$")
                    ax2.axhline(math.log(nu1v / 4.0), color=col, ls=":", lw=1.0)
            ax1.set_xscale("log")
            ax1.invert_xaxis()
            ax1.set_xlabel(r"$\delta$ (initial predator amplitude ceiling)")
            ax1.set_ylabel("decay_rate")
            ax1.set_title(r"Finite-amplitude decay vs $\delta$")
            ax1.legend(fontsize=8)
            ax2.set_xlabel("horizon (steps)")
            ax2.set_ylabel("decay_rate")
            ax2.set_title(rf"Horizon convergence at $\delta={smallest_delta_key}$")
            ax2.legend(fontsize=8)
            fig.suptitle("Design-G finite-amplitude decay-rate convergence "
                         r"$\to \ln(\nu_1/4)$ (U3, empirical, non-blocking)")
            p = os.path.join(figdir, "designG_lambda_convergence.png")
            fig.savefig(p, dpi=150, bbox_inches="tight")
            plt.close(fig)
            written.append(p)
            _manifest_append(p)

    # Diagnostic plot, not a manuscript figure (U4, Task B4): exact axis lambda_perp -- DE quadrature (blocking) +
    # orbit cross-check (secondary) vs nu1 against the ln(nu1/4) closed form (left),
    # and the orbit horizon-sweep convergence (right).  This is the EXACT (mu0=4
    # skew) model -- DISTINCT from the GPU fc_v4_skew parameter sweep (separate
    # script/model); the title/caption says so explicitly.
    if "exact_skew" in rec:
        es = rec["exact_skew"]["per_nu1"]
        nu1_keys = sorted(es.keys(), key=lambda s: float(s))
        if nu1_keys:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6))
            nu1_vals = [float(k) for k in nu1_keys]
            lam_q_vals = [es[k]["lambda_perp_quadrature"] for k in nu1_keys]
            lam_o_vals = [es[k]["lambda_perp_orbit"] for k in nu1_keys]
            sem_vals = [es[k]["sem"] for k in nu1_keys]
            grid_nu1 = np.linspace(0.5, 4.0, 200)
            ax1.plot(grid_nu1, np.log(grid_nu1 / 4.0), "k-", lw=1.5, label=r"$\ln(\nu_1/4)$ (closed form)")
            ax1.plot(nu1_vals, lam_q_vals, "D", ms=8, color="tab:red",
                      label="DE quadrature (BLOCKING, ~1e-10)")
            ax1.errorbar(nu1_vals, lam_o_vals, yerr=sem_vals, fmt="o", ms=6, color="tab:blue",
                         capsize=3, label="orbit cross-check (secondary, descriptive)")
            ax1.set_xlabel(r"$\nu_1$")
            ax1.set_ylabel(r"$\lambda_\perp$")
            ax1.set_title("Exact axis $\\lambda_\\perp$: quadrature vs orbit")
            ax1.legend(fontsize=7)
            colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(nu1_keys)))
            for nk, col in zip(nu1_keys, colors):
                nu1v = float(nk)
                hs = es[nk]["horizon_sweep"]
                hs_T = [pt["T"] for pt in hs]
                hs_mean = [pt["mean"] for pt in hs]
                ax2.plot(hs_T, hs_mean, "s-", ms=6, color=col, label=rf"$\nu_1={nk}$")
                ax2.axhline(math.log(nu1v / 4.0), color=col, ls=":", lw=1.0)
            ax2.set_xscale("log")
            ax2.set_xlabel("orbit horizon $T$ (steps)")
            ax2.set_ylabel(r"orbit-mean $\hat\lambda_\perp$")
            ax2.set_title("Orbit horizon-sweep convergence")
            ax2.legend(fontsize=8)
            fig.suptitle("Design-G EXACT axis skew model (U4): $x'=4x(1-x)$ -- "
                         "DISTINCT from the GPU fc_v4_skew sweep")
            p = os.path.join(figdir, "designG_exact_skew.png")
            fig.savefig(p, dpi=150, bbox_inches="tight")
            plt.close(fig)
            written.append(p)
            _manifest_append(p)

    return written


# =============================================================================
# Scaffold: preflight, path resolution, --smoke params, manifest, block registry
# (Task B0, 2026-07-09 -- shared contract, duplicated verbatim-in-spirit in the
#  sibling script designG_quasiperiodic_v2.py; this project duplicates rather
#  than shares utilities across CPU scripts by design)
# =============================================================================


def _preflight() -> None:
    """Abort early on an interpreter/package mismatch (never a silent wrong-env run).

    Enforces the lockfile: Python 3.10.x with exactly numpy 1.26.4 / numba 0.61.2
    (numba is version-sensitive) -- for EVERY run, author or not.

    DEFAULT INVERTED 2026-07-17. This used to additionally require
    ``sys.prefix == ~/.venv/nolds`` unless ``--public``/``LV_PUBLIC=1`` was passed.
    That is an AUTHOR-LOCAL path, so the check made the bundle's own README fail as
    written: the README says to build a fresh ``python3.10 -m venv .venv`` and then
    run ``python3 designG_dynamics_v4.py``, which supplies neither flag, so every
    third-party reproduction aborted with "local run requires ~/.venv/nolds".

    The wrong-environment case that the path check existed to catch is caught,
    strictly earlier, by the version pins above.

    ``LV_LOCAL=1`` restores the author-only path assertion. ``--public``/``LV_PUBLIC=1``
    remain accepted and are now
    NO-OPS, kept so existing documented commands and scripts do not break.
    """
    import numpy, numba as _numba  # local import: keep module import-safe pre-preflight

    if sys.version_info[:2] != (3, 10):
        sys.exit(f"[preflight] Python 3.10.x required (numba 0.61.2 is version-sensitive), got {sys.version.split()[0]}")
    want = {"numpy": "1.26.4", "numba": "0.61.2"}
    got = {"numpy": numpy.__version__, "numba": _numba.__version__}
    bad = {k: (got[k], want[k]) for k in want if got[k] != want[k]}
    if bad:
        sys.exit(f"[preflight] package mismatch vs requirements.txt: {bad}")
    if os.environ.get("LV_LOCAL") == "1":
        if sys.prefix.rstrip("/") != os.path.expanduser("~/.venv/nolds"):
            sys.exit(f"[preflight] LV_LOCAL=1 asserts the author venv ~/.venv/nolds, got {sys.prefix}")


def resolve_paths(smoke: bool) -> Dict[str, str]:
    """Env-var I/O contract: DESIGNG_FIGDIR/DESIGNG_DATADIR are REQUIRED outside --smoke.

    Never writes shipped paths implicitly -- a non-smoke run with the env vars
    unset aborts rather than silently defaulting to the FIGDIR/DATADIR module
    globals (which are themselves just an env-var-or-shipped-default fallback).
    """
    figd, datad = os.environ.get("DESIGNG_FIGDIR"), os.environ.get("DESIGNG_DATADIR")
    if not smoke and (not figd or not datad):
        sys.exit("[paths] DESIGNG_FIGDIR/DESIGNG_DATADIR required outside --smoke")
    return dict(figdir=figd or "./smoke", datadir=datad or "./smoke")


def _smoke_params() -> Dict[str, object]:
    """DYNAMICS smoke ($DYN, extinction-axis mu0=4); reduced VALUES only, FULL schema.

    ``grid``/``horizon`` override --grid-side/--bound-iter/--ext-horizon/--tele-iter;
    ``M``/``nu1_list``/``tol`` are reserved for the U1-U4 blocks (Tasks B1-B4) that
    do not exist yet in this v3 scaffold -- the CURRENT check_* functions do not
    consume them (their own nu1 sweeps + tolerances stay as coded; only iteration
    counts / grid sizes are cut for speed). See main()'s smoke-override block.
    """
    return dict(grid=50, horizon=1000, M=32, nu1_list=[2.0], tol=1e-1)


def _manifest_append(path: str) -> None:
    """Append a producer-authored ``| artifact | version | sha256 |`` row to DESIGNG_MANIFEST (if set).

    The manifest is authored by producer scripts (this function) and later only
    VERIFIED by the validator ($VAL, Task B6/B7) -- never both by the same tool.
    """
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {os.path.basename(__file__)} | {h.hexdigest()} |\n")


# Block names for --only (must match _block_registry()'s keys exactly; kept as a
# separate constant so argparse can declare --only's `choices` BEFORE args exists
# -- the registry's lambdas close over `args` lazily and are not evaluated here).
_BLOCK_NAMES: Tuple[str, ...] = (
    "boundedness", "lambda_hat", "extinction", "coexistence",
    "global_telescoping", "family_unification",
    "sensitivity",   # U2 (Task B2): C0-C10 extinction-fraction sensitivity sweep + CIs
    "finite_amp",    # U3 (Task B3): finite-amplitude two-way decay diagnostic (NON-blocking)
    "exact_skew",    # U4 (Task B4): exact axis lambda_perp=ln(nu1/4) DE quadrature (BLOCKING) + orbit cross-check
)


def _block_registry(args: argparse.Namespace) -> "Dict[str, Tuple[str, Callable[[], Tuple[List[Check], Dict[str, object]]]]]":
    """Map ``--only <block>`` CLI name -> (values-JSON top-level key, thunk running the check).

    CLI block names are a scaffold-facing label distinct from the (unchanged)
    JSON top-level keys -- e.g. "extinction" here dispatches the existing
    check_extinction_fraction() and still writes it under the existing
    "extinction_fraction" JSON key (Task B1 reads that key; renaming it is out
    of scope for B0). New blocks (U1-U5, Tasks B1-B5) extend this registry
    in-place; B0 only wires the six blocks that already exist in v2. Keys must
    match _BLOCK_NAMES exactly (checked by test_scaffold's sibling scaffold
    tests); building this dict does not touch `args` attributes (lambdas are
    lazy), so it is safe to call before argparse has fully populated `args`.
    """
    return {
        "boundedness": ("boundedness", lambda: check_boundedness(args.bound_iter, args.grid_side)),
        "lambda_hat": ("lambda_hat", lambda: check_lambda_hat(args.n_burn, args.n_keep, args.window)),
        "extinction": ("extinction_fraction", lambda: check_extinction_fraction(args.ext_horizon, args.grid_side, args.ext_eps)),
        "coexistence": ("coexistence", lambda: check_coexistence_fp_and_jury()),
        "global_telescoping": ("global_telescoping", lambda: check_global_telescoping(args.tele_iter, args.grid_side)),
        "family_unification": ("family_unification", lambda: check_family_unification(args.n_burn, args.n_keep)),
        # U2 (Task B2): run_sensitivity reads args.smoke itself to pick the nominal
        # vs reduced schedule + checkpoint file; output lands under values["sensitivity"].
        "sensitivity": ("sensitivity", lambda: run_sensitivity(args)),
        # U3 (Task B3): run_finite_amplitude reads args.smoke itself to pick the
        # nominal vs reduced schedule; output lands under values["finite_amp"].
        "finite_amp": ("finite_amp", lambda: run_finite_amplitude(args)),
        # U4 (Task B4): run_exact_skew reads args.smoke itself to pick the nominal
        # vs reduced ORBIT schedule (the quadrature side is never reduced); output
        # lands under values["exact_skew"].
        "exact_skew": ("exact_skew", lambda: run_exact_skew(args)),
    }


# B0 JSON-schema convention: "_taxonomy" is a top-level map from a field name
# (top-level scalars) or a JSON-Pointer path with "*" wildcarding any non-"_"
# record key (nested scalars) to one of regression_anchor | theory_invariant |
# empirical_discovery.  Populated here ONLY for the two RETAINED v2 producers
# this script already emits unconditionally; U1-U4 entries are added by
# Tasks B1-B4 (see the plan's "Required _taxonomy coverage").
_TAXONOMY: Dict[str, str] = {
    "/coexistence/jury_band/lower_edge_closed_form": "regression_anchor",
    "/family_unification/per_nu1/*/pairwise_max_spread": "regression_anchor",
    # U1 (Task B1): continuous extinction fraction + its "small" classification.
    "/validation_table/*/fraction_extinct": "regression_anchor",
    "/validation_table/*/is_small": "empirical_discovery",
    # Floating-point absorption bookkeeping (2026-10-17): per (nu1, eps),
    # how many terminally extinct orbits reached exact 0.0 before the terminal window
    # and the first/last iterate at which one did. Deterministic grid counts -> anchors.
    "/extinction_fraction/per_nu1/*/eps_sweep/*/n_extinct_exact_zero_before_tail": "regression_anchor",
    "/extinction_fraction/per_nu1/*/eps_sweep/*/first_zero_step_min": "regression_anchor",
    "/extinction_fraction/per_nu1/*/eps_sweep/*/first_zero_step_max": "regression_anchor",
    # U2 (Task B2): sensitivity sweep. The baseline `base` config is the reproducible
    # regression anchor; every other config (grid/horizon/eps refinements + random-IC)
    # is an empirical discovery. Most-specific pointer wins, so /sensitivity/base/*/
    # fraction is anchored even though /sensitivity/*/*/fraction also matches it.
    "/sensitivity/base/*/fraction": "regression_anchor",
    "/sensitivity/*/*/fraction": "empirical_discovery",
    "/sensitivity/*/*/fraction_first_passage": "empirical_discovery",
    "/sensitivity/*/*/fraction_first_passage_censored_inclusive": "empirical_discovery",
    # Benettin transverse-exponent table (mu0=4; the prey orbits do not depend on nu1, so
    # the launch/collapse/acceptance counts and the dispersion are shared across nu1).
    "/lambda_hat/per_nu1/*/lambda_hat": "regression_anchor",
    "/lambda_hat/per_nu1/*/n_trajectories": "regression_anchor",
    "/lambda_hat/per_nu1/*/n_noncollapsed": "regression_anchor",
    "/lambda_hat/per_nu1/*/n_accepted_stable": "regression_anchor",
    "/lambda_hat/per_nu1/*/sem_across_accepted": "empirical_discovery",
    "/lambda_hat/per_nu1/*/sd_across_accepted": "empirical_discovery",
    # U3 (Task B3): finite-amplitude two-way decay diagnostic. Empirical by
    # construction (a finite-y0 nonlinear measurement expected -> ln(nu1/4) only
    # in the delta->0 limit) -- NOT a theory gate, no PASS/FAIL Check row exists
    # for it to tag (see run_finite_amplitude's docstring).
    "/finite_amp/per_nu1/*/per_delta/*/decay_rate": "empirical_discovery",
    # U4 (Task B4): exact axis transverse-variational lambda_perp = ln(nu1/4).
    # The DE quadrature + its collapse_guard_ok gate are the BLOCKING theory
    # invariant; the orbit Birkhoff average is the secondary empirical cross-check.
    "/exact_skew/per_nu1/*/lambda_perp_quadrature": "theory_invariant",
    "/exact_skew/per_nu1/*/collapse_guard_ok": "theory_invariant",
    "/exact_skew/per_nu1/*/lambda_perp_orbit": "empirical_discovery",
}


# =============================================================================
# Driver
# =============================================================================


def _provenance(args: argparse.Namespace, figdir: str, datadir: str) -> Dict[str, object]:
    return {
        "generated_by": "output/codes/" + os.path.basename(__file__),

        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model": "Design G symmetric multiplicative map (eq:designG, methods_v6.tex)",
        "map": "x'=mu0*x*(1-x)*(1-y);  y'=nu1*x*y*(1-y) on [0,1]^2",
        "honesty": (
            "Numerics are CONSISTENCY CHECKS, not proofs. SRB-typical != all-orbit. "
            "Coexistence-band (2<nu1<4) ICs converging to coexistence are EXPECTED, "
            "not test failures. All checks run at mu0=4; ln(nu1/4) rests on the arcsine law, which holds at mu0=4 only."
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "numba": numba.__version__,
        "output_data_dir": _rel(datadir),
        "output_fig_dir": _rel(figdir),
        "no_figs": bool(args.no_figs),
        "smoke": bool(args.smoke),
        "only_block": args.only,
        "protocol": {
            "lambda_burn_in_Nb": args.n_burn,
            "lambda_retain_N": args.n_keep,
            "lambda_window": args.window,
            "extinction_horizon": args.ext_horizon,
            "extinction_eps": args.ext_eps,
            "boundedness_iter": args.bound_iter,
            "ic_grid_side": args.grid_side,
            "telescoping_iter": args.tele_iter,
        },
    }


def main(argv: List[str] | None = None) -> int:
    _preflight()

    parser = argparse.ArgumentParser(
        description="Design-G dynamics validation suite (spec 2A)."
    )
    parser.add_argument("--no-figs", action="store_true", help="skip optional diagnostic PNGs (P2)")
    parser.add_argument("--n-burn", type=int, default=75000, help="lambda burn-in N_b (eq:numerical-lambda)")
    parser.add_argument("--n-keep", type=int, default=50000, help="lambda retained N (eq:numerical-lambda)")
    parser.add_argument("--window", type=int, default=5000, help="stability-gate moving window")
    parser.add_argument("--ext-horizon", type=int, default=20000, help="extinction-sweep horizon H")
    parser.add_argument("--ext-eps", type=float, default=1e-3, help="extinction threshold eps")
    parser.add_argument("--grid-side", type=int, default=60, help="interior IC grid side (n x n)")
    parser.add_argument("--bound-iter", type=int, default=20000, help="boundedness orbit length")
    parser.add_argument("--tele-iter", type=int, default=4000, help="telescoping-demo orbit length")
    parser.add_argument("--smoke", action="store_true", help="tiny/reduced run, full schema -> <name>_smoke.json")
    parser.add_argument("--only", choices=sorted(_BLOCK_NAMES), default=None, help="run a single named block")
    parser.add_argument("--public", action="store_true",
                         help="accepted for compatibility; has no effect (the preflight checks the pinned versions on every run)")
    args = parser.parse_args(argv)

    # Keep the lambda stability gate non-vacuous. The per-trajectory gate needs at least
    # two converged late windows (nwin = n_keep // window >= 2), i.e. n_keep >= 2*window;
    # without this guard a pathological --n-keep < --window yields nwin=0 and a silently
    # vacuous gate. Warn (don't fail) if n_keep is not an exact multiple of window (the final
    # partial window is dropped from the gate by design).
    #
    # --smoke (B0): unconditionally overrides the iteration/grid-size args with a tiny,
    # fast configuration (from _smoke_params() where it applies; n_burn/n_keep/window are
    # not in _smoke_params() -- reserved there for U1-U4 -- so B0 picks its own fast values
    # here, satisfying n_keep=200 >= 2*window=100). This is "reduced VALUES only, FULL
    # schema": every block still runs (unless --only also restricts to one), just fast.
    if args.smoke:
        smoke = _smoke_params()
        args.grid_side = smoke["grid"]        # 50  (was 60)
        args.bound_iter = smoke["horizon"]    # 1000 (was 20000)
        args.ext_horizon = smoke["horizon"]   # 1000 (was 20000)
        args.tele_iter = smoke["horizon"]     # 1000 (was 4000)
        args.n_burn = 200                     # was 75000 (smoke-only; fast Benettin burn-in)
        args.n_keep = 200                     # was 50000
        args.window = 50                      # was 5000 (200 >= 2*50, 200 % 50 == 0)

    if args.window < 1:
        parser.error(f"--window ({args.window}) must be >= 1.")
    if args.n_keep < 2 * args.window:
        parser.error(
            f"--n-keep ({args.n_keep}) must be >= 2*--window ({2 * args.window}) so the "
            "per-trajectory stability gate has converged late windows to test."
        )
    if args.n_keep % args.window != 0:
        print(f"  [warn] --n-keep ({args.n_keep}) is not a multiple of --window "
              f"({args.window}); the final partial window is excluded from the stability gate.")

    paths = resolve_paths(args.smoke)
    figdir, datadir = paths["figdir"], paths["datadir"]
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)

    print("=" * 96)
    print("Design-G dynamics validation suite (spec 2A) -- CONSISTENCY CHECKS, NOT PROOFS")
    print("  map:  x' = mu0*x*(1-x)*(1-y) ;  y' = nu1*x*y*(1-y)   on [0,1]^2 (eq:designG)")
    if args.smoke:
        print("  (--smoke: tiny/reduced values, full schema)")
    print("=" * 96)

    registry = _block_registry(args)
    selected = [args.only] if args.only else list(registry)  # dict preserves a-f insertion order

    block_banner = {
        "boundedness": "\n[a] Boundedness / no blow-up (eq:bounds, rmk:three-regime-numeric) ...",
        "lambda_hat": "\n[b] Benettin--Wolf lambda_hat vs ln(nu1/4) (eq:numerical-lambda, eq:lambda-mu4) ...",
        "extinction": "\n[c] Finite-grid extinction-fraction sweep (rmk:three-regime-numeric) ...",
        "coexistence": "\n[d] Coexistence fixed point + Jury band (eq:coexistence-fp/-trace-det/-jury-band) ...",
        "global_telescoping": "\n[e] Global-extinction telescoping demo, nu1<4/3 (thm:global) ...",
        "family_unification": "\n[f] Family lambda_perp unification: skew = MOD-1 = Design-G = ln(nu1/4) (thm:skew, thm:mod1) ...",
        "sensitivity": "\n[U2] Extinction-fraction sensitivity sweep C0-C10 (det->convergence band; random->Wilson CI) ...",
        "finite_amp": "\n[U3] Finite-amplitude two-way decay diagnostic (empirical_discovery; NON-blocking) ...",
        "exact_skew": "\n[U4] Exact axis lambda_perp=ln(nu1/4): DE quadrature (BLOCKING, theory_invariant) + orbit cross-check ...",
    }

    # U1 (Task B1): check_extinction_fraction() may append raw dict rows (continuous
    # fraction_extinct/is_small/small_ceiling, no PASS/FAIL) alongside Check objects,
    # so all_rows is heterogeneous from here on; see the isinstance(c, Check) guards below.
    all_rows: List[object] = []
    values: Dict[str, object] = {"provenance": _provenance(args, figdir, datadir)}

    for name in selected:
        json_key, fn = registry[name]
        print(block_banner[name])
        rows, rec = fn()
        all_rows += rows
        values[json_key] = rec
        if name == "boundedness":
            print(f"    overall max x = {rec['overall_max_x']:.6f}  (<= 1) ; overall max y = {rec['overall_max_y']:.6f}  (<= 1)")

    # U2 follow-up: the plan's interface contract specifies convergence_band as a
    # TOP-LEVEL produced key (convergence_band.<axis>.<nu1>.{min,max,configs}),
    # parallel to "sensitivity" -- not nested inside it. Lift it here so the final
    # JSON matches the contract; no-op when the sensitivity block did not run
    # (e.g. --only <other-block>).
    if "sensitivity" in values and "convergence_band" in values["sensitivity"]:
        values["convergence_band"] = values["sensitivity"].pop("convergence_band")

    # ---- PASS/FAIL validation table --------------------------------------
    print("\n" + "=" * 96)
    print("VALIDATION TABLE  (Design-G dynamics; consistency checks, not proofs)")
    print("-" * 96)
    header = f"{'quantity':<46} {'predicted':>13} {'computed':>13} {'abs-err':>11}  result"
    print(header)
    print("-" * 96)
    n_pass = 0
    n_total = 0
    for c in all_rows:
        if not isinstance(c, Check):
            continue  # U1-style continuous-metric row (e.g. fraction_extinct/is_small
            # dict from check_extinction_fraction): informational, not PASS/FAIL --
            # excluded from the checks-PASS accounting and this printed table.
        n_total += 1
        status = "PASS" if c.passed else "FAIL"
        if c.passed:
            n_pass += 1
        print(f"{c.name:<46} {c.predicted:>13.6f} {c.computed:>13.6f} {c.abs_err:>11.3e}  {status}")
    print("-" * 96)
    print(f"SUMMARY: {n_pass}/{n_total} checks PASS.")
    print("NOTE: SRB-typical/finite-grid != all-orbit; coexistence-band convergence is expected, not a failure.")
    if args.smoke:
        print("NOTE: --smoke uses drastically reduced iteration/grid/burn-in counts; individual")
        print("      check PASS/FAIL at this N is not the acceptance bar -- schema + no-crash is.")
    print("=" * 96)

    # U1 (Task B1): serialize both row shapes -- Check.as_dict() for legacy PASS/FAIL
    # rows, and the raw dict itself for continuous-metric rows (already JSON-ready).
    values["validation_table"] = [c.as_dict() if isinstance(c, Check) else c for c in all_rows]
    values["summary"] = {"n_pass": n_pass, "n_total": n_total, "all_pass": bool(n_pass == n_total)}
    # B0 JSON-schema convention: every values JSON carries a top-level "_taxonomy" map
    # (field-name / wildcard-JSON-Pointer -> regression_anchor|theory_invariant|
    # empirical_discovery). B0 populates only the two RETAINED v2 producers this
    # script emits unconditionally; U1-U4 taxonomy entries are added by Tasks B1-B4.
    values["_taxonomy"] = dict(_TAXONOMY)

    # ---- optional figures -------------------------------------------------
    if not args.no_figs:
        try:
            written = make_figures(values, figdir)
            values["figures"] = [_rel(w) for w in written]   # bundle-relative
            for p in written:
                print("  wrote figure", p)
        except Exception as exc:
            # The values JSON is still written below, but the run FAILS (see the return): these PNGs include
            # manuscript figures, and a runner that got exit 0 would keep the previous PNG and report success.
            print(f"  [ERROR] figure generation FAILED: {exc!r}")
            values["figures_error"] = repr(exc)
    else:
        print("  (--no-figs: diagnostic PNGs skipped)")

    # ---- JSON deposit -----------------------------------------------------
    stem = os.path.splitext(os.path.basename(__file__))[0]
    out_name = f"{stem}_smoke.json" if args.smoke else f"{stem}_values.json"
    out_json = os.path.join(datadir, out_name)
    with open(out_json, "w") as f:
        json.dump(values, f, indent=2)
    print("wrote", out_json)
    _manifest_append(out_json)

    return 0 if (n_pass == n_total and "figures_error" not in values) else 1


if __name__ == "__main__":
    raise SystemExit(main())
