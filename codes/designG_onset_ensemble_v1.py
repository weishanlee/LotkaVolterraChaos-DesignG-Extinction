#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G onset ensemble + resolution study (P3 item 4).

WHY THIS EXISTS
---------------
The single-seed sweep of `designG_bifurcation_v2.py` reports "chaos onset = predator
extinction onset" for one orbit, (x0,y0) = (0.2, 1e-3). The onset moves
under initial-condition changes, so a coincidence asserted from one seed is not
established. This script therefore measures both onsets over 24 seeds x 3 grid
resolutions and reports the spread.

RULES. The CURRENT detector, shared with designG_bifurcation_v2.py::detect_onsets,
is the "sustained" rule: chaos onset = the least mu0 beyond which lambda_max > 0 at
EVERY later swept point; extinction onset = the least mu0 beyond which the retained
predator stays below Y_EXT at every later swept point. Two other readouts are kept
for the bracket they define: "first positive anywhere" (fires at the first transient
sliver) and the HISTORICAL ">60% of later points" rule once used by the bifurcation
script (field `chaos_legacy_60pct`; no longer used anywhere else). Inside a periodic
window lambda_max dips back below zero, so a single onset number is genuinely
ill-defined; the point of this file is to report the bracket and how far the
answer moves across seeds and resolutions. The mu0 grid is rounded to the nearest
double of its decimal value (shared with the bifurcation script), and in both scripts the
largest Lyapunov exponent is the larger of two finite-time QR (Gram-Schmidt) estimates
propagated from the orthonormal frame ((1, 0), (0, 1)); at finite length the larger of the
two accumulated rates is not in general the growth rate of the largest singular value of the
whole product.

WHAT THIS COMPUTES
------------------
Over an ensemble of interior seeds and three sweep resolutions, per seed:

  chaos_first_positive : min{ mu0 : lambda_max(mu0) > 0 }
                         earliest evidence of chaos; sensitive to slivers.
  chaos_sustained      : min{ mu0 : lambda_max(mu) > 0 for ALL grid mu >= mu0 }
                         strict; sits above the last periodic window.
                         [chaos_first_positive, chaos_sustained] BRACKETS the onset.
  chaos_legacy_60pct   : the v2 rule, reported ONLY so the old number is auditable.
  ext_onset            : min{ mu0 : max(y over the kept window) < Y_EXT for ALL
                         grid mu >= mu0 }. Already a terminal-window criterion
                         (the whole kept window must stay below Y_EXT), so unlike
                         the extinction FRACTION of designG_dynamics_v3 this one
                         did not carry the first-passage defect.

and reports, across the ensemble, the median and full range of each, plus the
signed gap ext_onset - chaos_sustained.

HONESTY
-------
This script settles ONE question: how stable is the reported onset, and is the
"coincident collapse" claim robust to the seed? It does NOT establish causation
in either direction.

CPU only, numba JIT. Run (~2 min) in the bundle's codes directory:
    DESIGNG_FIGDIR=<dir> DESIGNG_DATADIR=<dir> \
      python designG_onset_ensemble_v1.py
(run_designG_figs.sh does this).
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

try:
    from numba import njit
except ModuleNotFoundError:
    def njit(*a, **k):
        if len(a) == 1 and callable(a[0]) and not k:
            return a[0]
        return lambda f: f

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.dirname(_HERE)
DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))
FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))

NU1 = float(os.environ.get("DESIGNG_NU1", "2.1"))   # same representative as v2
N_BURN = 75_000
N_KEEP = 400
Y_EXT = 1e-3
LYAP_BURN = 20_000
LYAP_ITER = 20_000
SEED = 20260717
N_SEEDS = int(os.environ.get("DESIGNG_ONSET_SEEDS", "24"))
MU0_LO, MU0_HI = 3.0, 4.0
RESOLUTIONS = (0.01, 0.005, 0.002)


@njit(cache=True)
def _lyap_spectrum_2d(mu0, nu1, x0, y0, n_burn, n_iter):
    """Both Lyapunov exponents of the FULL 2-D Design-G map (Benettin two-vector Gram-Schmidt).

    The first vector starts at (1, 0) with the arithmetic of the previous single-vector helper,
    operation for operation (so its exponent is bit-identical to it on every orbit on which no
    norm is exactly zero); the second starts at (0, 1) and is orthogonalised against the
    renormalised first vector after each step. On the invariant axis y = 0 the first vector
    stays tangential and the second recovers the transverse rate.
    Exact-zero convention as in designG_bifurcation_v2.m0_lyapunov: a
    first norm exactly 0 gives lam1 = -inf and the run continues with the second vector alone;
    a second norm exactly 0 gives lam2 = -inf. Returns (nan, nan) if the orbit escapes the
    square during burn-in; a later escape ends the average over the steps taken (as before)."""
    x, y = x0, y0
    for _ in range(n_burn):
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x, y = xn, yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            return np.nan, np.nan
    vx, vy = 1.0, 0.0
    ux, uy = 0.0, 1.0
    s1 = 0.0; s2 = 0.0
    m = 0
    alive1 = True; alive2 = True
    for _ in range(n_iter):
        # J = [[mu0(1-2x)(1-y), -mu0 x(1-x)], [nu1 y(1-y), nu1 x(1-2y)]]
        a = mu0 * (1.0 - 2.0 * x) * (1.0 - y)
        b = -mu0 * x * (1.0 - x)
        c = nu1 * y * (1.0 - y)
        d = nu1 * x * (1.0 - 2.0 * y)
        if alive1:
            nvx = a * vx + b * vy
            nvy = c * vx + d * vy
            nrm = np.sqrt(nvx * nvx + nvy * nvy)
            if not np.isfinite(nrm):
                return np.nan, np.nan          # invalid evolution: no exponent is defined
            if nrm == 0.0:
                alive1 = False                 # direction annihilated exactly: lam1 = -inf
            else:
                vx, vy = nvx / nrm, nvy / nrm
                s1 += np.log(nrm)
        if alive2:
            nux = a * ux + b * uy
            nuy = c * ux + d * uy
            if alive1:
                dot = nux * vx + nuy * vy
                nux = nux - dot * vx
                nuy = nuy - dot * vy
            nrm2 = np.sqrt(nux * nux + nuy * nuy)
            if not np.isfinite(nrm2):
                return np.nan, np.nan
            if nrm2 == 0.0:
                alive2 = False                 # direction annihilated exactly: lam2 = -inf
            else:
                ux, uy = nux / nrm2, nuy / nrm2
                s2 += np.log(nrm2)
        m += 1                                 # this step counts, including the one that annihilates the last direction
        if not alive1 and not alive2:
            break                              # both exponents are -inf (zero Jacobian at (0, 1))
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x, y = xn, yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            break
    if m == 0:
        return np.nan, np.nan
    lam1 = s1 / m if alive1 else -np.inf
    lam2 = s2 / m if alive2 else -np.inf
    return lam1, lam2


@njit(cache=True)
def _lyap_max_2d(mu0, nu1, x0, y0, n_burn, n_iter):
    """Largest Lyapunov exponent of the FULL 2-D Design-G map: the max of the two exponents of
    _lyap_spectrum_2d (largest also on the invariant axis y = 0). NaN if the orbit escapes the
    square during burn-in."""
    l1, l2 = _lyap_spectrum_2d(mu0, nu1, x0, y0, n_burn, n_iter)
    if np.isnan(l1) and np.isnan(l2):
        return np.nan
    return max(l1, l2)


@njit(cache=True)
def _ymax_kept(mu0, nu1, x0, y0, n_burn, n_keep):
    """Max predator value over the post-burn-in kept window (-1.0 if the orbit escapes)."""
    x, y = x0, y0
    for _ in range(n_burn):
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x, y = xn, yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            return -1.0
    ym = 0.0
    for _ in range(n_keep):
        xn = mu0 * x * (1.0 - x) * (1.0 - y)
        yn = nu1 * x * y * (1.0 - y)
        x, y = xn, yn
        if x < 0.0 or x > 1.0 or y < 0.0 or y > 1.0:
            return -1.0
        if y > ym:
            ym = y
    return ym


@njit(cache=True)
def _scan(grid, nu1, x0, y0):
    """lambda_max and ymax over the mu0 grid for one seed."""
    n = grid.size
    lam = np.empty(n)
    ym = np.empty(n)
    for i in range(n):
        lam[i] = _lyap_max_2d(grid[i], nu1, x0, y0, LYAP_BURN, LYAP_ITER)
        ym[i] = _ymax_kept(grid[i], nu1, x0, y0, N_BURN, N_KEEP)
    return lam, ym


def _onsets(grid, lam, ym):
    """Bracket the chaos onset; locate the extinction onset. All rules explicit."""
    pos = np.zeros(grid.size, dtype=bool)
    for i in range(grid.size):
        pos[i] = bool(np.isfinite(lam[i]) and lam[i] > 0.0)

    first_positive = float(grid[np.argmax(pos)]) if pos.any() else None

    sustained = None
    for i in range(grid.size):
        if pos[i:].all():
            sustained = float(grid[i])
            break

    legacy = None                                   # the v2 >60% rule, for audit only
    for i in range(grid.size):
        tail = lam[i:]
        finite = np.isfinite(tail)
        if pos[i] and finite.any() and float(np.mean(tail[finite] > 0.0)) > 0.6:
            legacy = float(grid[i])
            break

    below = ym < Y_EXT                              # note: escaped orbits (-1.0) count as below
    ext = None
    for i in range(grid.size):
        if below[i:].all():
            ext = float(grid[i])
            break
    return first_positive, sustained, legacy, ext


def main() -> int:
    os.makedirs(DATADIR, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # Seed 0 is v2's exact seed so the published number is reproduced inside the
    # ensemble rather than merely compared against from memory.
    seeds = [(0.2, 1e-3)]
    while len(seeds) < N_SEEDS:
        seeds.append((float(rng.uniform(0.05, 0.95)), float(10.0 ** rng.uniform(-4, -2))))

    out = {
        "nu1": NU1, "mu0_range": [MU0_LO, MU0_HI], "y_ext": Y_EXT,
        "n_burn": N_BURN, "n_keep": N_KEEP,
        "lyap_burn": LYAP_BURN, "lyap_iter": LYAP_ITER,
        "seed_rng": SEED, "n_seeds": len(seeds),
        "seed_0_is_v2_published_seed": [0.2, 1e-3],
        "resolutions": list(RESOLUTIONS),
        "per_resolution": {},
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }

    for dm in RESOLUTIONS:
        # nearest double of each DECIMAL grid value: identical float64 parameters in every script
        # that sweeps mu0 (designG_bifurcation_v2.py uses the same rounding)
        grid = np.round(np.arange(MU0_LO, MU0_HI + 1e-9, dm), 10)
        rows = []
        for (x0, y0) in seeds:
            lam, ym = _scan(grid, NU1, x0, y0)
            fp, sus, leg, ext = _onsets(grid, lam, ym)
            rows.append({"x0": x0, "y0": y0, "chaos_first_positive": fp,
                         "chaos_sustained": sus, "chaos_legacy_60pct": leg,
                         "ext_onset": ext})
        out["per_resolution"][f"{dm:g}"] = {"n_grid": int(grid.size), "per_seed": rows}

        def stat(key):
            v = [r[key] for r in rows if r[key] is not None]
            if not v:
                return None
            return {"median": float(np.median(v)), "min": float(min(v)),
                    "max": float(max(v)), "n": len(v)}

        gaps = [r["ext_onset"] - r["chaos_sustained"] for r in rows
                if r["ext_onset"] is not None and r["chaos_sustained"] is not None]
        out["per_resolution"][f"{dm:g}"]["summary"] = {
            "chaos_first_positive": stat("chaos_first_positive"),
            "chaos_sustained": stat("chaos_sustained"),
            "chaos_legacy_60pct": stat("chaos_legacy_60pct"),
            "ext_onset": stat("ext_onset"),
            "gap_ext_minus_chaos_sustained": (
                {"median": float(np.median(gaps)), "min": float(min(gaps)),
                 "max": float(max(gaps)), "n": len(gaps)} if gaps else None),
        }

        s = out["per_resolution"][f"{dm:g}"]["summary"]
        print(f"\n=== d(mu0) = {dm:g}   ({grid.size} grid points, {len(seeds)} seeds) ===")
        for k in ("chaos_first_positive", "chaos_sustained", "chaos_legacy_60pct",
                  "ext_onset", "gap_ext_minus_chaos_sustained"):
            v = s[k]
            print(f"  {k:32s} " + ("n/a" if v is None else
                  f"median {v['median']:+.4f}   range [{v['min']:+.4f}, {v['max']:+.4f}]   n={v['n']}"))

    p = os.path.join(DATADIR, "designG_onset_ensemble_v1.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
