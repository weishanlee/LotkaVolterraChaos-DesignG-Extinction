#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G stable-leaf illustration: the objects the positive-area proof constructs (mu0 = 4).

For the Design-G map Phi(x, y) = (4x(1-x)(1-y), nu1 x y (1-y)) this script computes, at one
chosen axis point p (its Birkhoff genericity is not established), two objects that the proof of
the positive-area extinction theorem
(Section 3, thm:basin-positive, via thm:transverse-criterion and cor:rare-predator) builds:

  (a) the local stable leaf through (p, 0), as the fixed point of the paper's forward
      Lyapunov--Perron operator (eq:LP-stable, eq:LP-unstable) on a truncated horizon,
      z(t) = (xi, zeta), leaf point (p + s_0 t + xi_0(t), t), with the slope s_0 of the slope
      series eq:slope-series as the tangent;
  (b) the matched reference product of cor:rare-predator: along the orbit started at the leaf
      point of height t, the ratio y_n / (t V_n(p)), V_n(p) = prod_{j<n} nu1 f^j(p), which the
      corollary bounds by |ratio - 1| < 1/107 on the theorem's set (40 steps, float64);
  (c) the long-run rate of that reference product, (1/N) ln V_N(p), along the EXACT orbit of p up
      to N = 1e5: the base point is p = sin^2(pi theta0 / 2) with theta0 = sqrt(2) - 1, whose
      binary expansion is known to every bit (integer square root), and the logistic map at 4 is
      conjugate to the tent map through x = sin^2(pi theta / 2), so theta_n is computed exactly in
      integer arithmetic and x_n from it in float64. The closed form ln(nu1/4) is the limit of
      this Birkhoff average for almost every point under the arcsine law; the finite averages at
      this particular p are compared with it, but the genericity of p, and the convergence of its
      averages to ln(nu1/4), are not established.

This is an ILLUSTRATION, not an instance of the theorem: the base point p and the heights shown
are not certified to lie in the theorem's block B x [0, c_0(nu1)], no constant of the proof is
computed here, and nothing below is used in any proof. The Lyapunov--Perron iteration
is run in floating point on a horizon of N steps; its convergence is checked and reported (with
the observed ratios of successive changes, which are not a certificate that the operator
contracts on a neighbourhood), and a height at which it does not converge is refused
(RuntimeError), never plotted.

Outputs (portable, mirror the sibling designG_*.py scripts):
    figures/designG_stable_leaf.png                       (DESIGNG_FIGDIR; Figure 5: leaf + tracking)
    figures/designG_stable_leaf_rates.png                 (DESIGNG_FIGDIR; Figure 12: the two decay rates)
    data/designG/designG_stable_leaf_v1_values.json       (DESIGNG_DATADIR)
and, if DESIGNG_MANIFEST is set, one producer row per artifact (| artifact | version | sha256 |).

Run (CPU only, NumPy only, no random numbers; a few seconds), in the bundle's codes directory:
    python designG_stable_leaf_v1.py [--no-figs]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import struct
import sys
from datetime import datetime, timezone

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))   # .../codes
_OUT = os.path.dirname(_HERE)                         # .../output  (or the bundle root)
DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))
FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))

VERSION = "designG_stable_leaf_v1"

MU0 = 4.0
THETA0_LABEL = "sqrt(2) - 1"      # conjugate angle of the base point, known to every bit
# Frozen binary64 inputs (2026-10-21). The base point and the height grid are
# exact hexadecimal literals, not recomputed at run time: on an AVX-512 machine NumPy's np.logspace
# (SIMD power) returned 0x1.2a5884e341e82p-16 instead of ...83p-16 at index 17, and the observed
# contraction ratio of lp_fixed_point amplified that last bit beyond the reproducer's tolerance; libm's
# sin may likewise differ by one ulp across platforms. The literals are the values that produced every
# shipped number -- sin(pi (sqrt2 - 1)/2)**2 and np.logspace(-9, -3, 25) as evaluated on the reference
# machine (2026-10-21, NumPy 1.26.4) -- and _check_frozen_inputs() below refuses an edited literal
# (SHA-256 of the packed array) and confirms the provenance to within a few ulp of those formulas
# (4 ulp for the heights; 8 ulp for P_BASE, since squaring roughly doubles a one-ulp difference in sin).
P_BASE_HEX = "0x1.77ad613b17f9cp-2"
P_BASE = float.fromhex(P_BASE_HEX)  # = sin^2(pi (sqrt2 - 1)/2) ~ 0.3669; irrational angle, never pre-critical
NU1_LIST = (2.1, 3.0)             # the manuscript's two showcase conversion rates
N_HORIZON = 60                    # truncation horizon of the Lyapunov--Perron sums
N_MATCH = 40                      # steps of the matched-product window (cor:rare-predator)
N_LONG = 100_000                  # exact-orbit horizon of the long-run reference rate (Figure 12)
EXTRA_BITS = 80                   # guard bits of the integer tent-map orbit beyond N_LONG
T_MATCH = 1e-6                    # height of the orbit used for the matched product
HEIGHTS_HEX = (                   # heights at which the leaf point is computed: 10**(-9 + k/4), k = 0..24
    "0x1.12e0be826d695p-30", "0x1.e8cf49f422410p-30", "0x1.b29ec0fbe4534p-29",
    "0x1.827021d36920bp-28", "0x1.5798ee2308c3ap-27", "0x1.31818e389568ap-26",
    "0x1.0fa3389d6eb40p-25", "0x1.e30c2a484368dp-25", "0x1.ad7f29abcaf48p-24",
    "0x1.7de1f1c6bac2cp-23", "0x1.538c06c4ca610p-22", "0x1.2de79a6d2a218p-21",
    "0x1.0c6f7a0b5ed8dp-20", "0x1.dd5a6e3869738p-20", "0x1.a86f0875fcf94p-19",
    "0x1.7961810874a9fp-18", "0x1.4f8b588e368f1p-17", "0x1.2a5884e341e83p-16",
    "0x1.09456549be1bdp-15", "0x1.d7b9e14a91d46p-15", "0x1.a36e2eb1c432dp-14",
    "0x1.74eea61c12623p-13", "0x1.4b96be9c2da2cp-12", "0x1.26d42cce9b24cp-11",
    "0x1.0624dd2f1a9fcp-10",
)
HEIGHTS = np.array([float.fromhex(h) for h in HEIGHTS_HEX])
FROZEN_INPUTS_SHA256 = "7158ffc646f0a3b7c7cca32f01a88cb06729fe1890bf6a5993c7a6a42d49615c"  # sha256(pack('<26d', *HEIGHTS, P_BASE))
BAND = 1.0 / 107.0                # the corollary's bound on |ratio - 1|


def _ulp_distance(a: float, b: float) -> int:
    """Number of binary64 steps between two positive finite floats."""
    ia = struct.unpack("<q", struct.pack("<d", a))[0]
    ib = struct.unpack("<q", struct.pack("<d", b))[0]
    return abs(ia - ib)


def _check_frozen_inputs() -> None:
    """Refuse edited input literals; confirm their provenance to within a few ulp (run at import).

    Identity is enforced by the SHA-256; the ulp bounds only confirm that the literals are the intended
    formula values and must accept every platform whose libm/SIMD rounding differs in the last bits."""
    packed = struct.pack("<26d", *(float(h) for h in HEIGHTS), P_BASE)
    if hashlib.sha256(packed).hexdigest() != FROZEN_INPUTS_SHA256:
        raise RuntimeError("frozen HEIGHTS/P_BASE literals do not match FROZEN_INPUTS_SHA256 (edited literal?)")
    ref = 10.0 ** np.linspace(-9.0, -3.0, 25)
    for k, (h, r) in enumerate(zip(HEIGHTS, ref)):
        if _ulp_distance(float(h), float(r)) > 4:
            raise RuntimeError(f"frozen height {k} is more than 4 ulp from 10**(-9 + {k}/4)")
    if _ulp_distance(P_BASE, math.sin(0.5 * math.pi * (math.sqrt(2.0) - 1.0)) ** 2) > 8:
        raise RuntimeError("frozen P_BASE is more than 8 ulp from sin^2(pi (sqrt2 - 1)/2)")


_check_frozen_inputs()


# ----------------------------------------------------------------------------
# The axis orbit, its cocycle and the slope series (eq:axis-cocycle, eq:slope-series)
# ----------------------------------------------------------------------------
def axis_orbit(p: float, n: int) -> np.ndarray:
    """x[0] = p, x[k+1] = 4 x[k] (1 - x[k]); length n + 1."""
    x = np.empty(n + 1)
    x[0] = p
    for k in range(n):
        x[k + 1] = MU0 * x[k] * (1.0 - x[k])
    return x


def axis_cocycle(x: np.ndarray, nu1: float):
    """a[n] = 4(1 - 2x[n]), b[n] = -x[n+1], c[n] = nu1 x[n] for n < len(x) - 1."""
    a = MU0 * (1.0 - 2.0 * x[:-1])
    b = -x[1:]
    c = nu1 * x[:-1]
    return a, b, c


def slope_series(p: float, nu1: float, n: int, tail: int = 80) -> np.ndarray:
    """s[m] = - sum_{j<tail} b[m+j] V_{m,j} / U_{m,j+1}, m < n (eq:slope-series, tail terms
    decay like (nu1/4)^j / 2^(j+1))."""
    x = axis_orbit(p, n + tail + 2)
    a, b, c = axis_cocycle(x, nu1)
    s = np.empty(n)
    for m in range(n):
        acc = 0.0
        U = 1.0   # U_{m,0}
        V = 1.0   # V_{m,0}
        for j in range(tail):
            U_next = U * a[m + j]          # U_{m,j+1}
            acc += b[m + j] * V / U_next
            V = V * c[m + j]               # V_{m,j+1}
            U = U_next
        s[m] = -acc
    return s


def _theta0_scaled(K: int) -> int:
    """floor((sqrt(2) - 1) * 2**K), exactly (integer square root)."""
    return math.isqrt(2 << (2 * K)) - (1 << K)


def exact_axis_orbit(n: int, extra_bits: int = EXTRA_BITS) -> np.ndarray:
    """x_0..x_n of the EXACT logistic orbit of p* = sin^2(pi (sqrt2 - 1)/2) under f(x) = 4x(1-x).

    f(sin^2(pi theta/2)) = sin^2(pi T(theta)/2) with the tent map T(theta) = 1 - |1 - 2 theta|, so
    theta_j = T^j(theta0) is computed on the integer A = floor(theta * 2^K), K = n + extra_bits:
    A -> 2A if A <= 2^(K-1), else 2^(K+1) - 2A, which is exact; after j steps the error in theta_j
    is below 2^(j-K) <= 2^-extra_bits. x_j = sin^2(pi theta_j / 2) is then evaluated in float64."""
    K = n + extra_bits
    A = _theta0_scaled(K)
    half, full2 = 1 << (K - 1), 1 << (K + 1)
    shift = K - 60
    scale = float(1 << 60)
    out = np.empty(n + 1)
    for j in range(n + 1):
        theta = (A >> shift) / scale
        out[j] = math.sin(0.5 * math.pi * theta) ** 2
        A = (A << 1) if A <= half else (full2 - (A << 1))
    return out


def reference_running_rate(nu1: float, N: int) -> np.ndarray:
    """(1/n) ln V_n(p*) = (1/n) sum_{j<n} ln(nu1 x_j) for n = 1..N, along the exact orbit."""
    x = exact_axis_orbit(N - 1)
    return np.cumsum(np.log(nu1 * x)) / np.arange(1, N + 1)


# ----------------------------------------------------------------------------
# The map in adapted coordinates and the Lyapunov--Perron operator
# ----------------------------------------------------------------------------
def _phi(x: float, y: float, nu1: float):
    return MU0 * x * (1.0 - x) * (1.0 - y), nu1 * x * y * (1.0 - y)


def lp_fixed_point(p: float, nu1: float, t: float, N: int = N_HORIZON, tol: float = 1e-14,
                   max_sweeps: int = 200, floor: float = 1e-15, trace: list | None = None) -> dict:
    """Picard iteration of the Lyapunov--Perron operator (eq:LP-stable, eq:LP-unstable) on the
    horizon n < N for the stable parameter t (the predator height of the leaf point).

    Adapted coordinates along the reference orbit p_n = (x_n, 0): w = P_n z, P_n = [[1, s_n],[0, 1]],
    z = (xi, zeta); the nonlinear remainder H_n(z) = P_{n+1}^{-1}(Phi(p_n + P_n z) - p_{n+1}) -
    diag(a_n, c_n) z is evaluated through the map itself. The operator
        zeta_n <- V_{0,n} t + sum_{j<n} (V_{0,n}/V_{0,j+1}) H^s_j(z_j),
        xi_n   <- - sum_{n<=j<N} (U_{0,n}/U_{0,j+1}) H^u_j(z_j)
    is iterated from z = 0 until the sup-norm change is <= max(tol * t, floor); the floor (1e-15)
    is the rounding noise of the remainder H, which is a difference of O(1) map values and
    cannot be resolved below about 1e-16. Returns a dict with the fixed point (xi, zeta), the leaf
    point p + s_0 t + xi_0, the sweep count, the last change, the largest ratio of successive
    changes measured above the noise floor ("contraction": an observed ratio, not a certificate
    that the operator contracts on a neighbourhood) and the flag "converged"; if
    `trace` is a list, every sweep's change is appended to it.
    Raises RuntimeError if the iteration does not converge (including any non-finite value):
    such a height is refused, not plotted."""
    x = axis_orbit(p, N + 1)
    a, b, c = axis_cocycle(x, nu1)
    s = slope_series(p, nu1, N + 1)
    xi = np.zeros(N)
    zeta = np.zeros(N)
    changes = []
    converged = False
    for sweep in range(1, max_sweeps + 1):
      with np.errstate(over="ignore", invalid="ignore"):   # non-finite values are detected explicitly below
          # H_j(z_j) for every j < N through the map in adapted coordinates
          Hu = np.empty(N)
          Hs = np.empty(N)
          for j in range(N):
              X = x[j] + xi[j] + s[j] * zeta[j]
              Y = zeta[j]
              Xn, Yn = _phi(X, Y, nu1)
              wx = Xn - x[j + 1]
              wy = Yn
              zu_full = wx - s[j + 1] * wy      # P_{j+1}^{-1} w
              zs_full = wy
              Hu[j] = zu_full - a[j] * xi[j]
              Hs[j] = zs_full - c[j] * zeta[j]
          if not (np.all(np.isfinite(Hu)) and np.all(np.isfinite(Hs))):
              raise RuntimeError("LP iteration did not converge (non-finite remainder)")
          zeta_new = np.empty(N)
          xi_new = np.empty(N)
          # contracting coordinate: forward variation of constants
          for n in range(N):
              V0n = 1.0
              for i in range(n):
                  V0n *= c[i]
              acc = V0n * t
              ratio = 1.0                        # V_{0,n}/V_{0,j+1} for j = n-1 is the empty product
              for j in range(n - 1, -1, -1):
                  acc += ratio * Hs[j]
                  ratio *= c[j]                  # now V_{0,n}/V_{0,j}
              zeta_new[n] = acc
          # expanding coordinate: the unique future-bounded solution, tail truncated at N
          for n in range(N):
              acc = 0.0
              denom = 1.0
              for j in range(n, N):
                  denom *= a[j]                  # prod_{i=n..j} a_i = U_{0,j+1}/U_{0,n}
                  acc += Hu[j] / denom
              xi_new[n] = -acc
          if not (np.all(np.isfinite(zeta_new)) and np.all(np.isfinite(xi_new))):
              raise RuntimeError("LP iteration did not converge (non-finite iterate)")
          change = max(float(np.max(np.abs(zeta_new - zeta))), float(np.max(np.abs(xi_new - xi))))
          changes.append(change)
          if trace is not None:
              trace.append(change)
          xi, zeta = xi_new, zeta_new
          if change <= max(tol * t, floor):
              converged = True
              break
    if not converged:
        raise RuntimeError("LP iteration did not converge")
    ratios = [changes[k] / changes[k - 1] for k in range(1, len(changes)) if changes[k - 1] > 100.0 * floor]
    return {
        "xi": xi, "zeta": zeta, "s0": float(s[0]), "sweeps": sweep, "last_change": float(changes[-1]),
        "contraction": float(max(ratios)) if ratios else 0.0,
        "leaf_x0": float(p + s[0] * t + xi[0]), "converged": True,
    }


# ----------------------------------------------------------------------------
# The matched reference product (cor:rare-predator)
# ----------------------------------------------------------------------------
def matched_product(p: float, nu1: float, t: float, n_max: int = N_MATCH) -> dict:
    """Iterate the full map from the leaf point of height t (float64, n_max steps, within the
    horizon on which float64 rounding, amplified about twofold per step, stays small) and compare
    the predator with the matched reference product t V_n(p), V_n(p) = prod_{j<n} nu1 f^j(p):
      ratio_n = y_n / (t V_n(p)) (n <= n_max), max |ratio_n - 1|,
      rate_actual = (ln y_{n_max} - ln t)/n_max, rate_reference = (1/n_max) sum_{j<n_max} ln(nu1 f^j p),
      running_rate_actual = (1/n) ln(y_n/t), n = 1..n_max,
      itinerary_agreement = #{n <= n_max : sign(x_n - 1/2) = sign(f^n p - 1/2)}.
    The finite-window rates are finite Birkhoff averages, compared with the almost-everywhere limit
    ln(nu1/4) (the genericity of p is not established); the long-run rate is computed along the
    exact orbit by reference_running_rate."""
    x0 = lp_fixed_point(p, nu1, t)["leaf_x0"]
    xref = axis_orbit(p, n_max)
    xs = np.empty(n_max + 1)
    ys = np.empty(n_max + 1)
    xs[0], ys[0] = x0, t
    for n in range(n_max):
        xs[n + 1], ys[n + 1] = _phi(xs[n], ys[n], nu1)
    V = np.empty(n_max + 1)
    V[0] = 1.0
    for n in range(n_max):
        V[n + 1] = V[n] * nu1 * xref[n]
    ratio = ys / (t * V)
    agree = int(np.sum(np.sign(xs - 0.5) == np.sign(xref - 0.5)))
    rate_actual = (math.log(ys[n_max]) - math.log(t)) / n_max
    rate_reference = float(np.mean(np.log(nu1 * xref[:n_max])))
    running_actual = [(math.log(ys[n]) - math.log(t)) / n for n in range(1, n_max + 1)]
    return {
        "t": t, "n_max": n_max, "leaf_x0": x0, "ratio": ratio.tolist(),
        "max_abs_ratio_minus_1": float(np.max(np.abs(ratio - 1.0))),
        "rate_actual": float(rate_actual), "rate_reference": rate_reference,
        "running_rate_actual": running_actual,
        "closed_form_rate": math.log(nu1 / 4.0), "itinerary_agreement": agree,
        "y": ys.tolist(), "reference_product": (t * V).tolist(),
    }


# ----------------------------------------------------------------------------
# Figure and values
# ----------------------------------------------------------------------------
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _manifest_append(path: str) -> None:
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {VERSION}.py | {_sha256(path)} |\n")


def compute_all() -> dict:
    out = {}
    for nu1 in NU1_LIST:
        heights, leaf_x0, sweeps, contr = [], [], [], []
        t_fail = None
        for t in HEIGHTS:
            try:
                r = lp_fixed_point(P_BASE, nu1, float(t))
            except RuntimeError:
                t_fail = float(t)
                break
            heights.append(float(t)); leaf_x0.append(r["leaf_x0"]); sweeps.append(r["sweeps"]); contr.append(r["contraction"])
        s0 = float(slope_series(P_BASE, nu1, 1)[0])
        mp = matched_product(P_BASE, nu1, T_MATCH)
        psi = [lx - P_BASE - s0 * hh for lx, hh in zip(leaf_x0, heights)]
        half = max(1, len(heights) // 2)
        psi2 = float(np.median([ps / hh ** 2 for ps, hh in zip(psi[:half], heights[:half])]))
        rr = reference_running_rate(nu1, N_LONG)
        idx = np.unique(np.round(np.logspace(0, math.log10(N_LONG), 400)).astype(int)) - 1
        out[f"nu1_{nu1:g}"] = {
            "s0": s0, "heights": heights, "leaf_x0": leaf_x0, "sweeps_max": int(max(sweeps)),
            "contraction_max": float(max(contr)), "t_max_converged": float(max(heights)),
            "first_refused_height": t_fail, "matched": mp,
            "psi": psi, "psi_over_t2_small_t": psi2,
            "reference_rate_exact": {"n": (idx + 1).tolist(), "rate": rr[idx].tolist()},
            "reference_rate_at_40": float(rr[N_MATCH - 1]),
            "reference_rate_N_long": float(rr[-1]),
            "closed_form_rate": math.log(nu1 / 4.0),
        }
    return out


# Figure style: both figures are drawn on a 6.5-inch-wide canvas; the manuscript's text width is 6.27 in, and it
# includes Figure 5 at 0.86 of that and Figure 12 at its own width, so the point sizes below print slightly smaller
# (the v18 figure was a 9-inch canvas shown at 5.2 in, which printed its legends at about 4 pt).
_STYLE = {"font.size": 8.5, "axes.titlesize": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
          "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "lines.markersize": 3.0}


def _tex_number(v: float) -> str:
    """Three significant digits in mathtext, with a power of ten written out (-3.03\times10^{4}, not -3.03e+04)."""
    if v == 0 or 1e-3 <= abs(v) < 1e4:
        return f"{v:.3g}"
    e = int(math.floor(math.log10(abs(v)))); m = v / 10 ** e
    return rf"{m:.2f}\times10^{{{e}}}"


def make_figure_leaf(res: dict, out_png: str) -> None:
    """Figure 5: (a) the computed leaf's deviation from its tangent, (b) tracking by the matched reference product;
    columns nu1 = 2.1 and 3."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with matplotlib.rc_context(_STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.75))
        for col, nu1 in enumerate(NU1_LIST):
            r = res[f"nu1_{nu1:g}"]
            h = np.array(r["heights"])
            # (a) the leaf: deviation psi(t) = X(t) - p - s0 t from the slope-series tangent
            ax = axes[0, col]
            psi = np.abs(np.array(r["psi"]))
            ax.loglog(h, psi, "o", color="tab:blue", label=r"$|\psi(t)|$ (computed)")
            ax.loglog(h, abs(r["psi_over_t2_small_t"]) * h ** 2, "--", color="tab:red", lw=1.1,
                      label=rf"$|\psi_2|\,t^2$, $\psi_2={_tex_number(r['psi_over_t2_small_t'])}$")
            ax.set_xlabel(r"height $t=y_0$")
            if col == 0:
                ax.set_ylabel(r"$|X(t)-p-s_0t|$")
            ax.set_title(rf"(a) leaf through $(p,0)$, $\nu_1={nu1:g}$")
            ax.legend(loc="upper left", handlelength=1.6)
            # (b) tracking by the matched reference product (40 steps)
            mp = r["matched"]
            n = np.arange(mp["n_max"] + 1)
            ax = axes[1, col]
            l1, = ax.plot(n, np.log(mp["y"]), "o", color="tab:blue", label=r"$\ln y_n$ (leaf point, $t=10^{-6}$)")
            l2, = ax.plot(n, np.log(mp["reference_product"]), "-", color="tab:orange", lw=1.2,
                          label=r"$\ln(t\,V_n(p))$ (matched reference)")
            ax.set_xlabel("$n$")
            if col == 0:
                ax.set_ylabel(r"$\ln y_n$")
            ax.set_title(rf"(b) tracking, $\nu_1={nu1:g}$")
            ax2 = ax.twinx()
            l3, = ax2.plot(n, mp["ratio"], "-", color="tab:green", lw=1.0, label=r"ratio $y_n/(t\,V_n(p))$")
            band = ax2.axhspan(1.0 - BAND, 1.0 + BAND, color="tab:green", alpha=0.12, lw=0,
                               label=r"$|\mathrm{ratio}-1|<1/107$")
            ax2.set_ylim(1.0 - 2.0 * BAND, 1.0 + 2.0 * BAND)
            if col == 1:
                ax2.set_ylabel("ratio", color="tab:green")
            if col == 0:
                handles = [l1, l2, l3, band]
        fig.legend(handles, [hd.get_label() for hd in handles], loc="lower center", ncol=2,
                   bbox_to_anchor=(0.5, 0.0), handlelength=1.8, frameon=False)
        fig.tight_layout(pad=0.4, h_pad=0.8, w_pad=0.8, rect=(0, 0.075, 1, 1))
        fig.savefig(out_png, dpi=300)
        plt.close(fig)


def make_figure_rates(res: dict, out_png: str) -> None:
    """Figure 12: the reference rate N^-1 ln V_N(p) along the exact orbit of p (up to 1e5 steps), the actual predator
    rate n^-1 ln(y_n/t) for n <= 40, and the almost-everywhere value ln(nu1/4); columns nu1 = 2.1 and 3."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with matplotlib.rc_context(_STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.0))
        for col, nu1 in enumerate(NU1_LIST):
            r = res[f"nu1_{nu1:g}"]; mp = r["matched"]; rr = r["reference_rate_exact"]
            ax = axes[col]
            ax.semilogx(rr["n"], rr["rate"], "-", color="tab:orange", lw=1.2,
                        label=r"reference rate $\frac{1}{N}\ln V_N(p)$, exact orbit")
            ax.semilogx(np.arange(1, mp["n_max"] + 1), mp["running_rate_actual"], "o", ms=2.5, color="tab:blue",
                        label=r"predator rate $\frac{1}{n}\ln(y_n/t)$, $n\leq40$")
            ax.axhline(r["closed_form_rate"], color="k", ls="--", lw=1.0, label=r"$\ln(\nu_1/4)$ (almost-everywhere value)")
            ax.set_xlabel("$N$, $n$ (steps)")
            if col == 0:
                ax.set_ylabel("running decay rate")
            ax.set_title(rf"$\nu_1={nu1:g}$")
            lo = min(min(rr["rate"][5:]), r["closed_form_rate"]) - 0.15
            hi = max(max(rr["rate"][5:]), r["closed_form_rate"]) + 0.15
            ax.set_ylim(lo, hi)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), handlelength=1.8,
                   frameon=False)
        fig.tight_layout(pad=0.4, w_pad=0.8, rect=(0, 0.11, 1, 1))
        fig.savefig(out_png, dpi=300)
        plt.close(fig)


def values_dict(res: dict) -> dict:
    return {
        "provenance": {
            "script": f"{VERSION}.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "status": ("illustration of the objects the proof of thm:basin-positive constructs; "
                       "not an instance of the theorem (base point and heights uncertified); used in no proof"),
        },
        "parameters": {"mu0": MU0, "p": P_BASE, "theta0": THETA0_LABEL, "nu1": list(NU1_LIST), "N_horizon": N_HORIZON,
                       "N_match": N_MATCH, "N_long": N_LONG, "extra_bits": EXTRA_BITS, "t_match": T_MATCH,
                       "band_1_over_107": BAND, "heights": [float(t) for t in HEIGHTS]},
        **res,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-figs", action="store_true", help="write only the values JSON")
    args = ap.parse_args(argv)
    os.makedirs(DATADIR, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    res = compute_all()
    out_json = os.path.join(DATADIR, f"{VERSION}_values.json")
    with open(out_json, "w") as fh:
        json.dump(values_dict(res), fh, indent=2)
    _manifest_append(out_json)
    print(f"[{VERSION}] values -> {out_json}")
    if not args.no_figs:
        for name, maker in (("designG_stable_leaf.png", make_figure_leaf),
                            ("designG_stable_leaf_rates.png", make_figure_rates)):
            out_png = os.path.join(FIGDIR, name)
            maker(res, out_png)
            _manifest_append(out_png)
            print(f"[{VERSION}] figure -> {out_png}")
    for nu1 in NU1_LIST:
        r = res[f"nu1_{nu1:g}"]; mp = r["matched"]
        print(f"  nu1={nu1:g}: p={P_BASE:.6f} s0={r['s0']:+.6f}  t_max_converged={r['t_max_converged']:.1e}  "
              f"max|ratio-1|={mp['max_abs_ratio_minus_1']:.2e}  rate40={r['reference_rate_at_40']:+.4f}  "
              f"rate(N=1e5)={r['reference_rate_N_long']:+.4f} (closed form {r['closed_form_rate']:+.4f})  "
              f"itinerary {mp['itinerary_agreement']}/{mp['n_max'] + 1}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
