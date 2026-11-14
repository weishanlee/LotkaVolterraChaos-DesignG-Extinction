"""Tests for the stable-leaf illustration producer designG_stable_leaf_v1.py (added 2026-10-21).

The producer computes, for the Design-G map at mu0 = 4, the objects the proof of the
positive-area theorem constructs -- the local stable leaf through an axis point as the fixed
point of the paper's Lyapunov--Perron operator, and the matched reference product of the
rare-predator corollary -- and plots them as an ILLUSTRATION (not an instance of the theorem:
the base point and the heights are not certified to lie in the theorem's block). The base point
is p = sin^2(pi theta0 / 2) with theta0 = sqrt(2) - 1, whose binary expansion is known to every
bit, so its EXACT logistic orbit is computed through the tent-map conjugacy.

Run from the folder that holds codes/: set -a; source codes/.designg_bindings.sh; set +a;
     python -m pytest codes/test_designG_stable_leaf.py -q
"""
import math
import pytest
import designG_stable_leaf_v1 as leaf


def test_lp_iteration_converges_with_observed_ratios_below_one():
    # "contraction" is the largest OBSERVED ratio of successive sup-norm changes, not a certificate
    # that the operator contracts on a neighbourhood.
    r = leaf.lp_fixed_point(leaf.P_BASE, 2.1, 1e-6)
    assert r["converged"] and r["sweeps"] <= 50 and r["contraction"] < 1.0


def test_leaf_slope_matches_the_slope_series():
    # Tangency of the computed leaf to (s0, 1) at t = 0. The secant slope q(t) = (X(t) - p)/t carries
    # the leaf's curvature psi(t) = O(t^2) as a term linear in t (here about -424 t: the orbit passes
    # near the critical point at step 6) and rounding noise of order 1e-15/t; Richardson extrapolation
    # 2 q(t) - q(2t) removes the linear term, leaving O(t^2) plus noise.
    p = leaf.P_BASE
    s0 = leaf.slope_series(p, 2.1, 1)[0]
    q = lambda t: (leaf.lp_fixed_point(p, 2.1, t)["leaf_x0"] - p) / t
    assert abs(2.0 * q(1e-8) - q(2e-8) - s0) < 1e-6


def test_matched_product_stays_within_the_corollary_band():
    for nu1 in (2.1, 3.0):
        assert leaf.matched_product(leaf.P_BASE, nu1, 1e-6)["max_abs_ratio_minus_1"] < 1 / 107
    # The caption of Figure 5 also quotes, at nu1 = 3, the maximum over the first 30 steps (5.0e-4): up to
    # there the ratio does not depend on the last digits of the starting point, while the overall maximum
    # (1.3e-3) is reached at the last step, where it does.
    dev = [abs(r - 1.0) for r in leaf.matched_product(leaf.P_BASE, 3.0, 1e-6)["ratio"]]
    assert 4.95e-4 <= max(dev[:31]) < 5.05e-4 and dev.index(max(dev)) == 40


def test_lp_iteration_refuses_a_height_where_it_does_not_converge():
    with pytest.raises(RuntimeError):
        leaf.lp_fixed_point(leaf.P_BASE, 2.1, 0.3, max_sweeps=20)


def test_exact_axis_orbit_agrees_with_the_double_precision_orbit():
    # Two independent computations of the same orbit: the tent-map conjugacy applied to the exact
    # binary expansion of theta0 = sqrt(2) - 1, and the logistic map iterated in float64 from
    # P_BASE. They agree while the float64 rounding, amplified about twofold per step, stays small.
    xs_exact = leaf.exact_axis_orbit(20)
    xs_float = leaf.axis_orbit(leaf.P_BASE, 20)
    assert len(xs_exact) == len(xs_float) == 21
    assert max(abs(a - b) for a, b in zip(xs_exact, xs_float)) < 1e-8


def test_reference_running_rate_at_N_long_is_within_0p05_of_the_closed_form():
    # ln(nu1/4) is the limit of the Birkhoff average of ln(nu1 x) for almost every point under the
    # arcsine law. Along the exact orbit of the chosen base point the running rate at N = 1e5 is
    # compared with it at the loose tolerance 0.05. This is a finite numerical check that pins the
    # computed value against regressions; it is not an error bound (no convergence rate is
    # established for this deterministic orbit) and not a certificate that this base point is
    # generic.
    for nu1 in (2.1, 3.0):
        rate = leaf.reference_running_rate(nu1, leaf.N_LONG)[-1]
        assert abs(rate - math.log(nu1 / 4)) < 0.05


# ---------------------------------------------------------------------------------------------
# 2026-10-21: the height grid and the base point are FROZEN binary64
# inputs. np.logspace (NumPy's SIMD power) gave a different last bit at height index 17 on an
# AVX-512 machine, and the observed contraction ratio amplified it beyond the reproducer's
# tolerance. These tests run a fresh copy of the producer under a simulated different CPU path.
# ---------------------------------------------------------------------------------------------
import types
import numpy as np

CANONICAL_HEX = (
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


def _fresh_leaf_module(source_edit=None):
    """Execute a fresh copy of the producer module (optionally with one exact source substitution)."""
    src = open(leaf.__file__, encoding="utf-8").read()
    if source_edit is not None:
        old, new = source_edit
        assert src.count(old) == 1, f"expected exactly one occurrence of {old!r} in the producer"
        src = src.replace(old, new)
    mod = types.ModuleType("designG_stable_leaf_v1_fresh")
    mod.__file__ = leaf.__file__
    exec(compile(src, leaf.__file__, "exec"), mod.__dict__)
    return mod


def test_stable_leaf_inputs_do_not_depend_on_cpu_dispatched_functions(monkeypatch):
    orig_logspace = np.logspace
    monkeypatch.setattr(np, "logspace", lambda *a, **k: np.nextafter(orig_logspace(*a, **k), np.inf))
    orig_sin = math.sin
    monkeypatch.setattr(math, "sin", lambda v: math.nextafter(orig_sin(v), math.inf))
    m = _fresh_leaf_module()
    assert [float(h).hex() for h in m.HEIGHTS] == list(CANONICAL_HEX)
    assert float(m.P_BASE).hex() == "0x1.77ad613b17f9cp-2"


def test_stable_leaf_frozen_inputs_reject_an_edited_literal():
    with pytest.raises(RuntimeError):
        _fresh_leaf_module(source_edit=('"0x1.2a5884e341e83p-16"', '"0x1.2a5884e341e82p-16"'))
