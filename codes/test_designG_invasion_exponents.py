# -*- coding: utf-8 -*-
"""Exact checks for designG_invasion_exponents_v1.py (prop:controls, lem:meanvar).

Run from the bundle root or output/codes:
    python -m pytest test_designG_invasion_exponents.py -q
"""
import json
import math
import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import designG_invasion_exponents_v1 as m  # noqa: E402


def test_exact_rational_multipliers_at_showcase_values():
    nu = Fraction(21, 10)
    assert Fraction(3, 4) * nu == Fraction(63, 40) and Fraction(63, 40) > 1          # fixed point 3/4
    assert nu * nu * Fraction(5, 16) == Fraction(2205, 1600) and Fraction(2205, 1600) > 1  # period 2 (squared)
    assert 3 * (nu / 4) ** 3 == Fraction(27783, 64000) and Fraction(27783, 64000) < 1    # period 3 (cubed): no contrast at 2.1
    assert nu / 2 == Fraction(21, 20) and Fraction(21, 20) > 1                        # mu0=2 control
    assert 3 * Fraction(3, 4) ** 3 == Fraction(81, 64) and Fraction(81, 64) > 1        # period 3 (cubed) at nu1=3
    assert nu ** 3 == Fraction(9261, 1000) and nu ** 3 < 16                           # (2.1)^3 < 16 = (4^(2/3))^3
    checks = m.exact_sign_checks()
    assert all(c.get("gt_one", c.get("holds")) is True or c.get("gt_one") is False for c in checks.values())
    assert checks["period3_cubed_multiplier_at_2.1"]["gt_one"] is False
    assert checks["period3_cubed_multiplier_at_3"]["gt_one"] is True


def test_period3_orbit_moments_and_derivative():
    a = m.period3_points()
    f = lambda x: 4 * x * (1 - x)
    assert abs(f(a[0]) - a[1]) < 1e-15 and abs(f(a[1]) - a[2]) < 1e-15 and abs(f(a[2]) - a[0]) < 1e-15
    assert abs(sum(a) - 1.5) < 1e-14                     # mean 1/2
    assert abs(sum(x * x for x in a) / 3 - 0.375) < 1e-14  # second moment 3/8 -> variance 1/8
    assert abs(a[0] * a[1] * a[2] - 3 / 64) < 1e-15
    assert abs(abs(math.prod(4 * (1 - 2 * x) for x in a)) - 8.0) < 1e-12   # |(f^3)'| = 8: repelling


def test_thresholds_and_closed_forms():
    assert abs(m.NU_PERIOD3 - 2.773445097402539) < 1e-12
    assert abs(m.NU_MEANVAR - 2.5198420997897464) < 1e-12
    assert m.lam_period3(m.NU_PERIOD3) == pytest.approx(0.0, abs=1e-12)
    assert m.meanvar_bound(m.NU_MEANVAR) == pytest.approx(0.0, abs=1e-12)
    assert m.lam_arcsine(4.0) == pytest.approx(0.0, abs=1e-15)
    assert m.lam_arcsine(2.1) == pytest.approx(math.log(2.1 / 4), abs=1e-15)
    assert m.lam_period3(3.0) == pytest.approx(math.log(3 ** (4 / 3) / 4), abs=1e-12) and m.lam_period3(3.0) > 0
    assert m.meanvar_bound(2.1) == pytest.approx(-0.18225889601721634, abs=1e-12) and m.meanvar_bound(2.1) < 0
    # period-5 law is negative for every nu1 <= 4 and below the proved bound ln nu + 2 ln pi - 6 ln 2
    for nu in (1.0, 2.1, 3.0, 4.0):
        assert m.lam_period_p(nu, 5) < 0
        assert m.lam_period_p(nu, 5) <= math.log(nu) + 2 * math.log(math.pi) - 6 * math.log(2) + 1e-12
    pts = m.period_p_points(5)
    f = lambda x: 4 * x * (1 - x)
    x = pts[0]
    for _ in range(5):
        x = f(x)
    assert abs(x - pts[0]) < 1e-12 and all(0 < p < 1 and abs(p - 0.5) > 1e-6 for p in pts)


def test_script_outputs(tmp_path, monkeypatch):
    monkeypatch.setenv("DESIGNG_DATADIR", str(tmp_path / "data"))
    monkeypatch.setenv("DESIGNG_FIGDIR", str(tmp_path / "figs"))
    monkeypatch.setenv("DESIGNG_MANIFEST", str(tmp_path / "manifest.md"))
    monkeypatch.setattr(m, "DATADIR", str(tmp_path / "data"))
    monkeypatch.setattr(m, "FIGDIR", str(tmp_path / "figs"))
    assert m.main(["--grid", "200"]) == 0
    js = json.load(open(tmp_path / "data" / "designG_invasion_exponents_v1_values.json"))
    r = js["exponents"]["nu1=2.1"]
    assert r["arcsine_ln_nu_over_4"] == pytest.approx(math.log(2.1 / 4))
    assert r["period3_same_mean_variance"] == pytest.approx(math.log(2.1 / 4) + math.log(3) / 3)
    assert r["meanvar_upper_bound"] < 0 and js["exponents"]["nu1=3"]["period3_same_mean_variance"] > 0
    assert (tmp_path / "figs" / "designG_invasion_exponents.png").stat().st_size > 10_000
    rows = (tmp_path / "manifest.md").read_text().strip().splitlines()
    assert len(rows) == 2 and all(r.startswith("| ") and "designG_invasion_exponents_v1.py" in r for r in rows)
