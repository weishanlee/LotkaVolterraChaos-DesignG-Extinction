#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G controlled-comparison figure: on-axis invasion exponents versus nu1 (mu0 = 4).

Plots the CLOSED-FORM invasion exponent lambda(m) = ln nu1 + int ln x dm of a rare
predator along several prey laws m of the full-chaos logistic map f(x) = 4x(1-x)
(Section 4, prop:controls and lem:meanvar), as a function of the conversion rate nu1:

    arcsine (physical) law  mu4 :  ln(nu1/4)                      -- eq:lambda-mu4
    period-3 law (same mean 1/2 and variance 1/8 as mu4) :
                                   ln(nu1/4) + (ln 3)/3           -- eq:period3
    mu0 = 2 physical control delta_{1/2} (same mean, zero variance) :  ln(nu1/2)
    fixed point 3/4 :              ln(3 nu1/4)
    period-2 orbit (5 +- sqrt5)/8 (mean 5/8, variance 5/64) :  ln(nu1 sqrt5 / 4)
    mean--variance upper bound (every law with mean 1/2, variance 1/8) :
                                   ln(nu1 4^(-2/3))               -- eq:meanvar-bound
    period-5 law on sin^2(2^j pi/33), j<5 (the fixed example of prop:controls(iv): negative for every 0 < nu1 <= 4;
    the arbitrarily negative exponents of that proposition come from the growing-period family, not from this one law)

Every curve is an exact formula evaluated in floating point; the figure is an
illustration of proved statements, not evidence for them, and it compares
invasion EXPONENTS (growth rates of the normal cocycle), never outcomes: the
period-3 and period-5 orbits are repelling on the axis and are not physical
attractors.  Vertical markers: the showcase values nu1 = 2.1 and nu1 = 3; the
thresholds 4/3 (all-orbit extinction), 2 (same-mean physical control turns
positive), 4^(2/3) ~ 2.5198 (mean--variance obstruction ends), 4 3^(-1/3) ~ 2.7734
(same-mean, same-variance period-3 control turns positive).

Outputs (portable, mirror the sibling designG_*.py scripts):
    figures/designG_invasion_exponents.png                    (DESIGNG_FIGDIR)
    data/designG/designG_invasion_exponents_v1_values.json    (DESIGNG_DATADIR)
and, if DESIGNG_MANIFEST is set, one producer row per artifact (| artifact | version | sha256 |).

Run (CPU only, a few seconds), in the bundle's codes directory:
    python designG_invasion_exponents_v1.py [--no-figs] [--grid N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
from datetime import datetime, timezone
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))   # .../codes
_OUT = os.path.dirname(_HERE)                         # .../output  (or the bundle root)
DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))
FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))

VERSION = "designG_invasion_exponents_v1"

LN2 = math.log(2.0)
LN3 = math.log(3.0)
NU_ALLORBIT = 4.0 / 3.0            # all-orbit extinction threshold (thm:global)
NU_MEANVAR = 4.0 ** (2.0 / 3.0)     # mean--variance obstruction ends (lem:meanvar)
NU_PERIOD3 = 4.0 * 3.0 ** (-1.0 / 3.0)   # period-3 control turns positive (eq:period3)
NU_SHOWCASE = (2.1, 3.0)


# ----------------------------------------------------------------------------
# Closed-form exponents (prop:controls, lem:meanvar)
# ----------------------------------------------------------------------------
def lam_arcsine(nu):
    return np.log(nu / 4.0)


def lam_period3(nu):
    return np.log(nu / 4.0) + LN3 / 3.0


def lam_mu2_control(nu):
    return np.log(nu / 2.0)


def lam_fixed_point(nu):
    return np.log(3.0 * nu / 4.0)


def lam_period2(nu):
    return np.log(nu * math.sqrt(5.0) / 4.0)


def meanvar_bound(nu):
    return np.log(nu * 4.0 ** (-2.0 / 3.0))


def period_p_points(p: int):
    """The periodic orbit sin^2(2^j pi/(2^p+1)), j < p, of f(x)=4x(1-x) (prop:controls(iv))."""
    theta = math.pi / (2 ** p + 1)
    return [math.sin(2 ** j * theta) ** 2 for j in range(p)]


def lam_period_p(nu, p: int = 5):
    pts = period_p_points(p)
    return np.log(nu) + sum(math.log(x) for x in pts) / p


def period3_points():
    return [math.sin(2 ** k * math.pi / 9) ** 2 for k in range(3)]


# ----------------------------------------------------------------------------
# Exact rational sign checks quoted in the manuscript (prop:controls)
# ----------------------------------------------------------------------------
def exact_sign_checks():
    nu = Fraction(21, 10)
    checks = {
        "fixed_point_multiplier_at_2.1": {"value": str(Fraction(3, 4) * nu), "gt_one": Fraction(3, 4) * nu > 1},
        "period2_squared_multiplier_at_2.1": {"value": str(nu * nu * Fraction(5, 16)), "gt_one": nu * nu * Fraction(5, 16) > 1},
        "period3_cubed_multiplier_at_2.1": {"value": str(3 * (nu / 4) ** 3), "gt_one": 3 * (nu / 4) ** 3 > 1},
        "mu2_control_multiplier_at_2.1": {"value": str(nu / 2), "gt_one": nu / 2 > 1},
        "period3_cubed_multiplier_at_3": {"value": str(3 * Fraction(3, 4) ** 3), "gt_one": 3 * Fraction(3, 4) ** 3 > 1},
        "meanvar_cube_check_at_2.1": {"value": f"{nu**3} < 16", "holds": nu ** 3 < 16},
    }
    return checks


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _manifest_append(path: str) -> None:
    """Append a producer-authored ``| artifact | version | sha256 |`` row to DESIGNG_MANIFEST (if set)."""
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {VERSION}.py | {_sha256(path)} |\n")


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
def make_figure(out_png: str, grid: int = 2000) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nu = np.linspace(0.5, 4.0, grid)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axvspan(NU_PERIOD3, 4.0, color="tab:orange", alpha=0.12, lw=0,
               label=r"period-3 vs arcsine sign contrast: $4\cdot3^{-1/3}<\nu_1<4$")
    ax.plot(nu, lam_arcsine(nu), color="tab:blue", lw=2.2, label=r"arcsine $\mu_4$ (chaotic): $\ln(\nu_1/4)$")
    ax.plot(nu, lam_period3(nu), color="tab:orange", lw=2.0, label=r"period 3, same mean and variance: $\ln(\nu_1/4)+\frac{1}{3}\ln 3$")
    ax.plot(nu, lam_mu2_control(nu), color="tab:green", lw=1.6, ls="--", label=r"$\mu_0=2$ physical control $\delta_{1/2}$: $\ln(\nu_1/2)$")
    ax.plot(nu, lam_fixed_point(nu), color="tab:gray", lw=1.2, ls=":", label=r"fixed point $3/4$: $\ln(3\nu_1/4)$")
    ax.plot(nu, lam_period2(nu), color="tab:purple", lw=1.2, ls=":", label=r"period 2: $\ln(\nu_1\sqrt{5}/4)$")
    ax.plot(nu, meanvar_bound(nu), color="tab:red", lw=1.6, ls="-.", label=r"mean–variance upper bound: $\ln(\nu_1 4^{-2/3})$")
    ax.plot(nu, lam_period_p(nu, 5), color="tab:brown", lw=1.0, ls=(0, (1, 1)), label=r"period 5 on $\sin^2(2^j\pi/33)$ (nonchaotic, very negative)")
    for v, txt in ((NU_ALLORBIT, r"$4/3$"), (2.0, r"$2$"), (NU_MEANVAR, r"$4^{2/3}$"), (NU_PERIOD3, r"$4\cdot3^{-1/3}$")):
        ax.axvline(v, color="k", lw=0.6, ls=":", alpha=0.7)
        ax.text(v, 0.965, txt, ha="center", va="top", fontsize=7.5, transform=ax.get_xaxis_transform(),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    for v in NU_SHOWCASE:
        ax.axvline(v, color="tab:red", lw=1.0, alpha=0.8)
        ax.text(v, 0.80, rf"$\nu_1={v:g}$", ha="center", va="bottom", fontsize=8, color="tab:red",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(0.5, 4.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel(r"predator conversion rate $\nu_1$")
    ax.set_ylabel(r"on-axis invasion exponent $\lambda(m)=\ln\nu_1+\int\ln x\,dm$")
    ax.legend(fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    ax.set_title(r"Invasion exponents along prey laws at $\mu_0=4$ (closed forms)", fontsize=9.5, pad=10)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.34)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Values JSON
# ----------------------------------------------------------------------------
def values_dict() -> dict:
    a = period3_points()
    f3prime = abs(np.prod([4.0 * (1.0 - 2.0 * x) for x in a]))
    rows = {}
    for v in NU_SHOWCASE:
        rows[f"nu1={v:g}"] = {
            "arcsine_ln_nu_over_4": float(lam_arcsine(v)),
            "period3_same_mean_variance": float(lam_period3(v)),
            "mu0_2_physical_control": float(lam_mu2_control(v)),
            "fixed_point_3_over_4": float(lam_fixed_point(v)),
            "period2": float(lam_period2(v)),
            "meanvar_upper_bound": float(meanvar_bound(v)),
            "period5_example": float(lam_period_p(v, 5)),
        }
    return {
        "provenance": {
            "script": f"{VERSION}.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "status": "closed-form evaluation; illustration of prop:controls and lem:meanvar, not evidence for them",
        },
        "thresholds": {
            "all_orbit_extinction_4_over_3": NU_ALLORBIT,
            "same_mean_physical_control_positive_above": 2.0,
            "meanvar_obstruction_ends_4^(2/3)": NU_MEANVAR,
            "period3_contrast_begins_4*3^(-1/3)": NU_PERIOD3,
            "arcsine_exponent_zero": 4.0,
            "3^(4/3)/4_at_nu1_3": 3.0 ** (4.0 / 3.0) / 4.0,
        },
        "period3_orbit": {
            "points": a, "sum": float(sum(a)), "product": float(a[0] * a[1] * a[2]),
            "expected_sum": 1.5, "expected_product": 3.0 / 64.0, "abs_f3_prime": float(f3prime),
            "mean": float(sum(a) / 3.0), "second_moment": float(sum(x * x for x in a) / 3.0), "expected_second_moment": 0.375,
        },
        "exponents": rows,
        "exact_rational_sign_checks": exact_sign_checks(),
        "meanvar_equality_law": {"atoms": [0.25, 1.0], "weights": [2.0 / 3.0, 1.0 / 3.0], "sup_E_ln_x": -2.0 * math.log(4.0) / 3.0},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-figs", action="store_true", help="write only the values JSON")
    ap.add_argument("--grid", type=int, default=2000, help="number of nu1 grid points for the curves")
    args = ap.parse_args(argv)
    os.makedirs(DATADIR, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    out_json = os.path.join(DATADIR, f"{VERSION}_values.json")
    with open(out_json, "w") as fh:
        json.dump(values_dict(), fh, indent=2)
    _manifest_append(out_json)
    print(f"[{VERSION}] values -> {out_json}")
    if not args.no_figs:
        out_png = os.path.join(FIGDIR, "designG_invasion_exponents.png")
        make_figure(out_png, grid=args.grid)
        _manifest_append(out_png)
        print(f"[{VERSION}] figure -> {out_png}")
    v = values_dict()["exponents"]
    for k, r in v.items():
        print(f"  {k}: arcsine {r['arcsine_ln_nu_over_4']:+.4f}  period3 {r['period3_same_mean_variance']:+.4f}  "
              f"mu0=2 {r['mu0_2_physical_control']:+.4f}  bound {r['meanvar_upper_bound']:+.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
