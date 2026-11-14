#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact rational interval certificate of an attracting period-two orbit of the Design-G map
below the fixed-point Jury edge (prop:two-cycle), plus two numerical observations.

The map is Phi(x, y) = (mu0 x (1 - x)(1 - y), nu1 x y (1 - y)) at mu0 = 4, nu1 = 46/25 = 1.84,
a conversion rate strictly between the all-orbit threshold 4/3 and the Jury edge
nu1_J = 8 sqrt(3) - 12 = 1.8564..., where the interior fixed point is NOT a hyperbolic sink.

CERTIFICATE (standard library only; every logical decision uses fractions.Fraction).
  Let B be the closed coordinate box of radius r about the rational centre c, and let
  C be a fixed invertible rational 2x2 matrix (a floating approximation of the inverse of
  D(Phi^2 - I)(c), rationalised; its value is irrelevant to the logic).  Put
      T(z) = z - C (Phi^2(z) - z).
  (i)   Interval evaluation gives an enclosure M of {D(Phi^2)(z) : z in B} as
        DPhi(Phi(B)) . DPhi(B), hence an enclosure E = I - C (M - I) of DT on B.
  (ii)  Each component of T is C^1, so the mean-value theorem gives
        T(B) subset K := T(c) + E . [-r, r]^2 (componentwise, B being convex).
        We check K subset int(B) and q := max row sum of |E| < 1.
  (iii) By (ii) T maps the complete metric space B into itself and is a contraction
        (norm infinity), so Banach's theorem gives a unique fixed point p of T in B;
        C being invertible, Phi^2(p) = p.
  (iv)  The first coordinates of B and Phi(B) are disjoint, so Phi(p) != p: the least
        period of p is two.  B and Phi(B) lie in the open square: the orbit is interior.
  (v)   D(Phi^2)(p) lies in M.  The strict Jury (Schur-Cohn) conditions |d| < 1,
        1 - t + d > 0, 1 + t + d > 0 for the trace t and determinant d, verified on the
        interval enclosure, place both multipliers inside the open unit disc: the orbit is
        a hyperbolic sink with an open basin.
  The certificate proves existence, interior location, least period two and attraction at
  this one parameter; hyperbolicity makes the sink persist on an open parameter
  neighbourhood (implicit function theorem), but no parameter interval is certified here.

OBSERVATIONS (numpy, double precision, NOT certified):
  * the fraction of a 200 x 200 interior grid attracted to the two-cycle at nu1 = 46/25 after
    20,000 iterations (extinct / two-cycle / fixed point / other, by terminal classification);
  * a scan of nu1 over (4/3, nu1_J): the least period (<= 64) of the terminal state of each
    seed of a 20 x 20 grid, reported as the set of periods seen and the interior fraction.

All four strict Jury inequalities are asserted;
the fourth, 1 + d > 0,
is redundant (it is half the sum of 1 - t + d > 0 and 1 + t + d > 0) and is kept only to make
the unit-disc condition |d| < 1 explicit. The other three assertions are sufficient.

Run:
    python certify_period_two_v1.py [--out data/designG] [--no-scan]
Writes <out>/designG_period_two_certificate_v1.json and appends a manifest row when
DESIGNG_MANIFEST is set.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime, timezone
from fractions import Fraction as Q
from typing import Dict, List, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.dirname(_HERE)
MU0 = Q(4)
NU1 = Q(46, 25)
NU1_J = 8.0 * math.sqrt(3.0) - 12.0


# ----------------------------------------------------------------------------------------
# Exact interval arithmetic on Fraction endpoints
# ----------------------------------------------------------------------------------------
class I:
    """Closed interval [lo, hi] with Fraction endpoints; every operation is exact."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @staticmethod
    def of(a):
        return a if isinstance(a, I) else I(a)

    def __add__(self, b):
        b = I.of(b); return I(self.lo + b.lo, self.hi + b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, b):
        return self + (-I.of(b))

    def __rsub__(self, b):
        return I.of(b) + (-self)

    def __mul__(self, b):
        b = I.of(b)
        c = (self.lo * b.lo, self.lo * b.hi, self.hi * b.lo, self.hi * b.hi)
        return I(min(c), max(c))
    __rmul__ = __mul__

    def absmax(self) -> Q:
        return max(abs(self.lo), abs(self.hi))

    def inside_open(self, other: "I") -> bool:
        return other.lo < self.lo and self.hi < other.hi

    def floats(self) -> List[float]:
        return [float(self.lo), float(self.hi)]

    def exact(self) -> List[str]:
        return [str(self.lo), str(self.hi)]


def phi(z, mu0: Q, nu1: Q):
    x, y = z
    return [mu0 * x * (1 - x) * (1 - y), nu1 * x * y * (1 - y)]


def dphi(z, mu0: Q, nu1: Q):
    x, y = z
    return [[mu0 * (1 - 2 * x) * (1 - y), -mu0 * x * (1 - x)],
            [nu1 * y * (1 - y), nu1 * x * (1 - 2 * y)]]


def mm(A, B):
    return [[A[i][0] * B[0][j] + A[i][1] * B[1][j] for j in range(2)] for i in range(2)]


def mv(A, v):
    return [A[i][0] * v[0] + A[i][1] * v[1] for i in range(2)]


# ----------------------------------------------------------------------------------------
# Floating preparation: locate the two-cycle and an approximate inverse (values only)
# ----------------------------------------------------------------------------------------
def _float_two_cycle(mu0: float, nu1: float) -> Tuple[Tuple[float, float], List[List[float]]]:
    """Newton on Phi^2(z) - z = 0 from a point of the numerically observed two-cycle.
    Returns the point and the floating inverse of D(Phi^2 - I) there (used only to build C)."""
    def f(z):
        x, y = z
        return (mu0 * x * (1 - x) * (1 - y), nu1 * x * y * (1 - y))

    def df(z):
        x, y = z
        return [[mu0 * (1 - 2 * x) * (1 - y), -mu0 * x * (1 - x)], [nu1 * y * (1 - y), nu1 * x * (1 - 2 * y)]]

    z = (0.5, 0.1)
    for _ in range(20000):           # settle on the attractor numerically
        z = f(z)
    for _ in range(60):              # Newton for Phi^2(z) = z
        z1 = f(z)
        J = [[a - (1.0 if i == j else 0.0) for j, a in enumerate(row)] for i, row in
             enumerate([[sum(df(z1)[i][k] * df(z)[k][j] for k in range(2)) for j in range(2)] for i in range(2)])]
        g = [f(z1)[0] - z[0], f(z1)[1] - z[1]]
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        dx = (J[1][1] * g[0] - J[0][1] * g[1]) / det
        dy = (-J[1][0] * g[0] + J[0][0] * g[1]) / det
        z = (z[0] - dx, z[1] - dy)
    z1 = f(z)
    M = [[sum(df(z1)[i][k] * df(z)[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    J = [[M[i][j] - (1.0 if i == j else 0.0) for j in range(2)] for i in range(2)]
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    Cinv = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]]
    return z, Cinv


# ----------------------------------------------------------------------------------------
# The certificate
# ----------------------------------------------------------------------------------------
def certify(mu0: Q = MU0, nu1: Q = NU1, centre=None, radius: Q = Q(1, 10 ** 10), C=None) -> Dict[str, object]:
    """Run the exact certificate; returns a record with every check and its outcome.
    ``certified`` is True only if EVERY check passed. No float enters a decision."""
    if centre is None or C is None:
        zf, Cf = _float_two_cycle(float(mu0), float(nu1))
        if centre is None:
            centre = [Q(zf[0]), Q(zf[1])]              # exact binary rationals of the floats
        if C is None:
            C = [[I(Q(v)) for v in row] for row in Cf]
    c = [Q(centre[0]), Q(centre[1])]
    r = Q(radius)
    checks: List[Dict[str, object]] = []

    def chk(name: str, ok: bool, detail: str) -> bool:
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        return bool(ok)

    detC = C[0][0] * C[1][1] - C[0][1] * C[1][0]
    chk("C invertible (exact)", detC.lo > 0 or detC.hi < 0, f"det C in {detC.exact()}")

    B = [I(c[0] - r, c[0] + r), I(c[1] - r, c[1] + r)]
    PB = phi(B, mu0, nu1)                                     # enclosure of Phi(B)
    M = mm(dphi(PB, mu0, nu1), dphi(B, mu0, nu1))             # encloses D(Phi^2) on B
    ident = [[I(1), I(0)], [I(0), I(1)]]
    Dh = [[M[i][j] - ident[i][j] for j in range(2)] for i in range(2)]
    E = [[ident[i][j] - (C[i][0] * Dh[0][j] + C[i][1] * Dh[1][j]) for j in range(2)] for i in range(2)]
    q = max(E[i][0].absmax() + E[i][1].absmax() for i in range(2))
    chk("contraction: max row sum of |E| = |I - C(D(Phi^2) - I)| < 1 on B", q < 1, f"q <= {float(q):.3e}")

    pc = [I(v) for v in c]
    h = [a - b for a, b in zip(phi(phi(pc, mu0, nu1), mu0, nu1), pc)]    # Phi^2(c) - c, exact
    Ch = mv(C, h)
    ER = mv(E, [I(-r, r), I(-r, r)])
    K = [pc[i] - Ch[i] + ER[i] for i in range(2)]
    chk("mean-value enclosure T(B) subset K subset int B", all(K[i].inside_open(B[i]) for i in range(2)),
        f"K = {[k.exact() for k in K]}")

    chk("B and Phi(B) interior to the open unit square",
        all(0 < v.lo and v.hi < 1 for v in B + PB), f"B = {[b.exact() for b in B]}, Phi(B) = {[b.exact() for b in PB]}")
    chk("least period two: first coordinates of Phi(B) and B disjoint",
        PB[0].hi < B[0].lo or B[0].hi < PB[0].lo, f"x(Phi(B)) = {PB[0].exact()}, x(B) = {B[0].exact()}")

    t = M[0][0] + M[1][1]
    d = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    jury = {"1-d": 1 - d, "1+d": 1 + d, "1-t+d": 1 - t + d, "1+t+d": 1 + t + d}
    for name, val in jury.items():
        chk(f"strict Jury condition {name} > 0 on the enclosure of D(Phi^2)(p)", val.lo > 0, f"{name} in {val.exact()}")

    ok = all(ch["ok"] for ch in checks)
    # Floating multipliers of the midpoint matrix (display only; the decisions are above).
    tm, dm = float((t.lo + t.hi) / 2), float((d.lo + d.hi) / 2)
    disc = tm * tm - 4 * dm
    mult = ([(tm - math.sqrt(disc)) / 2, (tm + math.sqrt(disc)) / 2] if disc >= 0
            else [tm / 2, math.sqrt(-disc) / 2])
    return {
        "mu0": str(mu0), "nu1": str(nu1), "nu1_float": float(nu1), "nu1_J": NU1_J,
        "all_orbit_threshold": 4.0 / 3.0,
        "centre": [str(v) for v in c], "centre_float": [float(v) for v in c], "radius": str(r),
        "C": [[str(v.lo) for v in row] for row in C],
        "contraction_bound_q": float(q),
        "box_B": [b.exact() for b in B], "box_Phi_B": [b.exact() for b in PB],
        "point_p_float": [float(v) for v in c], "point_Phi_p_float": [float((v.lo + v.hi) / 2) for v in PB],
        "trace_D_Phi2_enclosure": t.exact(), "det_D_Phi2_enclosure": d.exact(),
        "trace_D_Phi2_float": t.floats(), "det_D_Phi2_float": d.floats(),
        "jury_enclosures": {k: v.exact() for k, v in jury.items()},
        "multipliers_midpoint_float": mult if disc >= 0 else {"real_part": mult[0], "imag_part": mult[1]},
        "multipliers_real": bool(disc >= 0),
        "checks": checks, "certified": ok,
        "statement": ("interior hyperbolic attracting period-two orbit of the Design-G map at mu0 = 4, "
                      "nu1 = 46/25 (certified by exact rational interval arithmetic)" if ok else "NOT certified"),
    }


# ----------------------------------------------------------------------------------------
# Numerical observations (double precision; not certified)
# ----------------------------------------------------------------------------------------
def basin_fraction(nu1: float, n_side: int = 200, n_iter: int = 20000, mu0: float = 4.0) -> Dict[str, object]:
    import numpy as np
    g = (np.arange(n_side) + 0.5) / n_side
    X, Y = np.meshgrid(g, g, indexing="ij")
    x, y = X.ravel().copy(), Y.ravel().copy()
    for _ in range(n_iter):
        xn = mu0 * x * (1 - x) * (1 - y)
        y = nu1 * x * y * (1 - y)
        x = xn
    x1 = mu0 * x * (1 - x) * (1 - y); y1 = nu1 * x * y * (1 - y)
    x2 = mu0 * x1 * (1 - x1) * (1 - y1); y2 = nu1 * x1 * y1 * (1 - y1)
    extinct = y < 1e-12
    fixed = (~extinct) & (np.hypot(x1 - x, y1 - y) < 1e-8)
    two = (~extinct) & (~fixed) & (np.hypot(x2 - x, y2 - y) < 1e-8)
    other = ~(extinct | fixed | two)
    n = x.size
    return {"nu1": nu1, "grid": f"{n_side}x{n_side} midpoint lattice", "n_iter": n_iter,
            "fraction_extinct": float(extinct.sum() / n), "fraction_two_cycle": float(two.sum() / n),
            "fraction_fixed_point": float(fixed.sum() / n), "fraction_other": float(other.sum() / n),
            "classification": "terminal state after n_iter steps: y < 1e-12 extinct; |Phi(z)-z| < 1e-8 fixed; "
                              "|Phi^2(z)-z| < 1e-8 two-cycle; else other (double precision, not certified)"}


def nu1_scan(nu1_lo: float = 1.34, nu1_hi: float = 1.856, step: float = 0.004, n_side: int = 20,
             n_iter: int = 20000, max_period: int = 64, mu0: float = 4.0) -> List[Dict[str, object]]:
    import numpy as np
    g = (np.arange(n_side) + 0.5) / n_side
    X, Y = np.meshgrid(g, g, indexing="ij")
    out: List[Dict[str, object]] = []
    k = 0
    while True:
        nu1 = round(nu1_lo + k * step, 6)
        if nu1 > nu1_hi + 1e-12:
            break
        k += 1
        x, y = X.ravel().copy(), Y.ravel().copy()
        for _ in range(n_iter):
            xn = mu0 * x * (1 - x) * (1 - y)
            y = nu1 * x * y * (1 - y)
            x = xn
        extinct = y < 1e-12
        period = np.zeros(x.size, dtype=np.int64)       # 0 = none detected within max_period
        xs, ys = x.copy(), y.copy()
        for p in range(1, max_period + 1):
            xn = mu0 * xs * (1 - xs) * (1 - ys)
            ys = nu1 * xs * ys * (1 - ys)
            xs = xn
            hit = (period == 0) & (~extinct) & (np.hypot(xs - x, ys - y) < 1e-7)
            period[hit] = p
        interior = ~extinct
        periods = sorted(int(v) for v in set(period[interior].tolist()) if v > 0)
        n_unresolved = int(((period == 0) & interior).sum())
        out.append({"nu1": nu1, "fraction_extinct": float(extinct.sum() / x.size),
                    "fraction_interior": float(interior.sum() / x.size),
                    "periods_observed": periods, "n_interior_no_period_le_64": n_unresolved})
    return out


def _manifest_append(path: str) -> None:
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    import hashlib
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {os.path.basename(__file__)} | {h} |\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exact certificate of the attracting two-cycle at (mu0, nu1) = (4, 46/25).")
    ap.add_argument("--out", default=os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG")))
    ap.add_argument("--no-scan", action="store_true", help="skip the numerical basin fraction and nu1 scan")
    args = ap.parse_args(argv)
    rec = certify()
    print("=" * 88)
    print("Exact rational interval certificate: attracting period-two orbit of Design-G at mu0=4, nu1=46/25")
    print("=" * 88)
    for ch in rec["checks"]:
        print(f"  [{'OK ' if ch['ok'] else 'FAIL'}] {ch['check']}")
    print(f"  p ~ {rec['point_p_float']}, Phi(p) ~ {rec['point_Phi_p_float']}")
    print(f"  trace in {rec['trace_D_Phi2_float']}, det in {rec['det_D_Phi2_float']}, multipliers ~ {rec['multipliers_midpoint_float']}")
    print(f"  CERTIFIED = {rec['certified']}")
    values: Dict[str, object] = {
        "provenance": {"generated_by": "codes/" + os.path.basename(__file__), "python": sys.version.split()[0],
                       "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                       "note": "certificate: exact fractions.Fraction arithmetic; observations: numpy double precision"},
        "certificate": rec,
    }
    if not args.no_scan:
        print("  numerical observations (double precision, not certified) ...")
        values["basin_fraction_nu1_46_25"] = basin_fraction(float(NU1))
        bf = values["basin_fraction_nu1_46_25"]
        print(f"  basin fraction at nu1=46/25 ({bf['grid']}): two-cycle {bf['fraction_two_cycle']:.4f}, "
              f"extinct {bf['fraction_extinct']:.4f}, fixed {bf['fraction_fixed_point']:.4f}, other {bf['fraction_other']:.4f}")
        values["nu1_scan"] = nu1_scan()
        rows = values["nu1_scan"]
        print("  nu1 scan (20x20 grid, 20000 steps): nu1 -> interior fraction, periods observed")
        for row in rows[::10]:
            print(f"    {row['nu1']:.3f}: interior {row['fraction_interior']:.3f}, periods {row['periods_observed']}")
    os.makedirs(args.out, exist_ok=True)
    out_json = os.path.join(args.out, "designG_period_two_certificate_v1.json")
    with open(out_json, "w") as fh:
        json.dump(values, fh, indent=2)
    print("wrote", out_json)
    _manifest_append(out_json)
    return 0 if rec["certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
