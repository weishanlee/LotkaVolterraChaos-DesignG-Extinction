#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pytest harness of the Design-G bundle (dynamics and quasiperiodic producers,
validator, reproduction check).

Run (from the directory that holds ``codes/``; README.md gives the command for all four test files):
    set -a; . codes/.designg_bindings.sh; set +a
    python -m pytest codes/test_designG_v4.py -q

$DYN/$QP/... are BOUND names (see ``codes/.designg_bindings.sh``, the
single source of truth) -- every test below uses
``_load("DYN")``/``_load("QP")``/``DYN_PATH``/``QP_PATH``, never a literal
version string or ``_latest()``-style discovery.

"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from typing import Any, Dict, Iterator, List

_PY = os.environ.get("DESIGNG_PYTHON", sys.executable)  # the interpreter running pytest; override with DESIGNG_PYTHON


def _bound(var):   # read a BOUND path from the exported .designg_bindings.sh (single source of truth)
    p = os.environ.get(var); assert p, f"{var} not exported -- from the tree root run: set -a; source codes/.designg_bindings.sh; set +a (output/codes/.designg_bindings.sh in the authoring tree)"; return p


def _load(var):
    p = _bound(var); spec = importlib.util.spec_from_file_location(os.path.basename(p)[:-3], p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


DYN_PATH = _bound("DYN"); QP_PATH = _bound("QP")


def test_scaffold():
    m = _load("DYN"); assert hasattr(m, "_preflight") and hasattr(m, "resolve_paths")


def test_scaffold_qp():
    """Sibling of test_scaffold for $QP: B0's acceptance bar ("each with _preflight,
    resolve_paths, ... --smoke/--only/default CLI dispatch") explicitly covers both
    scripts, so this checks $QP's scaffold surface the same way test_scaffold checks $DYN.
    """
    m = _load("QP")
    assert hasattr(m, "_preflight") and hasattr(m, "resolve_paths") and hasattr(m, "_smoke_params")


def test_smoke_cli_dispatch_present(tmp_path, monkeypatch):
    """Both scripts declare --smoke/--only/--public on their argparse parser (Task B0
    Step 3 "CLI dispatch contract"), and DYN's --only choices match its block registry.
    The exit status belongs to that contract: a run whose figure step fails does not exit 0
    (2026-11-10; the dynamics and the ergodic producer).
    """
    dyn = _load("DYN")
    qp = _load("QP")
    for m in (dyn, qp):
        assert hasattr(m, "main"), f"{m.__name__} missing main()"
    # DYN: --only's choices must equal the block registry's keys (single source of
    # truth for block names -- see _BLOCK_NAMES's docstring).
    assert hasattr(dyn, "_BLOCK_NAMES") and hasattr(dyn, "_block_registry")
    registry_keys = set(dyn._block_registry(None).keys())  # safe: lambdas are lazy, args unused at build time
    assert registry_keys == set(dyn._BLOCK_NAMES), (registry_keys, dyn._BLOCK_NAMES)
    # --public is accepted and does nothing (both preflights); its help must say so
    for path in (DYN_PATH, QP_PATH):
        r = subprocess.run([_PY, path, "--help"], capture_output=True, text=True)
        helptext = " ".join(r.stdout.split())
        assert r.returncode == 0 and "--public" in helptext and "has no effect" in helptext and "relax" not in helptext, (path, r.stdout[-600:], r.stderr[-300:])
    # The PNGs these two producers draw are manuscript figures (dynamics: Figures 2, 4, 11 and S2;
    # ergodic: Figures 7, 8 and S1) and run_designG_figs.sh relies on the exit status. A figure step that
    # raises was a warning and exit 0; it must fail the run, with the values JSON still written. Nothing is
    # computed here: the block registry and the checks are stubs, only main()'s bookkeeping runs.
    def no_figures(*_a, **_k):
        raise RuntimeError("figure backend unavailable (test)")
    for var in ("DESIGNG_FIGDIR", "DESIGNG_DATADIR"):
        monkeypatch.setenv(var, str(tmp_path))
    monkeypatch.delenv("DESIGNG_MANIFEST", raising=False)
    monkeypatch.setattr(dyn, "_block_registry", lambda args: {"coexistence": ("coexistence", lambda: ([], {}))})
    assert dyn.main(["--smoke", "--only", "coexistence", "--no-figs"]) == 0      # control: nothing drawn, nothing failed
    monkeypatch.setattr(dyn, "make_figures", no_figures)
    assert dyn.main(["--smoke", "--only", "coexistence"]) != 0, "a failed figure step must fail the dynamics producer"
    assert "figures_error" in json.load(open(tmp_path / "designG_dynamics_v4_smoke.json"))
    erg = _load_path("designG_ergodic_quadrature_v2.py")
    monkeypatch.setattr(erg, "DATADIR", str(tmp_path)); monkeypatch.setattr(erg, "FIGDIR", str(tmp_path))
    for name in ("check_arcsine_moments", "check_equilibrium_measure", "check_subaction", "check_I_family"):
        monkeypatch.setattr(erg, name, lambda *_a, **_k: ([], {}))
    assert erg.main(["--max-period", "1", "--no-figs"]) == 0                     # control (the period-one enumeration is real)
    monkeypatch.setattr(erg, "make_figures", no_figures)
    assert erg.main(["--max-period", "1"]) != 0, "a failed figure step must fail the ergodic producer"
    assert "figures_error" in json.load(open(tmp_path / "designG_ergodic_quadrature_v2_values.json"))


# =============================================================================
# _taxonomy JSON-schema-convention coverage (Task B0 gotcha 6): at B0, most of
# the coverage list below is not produced yet (Tasks B1-B5 add those blocks) --
# this test checks coverage of whatever check-bearing scalar IS present, so it
# passes (near-)trivially now and becomes enforcing as B1-B5 land. It also
# doubles as the --smoke end-to-end sanity check for both scripts (subprocess,
# real CLI invocation, real JSON on disk).
# =============================================================================

# Exact wildcard pointers from the plan's "Required _taxonomy coverage" (Task B0
# Step 3).
# Presence is checked per-JSON (DYN's smoke JSON, QP's smoke JSON); "rho" /
# "largest_lyap" are bare top-level field names (no leading "/"), the rest are
# JSON-Pointer paths with "*" wildcarding any non-"_" record key.
EXPECTED_TAXONOMY: Dict[str, str] = {
    "/validation_table/*/fraction_extinct": "regression_anchor",
    "/validation_table/*/is_small": "empirical_discovery",
    "/sensitivity/base/*/fraction": "regression_anchor",
    "/sensitivity/*/*/fraction": "empirical_discovery",
    "/finite_amp/per_nu1/*/per_delta/*/decay_rate": "empirical_discovery",
    "/exact_skew/per_nu1/*/lambda_perp_quadrature": "theory_invariant",
    "/exact_skew/per_nu1/*/collapse_guard_ok": "theory_invariant",
    "/exact_skew/per_nu1/*/lambda_perp_orbit": "empirical_discovery",
    "rho": "regression_anchor",
    "largest_lyap": "empirical_discovery",
    # Retained v2 producers ($DYN still emits these unconditionally; B0 populates
    # their _taxonomy entries directly since no later B-task "adds" them).
    "/coexistence/jury_band/lower_edge_closed_form": "regression_anchor",
    "/arcsine_moments/E_lnx": "regression_anchor",
    "/family_unification/per_nu1/*/pairwise_max_spread": "regression_anchor",
}


def _resolve_pointer(obj: Any, tokens: List[str]) -> Iterator[Any]:
    """Yield every value reachable via `tokens` ('*' = any non-'_' dict key, or any list index)."""
    if not tokens:
        yield obj
        return
    head, rest = tokens[0], tokens[1:]
    if head == "*":
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(k, str) and k.startswith("_"):
                    continue  # underscore-prefixed metadata keys are never wildcard-matched
                yield from _resolve_pointer(v, rest)
        elif isinstance(obj, list):
            for v in obj:
                yield from _resolve_pointer(v, rest)
        return
    if isinstance(obj, dict) and head in obj:
        yield from _resolve_pointer(obj[head], rest)


def _pointer_present(root: Dict[str, object], pointer: str) -> bool:
    if not pointer.startswith("/"):
        return isinstance(root, dict) and pointer in root
    tokens = pointer.strip("/").split("/")
    for _ in _resolve_pointer(root, tokens):
        return True
    return False


def test_taxonomy_coverage(tmp_path):
    """Runs in --smoke (both $DYN and $QP); fails if a PRESENT check-bearing
    scalar lacks a matching (or correctly-kinded) _taxonomy entry. Fields not
    yet produced (most of the U1-U5 list, pre Tasks B1-B5) are skipped, not
    violations: the test passes trivially at B0 and becomes enforcing as B1-B5
    add scalars.
    """
    env = dict(os.environ, DESIGNG_FIGDIR=str(tmp_path), DESIGNG_DATADIR=str(tmp_path))
    violations: List[str] = []
    for var, path in (("DYN", DYN_PATH), ("QP", QP_PATH)):
        # NOT check=True: under --smoke's drastically reduced sample sizes, DYN's
        # own PASS/FAIL summary may legitimately return exit code 1 (some physics
        # checks miss tolerance at tiny N -- see designG_dynamics_v3.py's smoke
        # NOTE) -- that is not a crash. The real acceptance bar is "the smoke JSON
        # got written", asserted explicitly below;
        # a genuine crash (preflight/path/traceback) leaves no file and the
        # assertion below fails with the captured stderr for diagnosis.
        result = subprocess.run([_PY, path, "--smoke"], env=env, capture_output=True, text=True)
        stem = os.path.basename(path)[:-3]
        smoke_file = tmp_path / f"{stem}_smoke.json"
        assert smoke_file.exists(), (
            f"{var}: --smoke did not write {smoke_file.name} (returncode={result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        v = json.load(open(smoke_file))
        assert "_taxonomy" in v, f"{var}: values JSON missing top-level _taxonomy map"
        tax = v["_taxonomy"]
        for ptr, kind in EXPECTED_TAXONOMY.items():
            if not _pointer_present(v, ptr):
                continue  # not produced yet (Tasks B1-B5) -- not a B0 violation
            if ptr not in tax:
                violations.append(f"{var}: present field {ptr!r} has no _taxonomy entry")
            elif tax[ptr] != kind:
                violations.append(f"{var}: _taxonomy[{ptr!r}] = {tax[ptr]!r}, expected {kind!r}")
    assert not violations, "taxonomy coverage violations:\n" + "\n".join(violations)


# =============================================================================
# Task B1 (U1): continuous extinction fraction + a separate `is_small` boolean
# in the validation table.
# `--only extinction --smoke` restricts the run to the block
# that builds `validation_table` from check_extinction_fraction().
# =============================================================================


def test_u1_continuous_fraction_and_flag(tmp_path):   # FAST/smoke: schema only
    import json, subprocess, os
    env = dict(os.environ, DESIGNG_FIGDIR=str(tmp_path), DESIGNG_DATADIR=str(tmp_path))
    subprocess.run([_PY, DYN_PATH, "--only", "extinction", "--smoke"], check=True, env=env)
    stem = os.path.basename(DYN_PATH)[:-3]
    v = json.load(open(tmp_path/f"{stem}_smoke.json"))
    row = {r["label"]: r for r in v["validation_table"]}["ext_frac_nu1_2"]
    assert type(row["fraction_extinct"]) is float            # continuous, not a bool recode to 0
    assert type(row["is_small"]) is bool and "small_ceiling" in row
    assert "fraction_extinct" in row and "is_small" in row   # distinct FIELDS (don't compare values: 0.0==False in Python)
    # the note of an "is small" row states the criterion it measures (threshold, terminal window,
    # horizon, as recorded by the same run) and names no destination of the other orbits
    ef = v["extinction_fraction"]
    small = [r for r in v["validation_table"] if str(r.get("label", "")).startswith("extinction fraction small")]
    assert small, "no 'extinction fraction small' row in the smoke record"
    crit = f"(terminal criterion: y < {ef['eps']:g} over the last {ef['tail']} of {ef['horizon']} iterates; "
    for r in small:
        assert crit in r["note"] and "coexistence" not in r["note"], r["note"]
    assert v["_taxonomy"]["/validation_table/*/fraction_extinct"] == "regression_anchor"


# =============================================================================
# Task B2 (U2): extinction-fraction sensitivity sweep (C0-C10) + confidence
# intervals.
# `--only sensitivity --smoke` restricts the run to run_sensitivity(). The point
# tested is the CI SEMANTICS BY DESIGN: deterministic grid
# fractions are exact-for-that-grid (no CI); only random-IC configs get Wilson.
# =============================================================================


def test_u2_ci_semantics_by_design(tmp_path):
    import json, subprocess, os
    env = dict(os.environ, DESIGNG_FIGDIR=str(tmp_path), DESIGNG_DATADIR=str(tmp_path))
    subprocess.run([_PY, DYN_PATH,
                    "--only","sensitivity","--smoke"], check=True, env=env)  # smoke = reduced
    stem = os.path.basename(DYN_PATH)[:-3]
    s = json.load(open(tmp_path/f"{stem}_smoke.json"))["sensitivity"]
    # random-IC configs get Wilson; deterministic configs get none + convergence band
    assert s["ic.runif"]["2"]["ci_method"] == "wilson"
    assert s["base"]["2"]["ci_method"] == "none" and s["base"]["2"]["ci_low"] is None


# =============================================================================
# Task B3 (U3): finite-amplitude two-way decay diagnostic, NON-blocking.
# `--only
# finite_amp --smoke` restricts the run to run_finite_amplitude(). This block is
# tagged empirical_discovery (an empirical diagnostic, NOT a theory gate) -- the
# crux checked here is BOTH the full delta-sweep schema AND that it introduces
# NO blocking PASS/FAIL Check row (validation_table stays empty under --only
# finite_amp, mirroring run_sensitivity's contribution of zero Check rows).
# =============================================================================


def test_u3_finite_amp_schema_and_nonblocking(tmp_path):
    import json, subprocess, os
    env = dict(os.environ, DESIGNG_FIGDIR=str(tmp_path), DESIGNG_DATADIR=str(tmp_path))
    subprocess.run([_PY, DYN_PATH,
                    "--only", "finite_amp", "--smoke"], check=True, env=env)
    stem = os.path.basename(DYN_PATH)[:-3]
    v = json.load(open(tmp_path/f"{stem}_smoke.json"))
    fa = v["finite_amp"]
    assert fa["interpretation"] == (
        "finite-amplitude; distinct from the axis variational exponent (U4); "
        "expected -> ln(nu1/4) as delta->0"
    )
    expected_deltas = {"0.01", "0.001", "0.0001"}
    seen_deltas = set()
    assert fa["per_nu1"], "finite_amp.per_nu1 must be non-empty"
    for nu1_key, nu1_rec in fa["per_nu1"].items():
        per_delta = nu1_rec["per_delta"]
        seen_deltas |= set(per_delta.keys())
        # FULL delta sweep for EVERY nu1 -- schema, not a cost knob (even under --smoke)
        assert set(per_delta.keys()) == expected_deltas, (nu1_key, set(per_delta.keys()))
        for delta_key, cell in per_delta.items():
            assert type(cell["decay_rate"]) is float
            assert type(cell["sem"]) is float
            assert isinstance(cell["horizon_sweep"], list) and len(cell["horizon_sweep"]) > 0
            for pt in cell["horizon_sweep"]:
                assert type(pt["decay_rate"]) is float
    assert seen_deltas == expected_deltas
    # tagged empirical_discovery (B0's EXPECTED_TAXONOMY / _TAXONOMY convention)
    assert v["_taxonomy"]["/finite_amp/per_nu1/*/per_delta/*/decay_rate"] == "empirical_discovery"
    # NON-blocking: no PASS/FAIL Check row at all under this --only run (the block
    # contributes zero rows to all_rows -- structurally there is no tolerance gate
    # on decay_rate to pass or fail).
    assert v["validation_table"] == []
    assert v["summary"] == {"n_pass": 0, "n_total": 0, "all_pass": True}


# =============================================================================
# Task B4 (U4): exact axis transverse-variational lambda_perp = ln(nu1/4) --
# quadrature primary (BLOCKING) + orbit cross-check.
# `--only exact_skew` restricts
# the run to run_exact_skew(). The double-exponential (tanh-sinh) quadrature
# is the BLOCKING theory invariant (tol 1e-6, reaches ~1e-10); the direct
# float64 orbit Birkhoff average is a SECONDARY, descriptive empirical
# cross-check (tol 2e-2).
# =============================================================================

import math
def test_u4_quadrature_is_exact():
    m = _load("DYN")
    for nu1 in (2.0, 3.0, 4.0):
        assert abs(m.lambda_perp_analytic(nu1) - math.log(nu1/4)) < 1e-15   # closed form, exact
        lam = m.lambda_perp_quadrature(nu1)                                 # DE quadrature
        assert abs(lam - math.log(nu1/4)) <= 1e-6                           # reaches ~1e-10


def test_u4_orbit_cross_check_agrees():
    m = _load("DYN")
    lam_q = m.lambda_perp_quadrature(3.0)
    lam_o, xs = m.lambda_perp_orbit(3.0, T=10000)
    assert m.collapse_guard(xs)                 # orbit did not degenerate
    assert abs(lam_o - lam_q) <= 2e-2           # descriptive agreement


# =============================================================================
# Defensive fix: lambda_perp_orbit's math.log(nu1*x) can
# raise ValueError ("math domain error") if a drawn orbit hits x=0.0 exactly
# (e.g. x0=0.5 -> 1.0 -> 0.0 in 2 logistic steps). Before this fix the ONLY
# crash guard was the _is_exceptional_x0 IC pre-filter (which DOES catch
# x0=0.5 as a low-order dyadic rational, p=1<=20) -- a numerically-degenerating
# x0 that ISN'T literal low-order-dyadic could still reach x=0.0 after enough
# iterations and abort the whole _run_one_exact_skew_nu1 horizon-sweep loop.
# The fix wraps the per-IC lambda_perp_orbit(...) call inside that loop in
# try/except ValueError, excluding the crashing IC (same excluded_E_ic_count
# bucket as _is_exceptional_x0) instead of propagating the crash. This test
# exercises that exact call-site pattern directly (lambda_perp_orbit wrapped
# the same way), decoupled from whether the
# pre-filter would normally intercept x0=0.5 first in the real pipeline.
# =============================================================================


def test_u4_orbit_crash_safety_excludes_not_aborts():
    import pytest
    m = _load("DYN")
    with pytest.raises(ValueError):
        m.lambda_perp_orbit(3.0, T=50, B=5, x0=0.5)   # confirms the crash premise

    vals: List[float] = []
    n_excluded = 0
    for x0 in (0.5, 0.401):                            # one crasher, one well-behaved IC
        try:
            mean_h, xs_h = m.lambda_perp_orbit(3.0, T=50, B=5, x0=x0)
        except ValueError:
            n_excluded += 1
            continue
        if m.collapse_guard(xs_h):
            vals.append(mean_h)
        else:
            n_excluded += 1
    assert n_excluded == 1   # the crasher was excluded, not left to abort the loop
    assert len(vals) == 1    # the well-behaved IC still made it into aggregation


# =============================================================================
# Task B5 (U5): quasiperiodic invariant-circle diagnostics.
# Checks the
# quasiperiodic (U5) quantities: rotation number rho, invariant-circle residual,
# FFT-peak SNR, and a Benettin Lyapunov spectrum on the mu0=3.5, nu1=4.9
# observed attracting-circle regime (nu1=4.9, above the unit-circle locus). `--smoke` KEEPS that regime (see $QP's
# _smoke_params()/runtime assert) -- this is the FAST/smoke schema+taxonomy
# check; the largest_lyap ~ 0 / normal_lyap < 0 signature is confirmed
# separately against the full (non-smoke) run, not gated in this fast test.
# =============================================================================


def test_u5_quasiperiodic_diagnostics(tmp_path):   # FAST/smoke: schema + taxonomy only
    import json, subprocess, os
    env = dict(os.environ, DESIGNG_FIGDIR=str(tmp_path), DESIGNG_DATADIR=str(tmp_path))
    subprocess.run([_PY, QP_PATH, "--smoke"], check=True, env=env)
    stem = os.path.basename(QP_PATH)[:-3]
    v = json.load(open(tmp_path/f"{stem}_smoke.json"))
    for k in ("rho","largest_lyap","normal_lyap","invariant_circle_residual","fft_peak_snr"):
        assert k in v and k in v["_taxonomy"]
    assert v["_taxonomy"]["largest_lyap"] == "empirical_discovery"   # finite-time diagnostic, not exact


# =============================================================================
# Task B6 (U6): machine-checkable validator ($VAL) + claim inventory ($CLAIMS).
# ONE
# comprehensive test: it exercises (a) the --smoke
# validator end-to-end over the REAL $CLAIMS + fresh smoke artifacts (exit 0;
# every per-check kind in the taxonomy set; every numeric anchor SKIPPED_SMOKE;
# non-blocking manuscript check skipped), (b) the one-scalar selector resolution
# incl. the ZERO/MULTIPLE hard-fail (SelectorError), multi-named-segment
# substitution, list field-value + underscore-key skipping, and most-specific
# _taxonomy kind lookup, and (c) the 6-env-var abort. $VAL/$CLAIMS are BOUND
# names (.designg_bindings.sh) -- read via _bound(), never a literal path.
# =============================================================================

VAL_PATH = _bound("VAL"); CLAIMS_PATH = _bound("CLAIMS")


def _b6_make_smoke_dir(tmp_path):
    """Fresh smoke artifacts (DYN+QP --smoke) + a producer-authored manifest.
    NOT check=True: full --smoke (all blocks) may exit 1 at tiny N (see
    test_taxonomy_coverage's NOTE); the real bar is 'the smoke JSON was written'."""
    figs = tmp_path / "figs"; data = tmp_path / "data"; figs.mkdir(); data.mkdir()
    manifest = tmp_path / "manifest.md"
    env = dict(os.environ, DESIGNG_FIGDIR=str(figs), DESIGNG_DATADIR=str(data),
               DESIGNG_MANIFEST=str(manifest))
    for path in (DYN_PATH, QP_PATH):
        subprocess.run([_PY, path, "--smoke"], env=env, capture_output=True, text=True)
        stem = os.path.basename(path)[:-3]
        assert (data / f"{stem}_smoke.json").exists(), f"{stem} --smoke wrote no JSON"
    return figs, data, manifest


def test_u6_validator_smoke_selectors_gates(tmp_path):
    import pytest
    m = _load("VAL")   # importable validator module

    # --- (b) one-scalar selector resolution + zero/multiple HARD-FAIL ---
    root = {"per_nu1": {"2": {"v": 1.0}, "3": {"v": 2.0}, "_meta": {"v": 9.0}}}
    assert m.resolve_one(root, "/per_nu1/*/v", {"nu1": "2"}) == 1.0        # segment-key match -> exactly one
    assert m.resolve_scalars(root, "/per_nu1/*/v", {}) == [1.0, 2.0]       # '*' skips underscore keys
    with pytest.raises(m.SelectorError):
        m.resolve_one(root, "/per_nu1/*/v", {"nu1": "99"})                 # ZERO -> hard fail
    with pytest.raises(m.SelectorError):
        m.resolve_one(root, "/per_nu1/*/v", {})                           # MULTIPLE -> hard fail
    lst = {"validation_table": [{"label": "a", "fraction_extinct": 0.2},
                                {"label": "b", "fraction_extinct": 0.5}]}
    assert m.resolve_one(lst, "/validation_table/*/fraction_extinct", {"label": "b"}) == 0.5
    with pytest.raises(m.SelectorError):
        m.resolve_one(lst, "/validation_table/*/fraction_extinct", {})     # 2 records -> hard fail
    multi = {"sensitivity": {"base": {"2": {"fraction": 0.33}}, "ic.runif": {"2": {"fraction": 0.31}}}}
    assert m.resolve_one(multi, "/sensitivity/{config}/{nu1}/fraction",
                         {"config": "base", "nu1": "2"}) == 0.33           # multi-named-segment -> one
    tax = {"/sensitivity/base/*/fraction": "regression_anchor",
           "/sensitivity/*/*/fraction": "empirical_discovery"}
    assert m.taxonomy_kind_for("/sensitivity/{config}/{nu1}/fraction",
                               {"config": "base", "nu1": "2"}, tax) == "regression_anchor"   # most-specific wins
    assert m.taxonomy_kind_for("/sensitivity/{config}/{nu1}/fraction",
                               {"config": "ic.runif", "nu1": "2"}, tax) == "empirical_discovery"

    # --- (a) --smoke validator end-to-end over the REAL $CLAIMS ---
    figs, data, manifest = _b6_make_smoke_dir(tmp_path)
    report = tmp_path / "r.md"; valjson = tmp_path / "v.json"
    env = dict(os.environ, DESIGNG_VALIDATE_MODE="smoke", DESIGNG_VALIDATE_IN=str(data),
               DESIGNG_VALIDATE_FIGDIR=str(figs), DESIGNG_CLAIMS_IN=CLAIMS_PATH,
               DESIGNG_REPORT_OUT=str(report), DESIGNG_VALJSON_OUT=str(valjson),
               DESIGNG_MANIFEST=str(manifest))
    r = subprocess.run([_PY, VAL_PATH], env=env, capture_output=True, text=True)
    assert r.returncode == 0, f"smoke validator exit {r.returncode}\nSTDOUT:{r.stdout}\nSTDERR:{r.stderr}"
    assert report.exists() and valjson.exists()
    v = json.load(open(valjson))
    assert v["overall"] == "PASS" and v["exit_code"] == 0
    checks = v["checks"]; assert checks, "no per-claim checks emitted"
    for c in checks:
        assert c["kind"] in ("regression_anchor", "theory_invariant", "empirical_discovery"), c
        assert c["status"] == "SKIPPED_SMOKE", c                          # never PASS/FAIL/hard-fail in smoke
        assert c.get("manuscript_change_required") is None, c             # non-blocking, skipped in smoke
    assert {"regression_anchor", "theory_invariant", "empirical_discovery"} <= {c["kind"] for c in checks}
    assert v["gates"]["audit_integrity"] == "OK"                          # every claim:<id> ref resolves
    assert v["summary"]["n_audit_rows"] > 0                               # coverage audit rendered

    # --- (c) 6-env-var abort: omit DESIGNG_MANIFEST -> exit 2 ---
    keep = {k: val for k, val in os.environ.items()
            if k not in ("DESIGNG_VALIDATE_IN", "DESIGNG_VALIDATE_FIGDIR", "DESIGNG_CLAIMS_IN",
                         "DESIGNG_REPORT_OUT", "DESIGNG_VALJSON_OUT", "DESIGNG_MANIFEST", "DESIGNG_VALIDATE_MODE")}
    keep.update(DESIGNG_VALIDATE_IN=str(data), DESIGNG_VALIDATE_FIGDIR=str(figs), DESIGNG_CLAIMS_IN=CLAIMS_PATH,
                DESIGNG_REPORT_OUT=str(tmp_path / "r2.md"), DESIGNG_VALJSON_OUT=str(tmp_path / "v2.json"))
    r2 = subprocess.run([_PY, VAL_PATH], env=keep, capture_output=True, text=True)
    assert r2.returncode == 2, (r2.returncode, r2.stdout, r2.stderr)
    assert "DESIGNG_MANIFEST" in r2.stderr


def test_quote_check_negative_controls(tmp_path):
    """Mutations kept as executable negative controls: a sign
    (adjacent or across TeX spacing), a decimal or scientific-notation continuation, and a
    corrupted Appendix copy must each defeat the quotation check; the clean text passes."""
    m = _load("VAL")
    claim = {"id": "ctl", "quoted_value": 3.96, "manuscript_tracked": True, "tex": "r_clean.tex",
             "manuscript_token": "3.96", "display_tol": 0.005, "quoted_in": {"r_clean.tex": 1, "a_clean.tex": 1}}
    def run(rtext, atext="row $3.96$ ok"):
        (tmp_path / "r_clean.tex").write_text(rtext); (tmp_path / "a_clean.tex").write_text(atext)
        return m.manuscript_quote_check(claim, str(tmp_path))["status"]
    assert run("gives $3.96$ here") == "FOUND"
    for bad in ("gives $-3.96$", "gives $- 3.96$", r"gives $-\,3.96$", r"gives $3.96\times10^{8}$",
                "gives $3.961$", "gives $13.96$", "gives $3.96e3$", r"gives $3.96\,\times10^{8}$",
                r"gives $3.96\cdot 10^{2}$", r"gives $3.96\times10^{+8}$", "gives $-\n3.96$", "gives $-~3.96$",
                # 2026-11-10: only printed text quotes a number -- not a comment, a label key or a graphics option
                "gives nothing % the old value was $3.96$", r"gives nothing \label{eq:3.96}",
                r"gives \includegraphics[width=3.96in]{figures/f.png}"):
        assert run(bad) != "FOUND", bad
    assert run("gives $3.96$", atext="row $0.7500$") == "TOO_FEW"      # the Appendix copy is checked too
    # The failure message names the file that fell short, with its counts -- not the claim's primary
    # file, which is complete here
    msg = m.quote_shortfall_message("ctl", 3.96, m.manuscript_quote_check(claim, str(tmp_path)))
    assert "a_clean.tex: 0 of 1 pinned" in msg and "r_clean.tex" not in msg, msg
    nclaim = {"id": "ctl2", "quoted_value": 0.0111, "manuscript_tracked": True, "tex": "r_clean.tex", "display_tol": 0.00005}
    for bad in ("value $- 0.0111$", r"value $-\,0.0111$", r"value $0.0111\times 10^{2}$", r"value $0.0111\cdot10^{8}$",
                r"value $0.0111\times10^{+8}$", r"value $0.0111\,\times10^{8}$", "value $-\n0.0111$"):
        (tmp_path / "r_clean.tex").write_text(bad)
        assert m.manuscript_quote_check(nclaim, str(tmp_path))["status"] != "FOUND", bad
    (tmp_path / "r_clean.tex").write_text("value $0.0111$ only")
    assert m.manuscript_quote_check(nclaim, str(tmp_path))["status"] == "FOUND"
    # A quoted value whose repr is in exponent form keeps the half unit of its own last digit (5e-06 for
    # 1e-05, where 0.5 was taken): the other small numbers of a file are not occurrences of it
    eclaim = {"id": "ctl3", "quoted_value": 1e-05, "manuscript_tracked": True, "tex": "r_clean.tex", "display_tol": 5e-07}
    (tmp_path / "r_clean.tex").write_text("values $0.3$, $0.49$ and $-0.2$, and the step $0.00002$")
    assert m.manuscript_quote_check(eclaim, str(tmp_path))["status"] == "NOT_FOUND"
    (tmp_path / "r_clean.tex").write_text(r"values $0.3$ and $1.0\times10^{-5}$")
    r5 = m.manuscript_quote_check(eclaim, str(tmp_path))
    assert (r5["status"], r5["n_occurrences"]) == ("FOUND", 1), r5
    assert [m._half_ulp(q, repr(q)) for q in (5.7e-05, 0.2269, 817, 570.0)] == [5e-07, 5e-05, 0.5, 0.05]
    # With the shipped inventory's own tolerances, the four-decimal number that one more extinct orbit of the
    # 3600 would print is not an occurrence of an extinction fraction (0.001 made 0.2272 an occurrence of 0.2269)
    shipped = {c["id"]: c for c in json.load(open(CLAIMS_PATH))["claims"]}
    for cid, near in (("ext_frac_nu1_2", "0.2272"), ("ext_frac_nu1_3", "0.0114")):
        sc = dict(shipped[cid], tex="r_clean.tex", quoted_in={"r_clean.tex": 1})
        (tmp_path / "r_clean.tex").write_text(f"fraction ${near}$")
        assert m.manuscript_quote_check(sc, str(tmp_path))["status"] == "NOT_FOUND", (cid, near)
        (tmp_path / "r_clean.tex").write_text(f"fraction ${sc['quoted_value']}$")
        assert m.manuscript_quote_check(sc, str(tmp_path))["status"] == "FOUND", cid


def test_quote_check_grouped_numerals_are_whole_numbers(tmp_path):
    """The TeX thousands separator ``{,}`` continues a number. The prefix of
    ``40{,}000`` is not an occurrence of 40 (exact-token path and numeric path), ``40{,}000`` is one
    occurrence of 40000, a tail group (``600`` in ``3{,}600``) is not a number of its own, and the
    shipped claim ``absorb_before_tail_nu1_3`` is quoted only in the file that states it, the Supplement,
    where it counts the occurrences of 40 (its pinned line, the table row and two step counts, n <= 40 and N = 40),
    none inside a grouped numeral; the two 40s of Section 6 are the step counts of a figure caption."""
    m = _load("VAL")
    def run(claim, text):
        (tmp_path / "a_clean.tex").write_text(text)
        return m.manuscript_quote_check(claim, str(tmp_path))
    def cl(cid, value, token=None, tol=0.5):
        c = {"id": cid, "quoted_value": value, "manuscript_tracked": True, "tex": "a_clean.tex",
             "display_tol": tol, "quoted_in": {"a_clean.tex": 1}}
        if token: c["manuscript_token"] = token
        return c
    for c in (cl("tok40", 40, "40"), cl("num40", 40)):
        assert run(c, "all $40$ at")["status"] == "FOUND", c["id"]
        assert run(c, "and $40{,}000$ retained iterations")["status"] == "NOT_FOUND", c["id"]
        r = run(c, "$40{,}000$ predator seeds; all $40$ at $\\nu_1=3$; entry $40/3600$")
        assert (r["status"], r["n_occurrences"]) == ("FOUND", 2), (c["id"], r)
    # the digits of an exponent are not a number of their own, on either matching path
    for c in (cl("tok40", 40, "40"), cl("num40", 40)):
        for good in ("all $40$ at", r"$n\le 40$", "$N=40$", "$x_{40}$"):
            assert run(c, good)["status"] == "FOUND", (c["id"], good)
        for bad in ("$1e40$", "$1E+40$", r"$1\times10^{40}$", r"$1{\times}10^{40}$", r"$10^{40}$", r"$10^{-40}$", r"$3\cdot 10^{ 40 }$"):
            assert run(c, bad)["status"] == "NOT_FOUND", (c["id"], bad)
    assert run(cl("num40000", 40000), "and $40{,}000$ retained")["status"] == "FOUND"
    assert run(cl("num40000", 40000), "and $40000$ retained")["status"] == "FOUND"
    assert run(cl("tok40000", 40000, "40{,}000"), "and $40{,}000$ retained")["status"] == "FOUND"
    assert run(cl("num3600", 3600), "entry $3{,}600$ points")["status"] == "FOUND"
    assert run(cl("num600", 600), "entry $3{,}600$ points")["status"] == "NOT_FOUND"
    assert run(cl("num140000", 140000), "a $140{,}000$-step run")["status"] == "FOUND"
    # a malformed grouping (a fourth digit after a group) is no number at all: neither its grouped head nor its tail matches
    for value, text in ((40000, "$40{,}0001$"), (3600, "$3{,}6000$"), (-30300, "$-30{,}3000$"), (1, "$40{,}0001$"), (0, "$3{,}6000$")):
        assert run(cl("digit_tail", value, tol=0), text)["status"] == "NOT_FOUND", (value, text)
    # the shipped inventory against the shipped text, file by file. The claim is stated in the Supplement (the Appendix until
    # v27), whose four 40s are the claim's two and two of the four step counts (n <= 40, N = 40) that v29 moved there from v27's
    # Results, none inside a grouped numeral. Section 6 prints 40 only as the other two step counts, which v30 moved there
    # with the caption of the decay-rate figure, so Section 6 is not among the claim's files (2026-11-08;
    # amendment 27 ruled the same for v18's Results). The coverage audit reads no two-digit number, so no audit row can name
    # them; this test pins them instead: a third 40 in Section 6, or a 40 there in another sentence, fails here.
    claim = next(c for c in json.load(open(CLAIMS_PATH))["claims"] if c["id"] == "absorb_before_tail_nu1_3")
    assert claim["quoted_in"] == {"supplementary_clean.tex": 4}, claim["quoted_in"]
    assert claim["tex"] == "supplementary_clean.tex" and claim["manuscript_min_occurrences"] == 4, claim
    res = os.path.join(os.path.dirname(_CODES_DIR), "resultsAndDiscussions_clean.tex")
    res_src = open(res).read().split("\n")
    res_hits = [res_src[ln - 1] for ln, _ in m._quote_hits_in_file(claim, res)]
    assert len(res_hits) == 2 and "n \\le 40" in res_hits[0] and "after $40$ steps" in res_hits[1], res_hits
    for name, n in claim["quoted_in"].items():
        path = os.path.join(os.path.dirname(_CODES_DIR), name)
        hits = m._quote_hits_in_file(claim, path)
        src = open(path).read().split("\n")
        assert len(hits) == n, (name, hits)
        assert all("{,}" not in src[ln - 1] for ln, _ in hits), (name, hits)
        if name == claim["tex"]:
            assert claim["line"] in [ln for ln, _ in hits], (name, hits)


# ---------------------------------------------------------------------------
# 2026-10-18: initial-condition archive, exceptional seeds,
# validator negative controls, exact enumeration coverage, certificate
# ---------------------------------------------------------------------------
_CODES_DIR = os.path.dirname(DYN_PATH)


def _load_path(name):
    p = os.path.join(_CODES_DIR, name)
    spec = importlib.util.spec_from_file_location(name[:-3], p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def test_r04_ic_archive_roundtrip_verify_and_tamper(tmp_path):
    """The archive writer/loader keeps the exact bytes (sha256 over C-order float64) and the
    loader REFUSES an edited array: a one-ULP change in a single predator seed -- the size of
    the CPU-dispatch effect -- is detected, not silently used."""
    import numpy as np, pytest
    dyn = _load("DYN")
    ics = dyn._random_ics("rlogy", 500, dyn._SENS_SEED)
    meta = dyn.write_ic_archive(ics, "rlogy", 500, dyn._SENS_SEED, str(tmp_path))
    assert meta["dtype"] == "float64" and meta["shape"] == [500, 2]
    assert meta["sha256_data"] == __import__("hashlib").sha256(np.ascontiguousarray(ics).tobytes()).hexdigest()
    assert "cpu_features_enabled" in meta["runtime"] and meta["runtime"]["numpy"] == np.__version__
    arr, meta2 = dyn.load_ic_archive("rlogy", 500, dyn._SENS_SEED, str(tmp_path))
    assert np.array_equal(arr, ics) and meta2["sha256_data"] == meta["sha256_data"]
    # tamper: nudge one predator seed by one ULP and rewrite the .npy without the sidecar
    bad = ics.copy(); bad[7, 1] = np.nextafter(bad[7, 1], 1.0)
    np.save(os.path.join(str(tmp_path), dyn._ic_archive_name("rlogy", 500, dyn._SENS_SEED) + ".npy"), bad)
    with pytest.raises(ValueError):
        dyn.load_ic_archive("rlogy", 500, dyn._SENS_SEED, str(tmp_path))
    # absent archive -> None (the producer then records a fresh draw, never an error)
    assert dyn.load_ic_archive("rlogy", 501, dyn._SENS_SEED, str(tmp_path)) is None


def test_r04_shipped_archive_matches_sidecar_and_reports_dispatch_agreement():
    """The SHIPPED archive of the C10 (log-uniform) seeds verifies against its sidecar and
    carries the documented checksum; the fresh draw on this host is compared to it with the
    ULP distance (informational: 0 on the authors' AVX2 host, 11,147 one-ULP changes on an
    AVX-512 host per an external measurement -- a platform fact, not a test failure)."""
    import numpy as np
    dyn = _load("DYN")
    d = dyn._ic_archive_dir()
    assert os.path.isdir(d), f"shipped IC archive directory missing: {d}"
    arr, meta = dyn.load_ic_archive("rlogy", 40000, dyn._SENS_SEED, d)
    assert arr.shape == (40000, 2) and meta["seed"] == dyn._SENS_SEED
    assert meta["sha256_data"] == "384d22a666e30c09903596fef3478d45292f9f134fe61270605b9d798674e1ef"
    fresh = dyn._random_ics("rlogy", 40000, dyn._SENS_SEED)
    dist = dyn._ulp_distance(arr, fresh)
    n_diff = int(np.count_nonzero(dist)); max_ulp = int(dist.max())
    print(f"fresh draw vs shipped archive on this host: {n_diff} value(s) differ, max {max_ulp} ULP; "
          f"numpy CPU features {dyn._numpy_dispatch_info().get('cpu_features_enabled')}")
    assert np.array_equal(arr[:, 0], fresh[:, 0]), "prey seeds involve no transcendental transform and must agree exactly"
    assert max_ulp <= 4, "a disagreement beyond a few ULP is not the documented dispatch effect; investigate"
    ics, prov = dyn.resolve_ics("rlogy", 40000, dyn._SENS_SEED)
    assert prov["source"] == "archive" and prov["used_sha256"] == meta["sha256_data"]
    assert prov["n_values_differ"] == n_diff and prov["max_ulp"] == max_ulp
    assert np.array_equal(ics, arr)


def test_r04_resample_mode_uses_the_fresh_draw(monkeypatch):
    import numpy as np
    dyn = _load("DYN")
    monkeypatch.setenv("DESIGNG_IC_MODE", "resample")
    ics, prov = dyn.resolve_ics("rlogy", 40000, dyn._SENS_SEED)
    assert prov["source"].startswith("resampled") and prov["archive_sha256"] is None
    assert np.array_equal(ics, dyn._random_ics("rlogy", 40000, dyn._SENS_SEED))


def test_r15_exceptional_seed_prefilter_uses_angles_not_dyadic_states():
    """Oracle in exact rationals: 1/4 -> 3/4 -> 3/4 under 4x(1-x), so the dyadic STATE 1/4 is
    NOT exceptional (it never reaches 0); 1/2 -> 1 -> 0 is; 0 and 1 are."""
    from fractions import Fraction
    dyn = _load("DYN")
    x = Fraction(1, 4); orbit = [x]
    for _ in range(5):
        x = 4 * x * (1 - x); orbit.append(x)
    assert orbit[1] == Fraction(3, 4) and all(v == Fraction(3, 4) for v in orbit[1:])   # exact oracle
    assert dyn._is_exceptional_x0(0.25) is False
    assert dyn._is_exceptional_x0(0.5) is True and dyn._is_exceptional_x0(0.0) is True and dyn._is_exceptional_x0(1.0) is True
    assert dyn._is_exceptional_x0(0.375) is False          # 3/8: dyadic state, orbit 3/8 -> 15/16 -> 15/64 ... never 0
    assert dyn._is_exceptional_x0(0.123456789) is False


def test_r05_validator_rejects_nonfinite_values_unless_declared():
    """A NaN (or an undeclared infinity) from a producer is a HARD failure, never CONSISTENT /
    OK_NONBLOCKING; a claim may declare an admissible infinity (allowed_nonfinite)."""
    m = _load("VAL")
    base = {"id": "ctl_nan", "producer": "p", "producer_field": "/v", "producer_selector": {},
            "producer_kind": "empirical_discovery", "manuscript_tracked": True, "quoted_value": -0.0495,
            "display_tol": 0.00005, "block": "quasiperiodic", "producer_emits_no_taxonomy": True, "tex": "x_clean.tex"}
    for bad in (float("nan"), float("inf"), float("-inf")):
        r = m.evaluate_claim(dict(base), {"p": ({"v": bad}, "p.json")}, "full")
        assert r["status"] == "NONFINITE" and r["hard_fail"] is True, r
    ok = m.evaluate_claim(dict(base, allowed_nonfinite=["-inf"]), {"p": ({"v": float("-inf")}, "p.json")}, "full")
    assert ok["status"] != "NONFINITE" and not ok["hard_fail"], ok
    still = m.evaluate_claim(dict(base, allowed_nonfinite=["-inf"]), {"p": ({"v": float("nan")}, "p.json")}, "full")
    assert still["status"] == "NONFINITE"
    # a wrong-sign empirical value is reported as a required manuscript change, never as consistent
    wrong = m.evaluate_claim(dict(base), {"p": ({"v": +0.5}, "p.json")}, "full")
    assert wrong["status"] == "OK_NONBLOCKING" and wrong["manuscript_change_required"] is True
    right = m.evaluate_claim(dict(base), {"p": ({"v": -0.04951}, "p.json")}, "full")
    assert right["manuscript_change_required"] is False


def _b6_full_env(tmp_path, data_dir, texdir, strict=False):
    """A FULL-mode validator environment over a given producer data dir. The figure dir is a
    temporary copy of the shipped manuscript figures plus a placeholder file for every
    producer-diagnostic PNG a block requires but the manuscript does not include (the
    completeness gate checks their existence only); the manifest is built from the copies."""
    import hashlib, glob, re, shutil
    root = os.path.dirname(_CODES_DIR)
    figdir = tmp_path / "figs"
    if not figdir.exists():
        figdir.mkdir()
        for png in glob.glob(os.path.join(root, "figures", "*.png")):
            shutil.copy(png, figdir)
        m = _load("VAL")
        for spec in m.BLOCK_SPEC.values():
            for fig in spec.get("figures", []):
                if not (figdir / fig).exists():
                    (figdir / fig).write_bytes(b"placeholder: producer diagnostic figure, not a manuscript figure")
    figdir = str(figdir)
    manifest = tmp_path / "manifest.md"
    with open(manifest, "w") as fh:
        for png in sorted(glob.glob(os.path.join(figdir, "*.png"))):
            fh.write(f"| {png} | designG_dynamics_v4.py | {hashlib.sha256(open(png, 'rb').read()).hexdigest()} |\n")
    env = dict(os.environ, DESIGNG_VALIDATE_MODE="full", DESIGNG_VALIDATE_IN=str(data_dir),
               DESIGNG_VALIDATE_FIGDIR=figdir, DESIGNG_CLAIMS_IN=CLAIMS_PATH,
               DESIGNG_REPORT_OUT=str(tmp_path / "r.md"), DESIGNG_VALJSON_OUT=str(tmp_path / "v.json"),
               DESIGNG_MANIFEST=str(manifest), DESIGNG_TEXDIR=str(texdir))
    if strict:
        env["DESIGNG_VALIDATE_STRICT"] = "1"
    return env


def test_r05_readiness_and_strict_mode_on_corrupted_copies(tmp_path):
    """Full-mode negative controls on COPIES of the shipped reference data and section files:
    (a) the clean copy is SUBMISSION_READY; (b) an empirical sign flip is reported as
    CHANGES_REQUIRED and, under DESIGNG_VALIDATE_STRICT=1, fails the run; (c) a NaN fails the
    run outright; (d) editing one of several occurrences of a duplicated onset number in
    the Supplement (Appendix A.2 until v27; 3.958 -> 3.956 on its own line) is caught (claim and audit), the full run FAILS."""
    import shutil, glob, re
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    data = tmp_path / "data"; shutil.copytree(shipped, data, ignore=shutil.ignore_patterns("superseded", "_superseded*", "ics"))
    texdir = tmp_path / "tex"; texdir.mkdir()
    for f in glob.glob(os.path.join(root, "*_clean.tex")):
        shutil.copy(f, texdir)

    def run(strict=False):
        r = subprocess.run([_PY, VAL_PATH], env=_b6_full_env(tmp_path, data, texdir, strict), capture_output=True, text=True)
        v = json.load(open(tmp_path / "v.json"))
        return r.returncode, v, r.stdout + r.stderr

    rc, v, out = run()
    assert rc == 0 and v["overall"] == "PASS" and v["readiness"] == "SUBMISSION_READY", out[-2000:]

    # (b) sign flip of the empirical QP normal exponent -> CHANGES_REQUIRED; strict -> exit 1
    qp = data / "designG_quasiperiodic_v2_values.json"; orig = qp.read_text()
    doc = json.loads(orig); doc["normal_lyap"] = +0.5; qp.write_text(json.dumps(doc))
    rc, v, out = run()
    assert rc == 0 and v["readiness"].startswith("CHANGES_REQUIRED"), (rc, v["readiness"], out[-1500:])
    rc, v, out = run(strict=True)
    assert rc == 1 and v["overall"] == "FAIL" and v["readiness"].startswith("CHANGES_REQUIRED"), (rc, v["readiness"])
    # 2026-11-10: a failed gate is a failure, also beside a required change. With the manifest emptied as
    # well, no claim row fails, and the strict run's readiness must read FAIL (it kept the label CHANGES_REQUIRED).
    empty = tmp_path / "empty_manifest.md"; empty.write_text("")
    env = dict(_b6_full_env(tmp_path, data, texdir, strict=True), DESIGNG_MANIFEST=str(empty))
    r = subprocess.run([_PY, VAL_PATH], env=env, capture_output=True, text=True)
    v = json.load(open(tmp_path / "v.json"))
    assert r.returncode == 1 and v["gates"]["manifest"] and v["manuscript_changes"], (r.returncode, v["gates"]["manifest_status"])
    assert not any(c.get("hard_fail") or c.get("blocking_fail") for c in v["checks"])
    assert v["overall"] == "FAIL" and v["readiness"] == "FAIL", v["readiness"]
    # (c) NaN -> NONFINITE hard failure
    doc["normal_lyap"] = float("nan"); qp.write_text(json.dumps(doc))
    rc, v, out = run()
    assert rc == 1 and any(c["status"] == "NONFINITE" for c in v["checks"]), out[-1500:]
    qp.write_text(orig)
    rc, v, out = run(); assert rc == 0 and v["readiness"] == "SUBMISSION_READY"

    # 2026-11-10: a deterministic count that moves by one orbit of its grid fails its anchor -- 41 instead of
    # 40 of the 3600 orbits at nu1 = 3, 3288 instead of 3287 of the 10^4 orbits of the base grid at nu1 = 2 (the
    # tolerances 0.001 and 0.005 let three and fifty orbits pass)
    dyn = data / "designG_dynamics_v4_values.json"; dyn_orig = dyn.read_text(); ddoc = json.loads(dyn_orig)
    row = next(r_ for r_ in ddoc["validation_table"] if isinstance(r_, dict) and r_.get("label") == "ext_frac_nu1_3")
    assert row["fraction_extinct"] == 40 / 3600 and ddoc["sensitivity"]["base"]["2"]["fraction"] == 0.3287
    row["fraction_extinct"] = 41 / 3600; ddoc["sensitivity"]["base"]["2"]["fraction"] = 0.3288
    dyn.write_text(json.dumps(ddoc))
    rc, v, out = run()
    failed = {c["name"] for c in v["checks"] if c.get("blocking_fail")}
    assert rc == 1 and failed == {"ext_frac_nu1_3", "sens_base_nu1_2"}, (rc, sorted(failed), out[-1500:])
    # ... and both quotations are reported as no longer agreeing with the record (display tolerance: half a unit of
    # the last printed digit)
    assert {c["name"] for c in v["manuscript_changes"]} == {"ext_frac_nu1_3", "sens_base_nu1_2"}, v["manuscript_changes"]
    dyn.write_text(dyn_orig)
    rc, v, out = run(); assert rc == 0 and v["readiness"] == "SUBMISSION_READY"

    # (d) a duplicated onset quotation edited at one site (the onset-ensemble paragraph moved from the
    # Results to Appendix A.2 in v18 and to the Supplement in v29; the fixture follows the text, the assertions are unchanged)
    res = texdir / "supplementary_clean.tex"; txt = res.read_text()
    lines = txt.split("\n")
    idx = next(i for i, ln in enumerate(lines) if ln.startswith("$3.958$."))
    lines[idx] = lines[idx].replace("$3.958$", "$3.956$", 1)
    # 2026-11-10: the same run also edits the Supplement's quotation of the two-cycle point, whose claim tc_p_x
    # has the Methods as its primary file; the failure must be reported for the Supplement (the primary file was named)
    assert any("0.502328" in ln for ln in lines)
    lines = [ln.replace("0.502328", "0.502329") for ln in lines]; res.write_text("\n".join(lines))
    rc, v, out = run()
    assert rc == 1, out[-2500:]
    assert any(c["name"] == "ons_seed0_sustained_0p002" and c.get("manuscript_quote", {}).get("status") != "FOUND" for c in v["checks"])
    assert "STALE" in str(v["gates"]["audit_integrity"])
    tc = [q for q in v["gates"]["manuscript_quotes"] if q.startswith("tc_p_x:")]
    assert len(tc) == 1 and "supplementary_clean.tex: 0 of " in tc[0] and "methods_clean.tex" not in tc[0], tc
    # (e) a claim bound to the text of its line. The closed-form exponent at nu1 = 4 is printed as the
    # cell "$0$" of Table S4, too coarse for an occurrence count; the claim pins the text of that row, and an
    # edited cell fails the run.
    sub = tmp_path / "site"; sub.mkdir()
    _, data2, tex2 = _val_copies(sub)
    claim = next(c for c in json.load(open(CLAIMS_PATH))["claims"] if c["id"] == "lambda_perp_nu1_4")
    sites = {s["line"]: s["text"] for s in claim["pinned_text"]}
    assert claim["line"] in sites and len(sites) == 2, claim
    sup = tex2 / "supplementary_clean.tex"; lines2 = sup.read_text().split("\n")
    row = min(sites)                                   # the row of Table S4
    assert sites[row] in lines2[row - 1] and "& $0$ &" in lines2[row - 1], lines2[row - 1]
    rc, v, out = _val_run(sub, data2, tex2)
    assert rc == 0 and not v["gates"]["manuscript_quotes"], out[-800:]
    lines2[row - 1] = lines2[row - 1].replace("& $0$ &", "& $1$ &", 1); sup.write_text("\n".join(lines2))
    rc, v, out = _val_run(sub, data2, tex2)
    assert rc == 1 and any("lambda_perp_nu1_4" in g and "SITE_CHANGED" in g for g in v["gates"]["manuscript_quotes"]), (rc, v["gates"]["manuscript_quotes"])


def test_r02_exact_enumeration_covers_every_orbit_independently_of_the_maximiser():
    """Coverage is MEASURED: the exact-angle enumeration visits 8,799 interior orbits /
    130,485 points through period 16, the per-period counts equal the Moebius/Lyndon
    formula, every orbit's ln x average is evaluated, and the period-2 orbit matches its
    closed form ln(5/16)/2 -- none of which depends on the maximising value ln(3/4)."""
    import math
    erg = _load_path("designG_ergodic_quadrature_v2.py")
    rec, rows = erg.enumerate_periodic_orbits(16)
    assert rec["coverage_complete"] is True and rec["least_period_count_mismatches"] == 0
    assert rec["interior_orbits_enumerated"] == 8799 == rec["interior_orbits_evaluated"]
    assert rec["interior_points_enumerated"] == 130485 and rec["boundary_orbits_excluded"] == 1
    ev = rec["evaluated_orbits_by_least_period"]
    assert ev[1] == 1 and ev[2] == 1 and ev[3] == 2 and ev[16] == 4080 and sum(ev.values()) == 8799
    assert abs(rec["max_orbit_avg_lnx_by_least_period"][2] - 0.5 * math.log(5.0 / 16.0)) < 1e-12
    assert rec["argmax_orbit_angle"] == "1/3" and abs(rec["max_orbit_avg_lnx"] - math.log(0.75)) < 1e-12
    assert rec["n_least_periods_2_to_max_with_avg_at_or_above_ln34"] == 0
    assert all(r.passed for r in rows), [r.name for r in rows if not r.passed]
    # a short enumeration is exact too, and its counts are the formula's (periods 1..4)
    rec4, _ = erg.enumerate_periodic_orbits(4)
    assert rec4["interior_orbits_enumerated"] == 1 + 1 + 2 + 3 and rec4["coverage_complete"] is True


def test_r01_two_cycle_certificate_passes_and_fails_on_a_wrong_centre():
    """The exact certificate certifies the two-cycle at (4, 46/25); a centre displaced by 1e-6
    (the box no longer contains the orbit) or the Jury-edge-plus parameter with the same
    box must NOT certify -- the certificate can fail, so its pass means something."""
    from fractions import Fraction as Q
    cert = _load_path("certify_period_two_v1.py")
    rec = cert.certify()
    assert rec["certified"] is True and all(c["ok"] for c in rec["checks"])
    assert rec["contraction_bound_q"] < 1e-6 and rec["multipliers_real"] is True
    lo, hi = rec["multipliers_midpoint_float"]
    assert -1 < lo < 0 < hi < 1
    c = [Q(v) for v in rec["centre"]]
    bad = cert.certify(centre=[c[0] + Q(1, 10 ** 6), c[1]])
    assert bad["certified"] is False
    off = cert.certify(nu1=Q(19, 10), centre=c)   # nu1 = 1.9 > nu1_J: no two-cycle near this point
    assert off["certified"] is False


def test_r16_common_horizon_continuation_matches_a_synchronous_oracle():
    """Independent oracle: while the predator stays representable, iterate
    the synchronous Design-G map directly in float64 and compare log10 y_N with the kernel's
    log-space continuation below the floor. (x0, y0) = (0.2, 1e-249) at mu0 = 4, nu1 = 1
    with floor 1e-250 crosses the floor at the first step and stays above 5e-324 for N = 8.
    A second case starts at the prey value 1/2, whose image is exactly 1 and then 0: the
    predator must survive one more step (y_2 = y_1 with x_1 = 1) before it becomes 0."""
    import math, numpy as np
    dyn = _load("DYN")
    mu0, nu1, floor = 4.0, 1.0, 1e-250

    def oracle_log10_yN(x0, y0, N):
        x, y = x0, y0
        for _ in range(N):
            x, y = mu0 * x * (1.0 - x) * (1.0 - y), nu1 * x * y * (1.0 - y)
        return math.log10(y) if y > 0.0 else -math.inf

    grid = np.array([[0.2, 1e-249]])
    out = dyn._designG_telescoping_rate(grid, mu0, nu1, 8, floor)
    max_log10_y_h = out[6]
    assert abs(max_log10_y_h - oracle_log10_yN(0.2, 1e-249, 8)) < 1e-9, (max_log10_y_h, oracle_log10_yN(0.2, 1e-249, 8))
    # boundary case: x0 = 1/2 -> x1 = 1 -> x2 = 0; y2 = nu1 * x1 * y1 = y1 (not 0), y3 = 0
    grid2 = np.array([[0.5, 1e-251]])
    out2 = dyn._designG_telescoping_rate(grid2, mu0, nu1, 2, floor)
    assert abs(out2[6] - oracle_log10_yN(0.5, 1e-251, 2)) < 1e-9, (out2[6], oracle_log10_yN(0.5, 1e-251, 2))
    out3 = dyn._designG_telescoping_rate(grid2, mu0, nu1, 3, floor)
    assert out3[6] < -299.0 and oracle_log10_yN(0.5, 1e-251, 3) == -math.inf   # predator exactly 0 at N = 3


def test_r05_validator_rejects_json_strings_and_malformed_required_fields(tmp_path):
    """Negative controls on COPIES of the shipped data: a JSON STRING "NaN" or
    "not a number" where a number is required, and a NaN in a completeness-required field
    that no direct claim compares, must each fail the full run; a malformed value must never
    leave the previous (successful) report and VALJSON in place -- the run writes a FAILED one."""
    import shutil, glob
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    data = tmp_path / "data"; shutil.copytree(shipped, data, ignore=shutil.ignore_patterns("superseded", "_superseded*", "ics"))
    texdir = tmp_path / "tex"; texdir.mkdir()
    for f in glob.glob(os.path.join(root, "*_clean.tex")):
        shutil.copy(f, texdir)

    def run():
        r = subprocess.run([_PY, VAL_PATH], env=_b6_full_env(tmp_path, data, texdir, strict=True), capture_output=True, text=True)
        v = json.load(open(tmp_path / "v.json"))
        return r.returncode, v, r.stdout + r.stderr

    rc, v, out = run(); assert rc == 0 and v["readiness"] == "SUBMISSION_READY", out[-1500:]
    clean_report = (tmp_path / "r.md").read_bytes(); clean_val = (tmp_path / "v.json").read_bytes()

    qp = data / "designG_quasiperiodic_v2_values.json"; qp_orig = qp.read_text()
    # (a) a JSON string where a number is required (a direct claim)
    doc = json.loads(qp_orig); doc["normal_lyap"] = "NaN"; qp.write_text(json.dumps(doc))
    rc, v, out = run()
    assert rc == 1 and any(c["status"] == "NONNUMERIC" and c["name"] == "qp_normal_lyap_nu1_4p9" for c in v["checks"]), out[-1500:]
    assert (tmp_path / "r.md").read_bytes() != clean_report and (tmp_path / "v.json").read_bytes() != clean_val   # a fresh FAILED record replaced the clean one
    # (b) a NaN in a completeness-required field that no direct claim compares
    doc = json.loads(qp_orig); doc["fft_peak_snr"] = float("nan"); qp.write_text(json.dumps(doc))
    rc, v, out = run()
    assert rc == 1 and any("fft_peak_snr" in m for m in v["gates"]["completeness"]), (rc, v["gates"]["completeness"][:3])
    qp.write_text(qp_orig)
    # (c) a string in a sensitivity fraction (an empirical claim with a Wilson interval)
    dyn = data / "designG_dynamics_v4_values.json"; dyn_orig = dyn.read_text()
    doc = json.loads(dyn_orig); doc["sensitivity"]["ic.runif"]["2"]["fraction"] = "not a number"; dyn.write_text(json.dumps(doc))
    rc, v, out = run()
    assert rc == 1 and any(c["status"] in ("NONNUMERIC", "SCHEMA_ERROR") for c in v["checks"]), out[-1500:]
    dyn.write_text(dyn_orig)
    rc, v, out = run(); assert rc == 0 and v["readiness"] == "SUBMISSION_READY"


def test_r05_validator_crash_replaces_stale_reports(tmp_path):
    """A producer JSON that cannot be parsed, or a container that is null or of
    the wrong type, must not leave the previous SUCCESSFUL report and VALJSON on disk -- the
    run exits 1 AND writes FAILED records in their place."""
    import shutil, glob
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    data = tmp_path / "data"; shutil.copytree(shipped, data, ignore=shutil.ignore_patterns("superseded", "_superseded*", "ics"))
    texdir = tmp_path / "tex"; texdir.mkdir()
    for f in glob.glob(os.path.join(root, "*_clean.tex")):
        shutil.copy(f, texdir)

    def run():
        r = subprocess.run([_PY, VAL_PATH], env=_b6_full_env(tmp_path, data, texdir, strict=True), capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    rc, out = run(); assert rc == 0, out[-1500:]
    clean_report = (tmp_path / "r.md").read_bytes(); clean_val = (tmp_path / "v.json").read_bytes()
    dyn = data / "designG_dynamics_v4_values.json"; orig = dyn.read_text()
    mutations = {
        "truncated JSON": orig[: len(orig) // 2],
        "_taxonomy null": json.dumps(dict(json.loads(orig), _taxonomy=None)),
        "sensitivity null": json.dumps(dict(json.loads(orig), sensitivity=None)),
        "sensitivity list": json.dumps(dict(json.loads(orig), sensitivity=[1, 2, 3])),
    }
    for name, text in mutations.items():
        dyn.write_text(text)
        rc, out = run()
        assert rc == 1, (name, rc, out[-800:])
        v = json.load(open(tmp_path / "v.json"))
        assert v["overall"] == "FAIL" and str(v["readiness"]).startswith("FAIL"), (name, v.get("readiness"))
        assert (tmp_path / "r.md").read_bytes() != clean_report and (tmp_path / "v.json").read_bytes() != clean_val, name
        assert b"FAIL" in (tmp_path / "r.md").read_bytes(), name
    dyn.write_text(orig)
    rc, out = run(); assert rc == 0


def test_r05_validator_config_aborts_replace_stale_reports(tmp_path):
    """An invalid
    DESIGNG_VALIDATE_MODE, an unknown DESIGNG_VALIDATE_SCOPE token and a missing required
    variable (DESIGNG_MANIFEST) end the run before evaluation. With the two output paths
    known, the previous PASS report and VALJSON must be replaced by FAILED records carrying
    the same exit status; each case starts from a fresh clean PASS."""
    import shutil, glob
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    data = tmp_path / "data"; shutil.copytree(shipped, data, ignore=shutil.ignore_patterns("superseded", "_superseded*", "ics"))
    texdir = tmp_path / "tex"; texdir.mkdir()
    for f in glob.glob(os.path.join(root, "*_clean.tex")):
        shutil.copy(f, texdir)

    def run(**overrides):
        env = _b6_full_env(tmp_path, data, texdir, strict=True)
        for k, v in overrides.items():
            if v is None:
                env.pop(k, None)
            else:
                env[k] = v
        r = subprocess.run([_PY, VAL_PATH], env=env, capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    for name, overrides, want_rc in (("invalid mode", {"DESIGNG_VALIDATE_MODE": "bogus"}, 2),
                                     ("unknown scope token", {"DESIGNG_VALIDATE_SCOPE": "bogus"}, 1),
                                     ("missing DESIGNG_MANIFEST", {"DESIGNG_MANIFEST": None}, 2)):
        rc, out = run(); assert rc == 0, (name, out[-1500:])
        v0 = json.load(open(tmp_path / "v.json")); assert v0["overall"] == "PASS", name
        clean_report = (tmp_path / "r.md").read_bytes(); clean_val = (tmp_path / "v.json").read_bytes()
        rc, out = run(**overrides)
        assert rc == want_rc, (name, rc, out[-600:])
        v = json.load(open(tmp_path / "v.json"))
        assert v["overall"] == "FAIL" and str(v["readiness"]).startswith("FAIL") and v["exit_code"] == want_rc, (name, v.get("readiness"), v.get("exit_code"))
        rep = (tmp_path / "r.md").read_bytes()
        assert rep != clean_report and b"**FAIL**" in rep and (tmp_path / "v.json").read_bytes() != clean_val, name



# =============================================================================
# Assurance negative controls (2026-10-18): success flags, the claims inventory, the QP regime
# parameters, Wilson records, checkpoint identity, an exact predator zero, the period-one enumeration
# =============================================================================


def _val_copies(tmp_path):
    """Copies of the shipped reference data and of the section files, as in the validator negative controls."""
    import shutil, glob
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    data = tmp_path / "data"; shutil.copytree(shipped, data, ignore=shutil.ignore_patterns("superseded", "_superseded*", "ics"))
    texdir = tmp_path / "tex"; texdir.mkdir()
    for f in glob.glob(os.path.join(root, "*_clean.tex")):
        shutil.copy(f, texdir)
    return shipped, data, texdir


def _val_run(tmp_path, data, texdir, **overrides):
    env = _b6_full_env(tmp_path, data, texdir, strict=True)
    for k, v in overrides.items():
        if v is None:
            env.pop(k, None)
        else:
            env[k] = v
    r = subprocess.run([_PY, VAL_PATH], env=env, capture_output=True, text=True)
    return r.returncode, json.load(open(tmp_path / "v.json")), r.stdout + r.stderr


def _flag_mutations():
    """Success-flag corruptions, each applied to a fresh clean record."""
    def cert_failed(d):
        d["certificate"]["certified"] = False
        for c in d["certificate"]["checks"]:
            c["ok"] = False
        d["certificate"]["statement"] = "NOT certified"

    def guard_failed(d):
        d["exact_skew"]["per_nu1"]["2"]["collapse_guard_ok"] = False

    def coverage_failed(d):
        d["periodic_enumeration"]["coverage_complete"] = False

    def count_mismatch(d):
        d["periodic_enumeration"]["least_period_count_mismatches"] = 1

    # Deletions: a missing member or a missing flag must fail, not shrink the set inspected.
    def check_flag_missing(d):
        del d["certificate"]["checks"][2]["ok"]

    def check_removed(d):
        del d["certificate"]["checks"][2]

    def guard_missing(d):
        del d["exact_skew"]["per_nu1"]["3"]["collapse_guard_ok"]

    def skew_member_removed(d):
        del d["exact_skew"]["per_nu1"]["3"]

    return [("certificate not certified", "designG_period_two_certificate_v1.json", cert_failed),
            ("collapse guard failed", "designG_dynamics_v4_values.json", guard_failed),
            ("enumeration coverage incomplete", "designG_ergodic_quadrature_v2_values.json", coverage_failed),
            ("least-period count mismatch", "designG_ergodic_quadrature_v2_values.json", count_mismatch),
            ("certificate check flag missing", "designG_period_two_certificate_v1.json", check_flag_missing),
            ("certificate check removed", "designG_period_two_certificate_v1.json", check_removed),
            ("collapse guard missing", "designG_dynamics_v4_values.json", guard_missing),
            ("exact-skew parameter removed", "designG_dynamics_v4_values.json", skew_member_removed)]


def test_r05_validator_rejects_failed_success_flags(tmp_path):
    """Records stating that the certificate, a collapse guard or the
    enumeration coverage FAILED were validated SUBMISSION_READY, because only the presence and
    the type of those flags was checked. Each corruption must now fail the strict full run."""
    _, data, texdir = _val_copies(tmp_path)
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and v["readiness"] == "SUBMISSION_READY", out[-1500:]
    for name, fn, mutate in _flag_mutations():
        p = data / fn; orig = p.read_text(); doc = json.loads(orig); mutate(doc); p.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and v["overall"] == "FAIL" and str(v["readiness"]).startswith("FAIL"), (name, rc, v.get("readiness"))
        assert v["gates"]["essential_flags"], (name, v["gates"].get("essential_flags_status"))
        p.write_text(orig)
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and v["readiness"] == "SUBMISSION_READY"


def test_r05_reproducer_rejects_failed_success_flags(tmp_path):
    """The reproducer compared numeric leaves only, so the same corruptions
    reproduced as PASS. A corrupted regenerated set against the clean shipped files must now fail
    ((A') boolean leaves and (D) success flags), and so must identical corruptions on BOTH sides
    ((D) alone)."""
    import shutil
    repro = os.path.join(_CODES_DIR, "reproduce_designG_v1.py")
    if not os.path.exists(repro):
        import pytest
        pytest.skip("reproduce_designG_v1.py ships only in the bundle (codes/); run the suite from the bundle root")
    shipped, data, _ = _val_copies(tmp_path)

    def run(regen, ship):
        r = subprocess.run([_PY, repro, "--regen", str(regen), "--shipped", str(ship)], capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    rc, out = run(data, shipped); assert rc == 0 and "OVERALL: PASS" in out, out[-1500:]
    both = tmp_path / "both"; shutil.copytree(data, both)
    for name, fn, mutate in _flag_mutations():
        for d in (data, both):
            p = d / fn; doc = json.loads(p.read_text()); mutate(doc); p.write_text(json.dumps(doc))
        rc, out = run(data, shipped)
        assert rc == 1 and "OVERALL: FAIL" in out and "(D) essential success flags" in out, (name, out[-1200:])
        rc, out = run(both, both)
        assert rc == 1 and "OVERALL: FAIL" in out, (name, out[-1200:])
        for d in (data, both):   # restore from the clean shipped file
            shutil.copy(os.path.join(shipped, fn), d / fn)
    rc, out = run(data, shipped); assert rc == 0, out[-800:]


def test_r05_validator_rejects_a_vacuous_inventory_and_never_calls_a_scope_ready(tmp_path):
    """A claims inventory replaced by {} disabled every scientific check and
    still reported SUBMISSION_READY over zero checks. The inventory is now schema-checked against a
    block list fixed in the validator, and a run restricted by DESIGNG_VALIDATE_SCOPE reports
    SCOPED_PASS, never whole-submission readiness."""
    _, data, texdir = _val_copies(tmp_path)
    good = json.load(open(CLAIMS_PATH))
    cases = {
        "empty object": {},
        "no claims": dict(good, claims=[]),
        "no coverage audit": dict(good, coverage_audit=[]),
        "duplicate id": dict(good, claims=good["claims"] + [good["claims"][0]]),
        "required block without claims": dict(good, claims=[c for c in good["claims"] if c.get("block") != "two_cycle"]),
        # 2026-11-10: the stable-leaf block (eight claims, Figures 5 and 12) is required like the others; the
        # validator's block list had no entry for it, so its loss was no inventory error
        "stable-leaf block without claims": dict(good, claims=[c for c in good["claims"] if c.get("block") != "stable_leaf"]),
        "audit row null": dict(good, coverage_audit=[None]),
        "audit row empty": dict(good, coverage_audit=[{}]),
        "audit row without lines": dict(good, coverage_audit=[dict(good["coverage_audit"][0], lines=[])] + good["coverage_audit"][1:]),
        "tracked claim without quoted_value": dict(good, claims=[{k: v for k, v in c.items() if not (c["id"] == "qp_normal_lyap_nu1_4p9" and k == "quoted_value")} for c in good["claims"]]),
        "claim without manuscript_tracked": dict(good, claims=[{k: v for k, v in c.items() if not (c["id"] == "qp_normal_lyap_nu1_4p9" and k == "manuscript_tracked")} for c in good["claims"]]),
    }
    for name, doc in cases.items():
        cp = tmp_path / "claims.json"; cp.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
        assert rc == 1 and str(v["readiness"]).startswith("FAIL") and v["gates"]["inventory"], (name, rc, v.get("readiness"))
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_VALIDATE_SCOPE="two_cycle")
    assert rc == 0 and v["readiness"].startswith("SCOPED_PASS"), (rc, v.get("readiness"), out[-800:])
    # The block can be named as a scope (it was an unknown token), its eight claims are then the ones
    # evaluated, and its completeness rule asks for the record of the illustration
    n_leaf = sum(1 for c in good["claims"] if c.get("block") == "stable_leaf")
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_VALIDATE_SCOPE="stable_leaf")
    assert rc == 0 and v["readiness"].startswith("SCOPED_PASS") and v["summary"]["n_evaluated"] == n_leaf >= 8, (rc, v.get("readiness"), out[-800:])
    leaf = data / "designG_stable_leaf_v1_values.json"; leaf_orig = leaf.read_text(); leaf.unlink()
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_VALIDATE_SCOPE="stable_leaf")
    assert rc == 1 and any("stable_leaf" in g for g in v["gates"]["completeness"]), (rc, v["gates"].get("completeness"))
    leaf.write_text(leaf_orig)


def test_r05_validator_qp_parameters_must_be_finite_and_match_the_regime(tmp_path, monkeypatch):
    """A string or a NaN in the QP output's parameters bypassed the regime
    check, and nu1 was not compared at all when the source regime is the documented default.
    Each of the record's two smoke indicators fails a full validation on its own, and the
    source regime that the output is compared with must itself be finite."""
    _, data, texdir = _val_copies(tmp_path)
    qp = data / "designG_quasiperiodic_v2_values.json"; orig = qp.read_text()
    for name, field, value in (("mu0 NaN", "mu0", float("nan")), ("mu0 string", "mu0", "invalid"),
                               ("nu1 wrong", "nu1", 2.1), ("nu1 string", "nu1", "4.9")):
        doc = json.loads(orig); doc[field] = value; qp.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and v["gates"]["qp_source"], (name, rc, v["gates"].get("qp_source_status"))
    # producer metadata cannot exempt a full validation: a wrong nu1 labelled as a smoke run fails, for both reasons
    doc = json.loads(orig); doc["nu1"] = 2.1; doc.setdefault("_provenance", {})["qp_source_status"] = "smoke"
    qp.write_text(json.dumps(doc))
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 1 and len(v["gates"]["qp_source"]) == 2, (rc, v["gates"].get("qp_source"))
    # each smoke indicator on its own, with the regime parameters left correct: the status, the boolean and the
    # boolean written as text
    assert json.loads(orig)["_provenance"]["smoke"] is False
    for name, key, value in (("status", "qp_source_status", "smoke"), ("boolean", "smoke", True), ("boolean as text", "smoke", "true")):
        doc = json.loads(orig); doc["_provenance"][key] = value; qp.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and len(v["gates"]["qp_source"]) == 1 and "smoke" in v["gates"]["qp_source"][0], (name, rc, v["gates"].get("qp_source"))
    qp.write_text(orig)
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and not v["gates"]["qp_source"], out[-800:]
    # the source regime must be finite too: a source file whose mu0 is the string "NaN" resolved to a NaN,
    # and the comparison with the output's 3.5 reported nothing
    m = _load("VAL")
    qp_root = json.loads(orig)
    monkeypatch.delenv("DESIGNG_QP_SOURCE", raising=False)
    assert m.qp_source_check(qp_root) == []
    src = tmp_path / "source.json"
    for field, value in (("mu0", "NaN"), ("mu0", "Infinity"), ("nu1", "NaN")):
        src.write_text(json.dumps({"mu0": 3.5, field: value}))
        monkeypatch.setenv("DESIGNG_QP_SOURCE", str(src))
        issues = m.qp_source_check(qp_root)
        assert issues and "finite number is required" in issues[0], (field, value, issues)
    src.write_text(json.dumps({"mu0": 3.5})); assert m.qp_source_check(qp_root) == []


def test_r05_validator_rejects_an_inconsistent_wilson_record(tmp_path):
    """A quotation inside the recorded Wilson interval was called consistent
    without checking that the interval belongs to the produced fraction and its count."""
    m = _load("VAL")
    _, data, texdir = _val_copies(tmp_path)
    dyn = data / "designG_dynamics_v4_values.json"; orig = dyn.read_text(); doc = json.loads(orig)
    cell = doc["sensitivity"]["ic.rlogy"]["3"]
    lo, hi, n = cell["ci_low"], cell["ci_high"], cell["n"]
    nt, nc = cell["n_total"], cell["n_censored"]
    assert m.wilson_record_error(cell["fraction"], lo, hi, n, nt, nc) is None
    assert m.wilson_record_error(0.9, lo, hi, n, nt, nc)              # outside its own interval
    assert m.wilson_record_error(cell["fraction"], hi, lo, n, nt, nc)  # reversed bounds
    assert m.wilson_record_error(cell["fraction"], True, hi, n, nt, nc)  # a bool is not a bound
    assert m.wilson_record_error(cell["fraction"], lo - 0.01, hi, n, nt, nc)  # not the Wilson interval of the count
    assert m.wilson_record_error(cell["fraction"], lo, hi, 19021.5, nt, nc)    # a non-integer denominator
    # the sample identity n = n_total - n_censored is part of the record, and both counts are required
    assert m.wilson_record_error(cell["fraction"], lo, hi, n, nt, nc + 1)
    assert m.wilson_record_error(cell["fraction"], lo, hi, n, None, None)
    assert m.wilson_record_error(cell["fraction"], lo, hi, n, nt, None)
    # endpoint roundoff of the producer's own formula is not an inconsistency: wilson(n, n) has an
    # upper bound of 1 - 2e-16 and wilson(0, n) a lower bound within 1e-16 of 0
    dmod = _load("DYN")
    for k_, n_ in ((40000, 40000), (0, 40000), (1, 1), (0, 1), (19022, 19022)):
        wlo, whi = dmod.wilson(k_, n_)
        assert m.wilson_record_error(k_ / n_, wlo, whi, n_, n_, 0) is None, (k_, n_, wlo, whi)
    cell["fraction"] = 0.9; dyn.write_text(json.dumps(doc))
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 1 and any(c.get("status") == "CI_INCONSISTENT" for c in v["checks"]), (rc, out[-1200:])
    # a valid Wilson record of 1/1 grafted onto the 40000-seed sample fails the sample identity, also
    # when both sample counts are deleted or nulled (missing counts must not switch the check off)
    for counts in ("kept", "deleted", "nulled"):
        doc = json.loads(orig); cell = doc["sensitivity"]["ic.rlogy"]["3"]
        cell["fraction"], cell["n"] = 1.0, 1; cell["ci_low"], cell["ci_high"] = dmod.wilson(1, 1)
        if counts == "deleted":
            del cell["n_total"], cell["n_censored"]
        elif counts == "nulled":
            cell["n_total"] = cell["n_censored"] = None
        dyn.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and any(c.get("status") == "CI_INCONSISTENT" for c in v["checks"]), (counts, rc, out[-1200:])
    # relabelling the method switches nothing off: the claim expects a Wilson record
    doc = json.loads(orig); cell = doc["sensitivity"]["ic.rlogy"]["3"]
    cell["ci_method"], cell["ci_low"], cell["ci_high"] = "none", float("nan"), float("nan")
    dyn.write_text(json.dumps(doc))
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 1 and any(c.get("status") == "CI_METHOD_MISMATCH" for c in v["checks"]), (rc, out[-1200:])
    dyn.write_text(orig)


def test_r05_validator_audit_coverage_cannot_be_vacuous(tmp_path):
    """A coverage audit could cite a claim whose manuscript check was switched off, and an audit of
    meaningless rows passed. A claim that provides coverage must be tracked; a row must name a known
    claim or exclusion category and existing lines; a row without a checkable number may only document
    a symbol; and in full mode every number the section files print must appear in some row."""
    import re
    _, data, texdir = _val_copies(tmp_path)
    good = json.load(open(CLAIMS_PATH))
    rows = good["coverage_audit"]
    cited = next(cid for a in rows for cid in re.findall(r"claim:([A-Za-z0-9_]+)", a["disposition"]))
    sym = next(i for i, a in enumerate(rows) if a["disposition"].startswith("exclusion:symbolic_mention"))
    def swap(i, **kw):
        return dict(good, coverage_audit=[dict(a, **kw) if j == i else a for j, a in enumerate(rows)])
    cases = {
        "cited claim untracked": dict(good, claims=[dict(c, manuscript_tracked=False) if c["id"] == cited else c
                                                    for c in good["claims"]]),
        "vacuous audit": dict(good, coverage_audit=[{"tex": "resultsAndDiscussions_clean.tex", "lines": [999999],
                                                     "literal": "no numeric evidence", "disposition": "claim:"}]),
        "unknown exclusion category": swap(0, disposition="exclusion:made_up (no reason)"),
        "symbolic row claiming another class": swap(sym, disposition="exclusion:closed_form_constant (relabelled)"),
        "symbolic row with a line past the end": swap(sym, lines=[10 ** 6]),
        "rows dropped": dict(good, coverage_audit=rows[:-20]),
    }
    for name, doc in cases.items():
        cp = tmp_path / "claims.json"; cp.write_text(json.dumps(doc))
        rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
        assert rc == 1 and str(v["readiness"]).startswith("FAIL"), (name, rc, v.get("readiness"), out[-600:])
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and v["readiness"] == "SUBMISSION_READY", out[-800:]


def test_r05_validator_coverage_sees_printed_numbers_after_commands(tmp_path):
    r"""The coverage mask hid the rest of a line after \ref or \caption, and every bracketed argument, so
    a printed number there needed no audit row, and coverage compared mantissas without sign or exponent.
    Only non-printed source is masked now and coverage compares complete signed values: deleting the row
    for the 1.6681e-13 printed after a \ref, or adding an unaudited number after a \ref, inside a caption,
    inside bracketed math, in e-notation, with a different exponent (in any supported spelling), as a bare
    power of ten (with TeX spacing, or too small for a binary float), with a minus sign, in a brace group
    that follows a layout command, after a star the command does not take, after a line break, or joined
    across a comment or a comment chain, with spacing inside its exponent, as an unreadable power, in
    verbatim text, or outside the supported precision or range, fails; an unaudited number inside a
    non-printed argument (a graphics width, \phantom) or a comment does not; a number belongs to the line
    of its first digit; and an audit row cannot be satisfied by exponent digits."""
    _, data, texdir = _val_copies(tmp_path)
    good = json.load(open(CLAIMS_PATH))
    res = texdir / "supplementary_clean.tex"; orig = res.read_text()
    rows = [a for a in good["coverage_audit"] if not a["literal"].startswith("1.6681e-13")]
    assert len(rows) == len(good["coverage_audit"]) - 1
    cp = tmp_path / "claims.json"; cp.write_text(json.dumps(dict(good, coverage_audit=rows)))
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
    assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (rc, out[-800:])
    lines = orig.split("\n")
    i_ref = next(i for i, l in enumerate(lines) if "\\ref{" in l and "\\caption" not in l)
    i_cap = next(i for i, l in enumerate(lines) if "\\caption{" in l)
    edits = {"after a ref": (i_ref, lambda s: s + " Table~\\ref{tab:appA-ci} gives a measured fraction of $0.123456789$."),
             "inside a caption": (i_cap, lambda s: s.replace("\\caption{", "\\caption{Measured fraction $0.123456789$. ", 1)),
             "inside bracketed math": (i_ref, lambda s: s + " $x\\in[0.123, 0.456]$"),
             # complete signed values: a covered mantissa or an unsigned value covers nothing else
             "e-notation exponent": (i_ref, lambda s: s + " The measured error is $7e-19$."),
             "times exponent": (i_ref, lambda s: s + " The measured value is $0.2727\\times10^{8}$."),
             "sign": (i_ref, lambda s: s + " The measured value is $-0.2727$.")}
    for name, (i, f) in edits.items():
        L = list(lines); L[i] = f(L[i]); res.write_text("\n".join(L))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (name, rc, out[-800:])
    # every supported exponent spelling keeps its exponent, and a bare power of ten is a number
    for name, extra in (("cdot exponent", "The measured value is $0.2727\\cdot10^{8}$."),
                        ("spaced exponent", "The measured value is $0.2727\\times  10^{8}$."),
                        ("thin-space exponent", "The measured value is $0.2727\\times\\,10^{8}$."),
                        ("line-broken exponent", "The measured value is $0.2727\\times\n10^{8}$."),
                        ("bare power", "The measured error is $10^{-19}$."),
                        ("group after hspace, same line", "\\hspace{1em}{The measured fraction is $0.123456789$.}"),
                        ("group after phantom, same line", "\\phantom{x}{The measured fraction is $0.123456789$.}")):
        res.write_text(orig + "\n" + extra + "\n")
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (name, rc, out[-800:])
    # a comment is not printed: a number in it needs no audit row
    res.write_text(orig + "\n% a note, not printed: 0.123456789\n")
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and v["readiness"] == "SUBMISSION_READY", ("comment", rc, out[-800:])
    # TeX tokens: a star the command does not take is printed material, "\\" is a line break (what
    # follows is text), a digit continued across a comment is one number, TeX spacing inside a power
    # changes nothing, and a tiny power keeps its magnitude and sign (no binary underflow)
    for name, extra in (("phantom star", "\\phantom*{0.123456789}"),
                        ("hphantom star", "\\hphantom*{0.123456789}"),
                        ("line break then hspace", "Line.\\\\hspace{0.123456789}"),
                        ("line break then phantom", "Line.\\\\phantom{0.123456789}"),
                        ("digits joined across a comment", "Measured fraction $0.2727% joins digits\n1$."),
                        ("spaced power 1", "$10^ {-19}$"), ("spaced power 2", "$10 ^{-19}$"),
                        ("spaced power 3", "$10^{ -19}$"), ("line-broken power", "$10\n^{-19}$"),
                        ("tiny power", "$10^{-325}$"), ("tiny negative power", "$-10^{-325}$")):
        res.write_text(orig + "\n" + extra + "\n")
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (name, rc, out[-800:])
    # comment chains, spacing inside an exponent, unreadable powers, verbatim text, precision and range
    for name, extra in (("digits joined across a comment chain", "Measured $0.2727% joins digits\n% another comment\n1$."),
                        ("spaced exponent sign", "$10^{- 19}$"), ("spaced exponent digits", "$10^{-1 9}$"),
                        ("spaced times exponent", "$0.2727\\times10^{- 19}$"), ("fractional power", "$10^{-1.5}$"),
                        ("verbatim text", "\\verb|\\phantom{0.123456789}|"),
                        ("thirty significant digits", "$0.27270000000000000000000000001$"),
                        ("power out of range", "$10^{-1000027}$"), ("negative power out of range", "$-10^{-1000027}$")):
        res.write_text(orig + "\n" + extra + "\n")
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (name, rc, out[-800:])
    for name, extra in (("phantom without star", "\\phantom{0.123456789}"),
                        ("line break then comment", "Line.\\\\% not printed 0.123456789")):
        res.write_text(orig + "\n" + extra + "\n")
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 0 and v["readiness"] == "SUBMISSION_READY", (name, rc, out[-800:])
    res.write_text(orig)
    # a number belongs to the line of its sign or first digit, not to spacing on the line before
    marked = orig + "\nMARKER.\n0.123456789\n"
    res.write_text(marked)
    n_num = marked.split("\n").index("0.123456789") + 1
    for lines_cited, want in (([n_num - 1], 1), ([n_num], 0)):
        row = {"tex": "supplementary_clean.tex", "lines": lines_cited, "literal": "0.123456789 (test row)",
               "disposition": "exclusion:protocol_parameter (test)"}
        cp.write_text(json.dumps(dict(good, coverage_audit=good["coverage_audit"] + [row])))
        rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
        assert rc == want, (lines_cited, rc, v["gates"].get("audit_integrity"), out[-600:])
    # a row keeps a line only through a number equal to one of its tokens at the finer of the two printed
    # precisions: beside a row that covers the number exactly, a row whose coarser token the number merely
    # rounds to (0.1235 for the printed 0.123456789) is a weak placement and fails
    exact_row = {"tex": "supplementary_clean.tex", "lines": [n_num], "literal": "0.123456789 (test row)",
                 "disposition": "exclusion:protocol_parameter (test)"}
    weak_row = dict(exact_row, literal="0.1235 (a coarser token that the printed number rounds to)")
    cp.write_text(json.dumps(dict(good, coverage_audit=good["coverage_audit"] + [exact_row, weak_row])))
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
    assert rc == 1 and "WEAK" in str(v["gates"].get("audit_integrity")), (rc, out[-800:])
    res.write_text(orig)
    # an audit row cannot be satisfied by the digits of an exponent: the Supplement prints 10^{-250}, not -250
    # (the line is located, not hard-coded: edits above it move it)
    n250 = next(i for i, ln in enumerate(orig.split("\n")) if "10^{-250}" in ln) + 1
    fake = dict(good, coverage_audit=good["coverage_audit"] + [{"tex": "supplementary_clean.tex", "lines": [n250],
                "literal": "-250 (fabricated)", "disposition": "exclusion:protocol_parameter (fabricated)"}])
    assert "10^{-250}" in orig.split("\n")[n250 - 1]
    cp.write_text(json.dumps(fake))
    rc, v, out = _val_run(tmp_path, data, texdir, DESIGNG_CLAIMS_IN=str(cp))
    assert rc == 1 and "STALE" in str(v["gates"].get("audit_integrity")), (rc, out[-800:])
    # a brace group that merely follows a layout command is printed text, not one of its arguments
    for name, extra in (("group after hspace", "\\hspace{0pt}{\nMeasured fraction $0.123456789$.\n}"),
                        ("group after setlength", "\\setlength{\\tabcolsep}{4pt}{Measured fraction $0.123456789$.}")):
        res.write_text(orig + "\n" + extra + "\n")
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and "UNCOVERED" in str(v["gates"].get("audit_integrity")), (name, rc, out[-800:])
    # control: a number in a non-printed argument needs no audit row
    assert "\\includegraphics[width=0.86\\linewidth]" in orig
    res.write_text(orig.replace("\\includegraphics[width=0.86\\linewidth]", "\\includegraphics[width=0.8765\\linewidth]", 1))
    rc, v, out = _val_run(tmp_path, data, texdir)
    assert rc == 0 and v["readiness"] == "SUBMISSION_READY", out[-800:]
    res.write_text(orig)
    # a multiplication sign in braces. 5{\times}10^{-3} is the number 0.005, one number and not 5 followed
    # by the power 10^{-3}, so an edited resolution of the onset ensemble or an edited horizon of Table S8 fails.
    m = _load("VAL")
    from decimal import Decimal
    for srctext, want in ((r"$5{\times}10^{-3}$", "0.005"), (r"$2{\cdot}10^{-3}$", "0.002"), (r"$7.5{\times}10^4$", "75000"),
                          (r"$7.5{ \times }10^{4}$", "75000"), (r"$5\times10^{-3}$", "0.005")):
        nums, problems = m._printed_numbers(srctext)
        assert [n[1] for n in nums] == [Decimal(want)] and not problems, (srctext, nums, problems)
    ctl = {"id": "ctl", "quoted_value": 7.5, "manuscript_tracked": True, "tex": "a_clean.tex", "display_tol": 0.05,
           "quoted_in": {"a_clean.tex": 1}}
    (tmp_path / "a_clean.tex").write_text(r"horizon $7.5{\times}10^4$")
    assert m.manuscript_quote_check(ctl, str(tmp_path))["status"] == "NOT_FOUND"      # 75000 is not 7.5
    assert m.manuscript_quote_check(dict(ctl, manuscript_token="7.5"), str(tmp_path))["status"] == "NOT_FOUND"
    sup = texdir / "supplementary_clean.tex"; sup_orig = sup.read_text()
    for old, new in ((r"5{\times}10^{-3}", r"9{\times}10^{-3}"), (r"7.5{\times}10^4", r"8.5{\times}10^4")):
        assert old in sup_orig, old
        sup.write_text(sup_orig.replace(old, new, 1))
        rc, v, out = _val_run(tmp_path, data, texdir)
        assert rc == 1 and v["gates"]["audit_integrity"] != "OK", (old, rc, v["gates"]["audit_integrity"])
    sup.write_text(sup_orig)


def test_r04_checkpoint_reuse_requires_matching_configuration_and_inputs(monkeypatch):
    """Any checkpointed config carrying `_params` used to be reused, so a
    run with DESIGNG_IC_MODE=resample returned archive-based results. Reuse now requires the same
    configuration and the same initial-condition identity (mode and SHA-256 of the array)."""
    import copy
    m = _load("DYN")
    root = os.path.dirname(_CODES_DIR)
    shipped = next(p for p in (os.path.join(_CODES_DIR, "data", "designG"), os.path.join(root, "data", "designG")) if os.path.isdir(p))
    ck = json.load(open(os.path.join(shipped, "designG_dynamics_v4_sensitivity_checkpoint.json")))
    sched, nu1 = m._sensitivity_schedule(False)
    cfg = {c["key"]: c for c in sched}
    budget = float(m._SENS_BUDGET_S_DEFAULT)
    monkeypatch.delenv("DESIGNG_IC_MODE", raising=False)
    for key in cfg:   # the shipped checkpoint matches a default (archive) request, config by config
        assert m._checkpoint_identity_mismatch(ck[key], cfg[key], nu1, 4.0, budget) is None, key
    monkeypatch.setenv("DESIGNG_IC_MODE", "resample")
    for key in ("ic.runif", "ic.rlogy"):
        why = m._checkpoint_identity_mismatch(ck[key], cfg[key], nu1, 4.0, budget)
        assert why and "ic_provenance.mode" in why, (key, why)
    assert m._checkpoint_identity_mismatch(ck["base"], cfg["base"], nu1, 4.0, budget) is None   # a grid config has no random inputs
    monkeypatch.delenv("DESIGNG_IC_MODE")
    edited = copy.deepcopy(ck["base"]); edited["_params"]["T"] = 5000
    assert "T:" in m._checkpoint_identity_mismatch(edited, cfg["base"], nu1, 4.0, budget)
    other = copy.deepcopy(ck["ic.rlogy"]); other["_params"]["ic_provenance"]["used_sha256"] = "0" * 64
    assert "used_sha256" in m._checkpoint_identity_mismatch(other, cfg["ic.rlogy"], nu1, 4.0, budget)
    bare = copy.deepcopy(ck["ic.rlogy"]); del bare["_params"]["ic_provenance"]
    assert m._checkpoint_identity_mismatch(bare, cfg["ic.rlogy"], nu1, 4.0, budget) == "no ic_provenance record"
    # completeness: a requested cell, a result field or a consistent sample missing -> recompute
    gone = copy.deepcopy(ck["base"]); del gone["3"]
    assert "cell nu1=3: missing" in m._checkpoint_identity_mismatch(gone, cfg["base"], nu1, 4.0, budget)
    thin = copy.deepcopy(ck["ic.rlogy"]); del thin["2"]["ci_high"]
    assert "ci_high" in m._checkpoint_identity_mismatch(thin, cfg["ic.rlogy"], nu1, 4.0, budget)
    odd = copy.deepcopy(ck["ic.rlogy"]); odd["4"]["n"] = odd["4"]["n"] + 1
    assert "inconsistent sample" in m._checkpoint_identity_mismatch(odd, cfg["ic.rlogy"], nu1, 4.0, budget)
    # present but unusable values are not reusable either
    for bad_val in (None, float("nan"), "uncomputed", 1.5):
        broken = copy.deepcopy(ck["base"]); broken["3"]["fraction"] = bad_val
        assert "unusable" in (m._checkpoint_identity_mismatch(broken, cfg["base"], nu1, 4.0, budget) or ""), bad_val
    noci = copy.deepcopy(ck["ic.rlogy"]); noci["3"]["ci_low"] = None
    assert "unusable Wilson interval" in m._checkpoint_identity_mismatch(noci, cfg["ic.rlogy"], nu1, 4.0, budget)
    detci = copy.deepcopy(ck["base"]); detci["2"]["ci_low"] = 0.1
    assert "carries an interval" in m._checkpoint_identity_mismatch(detci, cfg["base"], nu1, 4.0, budget)
    # consistency: the requested, non-empty sample; a count; the Wilson interval of that count
    for lo_hi in ((0.0, 0.0), (0.0, 1.0)):
        wrong = copy.deepcopy(ck["ic.rlogy"]); wrong["3"]["ci_low"], wrong["3"]["ci_high"] = lo_hi
        assert "not the Wilson interval" in (m._checkpoint_identity_mismatch(wrong, cfg["ic.rlogy"], nu1, 4.0, budget) or ""), lo_hi
    empty = copy.deepcopy(ck["base"]); empty["3"]["n"] = empty["3"]["n_total"] = empty["3"]["n_censored"] = 0
    assert "requested sample" in m._checkpoint_identity_mismatch(empty, cfg["base"], nu1, 4.0, budget)
    notcount = copy.deepcopy(ck["base"]); notcount["3"]["fraction"] = 0.04585
    assert "not a count" in m._checkpoint_identity_mismatch(notcount, cfg["base"], nu1, 4.0, budget)
    # coherence with the first-passage count, the censored-inclusive identity and the inputs
    fp0 = copy.deepcopy(ck["base"]); fp0["3"]["fraction_first_passage"] = 0.0
    assert "exceeds the first-passage count" in m._checkpoint_identity_mismatch(fp0, cfg["base"], nu1, 4.0, budget)
    fphalf = copy.deepcopy(ck["base"]); fphalf["3"]["fraction_first_passage"] = 0.06175
    assert "not a count" in m._checkpoint_identity_mismatch(fphalf, cfg["base"], nu1, 4.0, budget)
    inc0 = copy.deepcopy(ck["ic.rlogy"]); inc0["3"]["fraction_first_passage_censored_inclusive"] = 0.0
    assert "censored-inclusive" in m._checkpoint_identity_mismatch(inc0, cfg["ic.rlogy"], nu1, 4.0, budget)
    halfcens = copy.deepcopy(ck["base"])
    for c in ("2", "3", "4"):
        halfcens[c].update(n=5000, n_censored=5000, fraction=0.5, fraction_first_passage=0.5,
                           fraction_first_passage_censored_inclusive=0.75)
    assert "requested sample" in m._checkpoint_identity_mismatch(halfcens, cfg["base"], nu1, 4.0, budget)
    regraft = copy.deepcopy(ck["ic.rlogy"]); cell = regraft["3"]
    cell.update(n=20000, n_censored=20000, fraction=5454 / 20000, fraction_first_passage=5454 / 20000,
                fraction_first_passage_censored_inclusive=(5454 + 20000) / 40000)
    cell["ci_low"], cell["ci_high"] = m.wilson(5454, 20000)
    assert "requested sample" in m._checkpoint_identity_mismatch(regraft, cfg["ic.rlogy"], nu1, 4.0, budget)


def test_r16_telescoping_keeps_an_exact_predator_zero():
    """An exact predator zero was reported as the safe floor, and its logarithm as
    the finite sentinel -1e300. Independent oracle: from (x, y) = (1/2, 1e-249) at mu0 = 4, nu1 = 1
    the float map gives x: 1/2 -> 1 -> 0 and y1 = 5e-250, y2 = 5e-250, y3 = 0 exactly."""
    import numpy as np
    dyn = _load("DYN")
    r = dyn._designG_telescoping_rate(np.array([[0.5, 1e-249]]), 4.0, 1.0, 3, 1e-250)
    max_rate, stop_y, _pred, max_stop, min_stop, n_stopped, log10_h, excess = r
    assert stop_y == 0.0 and log10_h == -np.inf and excess == -np.inf, r
    # the first sub-floor value comes at step 3 = N, the horizon itself: not "before the horizon"
    assert (max_stop, min_stop, n_stopped) == (3, 3, 0), r
    r5 = dyn._designG_telescoping_rate(np.array([[0.5, 1e-249]]), 4.0, 1.0, 5, 1e-250)
    assert r5[5] == 1 and r5[3] == 3, r5   # with N = 5 the same crossing is strictly before it
    # ln(y2/y0)/2 with y2 = y0/2, from two logarithms near -573: double rounding allows ~1e-13
    assert abs(max_rate - 0.5 * np.log(0.5)) < 1e-12, r


def test_r02_period_one_enumeration_passes_its_own_checks():
    """The period-two closed-form oracle ran even for max_period=1, so a complete
    period-one enumeration necessarily failed it. Period one: the single interior fixed point 3/4."""
    import math
    erg = _load_path("designG_ergodic_quadrature_v2.py")
    rec, rows = erg.enumerate_periodic_orbits(1)
    assert rec["interior_orbits_enumerated"] == 1 and rec["coverage_complete"] is True
    assert abs(rec["max_orbit_avg_lnx"] - math.log(0.75)) < 1e-12
    assert all(r.passed for r in rows), [r.name for r in rows if not r.passed]
    assert not any("period-2" in r.name for r in rows)
    # The closed forms drawn in Figure 8 (the map of the parameter box): the coexistence point and its Jacobian,
    # the curve on which it has a multiplier -1, and the invasion threshold on the prey's attracting fixed point
    # and two-cycle, checked against the multipliers and against orbit averages of ln x.
    import numpy as np
    assert abs(float(erg.nu1_flip(3.0)) - 1.5) < 1e-14
    assert abs(float(erg.nu1_flip(4.0)) - (8.0 * math.sqrt(3.0) - 12.0)) < 1e-13
    assert abs(float(erg.nu1_all_orbit(4.0)) - 4.0 / 3.0) < 1e-15 and float(erg.nu1_star_closed(4.0)) == 4.0
    rng = np.random.default_rng(8)
    n_checked = 0
    for _ in range(600):
        m, v = rng.uniform(4.0 / 3.0, 4.0), rng.uniform(4.0 / 3.0, 4.0)
        xs, ys = erg.coexistence_point(m, v)
        if ys <= 1e-6 or abs(v - float(erg.nu1_flip(m))) < 1e-6:
            continue
        J = erg.coexistence_jacobian(m, v)
        assert abs(m * xs * (1 - xs) * (1 - ys) - xs) < 1e-13 and abs(v * xs * ys * (1 - ys) - ys) < 1e-13
        assert abs(1.0 - np.trace(J) + np.linalg.det(J) - m * ys) < 1e-12
        assert erg.coexistence_is_sink(m, v) == (v > float(erg.nu1_flip(m)) and ys < 0.5), (m, v)
        n_checked += 1
    assert n_checked > 300
    for mu0 in (2.0, 2.9, 3.2, 3.4):
        x = 0.3
        for _ in range(20000):
            x = mu0 * x * (1.0 - x)
        s = 0.0
        for _ in range(20000):
            s += math.log(x); x = mu0 * x * (1.0 - x)
        assert abs(math.exp(-s / 20000) - float(erg.nu1_star_closed(mu0))) < 1e-6, mu0
        assert float(erg.nu1_star_closed(mu0)) >= float(erg.nu1_all_orbit(mu0)) - 1e-12


def _load_reproducer():
    repro = os.path.join(_CODES_DIR, "reproduce_designG_v1.py")
    if not os.path.exists(repro):
        import pytest
        pytest.skip("reproduce_designG_v1.py ships only in the bundle (codes/); run the suite from the bundle root")
    spec = importlib.util.spec_from_file_location("reproduce_designG_v1_under_test", repro)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def test_reproducer_field_rule_for_the_noise_floor_diagnostic():
    """2026-10-21: the observed contraction ratio of the stable-leaf
    illustration is a quotient of sup-norm changes near 1e-13 and moved by 1.4e-5 on another CPU.
    Only that field has a stated rule (absolute 1e-3, finite, 0 <= c < 1 on both sides); every
    other leaf keeps the global tolerance."""
    r = _load_reproducer()
    fn = "designG_stable_leaf_v1_values.json"
    assert r.leaf_ok(fn, "/nu1_3/contraction_max", 0.3526105969288764, 0.3525962463526822, 1e-6, 1e-9)[0]
    assert not r.leaf_ok(fn, "/nu1_3/contraction_max", 1.0, 0.3525962463526822, 1e-6, 1e-9)[0]
    assert not r.leaf_ok(fn, "/nu1_3/contraction_max", float("nan"), 0.3525962463526822, 1e-6, 1e-9)[0]
    assert not r.leaf_ok(fn, "/nu1_3/contraction_max_x", 0.3526105969288764, 0.3525962463526822, 1e-6, 1e-9)[0]
    assert not r.leaf_ok(fn, "/nu1_3/leaf_x0/0", 0.36681, 0.36680, 1e-6, 1e-9)[0]
    # Equal values outside 0 <= c < 1, or not finite, pass the global
    # tolerance (_close accepts nan == nan and inf == inf), so only the rule's own range check rejects
    # them; a '*' in the rule never crosses a '/'; and on the shipped file the rule covers exactly the
    # two contraction_max leaves.
    for v in (1.0, -0.1, float("nan"), float("inf")):
        assert r._close(v, v, 1e-6, 1e-9)
        assert not r.leaf_ok(fn, "/nu1_3/contraction_max", v, v, 1e-6, 1e-9)[0], v
    assert r.field_rule(fn, "/nu1_3/nested/contraction_max") is None
    assert r.field_rule("designG_bifurcation_values.json", "/nu1_3/contraction_max") is None
    leaves = {}
    r._walk_numeric(json.load(open(os.path.join(_CODES_DIR, "data", "designG", fn))), "", leaves)
    assert sorted(k for k in leaves if r.field_rule(fn, k)) == ["/nu1_2.1/contraction_max", "/nu1_3/contraction_max"]
    # The range is required on BOTH sides -- one value just outside it fails the
    # rule although the pair passes the global tolerance.
    for a, b in ((0.9999999999, 1.0), (1.0, 0.9999999999), (0.0, -1e-10), (-1e-10, 0.0)):
        assert r._close(a, b, 1e-6, 1e-9) and not r.leaf_ok(fn, "/nu1_3/contraction_max", a, b, 1e-6, 1e-9)[0], (a, b)


def test_reproducer_overall_verdict_fails_when_any_producer_fails(tmp_path):
    """A past defect: the per-leaf result was assigned to the aggregate
    flag, so a matching producer later in PRODUCERS erased an earlier mismatch, missing file or
    field-bound failure, and the reproducer printed OVERALL: PASS and exited 0. The stable-leaf file
    is not the last producer, so the stable-leaf mismatch itself was masked. Each case must now fail the
    numeric section and the run; the value measured on an AVX-512 CPU must still pass."""
    import shutil
    repro = os.path.join(_CODES_DIR, "reproduce_designG_v1.py")
    if not os.path.exists(repro):
        import pytest
        pytest.skip("reproduce_designG_v1.py ships only in the bundle (codes/); run the suite from the bundle root")
    _shipped, data, _ = _val_copies(tmp_path)
    leaf, onset = "designG_stable_leaf_v1_values.json", "designG_onset_ensemble_v1.json"

    def dirs(case):
        regen, ship = tmp_path / case / "regen", tmp_path / case / "shipped"
        shutil.copytree(data, regen); shutil.copytree(data, ship)
        return regen, ship

    def edit(d, change):
        p = d / leaf; doc = json.loads(p.read_text()); change(doc["nu1_3"]); p.write_text(json.dumps(doc))

    def run(regen, ship):
        r = subprocess.run([_PY, repro, "--regen", str(regen), "--shipped", str(ship)], capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    regen, ship = dirs("near_miss")              # global tolerance, beside the field rule
    edit(regen, lambda b: b["leaf_x0"].__setitem__(0, b["leaf_x0"][0] + 1e-5))
    regen2, ship2 = dirs("field_bound")          # equal on both sides: only the range check rejects it
    for d in (regen2, ship2):
        edit(d, lambda b: b.__setitem__("contraction_max", 1.0))
    regen3, ship3 = dirs("missing")
    (regen3 / onset).unlink()
    for case, (rg, sh) in {"near_miss": (regen, ship), "field_bound": (regen2, ship2), "missing": (regen3, ship3)}.items():
        rc, out = run(rg, sh)
        assert rc == 1 and "OVERALL: FAIL" in out, (case, out[-1500:])
        assert "(A) numeric-leaf reproduction (regenerated vs shipped): FAIL" in out, (case, out[-1500:])
    regen4, ship4 = dirs("avx512_cpu")
    edit(regen4, lambda b: b.__setitem__("contraction_max", 0.3526105969288764))
    rc, out = run(regen4, ship4)
    assert rc == 0 and "OVERALL: PASS" in out and "/nu1_3/contraction_max: regen=0.3526105969288764" in out, out[-1500:]
    # 2026-11-10: nothing checked is not a pass. On identical data, part (B) passed with a claims file that
    # had no "claims" key, no entries, or only its first ten ("resolved=0 matched=0"); each must fail (B) and the run.
    # The PASS above already holds the frozen file to its 77 entries, 65 of them resolved.
    frozen = json.load(open(os.path.join(_CODES_DIR, "submitted_claims_designG_v1.json")))
    cp = tmp_path / "claims_case.json"
    def run_claims(doc, *extra):
        cp.write_text(json.dumps(doc))
        r = subprocess.run([_PY, repro, "--regen", str(data), "--shipped", str(data), "--claims", str(cp), *extra],
                           capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr
    for case, doc in {"no claims key": {}, "no entries": dict(frozen, claims=[]),
                      "first ten entries": dict(frozen, claims=frozen["claims"][:10])}.items():
        rc, out = run_claims(doc)
        assert rc == 1 and "(B) submitted_claims direct resolution: FAIL" in out and "OVERALL: FAIL" in out, (case, out[-800:])
    rc, out = run_claims(dict(frozen, claims=frozen["claims"][:10]), "--expect-claims", "10", "9")
    assert rc == 0 and "OVERALL: PASS" in out, out[-800:]     # another inventory states its own counts
