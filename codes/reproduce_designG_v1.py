#!/usr/bin/env python3
"""reproduce_designG_v1.py -- check that the Design-G producer numbers reproduce.

Reproducibility check for the Design-G numbers.  It takes a
directory of FRESHLY REGENERATED producer JSON (written by run_designG_figs.sh and
designG_quasiperiodic_v2.py: the nine result files of the eight producers listed in PRODUCERS below) and the
SHIPPED reference producer JSON. Four checks decide the exit status -- (A) the numeric leaves and (A') the
boolean leaves outside the volatile keys (VOLATILE_KEYS below: provenance blocks, times, host and path
fields, run times), (B) the frozen entries that resolve directly, (D) the essential success flags -- and a
fifth, (C), the comparison of the random initial conditions, is reported and never changes it:

  (A) numeric-leaf reproduction -- every numeric leaf of each regenerated
      producer JSON matches the shipped reference within tolerance (volatile
      metadata subtrees: provenance / timestamps / figure paths are ignored).
      Leaf-level agreement establishes that every number these producers emit
      reproduces; it does not verify derivations made by hand, the 12 frozen
      entries without a machine-resolvable pointer, or the manuscript text itself
      (agreement between text and producer JSON is the job of validate_designG_v1.py).

  (B) submitted-claims direct resolution -- for every submitted_claims entry
      whose json_source_field is a plain '<file>.json:/a/b/c' pointer (or a bare
      '/a/b/c' against the dynamics producer), resolve it against the REGENERATED
      JSON and confirm it equals the recorded json_value.  Entries with
      human-annotated pointers ( '{a,b}' sets, parentheticals, 'derived' ) are
      counted and reported; their recorded values are not checked (layer (A)
      compares the producer fields they derive from, not the derived values).  The file
      must hold the expected number of entries and of resolved entries (77 and
      65 for the frozen inventory; --expect-claims for another one), so that an
      empty or truncated claims file fails instead of passing with nothing checked.

  (C) initial-condition provenance (2026-10-18) -- for each random sensitivity
      configuration, the SHA-256 of the initial conditions the regenerated run used against
      that of the shipped run: STRICT when both records carry a valid digest and the two are
      equal (the archive loaded, or a fresh draw with the same bytes), PLATFORM when both are
      valid and differ (a fresh draw or a resample under another CPU dispatch), UNAVAILABLE
      when a record or a valid digest is missing on a side; the record's ``source`` field
      says how the inputs were obtained. Informational, never changes the exit code, but it
      tells the reader which kind of comparison (A) performed.

  (A') boolean-leaf reproduction and (D) essential success flags -- every
      boolean leaf outside the volatile subtrees (pass flags, Jury stability, chaotic labels,
      certificate checks) must be IDENTICAL in the regenerated and the shipped JSON, and on
      BOTH sides the result records must report success: the period-two certificate certified
      with every check ok, every exact-skew collapse guard ok, the periodic-orbit enumeration
      complete with zero least-period count mismatches. Numeric agreement alone used to pass
      records stating that certification or a guard had failed.

  Field rule (2026-10-21). Exactly one producer field has its own
      comparison rule: designG_stable_leaf_v1_values.json:/nu1_*/contraction_max, the largest
      OBSERVED ratio of successive sup-norm changes of the Lyapunov--Perron iteration (an
      illustration diagnostic quoted nowhere in the manuscript; not a contraction constant).
      At nu1 = 3 its maximizing sweep divides changes of about 1.4e-13, so rounding noise in
      the inputs moves it: perturbing each height by +-1..3 ulp moves it by up to 3.18e-5
      (by up to 5.6e-16 at nu1 = 2.1), and a one-ulp height difference on an AVX-512 CPU moved
      it by 1.435e-5 (40.6 times the global tolerance). The heights are now frozen literals, and
      this field is compared with an absolute tolerance of 1e-3 (about 30 times that noise
      floor) together with a strict check on BOTH sides that it is finite and 0 <= c < 1 --
      the property it records. Every other leaf keeps the global rtol/atol; the gate prints a
      [FIELD-TOL] line for each ruled leaf.

Exit 0 iff (A), (A'), (B) and (D) all pass.  Standard library only.

Usage:
  reproduce_designG_v1.py --regen <dir> [--shipped <dir>] [--claims <json>]
                          [--rtol 1e-6] [--atol 1e-9]
  Environment fallbacks: DESIGNG_REGEN_DIR, DESIGNG_SHIPPED_DIR,
  DESIGNG_SUBMITTED_CLAIMS.
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import math
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

PRODUCERS = (
    # v4 supersedes v3 (terminal persistence); v3's values JSON is archived under
    # data/designG/superseded/ and is deliberately NOT a producer any more.
    "designG_dynamics_v4_values.json",
    "designG_quasiperiodic_v2_values.json",
    "designG_ergodic_quadrature_v2_values.json",
    # added 2026-10-17: every producer of a manuscript figure has a shipped
    # reference, so the regression in (A) covers the onset ensemble, the closed-form invasion
    # exponents and the mu0-sweep onset values as well.
    "designG_onset_ensemble_v1.json",
    "designG_invasion_exponents_v1_values.json",
    # added 2026-10-21: the stable-leaf / matched-product illustration (deterministic, NumPy only).
    "designG_stable_leaf_v1_values.json",
    "designG_bifurcation_values.json",
    "designG_extinction_bif_values.json",
    # added 2026-10-18: the exact two-cycle certificate and its
    # numerical observations (basin fraction at nu1 = 46/25, nu1 scan).
    "designG_period_two_certificate_v1.json",
)

# (entries, entries resolved by a machine-readable pointer) of the frozen submitted_claims_designG_v1.json.
# Part (B) requires these counts: a claims file that was empty, truncated or without its "claims" key used to
# pass with "resolved=0 matched=0".
EXPECTED_CLAIM_COUNTS = (77, 65)

# Keys whose values are environment-dependent (wall-clock timings, timestamps,
# host, absolute figure paths, argv) and therefore excluded from numeric
# reproduction.  Everything else (grid sizes, node counts, cost budgets, and all
# scientific values) is deterministic and IS compared.
VOLATILE_KEYS = {
    "provenance", "_provenance", "figures", "utc", "generated", "generated_utc",
    "timestamp", "timestamp_utc", "host", "hostname", "output_data_dir",
    "output_fig_dir", "python", "platform", "argv", "cwd", "wall_s",
    # The initial-condition provenance block records the
    # CPU-dispatch state and the fresh draw's ULP agreement with the archive -- host
    # facts, compared separately in section (C) below, never as scientific leaves.
    "ic_provenance", "written_utc",
}


def ic_provenance_check(regen_dir: str, shipped_dir: str) -> Tuple[str, List[str]]:
    """(C) Initial-condition provenance of the random sensitivity configurations.

    Compares, for every random config of the dynamics producer, the checksum of the
    initial conditions the REGENERATED run used with the checksum the SHIPPED run used,
    and reports what the regenerated run did (archive loaded / fresh draw / resampled)
    together with its NumPy CPU-dispatch state. Three outcomes, kept apart:
      STRICT        -- identical inputs: any numeric difference in those blocks is a
                       computation difference, not a sampling one;
      PLATFORM      -- the run resampled or found no archive and its fresh draw differs
                       from the shipped inputs (typically one ULP in np.power under a
                       different CPU dispatch); differences in the random blocks are then a
                       finite-sample / platform effect to be judged against the Wilson
                       intervals;
      UNAVAILABLE   -- one side carries no provenance record (pre-2026-10-18 files).
    This section is informational: it never changes the exit code, it classifies.
    """
    fn = "designG_dynamics_v4_values.json"
    rp, sp = os.path.join(regen_dir, fn), os.path.join(shipped_dir, fn)
    if not (os.path.exists(rp) and os.path.exists(sp)):
        return "UNAVAILABLE", ["  dynamics producer JSON missing on one side"]
    r = json.load(open(rp)).get("sensitivity", {}); s_ = json.load(open(sp)).get("sensitivity", {})
    msgs: List[str] = []
    verdicts: List[str] = []
    for cfg in sorted(k for k in r if not str(k).startswith("_")):
        pr = ((r.get(cfg) or {}).get("_params") or {}).get("ic_provenance")
        ps = ((s_.get(cfg) or {}).get("_params") or {}).get("ic_provenance")
        if pr is None and ps is None:
            continue   # deterministic lattice config: no random inputs
        if pr is None or ps is None:
            verdicts.append("UNAVAILABLE")
            msgs.append(f"  [UNAVAILABLE] {cfg}: provenance missing on {'the regenerated' if pr is None else 'the shipped'} side")
            continue
        hr, hs = pr.get("used_sha256"), ps.get("used_sha256")
        import re as _re
        _hex64 = _re.compile(r"[0-9a-f]{64}")
        if not (isinstance(hr, str) and _hex64.fullmatch(hr) and isinstance(hs, str) and _hex64.fullmatch(hs)):
            # Two missing or malformed checksums (None, or a non-hex string) must not compare
            # equal and read as STRICT; only two valid 64-character lowercase hex digests are compared
            verdicts.append("UNAVAILABLE")
            msgs.append(f"  [UNAVAILABLE] {cfg}: a provenance block lacks a valid sha256 (regenerated {hr!r}, shipped {hs!r})")
            continue
        same = hr == hs
        rt = pr.get("runtime") or {}
        feats = rt.get("cpu_features_enabled")
        feats_s = ",".join(feats) if isinstance(feats, list) else str(feats)
        verdict = "STRICT" if same else "PLATFORM"
        verdicts.append(verdict)
        nd = pr.get("n_values_differ")
        diff_s = "" if nd is None else f", {nd} value(s) differ, max {pr.get('max_ulp')} ULP"
        agree_s = "==" if pr.get("identical_to_archive") else "!="
        dis = rt.get("NPY_DISABLE_CPU_FEATURES")
        dis_s = "" if not dis else f" (NPY_DISABLE_CPU_FEATURES={dis})"
        msgs.append(f"  [{verdict}] {cfg}: regenerated inputs {str(pr.get('used_sha256'))[:16]}... "
                    f"({pr.get('source')}; fresh draw {agree_s} archive{diff_s}) "
                    f"vs shipped {str(ps.get('used_sha256'))[:16]}... ({ps.get('source')}); "
                    f"numpy {rt.get('numpy')} on {rt.get('machine')} with {feats_s}{dis_s}")
    if not verdicts:
        return "UNAVAILABLE", ["  no random-config provenance on either side"]
    overall = "PLATFORM" if "PLATFORM" in verdicts else ("UNAVAILABLE" if "UNAVAILABLE" in verdicts else "STRICT")
    return overall, msgs


def _num(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _checkable(x: Any) -> bool:
    """A scalar this gate can verify: a real number OR a bool.

    Bools were previously excluded (``_num`` rejects them), so a claim whose json_value is
    True/False fell into the `annotated` bucket and was never checked -- e.g.
    ``nu1star_3p74_nonchaotic`` (/I_family/family/3.74/chaotic = False). Its correctness then
    rested on a prose note in the claims _meta, which no tool enforces. A certification only a
    human can see is not a gate, so bools are now first-class and compared exactly.
    """
    return isinstance(x, bool) or _num(x)


# Field-specific comparison rules: {producer file: ((JSON-pointer pattern, absolute tolerance,
# (lo, hi) with lo <= value < hi required on both sides, reason), ...)}. Matched with fnmatchcase
# component by component over the WHOLE pointer (a '*' never crosses a '/'), so near-miss and
# nested paths keep the global tolerance.
FIELD_RULES: Dict[str, Tuple[Tuple[str, float, Tuple[float, float], str], ...]] = {
    "designG_stable_leaf_v1_values.json": (
        ("/nu1_*/contraction_max", 1e-3, (0.0, 1.0),
         "observed finite-iteration ratio near the floating-point floor (noise <= 3.2e-5 under 1-3 ulp "
         "input perturbations); compared at 1e-3 with 0 <= c < 1 enforced on both sides"),
    ),
}


def field_rule(fn: str, path: str) -> Optional[Tuple[float, Tuple[float, float], str]]:
    """The field-specific rule for leaf `path` of producer file `fn`, or None. The pattern must
    have as many '/'-separated components as the pointer, each matched by fnmatchcase."""
    parts = path.split("/")
    for pattern, tol, bounds, why in FIELD_RULES.get(fn, ()):
        pparts = pattern.split("/")
        if len(pparts) == len(parts) and all(fnmatch.fnmatchcase(a, b) for a, b in zip(parts, pparts)):
            return tol, bounds, why
    return None


def leaf_ok(fn: str, path: str, regen: float, shipped: float, rtol: float, atol: float) -> Tuple[bool, str]:
    """Compare one numeric leaf: the field rule if one matches, else the global tolerance."""
    rule = field_rule(fn, path)
    if rule is None:
        return _close(regen, shipped, rtol, atol), ""
    tol, (lo, hi), _why = rule
    vals_ok = all(math.isfinite(v) and lo <= v < hi for v in (regen, shipped))
    diff = abs(regen - shipped) if vals_ok else float("nan")
    ok = vals_ok and diff <= tol
    note = (f"[FIELD-TOL] {fn}:{path}: regen={regen!r} shipped={shipped!r} |diff|={diff:.6g} "
            f"(field rule: abs <= {tol:g}, {lo:g} <= value < {hi:g} on both sides) -> {'OK' if ok else 'FAIL'}")
    return ok, note


def _close(a: float, b: float, rtol: float, atol: float) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= atol + rtol * abs(b)


def _walk_numeric(obj: Any, path: str, out: Dict[str, float]) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in VOLATILE_KEYS:
                continue
            _walk_numeric(v, path + "/" + str(k), out)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _walk_numeric(v, path + "/" + str(i), out)
    elif _num(obj):
        out[path] = float(obj)


def _walk_bool(obj: Any, path: str, out: Dict[str, bool]) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in VOLATILE_KEYS:
                continue
            _walk_bool(v, path + "/" + str(k), out)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _walk_bool(v, path + "/" + str(i), out)
    elif isinstance(obj, bool):
        out[path] = obj


def boolean_leaf_check(regen_dir: str, shipped_dir: str) -> Tuple[bool, List[str], int]:
    """(A') every non-volatile boolean leaf compared EXACTLY: `_num` excludes
    bools, so a flipped pass flag or certificate check used to be invisible to (A)."""
    msgs: List[str] = []
    ok = True
    total = 0
    for fn in PRODUCERS:
        rp, sp = os.path.join(regen_dir, fn), os.path.join(shipped_dir, fn)
        if not (os.path.exists(rp) and os.path.exists(sp)):
            continue   # a missing file already fails (A)
        rb: Dict[str, bool] = {}
        sb: Dict[str, bool] = {}
        _walk_bool(json.load(open(rp)), "", rb)
        _walk_bool(json.load(open(sp)), "", sb)
        common = set(rb) & set(sb)
        only_r, only_s = set(rb) - set(sb), set(sb) - set(rb)
        mism = [k for k in sorted(common) if rb[k] is not sb[k]]
        total += len(common)
        if mism or only_r or only_s:
            ok = False
            msgs.append(f"  [FAIL] {fn}: {len(common)} common boolean leaves, {len(mism)} differ, "
                        f"{len(only_r)} regen-only, {len(only_s)} shipped-only")
            for k in mism[:8]:
                msgs.append(f"      DIFFER {k}: regen={rb[k]!r} shipped={sb[k]!r}")
            for k in sorted(only_r)[:4]:
                msgs.append(f"      REGEN-ONLY {k}")
            for k in sorted(only_s)[:4]:
                msgs.append(f"      SHIPPED-ONLY {k}")
        elif common:
            msgs.append(f"  [OK] {fn}: {len(common)} boolean leaves identical")
    return ok, msgs, total


# (file, pointer with `*` over dict keys / list items, requirement) -- the success flags a result
# record must carry for its numbers to count as results at all.
ESSENTIAL_FLAGS = (
    ("designG_period_two_certificate_v1.json", "/certificate/certified", "True"),
    ("designG_period_two_certificate_v1.json", "/certificate/checks/*/ok", "True"),
    ("designG_dynamics_v4_values.json", "/exact_skew/per_nu1/*/collapse_guard_ok", "True"),
    ("designG_ergodic_quadrature_v2_values.json", "/periodic_enumeration/coverage_complete", "True"),
    ("designG_ergodic_quadrature_v2_values.json", "/periodic_enumeration/least_period_count_mismatches", "0"),
)
# The members those flags must cover: a deleted member or a deleted flag must fail, not shrink
# the set inspected -- the nine named checks of the two-cycle certificate and the three
# exact-skew parameters.
EXPECTED_CERT_CHECKS = (
    "C invertible (exact)",
    "contraction: max row sum of |E| = |I - C(D(Phi^2) - I)| < 1 on B",
    "mean-value enclosure T(B) subset K subset int B",
    "B and Phi(B) interior to the open unit square",
    "least period two: first coordinates of Phi(B) and B disjoint",
    "strict Jury condition 1-d > 0 on the enclosure of D(Phi^2)(p)",
    "strict Jury condition 1+d > 0 on the enclosure of D(Phi^2)(p)",
    "strict Jury condition 1-t+d > 0 on the enclosure of D(Phi^2)(p)",
    "strict Jury condition 1+t+d > 0 on the enclosure of D(Phi^2)(p)",
)
EXPECTED_SKEW_NU1 = ("2", "3", "4")
_MISSING = object()   # a member reached through `*` that lacks the requested field


def _flag_values(node: Any, toks: List[str], path: str, out: List[Tuple[str, Any]], wild: bool = False) -> None:
    if not toks:
        out.append((path, node))
        return
    tok, rest = toks[0], toks[1:]
    if tok == "*":
        if isinstance(node, dict):
            for k, v in node.items():
                if not str(k).startswith("_"):
                    _flag_values(v, rest, f"{path}/{k}", out, True)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                _flag_values(v, rest, f"{path}/{i}", out, True)
    elif isinstance(node, dict) and tok in node:
        _flag_values(node[tok], rest, f"{path}/{tok}", out, wild)
    elif wild:
        out.append((f"{path}/{tok}", _MISSING))   # a member without the flag fails, never skipped


def essential_flags_check(data_dir: str, label: str) -> Tuple[bool, List[str]]:
    """(D) every ESSENTIAL_FLAGS entry present and satisfied in one data directory: a flag must
    be the JSON bool true, a mismatch count the integer 0."""
    msgs: List[str] = []
    ok = True
    for fn, pointer, need in ESSENTIAL_FLAGS:
        fp = os.path.join(data_dir, fn)
        if not os.path.exists(fp):
            ok = False
            msgs.append(f"  [FAIL] {label} {fn}: missing (needs {pointer} {need})")
            continue
        hits: List[Tuple[str, Any]] = []
        _flag_values(json.load(open(fp)), [t for t in pointer.split("/") if t], "", hits)
        if not hits:
            ok = False
            msgs.append(f"  [FAIL] {label} {fn}{pointer}: absent (must be {need})")
            continue
        bad = [(pth, v) for pth, v in hits
               if not ((v is True) if need == "True" else (isinstance(v, int) and not isinstance(v, bool) and v == 0))]
        if bad:
            ok = False
            for pth, v in bad[:6]:
                shown = "<missing>" if v is _MISSING else repr(v)
                msgs.append(f"  [FAIL] {label} {fn}{pth} = {shown}: required {need}")
    # The members themselves.
    cp = os.path.join(data_dir, "designG_period_two_certificate_v1.json")
    if os.path.exists(cp):
        checks = (json.load(open(cp)).get("certificate") or {}).get("checks")
        names = [c.get("check") if isinstance(c, dict) else None for c in checks] if isinstance(checks, list) else None
        if names is None or len(names) != len(EXPECTED_CERT_CHECKS) or sorted(map(str, names)) != sorted(EXPECTED_CERT_CHECKS):
            ok = False
            msgs.append(f"  [FAIL] {label} certificate checks {names!r} are not the nine expected checks")
    dp = os.path.join(data_dir, "designG_dynamics_v4_values.json")
    if os.path.exists(dp):
        per = (json.load(open(dp)).get("exact_skew") or {}).get("per_nu1")
        keys = sorted(k for k in per if not str(k).startswith("_")) if isinstance(per, dict) else None
        if keys != sorted(EXPECTED_SKEW_NU1):
            ok = False
            msgs.append(f"  [FAIL] {label} exact_skew parameters {keys!r}, expected {list(EXPECTED_SKEW_NU1)}")
    if ok:
        msgs.append(f"  [OK] {label}: {len(ESSENTIAL_FLAGS)} success-flag requirements satisfied")
    return ok, msgs


def numeric_leaf_check(regen_dir: str, shipped_dir: str, rtol: float, atol: float
                       ) -> Tuple[bool, List[str], int]:
    msgs: List[str] = []
    ok = True
    total = 0
    for fn in PRODUCERS:
        rp, sp = os.path.join(regen_dir, fn), os.path.join(shipped_dir, fn)
        if not os.path.exists(rp):
            ok = False
            msgs.append(f"  MISSING regenerated {fn}")
            continue
        if not os.path.exists(sp):
            ok = False
            msgs.append(f"  MISSING shipped {fn}")
            continue
        rleaves: Dict[str, float] = {}
        sleaves: Dict[str, float] = {}
        _walk_numeric(json.load(open(rp)), "", rleaves)
        _walk_numeric(json.load(open(sp)), "", sleaves)
        common = set(rleaves) & set(sleaves)
        only_r = set(rleaves) - set(sleaves)
        only_s = set(sleaves) - set(rleaves)
        mism = []
        notes = []
        for k in sorted(common):
            leaf_matches, note = leaf_ok(fn, k, rleaves[k], sleaves[k], rtol, atol)
            if note:
                notes.append("      " + note)
            if not leaf_matches:   # never assign the aggregate `ok` here: a later match would erase a failure
                mism.append(k)
        total += len(common)
        status = "OK" if (not mism and not only_r and not only_s) else "FAIL"
        if status != "OK":
            ok = False
        msgs.append(f"  [{status}] {fn}: {len(common)} common numeric leaves, "
                    f"{len(mism)} mismatched, {len(only_r)} regen-only, {len(only_s)} shipped-only")
        msgs.extend(notes)
        for k in mism[:8]:
            msgs.append(f"      MISMATCH {k}: regen={rleaves[k]!r} shipped={sleaves[k]!r}")
        for k in list(only_r)[:4]:
            msgs.append(f"      REGEN-ONLY {k}")
        for k in list(only_s)[:4]:
            msgs.append(f"      SHIPPED-ONLY {k}")
    return ok, msgs, total


def _resolve_pointer(doc: Any, pointer: str) -> Tuple[bool, Any]:
    """Resolve a '/a/b/c' pointer where each token is a dict key or list index.
    Keys may themselves contain '/', so try greedy multi-token joins too.

    Also supports two forms that previously made a claim silently UNRESOLVABLE (and so,
    previously, silently unchecked):
      ``*``      -- every non-underscore child; resolves only if they ALL AGREE, which is
                    exactly what a claim like /exact_skew/per_nu1/*/abs_err asserts (one
                    value across the sweep). Disagreement is a real finding, not a match.
      ``key[i]`` -- index into a list-valued key, e.g. horizon_sweep[0].
    """
    import re as _re
    # expand "key[i]" -> "key/i" so the existing list branch handles it
    pointer = _re.sub(r"\[(\d+)\]", r"/\1", pointer)
    tokens = [t for t in pointer.split("/") if t != ""]

    def rec(node: Any, toks: List[str]) -> Tuple[bool, Any]:
        if not toks:
            return True, node
        if toks[0] == "*" and isinstance(node, dict):
            vals = []
            for k, v in node.items():
                if str(k).startswith("_"):
                    continue
                ok, val = rec(v, toks[1:])
                if ok:
                    vals.append(val)
            if not vals:
                return False, None
            first = vals[0]
            for v in vals[1:]:
                if v != first:
                    return False, None      # children disagree -> NOT a single value
            return True, first
        # try consuming 1..len tokens as a single (possibly slash-containing) key
        if isinstance(node, dict):
            for take in range(1, len(toks) + 1):
                key = "/".join(toks[:take])
                if key in node:
                    ok, val = rec(node[key], toks[take:])
                    if ok:
                        return True, val
            return False, None
        if isinstance(node, list):
            try:
                idx = int(toks[0])
            except ValueError:
                return False, None
            if 0 <= idx < len(node):
                return rec(node[idx], toks[1:])
            return False, None
        return False, None

    return rec(doc, tokens)


def submitted_claims_check(claims_path: str, regen_dir: str, rtol: float, atol: float,
                           expect: Optional[Tuple[int, int]] = EXPECTED_CLAIM_COUNTS
                           ) -> Tuple[bool, List[str], Dict[str, int]]:
    msgs: List[str] = []
    docs: Dict[str, Any] = {}
    for fn in PRODUCERS:
        p = os.path.join(regen_dir, fn)
        if os.path.exists(p):
            docs[fn] = json.load(open(p))
    claims_doc = json.load(open(claims_path))
    claims = claims_doc.get("claims") if isinstance(claims_doc, dict) else None
    claims = [c for c in claims if isinstance(c, dict)] if isinstance(claims, list) else []
    stats = {"resolved": 0, "matched": 0, "unresolvable": 0, "annotated": 0}
    ok = True
    for c in claims:
        jv = c.get("json_value")
        jsf = c.get("json_source_field")
        cid = c.get("claim_id", "?")
        if not (_checkable(jv) and isinstance(jsf, str)):
            stats["annotated"] += 1
            continue
        if any(ch in jsf for ch in "{}") or "(" in jsf or jsf.startswith("derived"):
            stats["annotated"] += 1
            continue
        if ":" in jsf:
            fname, ptr = jsf.split(":", 1)
            fname = fname.strip()
        else:
            fname, ptr = "designG_dynamics_v4_values.json", jsf  # bare pointer -> dynamics
        # strip any trailing " (...)" annotation already excluded above; strip spaces
        ptr = ptr.strip()
        doc = docs.get(os.path.basename(fname))
        if doc is None:
            stats["unresolvable"] += 1
            msgs.append(f"  UNRESOLVABLE {cid}: producer {fname} not regenerated")
            ok = False
            continue
        found, val = _resolve_pointer(doc, ptr)
        if not found or not _checkable(val):
            stats["unresolvable"] += 1
            msgs.append(f"  UNRESOLVABLE {cid}: pointer {ptr!r} in {fname}")
            ok = False   # an unresolved claim is an UNVERIFIED claim.
            continue
        stats["resolved"] += 1
        if isinstance(jv, bool) or isinstance(val, bool):
            # exact identity, and the TYPES must agree: float(True)==1.0, so routing a bool
            # through _close would let True silently "match" a producer 1 (or False a 0).
            hit = isinstance(jv, bool) and isinstance(val, bool) and (jv is val)
        else:
            hit = _close(float(val), float(jv), rtol, atol)
        if hit:
            stats["matched"] += 1
        else:
            ok = False
            msgs.append(f"  CLAIM-MISMATCH {cid}: regen={val!r} vs submitted_claims json_value={jv!r} ({jsf})")
    # Nothing checked is not a pass: the file must hold the expected entries, and as many of them must resolve.
    if expect is not None and (len(claims), stats["resolved"]) != tuple(expect):
        ok = False
        msgs.append(f"  COUNT-MISMATCH: {len(claims)} entries, {stats['resolved']} resolved; "
                    f"expected {expect[0]} entries, {expect[1]} resolved")
    return ok, msgs, stats


def main(argv: Optional[List[str]] = None) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Check that the Design-G producer numbers reproduce.")
    ap.add_argument("--regen", default=os.environ.get("DESIGNG_REGEN_DIR"),
                    help="dir of freshly regenerated producer *_values.json")
    ap.add_argument("--shipped", default=os.environ.get("DESIGNG_SHIPPED_DIR",
                    os.path.join(here, "data", "designG")),
                    help="dir of shipped reference producer *_values.json")
    ap.add_argument("--claims", default=os.environ.get("DESIGNG_SUBMITTED_CLAIMS",
                    os.path.join(here, "submitted_claims_designG_v1.json")),
                    help="submitted_claims_designG_v1.json")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--expect-claims", type=int, nargs=2, metavar=("ENTRIES", "RESOLVED"),
                    default=list(EXPECTED_CLAIM_COUNTS),
                    help="entries and machine-resolved entries the claims file must hold (default: the frozen inventory)")
    args = ap.parse_args(argv)

    if not args.regen:
        sys.stderr.write("[reproduce] ERROR: --regen (or DESIGNG_REGEN_DIR) is required\n")
        return 2

    print("== Design-G reproduction check ==")
    print(f"regen:   {args.regen}")
    print(f"shipped: {args.shipped}")
    print(f"claims:  {args.claims}")
    print(f"tol:     rtol={args.rtol} atol={args.atol}\n")

    a_ok, a_msgs, a_total = numeric_leaf_check(args.regen, args.shipped, args.rtol, args.atol)
    print(f"(A) numeric-leaf reproduction (regenerated vs shipped): "
          f"{'PASS' if a_ok else 'FAIL'}  [{a_total} numeric leaves compared]")
    for m in a_msgs:
        print(m)
    a2_ok, a2_msgs, a2_total = boolean_leaf_check(args.regen, args.shipped)
    print(f"\n(A') boolean-leaf reproduction, exact (regenerated vs shipped): "
          f"{'PASS' if a2_ok else 'FAIL'}  [{a2_total} boolean leaves compared]")
    for m in a2_msgs:
        print(m)

    b_ok, b_msgs, b_stats = submitted_claims_check(args.claims, args.regen, args.rtol, args.atol,
                                                   expect=tuple(args.expect_claims))
    print(f"\n(B) submitted_claims direct resolution: {'PASS' if b_ok else 'FAIL'}")
    print(f"    resolved={b_stats['resolved']} matched={b_stats['matched']} "
          f"unresolvable={b_stats['unresolvable']} annotated(not value-checked)={b_stats['annotated']}")
    for m in b_msgs:
        print(m)

    c_verdict, c_msgs = ic_provenance_check(args.regen, args.shipped)
    print(f"\n(C) initial-condition provenance of the random sensitivity configs: {c_verdict}")
    for m in c_msgs:
        print(m)
    if c_verdict == "PLATFORM":
        print("    (inputs differ at the floating-point level: differences confined to the random-IC blocks are a "
              "platform/sampling effect, to be read against the Wilson intervals)")

    d_ok_r, d_msgs_r = essential_flags_check(args.regen, "regenerated")
    d_ok_s, d_msgs_s = essential_flags_check(args.shipped, "shipped")
    d_ok = d_ok_r and d_ok_s
    print(f"\n(D) essential success flags (certificate, collapse guards, enumeration coverage): {'PASS' if d_ok else 'FAIL'}")
    for m in d_msgs_r + d_msgs_s:
        print(m)

    overall = a_ok and a2_ok and d_ok and b_ok and b_stats["resolved"] == b_stats["matched"]
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}  [inputs: {c_verdict}]")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
