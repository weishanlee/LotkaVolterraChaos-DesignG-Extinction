#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Machine-checkable Design-G manuscript-vs-code validator.

Maps every quantitative Design-G number the CURRENT manuscript quotes (the claim
inventory ``current_claims_designG_v1.json``, DESIGNG_CLAIMS_IN) to a code
producer field, resolves each to EXACTLY ONE scalar via a selector, and runs TWO
checks from that one scalar:

  (1) producer-kind check  -- BLOCKING for regression_anchor / theory_invariant
      (the new scalar must match ``anchor_expected`` within ``anchor_tol``);
      empirical_discovery is NON-blocking.
  (2) manuscript-consistency check -- if ``manuscript_tracked``, compare
      ``quoted_value`` to the new scalar at ``display_tol``; a mismatch sets
      ``manuscript_change_required`` (NON-blocking, recorded with the tex file:line).

Modes (``DESIGNG_VALIDATE_MODE``):
  * ``smoke`` -- schema / taxonomy / selector-resolvability / coverage-of-present
    only; every numeric anchor -> ``SKIPPED_SMOKE`` (never blocks). Used by F1's
    clone gate and by this task's own Step-5 run.
  * ``full`` (default; B7 promoted run) -- blocking kinds must PASS; a
    completeness gate (all U1-U5 blocks + fields + _taxonomy + C0-C10 count +
    U2/U3/U4/U5 figures + manifest rows + referenced retained-v2 blocks) and a
    fallback gate (any sensitivity config with ``_fallback_used=true`` HARD-blocks)
    apply. ``DESIGNG_VALIDATE_SCOPE`` (e.g. ``U5`` or ``exact_skew``) restricts the
    full-mode completeness + checks to the named block(s).

Env contract (ALL SIX are REQUIRED -- the validator ABORTS if any is unset):
  DESIGNG_VALIDATE_IN     -- dir of regenerated producer ``*_{smoke,values}.json``.
  DESIGNG_VALIDATE_FIGDIR -- dir of regenerated figures (figure-completeness gate).
  DESIGNG_CLAIMS_IN       -- the claim inventory JSON ($CLAIMS).
  DESIGNG_REPORT_OUT      -- human-readable report ($REPORT) to write.
  DESIGNG_VALJSON_OUT     -- machine-readable result JSON ($VALJSON) to write.
  DESIGNG_MANIFEST        -- phase manifest to VERIFY (rows authored by the producers).

Run (canonical preamble, from the repo root):
    set -a; source output/codes/.designg_bindings.sh; set +a
    S=$(mktemp -d); mkdir -p "$S/figs" "$S/data"
    DESIGNG_FIGDIR="$S/figs" DESIGNG_DATADIR="$S/data" DESIGNG_MANIFEST="$S/manifest.md" ~/.venv/nolds/bin/python "$DYN" --smoke
    DESIGNG_FIGDIR="$S/figs" DESIGNG_DATADIR="$S/data" DESIGNG_MANIFEST="$S/manifest.md" ~/.venv/nolds/bin/python "$QP" --smoke
    DESIGNG_VALIDATE_MODE=smoke DESIGNG_VALIDATE_IN="$S/data" DESIGNG_VALIDATE_FIGDIR="$S/figs" \
      DESIGNG_CLAIMS_IN="$CLAIMS" DESIGNG_REPORT_OUT="$S/r.md" DESIGNG_VALJSON_OUT="$S/v.json" \
      DESIGNG_MANIFEST="$S/manifest.md" ~/.venv/nolds/bin/python "$VAL"


"""
from __future__ import annotations

import functools
import glob
import hashlib
import re
import json
import math
import os
import sys
import bisect
import decimal
from datetime import datetime, timezone
from decimal import Decimal
from typing import Any, Dict, Iterator, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REQUIRED_ENV = (
    "DESIGNG_VALIDATE_IN",
    "DESIGNG_VALIDATE_FIGDIR",
    "DESIGNG_CLAIMS_IN",
    "DESIGNG_REPORT_OUT",
    "DESIGNG_VALJSON_OUT",
    "DESIGNG_MANIFEST",
)

KINDS = ("regression_anchor", "theory_invariant", "empirical_discovery")
BLOCKING_KINDS = ("regression_anchor", "theory_invariant")

# U-number aliases accepted by DESIGNG_VALIDATE_SCOPE (in addition to block names).
_U_ALIAS = {
    "U1": "extinction", "U2": "sensitivity", "U3": "finite_amp",
    "U4": "exact_skew", "U5": "quasiperiodic",
}

# Producer stems.
_DYN = "designG_dynamics_v4"
_QP = "designG_quasiperiodic_v2"
_ERG = "designG_ergodic_quadrature_v2"
_BIF = "designG_extinction_bif"  # designG_bifurcation_v2.py's mu0-sweep onsets at nu1=2.1 (emits no _taxonomy)
_ONS = "designG_onset_ensemble_v1"  # designG_onset_ensemble_v1.py's 24-seed x 3-resolution onset ensemble (file <stem>.json, no _taxonomy)
_CERT = "designG_period_two_certificate_v1"  # certify_period_two_v1.py: exact two-cycle certificate + numerical observations (file <stem>.json, no _taxonomy)

# Full-mode completeness spec: block -> what a COMPLETE phase must contain.
#   stem          -- producer whose JSON carries the block.
#   require        -- (producer_field, selector) pairs that must resolve to >=1 scalar.
#   taxonomy_keys  -- _taxonomy pointers the producer must declare ([] for the
#                     sibling, which emits no _taxonomy -- documented).
#   figures        -- figure basenames that must exist in DESIGNG_VALIDATE_FIGDIR.
#   sens_configs   -- (sensitivity only) required C0-C10 config count.
BLOCK_SPEC: Dict[str, Dict[str, Any]] = {
    "extinction": dict(stem=_DYN,
                       require=[("/validation_table/*/fraction_extinct", {"label": "ext_frac_nu1_2"})],
                       taxonomy_keys=["/validation_table/*/fraction_extinct"],
                       figures=["designG_extinction_fraction.png"]),
    "sensitivity": dict(stem=_DYN,
                        require=[("/sensitivity/base/{nu1}/fraction", {"nu1": "2"})],
                        taxonomy_keys=["/sensitivity/base/*/fraction", "/sensitivity/*/*/fraction"],
                        figures=["designG_extinction_sensitivity.png"], sens_configs=11),
    "finite_amp": dict(stem=_DYN,
                       require=[("/finite_amp/per_nu1/*/per_delta/*/decay_rate", {})],
                       taxonomy_keys=["/finite_amp/per_nu1/*/per_delta/*/decay_rate"],
                       figures=["designG_lambda_convergence.png"]),
    "exact_skew": dict(stem=_DYN,
                       require=[("/exact_skew/per_nu1/*/lambda_perp_quadrature", {"nu1": "2"}),
                                ("/exact_skew/per_nu1/*/collapse_guard_ok", {"nu1": "2"})],
                       taxonomy_keys=["/exact_skew/per_nu1/*/lambda_perp_quadrature",
                                      "/exact_skew/per_nu1/*/collapse_guard_ok"],
                       figures=["designG_exact_skew.png"]),
    "coexistence": dict(stem=_DYN,
                        require=[("/coexistence/jury_band/lower_edge_closed_form", {})],
                        taxonomy_keys=["/coexistence/jury_band/lower_edge_closed_form"],
                        figures=["designG_jury_band.png"]),
    "family_unification": dict(stem=_DYN,
                               require=[("/family_unification/per_nu1/*/pairwise_max_spread", {"nu1": "4.00"})],
                               taxonomy_keys=["/family_unification/per_nu1/*/pairwise_max_spread"],
                               figures=[]),
    "lambda_hat": dict(stem=_DYN,
                       require=[("/lambda_hat/per_nu1/{nu1}/lambda_hat", {"nu1": "2.00"})],
                       taxonomy_keys=["/lambda_hat/per_nu1/*/lambda_hat", "/lambda_hat/per_nu1/*/n_trajectories",
                                      "/lambda_hat/per_nu1/*/n_noncollapsed", "/lambda_hat/per_nu1/*/n_accepted_stable",
                                      "/lambda_hat/per_nu1/*/sem_across_accepted"],
                       figures=["designG_lambda_vs_nu1.png"]),
    "sweep": dict(stem=_BIF,
                  require=[("/chaos_onset_mu0", {}), ("/extinction_onset_mu0", {})],
                  taxonomy_keys=[], figures=["designG_bifurcation_mu0.png"]),
    # The ensemble numbers Section 6 and Section S12 quote (medians, ranges and
    # the plotted seed's onsets at the three resolutions) are value-checked by direct claims
    # instead of being covered only by the occurrence audit of the listed lines.
    # The certified two-cycle (Proposition prop:two-cycle) and the
    # numerical basin fraction quoted beside it are value-checked against the certificate JSON.
    "two_cycle": dict(stem=_CERT,
                      require=[("/certificate/contraction_bound_q", {}), ("/certificate/point_p_float/0", {}),
                               ("/basin_fraction_nu1_46_25/fraction_two_cycle", {})],
                      taxonomy_keys=[], figures=[]),
    "onset_ensemble": dict(stem=_ONS,
                           require=[("/per_resolution/0.01/summary/chaos_sustained/median", {}),
                                    ("/per_resolution/0.002/summary/chaos_sustained/median", {}),
                                    ("/per_resolution/0.01/per_seed/0/chaos_sustained", {})],
                           taxonomy_keys=[], figures=[]),
    "quasiperiodic": dict(stem=_QP,
                          require=[("rho", {}), ("largest_lyap", {}), ("normal_lyap", {}),
                                   ("invariant_circle_residual", {}), ("fft_peak_snr", {})],
                          taxonomy_keys=["rho", "largest_lyap", "normal_lyap",
                                         "invariant_circle_residual", "fft_peak_snr"],
                          figures=["designG_quasiperiodic_panel.png", "designG_quasiperiodic_fft.png"]),
    "arcsine": dict(stem=_ERG,
                    require=[("/arcsine_moments/E_lnx", {}), ("/arcsine_moments/jensen_gap", {})],
                    taxonomy_keys=[], figures=[]),
    "subaction": dict(stem=_ERG, require=[("/subaction/max_H", {})], taxonomy_keys=[], figures=[]),
    "I_family": dict(stem=_ERG, require=[("/I_family/nu1star_4", {})], taxonomy_keys=[], figures=[]),
    # The stable-leaf illustration (designG_stable_leaf_v1.py: manuscript Figures 5 and 12). Eight tracked claims
    # carry this block; without an entry here the inventory check could not miss them, the block could not be
    # named in DESIGNG_VALIDATE_SCOPE, and no completeness rule asked for its record or its two figures.
    "stable_leaf": dict(stem="designG_stable_leaf_v1",
                        require=[(f"/nu1_{nu}/{field}", {}) for nu in ("2.1", "3")
                                 for field in ("psi_over_t2_small_t", "matched/max_abs_ratio_minus_1",
                                               "reference_rate_at_40", "reference_rate_N_long")],
                        taxonomy_keys=[], figures=["designG_stable_leaf.png", "designG_stable_leaf_rates.png"]),
}

# The blocks a FULL, unscoped run must find in the claims inventory are
# fixed HERE, independently of the inventory under check. An inventory replaced by `{}` used to
# derive an empty block list from its own absent claims and report SUBMISSION_READY over zero
# checks. (`global_telescoping` has claims but no completeness spec in BLOCK_SPEC.)
REQUIRED_CLAIM_BLOCKS = sorted(set(BLOCK_SPEC) | {"global_telescoping"})

# A result record whose own success flags report a failure must never be
# validated as ready -- a certificate that is not certified or has a failed check, a failed
# collapse guard, an incomplete periodic-orbit enumeration. Each entry: (producer stem, the blocks
# that bring it into scope, pointer with the `*` wildcard, requirement).
ESSENTIAL_FLAGS = [
    (_CERT, ("two_cycle",), "/certificate/certified", "True"),
    (_CERT, ("two_cycle",), "/certificate/checks/*/ok", "True"),
    (_DYN, ("exact_skew",), "/exact_skew/per_nu1/*/collapse_guard_ok", "True"),
    (_ERG, ("arcsine", "subaction", "I_family"), "/periodic_enumeration/coverage_complete", "True"),
    (_ERG, ("arcsine", "subaction", "I_family"), "/periodic_enumeration/least_period_count_mismatches", "0"),
]
# The members those flags must cover, fixed here: a deleted member or a deleted flag must fail,
# not shrink the set the gate inspects. The nine named checks of the two-cycle certificate and the
# three exact-skew parameters (kept in both modes).
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
# The section files an audit row may point at.
SECTION_FILES = ("abstract_clean.tex", "introduction_clean.tex", "methods_clean.tex", "resultsAndDiscussions_clean.tex",
                 "conclusions_clean.tex", "supplementary_clean.tex", "reportOnAITools_clean.tex")
# The exclusion categories an audit row may use, fixed here: a disposition that names no known claim and
# no known category documents nothing. The second set may stand for a literal with no checkable number
# (a symbol written out, not a value).
AUDIT_EXCLUSION_CATEGORIES = frozenset({
    "arithmetic_demonstration", "closed_form_band_edge", "closed_form_constant", "closed_form_figure_value",
    "derived_quantity", "diagnostic_error_bar", "external_measurement", "figure_input_parameter",
    "frozen_verification_record", "input_parameter_or_tolerance", "literature_constant", "not_a_quantity",
    "protocol_constant", "protocol_observation", "protocol_parameter", "regression_covered", "rounded_ceiling",
    "symbolic_mention"})
AUDIT_UNCHECKABLE = frozenset({"symbolic_mention", "not_a_quantity"})
# The section files whose printed numbers the coverage audit must account for (the declarations are excluded).
AUDIT_FILES = ("introduction_clean.tex", "methods_clean.tex", "resultsAndDiscussions_clean.tex",
               "conclusions_clean.tex", "supplementary_clean.tex", "abstract_clean.tex")


class SelectorError(ValueError):
    """A producer_field + selector that does NOT resolve to exactly one scalar."""


# ---------------------------------------------------------------------------
# Selector / JSON-Pointer resolution
# ---------------------------------------------------------------------------

def _split_pointer(pointer: str) -> List[str]:
    """Split a pointer into tokens. A bare top-level name (no leading '/') like
    'rho' is a single-token pointer; otherwise split on '/'."""
    if pointer.startswith("/"):
        return pointer.strip("/").split("/")
    return [pointer]


def _eq(a: Any, b: Any) -> bool:
    """Selector value equality: string-equal, or numeric-close."""
    if str(a) == str(b):
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
    except (TypeError, ValueError):
        return False


def _substitute_named(tokens: List[str], selector: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """Replace every ``{name}`` token by ``selector[name]`` (consuming that key).
    Returns (new_tokens, remaining_selector). Missing name -> SelectorError."""
    remaining = dict(selector)
    out: List[str] = []
    for t in tokens:
        if len(t) >= 2 and t[0] == "{" and t[-1] == "}":
            name = t[1:-1]
            if name not in remaining:
                raise SelectorError(f"named segment {{{name}}} not supplied by selector {selector!r}")
            out.append(str(remaining.pop(name)))
        else:
            out.append(t)
    return out, remaining


def _walk(obj: Any, tokens: List[str], star_keys: List[Any], star_recs: List[Any]) -> Iterator[Tuple[Any, List[Any], List[Any]]]:
    """Yield (leaf_value, star_keys, star_recs) for every path matching ``tokens``.
    '*' = any non-'_' dict key or any list index (mirrors the test's _resolve_pointer);
    a literal digit token indexes a list; a literal string descends a dict key."""
    if not tokens:
        yield obj, star_keys, star_recs
        return
    head, rest = tokens[0], tokens[1:]
    if head == "*":
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(k, str) and k.startswith("_"):
                    continue
                yield from _walk(v, rest, star_keys + [k], star_recs + [v])
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                yield from _walk(v, rest, star_keys + [i], star_recs + [v])
        return
    if isinstance(obj, dict) and head in obj:
        yield from _walk(obj[head], rest, star_keys, star_recs)
    elif isinstance(obj, list) and head.lstrip("-").isdigit():
        idx = int(head)
        if -len(obj) <= idx < len(obj):
            yield from _walk(obj[idx], rest, star_keys, star_recs)


def resolve_scalars(root: Any, producer_field: str, selector: Dict[str, Any]) -> List[Any]:
    """All scalar values matching ``producer_field`` under ``selector``.

    Selector semantics (per Task B6): a key naming a ``{named}`` segment substitutes
    it; a remaining key is a field-value constraint on the '*'-selected record, or
    (if not a record field) the required value of the '*'-segment key."""
    tokens, remaining = _substitute_named(_split_pointer(producer_field), dict(selector))
    candidates = list(_walk(root, tokens, [], []))

    def matches(star_keys: List[Any], star_recs: List[Any]) -> bool:
        for k, v in remaining.items():
            hit = False
            for rec in star_recs:  # field-value constraint on a '*'-selected record
                if isinstance(rec, dict) and k in rec and _eq(rec[k], v):
                    hit = True
                    break
            if not hit and any(_eq(sk, v) for sk in star_keys):  # segment-key constraint
                hit = True
            if not hit:
                return False
        return True

    return [leaf for leaf, sk, sr in candidates if matches(sk, sr)]


def resolve_one(root: Any, producer_field: str, selector: Dict[str, Any]) -> Any:
    """Resolve to EXACTLY ONE scalar. HARD-FAIL (SelectorError) on zero, multiple,
    or a non-scalar (dict/list) match."""
    vals = resolve_scalars(root, producer_field, selector)
    if len(vals) == 0:
        raise SelectorError(f"selector {selector!r} on {producer_field!r} matched ZERO scalars")
    if len(vals) > 1:
        raise SelectorError(
            f"selector {selector!r} on {producer_field!r} matched {len(vals)} scalars (expected exactly 1): {vals[:6]!r}")
    v = vals[0]
    if isinstance(v, (dict, list)):
        raise SelectorError(f"selector {selector!r} on {producer_field!r} matched a non-scalar {type(v).__name__}")
    return v


# ---------------------------------------------------------------------------
# Taxonomy kind lookup (concrete pointer -> most-specific _taxonomy kind)
# ---------------------------------------------------------------------------

def _taxonomy_lookup_tokens(producer_field: str, selector: Dict[str, Any]) -> List[str]:
    """Concrete tokens for a _taxonomy lookup: substitute {named} segments; leave
    bare '*' as '*' (the taxonomy keys are themselves '*'-wildcarded)."""
    tokens, _ = _substitute_named(_split_pointer(producer_field), dict(selector))
    return tokens


def _tax_key_matches(concrete: List[str], taxkey: str) -> Optional[int]:
    """If ``taxkey`` (a _taxonomy pointer) matches ``concrete`` tokens, return its
    specificity (number of literal, non-'*' segments); else None."""
    tk = _split_pointer(taxkey)
    if len(tk) != len(concrete):
        return None
    literals = 0
    for t, c in zip(tk, concrete):
        if t == "*":
            continue
        if t != c:
            return None
        literals += 1
    return literals


def taxonomy_kind_for(producer_field: str, selector: Dict[str, Any], taxonomy: Dict[str, str]) -> Optional[str]:
    """Most-specific matching _taxonomy kind for the (substituted) field, or None
    if no taxonomy key matches (e.g. a horizon-sweep cell deeper than any key)."""
    concrete = _taxonomy_lookup_tokens(producer_field, selector)
    best_kind, best_spec = None, -1
    for key, kind in taxonomy.items():
        spec = _tax_key_matches(concrete, key)
        if spec is not None and spec > best_spec:
            best_kind, best_spec = kind, spec
    return best_kind


# ---------------------------------------------------------------------------
# Producer JSON loading
# ---------------------------------------------------------------------------

def producer_path(data_dir: str, stem: str, mode: str) -> Optional[str]:
    """Resolve <stem>_{smoke,values}.json in ``data_dir``, preferring the
    mode-appropriate suffix. None if neither exists."""
    smoke = os.path.join(data_dir, stem + "_smoke.json")
    values = os.path.join(data_dir, stem + "_values.json")
    bare = os.path.join(data_dir, stem + ".json")   # the onset ensemble writes <stem>.json (no --smoke variant)
    order = [smoke, values, bare] if mode == "smoke" else [values, smoke, bare]
    for p in order:
        if os.path.exists(p):
            return p
    return None


def load_producers(data_dir: str, stems: List[str], mode: str) -> Dict[str, Tuple[Dict[str, Any], str]]:
    out: Dict[str, Tuple[Dict[str, Any], str]] = {}
    for stem in stems:
        p = producer_path(data_dir, stem, mode)
        if p is None:
            continue
        with open(p) as fh:
            out[stem] = (json.load(fh), p)
    return out


# ---------------------------------------------------------------------------
# CI-method helper (sensitivity: Wilson vs none)
# ---------------------------------------------------------------------------

def _sibling_scalar(root: Any, producer_field: str, selector: Dict[str, Any], sibling: str) -> Any:
    """Resolve the sibling field (same record, terminal segment replaced by
    ``sibling``) -- e.g. 'fraction' -> 'ci_method'. None if unresolvable."""
    toks = _split_pointer(producer_field)
    if not toks:
        return None
    sib_field = "/" + "/".join(toks[:-1] + [sibling]) if producer_field.startswith("/") else sibling
    try:
        return resolve_one(root, sib_field, selector)
    except SelectorError:
        return None



# ---------------------------------------------------------------------------
# Manuscript figure inventory (2026-10-17): the manuscript
# itself is parsed, so an empty or incomplete manifest, or a quoted number that
# no longer appears in the text, is caught here rather than silently passing.
# ---------------------------------------------------------------------------
_BUNDLE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TEXDIR: Optional[str] = None  # set in main(); the directory holding the *_clean.tex files


def manuscript_figures(texdir: str) -> List[str]:
    """Basenames of every ``\\includegraphics{figures/...}`` in ``texdir/*_clean.tex``."""
    names = set()
    for tex in sorted(glob.glob(os.path.join(texdir, "*_clean.tex"))):
        with open(tex) as fh:
            for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}", fh.read()):
                names.add(os.path.basename(m.group(1).strip()))
    return sorted(names)


# One complete-number grammar for BOTH matching paths: optional
# minus, optional TeX spacing (whitespace incl. line breaks, thin spaces \, \; \! \:, an empty
# group, ~), mantissa (digits, or digits grouped in thousands by the TeX separator "{,}":
# 40{,}000 is the one number 40000 and never an occurrence of 40; a grouped
# mantissa must end at a digit boundary, so a malformed 40{,}0001 matches nothing),
# optional exponent written as "\times 10^{b}", "\cdot 10^{b}" (b with an optional sign, braces
# optional, spacing allowed around the operator) or e-notation.
# Groups: (sign, mantissa, tex_exponent, e_exponent).
_SP = r"(?:\s|\\[,;!:]|\{\}|~)"
_GRP = r"\{,\}"
_NUM_RE = re.compile(r"(-?)" + _SP + r"*(\d{1,3}(?:" + _GRP + r"\d{3})+(?:\.\d+)?(?!\d)|\d+(?:\.\d+)?)"
                     r"(?:" + _SP + r"*\\(?:times|cdot)" + _SP + r"*10\^\{?([+-]?\d+)\}?|[eE]([+-]?\d+))?")
# After an exact token nothing may continue the number: no digit, no decimal part, no thousands
# group "{,}ddd", no multiplication/exponent operator (even after spacing), no e-notation.
_TOKEN_TAIL = r"(?!\d|\.\d|" + _GRP + r"\d|" + _SP + r"*\\(?:times|cdot)|" + _SP + r"*[eE][+-]?\d)"
_CONTINUATION_RE = re.compile(r"\.\d|" + _GRP + r"\d|" + _SP + r"*\\(?:times|cdot)|" + _SP + r"*[eE][+-]?\d")


def _signed_before(text: str, pos: int) -> bool:
    r"""Is the token starting at ``pos`` preceded by a sign, possibly across whitespace (line
    breaks included), TeX thin spaces (\, \; \! \:), an empty group or ~?"""
    i = pos
    while i > 0:
        if text[i - 1] in " \t\n\r~":
            i -= 1
            continue
        for sp in ("\\,", "\\;", "\\!", "\\:", "{}"):
            if text[max(0, i - len(sp)):i] == sp:
                i -= len(sp)
                break
        else:
            break
    return i > 0 and text[i - 1] in "+-"


def _half_ulp(x: float, text: str) -> float:
    """Half a unit of the last digit of ``text``, the repr of a quoted value. The repr is read as a decimal, so
    an exponent form keeps its own last digit: '1e-05' gives 5e-06 and '5.7e-05' gives 5e-07. (Counting the
    characters after the point gave 0.5 and 5e-06 for them, and a tolerance of 0.5 made every small number of
    the file an occurrence of the quoted value.)"""
    try:
        return float(Decimal(5).scaleb(Decimal(text).as_tuple().exponent - 1))
    except (decimal.InvalidOperation, TypeError, ValueError):   # not a finite decimal: the old reading
        dec = len(text.split(".")[1]) if "." in text else 0
        return 0.5 * 10.0 ** (-dec)


def _quote_hits_in_file(claim: Dict[str, Any], path: str) -> List[Tuple[int, str]]:
    """All matches of the claim's quoted number in one ``*_clean.tex`` file (line, token text)."""
    with open(path) as fh:
        text = _printed_text(fh.read())
    q = claim["quoted_value"]
    qf = float(q)
    base_tol = max(float(claim.get("display_tol") or 0.0), _half_ulp(qf, repr(q)))
    hits: List[Tuple[int, str]] = []
    token = claim.get("manuscript_token")
    if token:
        # Coarse quoted values (0, 1, 4) match too many incidental numbers under any tolerance,
        # so the inventory may pin the EXACT printed token (e.g. "0.0000", "1.0000", "4.000");
        # only that token, as a whole number token, counts.
        # A sign character before the token is excluded too: "-3.95" must not satisfy "3.95"; so is a
        # thousands separator: the "000" of "40{,}000" is not a token of its own.
        # The numeric path below captures the sign inside _NUM_RE.
        pat = re.compile(r"(?<![\d.+-])(?<!\{,\})" + re.escape(str(token)) + _TOKEN_TAIL)
        for m in pat.finditer(text):
            if _signed_before(text, m.start()):
                continue
            hits.append((text.count("\n", 0, m.start()) + 1, m.group(0)))
    for m in ([] if token else _NUM_RE.finditer(text)):
        if text[max(0, m.start(2) - 3):m.start(2)] == "{,}":
            continue  # the tail group of a grouped numeral (its head was not a complete number): not a number of its own
        try:
            val = float(m.group(2).replace("{,}", ""))
        except ValueError:
            continue
        if m.group(1):
            val = -val
        exp = m.group(3) or m.group(4)
        if exp:
            val *= 10.0 ** int(exp)
        if _CONTINUATION_RE.match(text[m.end():m.end() + 24] or ""):
            continue  # an unconsumed decimal / operator / e-notation continuation: not this number
        # The tolerance is the INVENTORY's: display_tol or half a unit of the quoted value's
        # last digit. The text token's own coarseness is deliberately not used (a bare "0"
        # or "1" must not match 0.2269), so a manuscript that prints more digits than the
        # inventory records is flagged -- that is an inventory defect, to be fixed there.
        if abs(val - qf) <= base_tol:
            hits.append((text.count("\n", 0, m.start()) + 1, m.group(0)))
    return hits


def manuscript_quote_check(claim: Dict[str, Any], texdir: Optional[str]) -> Optional[Dict[str, Any]]:
    """Does the number the inventory says the manuscript quotes actually appear in the manuscript?

    Searches each ``*_clean.tex`` file the claim names (resolved by basename inside
    ``texdir``) for a numeric token -- plain or ``a \\times 10^{b}`` -- within
    max(display_tol, half a unit of the last printed digit) of ``quoted_value``.
    Returns None for untracked claims. ``line`` in the inventory is informational;
    the lines actually found are reported, together with the number of occurrences.
    The inventory pins the minimum number of occurrences PER FILE: ``quoted_in`` maps
    each section file to its pinned count (a number repeated in a
    Supplement table is checked there too); without it, ``tex`` with
    ``manuscript_min_occurrences`` (default 1) is the single pinned file. Fewer matching
    tokens than pinned in any file is a TOO_FEW status, so a quotation edited at one of
    several sites is caught even while another site still carries the token.
    This remains a presence/occurrence check of the quoted NUMBER in the named
    section files; it does not bind the number to a particular sentence.
    """
    tex, q = claim.get("tex"), claim.get("quoted_value")
    if not texdir or not tex or q is None or not claim.get("manuscript_tracked"):
        return None
    files: Dict[str, int] = dict(claim.get("quoted_in") or {})
    primary = os.path.basename(str(tex))
    if not files:
        files = {primary: int(claim.get("manuscript_min_occurrences") or 1)}
    per_file: Dict[str, Dict[str, Any]] = {}
    statuses: List[str] = []
    for base, need in files.items():
        path = os.path.join(texdir, os.path.basename(base))
        if not os.path.exists(path):
            per_file[base] = {"status": "TEX_MISSING", "n_occurrences": 0, "min_occurrences": int(need), "lines": [], "tokens": []}
            statuses.append("TEX_MISSING")
            continue
        hits = _quote_hits_in_file(claim, path)
        st = "FOUND" if len(hits) >= int(need) else ("TOO_FEW" if hits else "NOT_FOUND")
        per_file[base] = {"status": st, "n_occurrences": len(hits), "min_occurrences": int(need),
                          "lines": [h[0] for h in hits[:8]], "tokens": sorted({h[1].strip() for h in hits})[:6]}
        statuses.append(st)
    if "TEX_MISSING" in statuses:
        status = "TEX_MISSING"
    elif all(st == "FOUND" for st in statuses):
        status = "FOUND"
    elif any(st in ("FOUND", "TOO_FEW") for st in statuses):
        status = "TOO_FEW"
    else:
        status = "NOT_FOUND"
    main = per_file.get(primary) or next(iter(per_file.values()))
    return {"status": status, "path": primary, "n_occurrences": main["n_occurrences"],
            "min_occurrences": main["min_occurrences"], "token": claim.get("manuscript_token"),
            "lines": main["lines"], "tokens": main["tokens"], "files": per_file}


def quote_shortfall_message(name: str, quoted: Any, mq: Dict[str, Any]) -> str:
    """The blocking message for a tracked number that is not found often enough. It names every pinned file
    that fell short, with the occurrences found and pinned there. (The claim's primary file used to be named,
    also when that file was complete and another pinned file was short.)"""
    short = [f"{base}: {pf.get('n_occurrences')} of {pf.get('min_occurrences')} pinned ({pf.get('status')})"
             for base, pf in (mq.get("files") or {}).items() if pf.get("status") != "FOUND"]
    where = "; ".join(short) if short else f"{mq.get('path')} ({mq.get('status')})"
    return f"{name}: quoted value {quoted!r} not found in {where}"

# ---------------------------------------------------------------------------
# Manifest verification (the validator VERIFIES; producers AUTHOR)
# ---------------------------------------------------------------------------

def verify_manifest(path: str, texdir: Optional[str] = None, figdir: Optional[str] = None) -> List[str]:
    """Check every ``| artifact | version | sha256 |`` row (file exists + sha matches), refuse an
    EMPTY manifest, and -- when ``texdir`` is given -- require a row for every figure the
    manuscript includes (``\\includegraphics{figures/...}`` in ``texdir/*_clean.tex``).
    When ``figdir`` is also given, that row must name the figure AT THE CONFIGURED FIGURE
    DIRECTORY (real-path equality), so a correctly hashed row for some OTHER copy of a
    same-named file cannot satisfy the coverage check.

    Relative artifact paths are resolved against the manifest's own directory. A manifest
    with no rows, or one that omits a manuscript figure, is an issue: the gate is not
    allowed to pass vacuously.
    """
    issues: List[str] = []
    if not os.path.exists(path):
        return [f"manifest missing: {path}"]
    rows = 0
    rowed: set = set()
    rowed_real: set = set()
    mdir = os.path.dirname(os.path.abspath(path))
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line.startswith("|"):
                continue
            parts = [p.strip() for p in line.strip("|").split("|")]
            if len(parts) < 3:
                continue
            fpath, producer, sha = parts[0], parts[1], parts[2]
            if fpath.lower() == "artifact":  # tolerate a header row if ever added
                continue
            rows += 1
            # The producer column must name a script shipped beside this validator. This is a
            # presence check on the producer NAME plus byte integrity of the artifact; it does
            # not authenticate source provenance by digest (stated in README).
            scripts = re.findall(r"[A-Za-z0-9_.-]+\.py", producer)
            here = os.path.dirname(os.path.abspath(__file__))
            if not scripts:
                issues.append(f"manifest row for {fpath}: producer column {producer!r} names no script")
            for sc in scripts:
                if not os.path.exists(os.path.join(here, sc)):
                    issues.append(f"manifest row for {fpath}: producer {sc} is not shipped beside the validator")
            if not os.path.isabs(fpath) and not os.path.exists(fpath):
                fpath = os.path.join(mdir, fpath)
            rowed.add(os.path.basename(fpath))
            rowed_real.add(os.path.realpath(fpath))
            if not os.path.exists(fpath):
                issues.append(f"manifest artifact missing on disk: {fpath}")
                continue
            got = hashlib.sha256(open(fpath, "rb").read()).hexdigest()
            if got != sha:
                issues.append(f"manifest sha256 mismatch: {fpath}")
    if rows == 0:
        issues.append("manifest has no artifact rows (an empty manifest is not a passing manifest)")
    if texdir:
        required = manuscript_figures(texdir)
        if not required:
            issues.append(f"no *_clean.tex with \\includegraphics found under {texdir}; cannot check figure coverage")
        for name in required:
            if figdir:
                want = os.path.realpath(os.path.join(figdir, name))
                if want not in rowed_real:
                    issues.append(f"manifest lacks a row for manuscript figure {name} at the configured figure directory ({want})")
            elif name not in rowed:
                issues.append(f"manifest lacks a row for manuscript figure {name}")
    return issues


# ---------------------------------------------------------------------------
# Per-claim evaluation (one scalar, two checks)
# ---------------------------------------------------------------------------

def evaluate_claim(claim: Dict[str, Any], producers: Dict[str, Tuple[Dict[str, Any], str]], mode: str) -> Dict[str, Any]:
    cid = claim["id"]
    kind = claim.get("producer_kind")
    tex = claim.get("tex")
    tex_ref = f"{tex}:{claim.get('line')}" if tex else "(untracked)"
    res: Dict[str, Any] = {
        "kind": kind, "name": cid, "producer": claim.get("producer"),
        "producer_field": claim.get("producer_field"), "selector": claim.get("producer_selector"),
        "block": claim.get("block"), "manuscript_tracked": bool(claim.get("manuscript_tracked")),
        "tex_ref": tex_ref, "actual": None, "expected": None, "tol": None,
        "status": None, "manuscript_change_required": None, "ci_method": None,
        "hard_fail": False, "blocking_fail": False, "detail": None,
    }
    if kind not in KINDS:
        res["status"] = "KIND_INVALID"
        res["detail"] = f"producer_kind {kind!r} not in {KINDS}"
        res["hard_fail"] = True
        return res

    prod = producers.get(claim.get("producer"))
    if prod is None:
        # Producer JSON absent. In smoke this is never a block (the sibling has no
        # smoke variant); in full it is caught by the completeness gate.
        res["status"] = "SKIPPED_SMOKE" if mode == "smoke" else "ABSENT"
        res["detail"] = f"producer JSON {claim.get('producer')!r} absent from DESIGNG_VALIDATE_IN"
        return res
    root = prod[0]

    # Resolve the ONE scalar (schema/selector check in both modes).
    try:
        val = resolve_one(root, claim["producer_field"], claim.get("producer_selector", {}))
    except SelectorError as exc:
        if mode == "smoke" and claim.get("smoke_optional"):
            # The reduced --smoke output legitimately lacks this field (fewer nu1 values,
            # fewer sensitivity configs, or no exact-zero event); the inventory says so.
            # Full mode is unaffected: zero/multiple matches remain a HARD FAIL there.
            res["status"] = "SKIPPED_SMOKE"
            res["detail"] = f"smoke_optional field absent from the --smoke output: {exc}"
            return res
        res["status"] = "SELECTOR_ERROR"
        res["detail"] = str(exc)
        res["hard_fail"] = True
        return res
    res["actual"] = val

    # A NON-FINITE producer value where the claim requires a finite
    # number is a HARD failure in both modes. NaN compares false against everything, so the
    # earlier code reported a NaN as CONSISTENT with any quoted number (and as OK_NONBLOCKING
    # for an empirical kind). A claim whose scientific contract admits an infinity lists it in
    # ``allowed_nonfinite`` (e.g. ["-inf"] for a tangential exponent at an exact critical hit).
    # The gate must act on the JSON TYPE, not only on Python floats. A JSON string
    # such as "NaN" or "not a number" converted silently via float(val) downstream (a string
    # "NaN" became a NaN that compared "consistent"); a bool is a number only if the claim says
    # so (value_type "bool"). Anything else that is not an int/float is NONNUMERIC.
    if isinstance(val, bool):
        if str(claim.get("value_type", "")).lower() != "bool":
            res["status"] = "NONNUMERIC"
            res["detail"] = "producer value is a bool but the claim does not declare value_type 'bool'"
            res["hard_fail"] = True
            return res
    elif not isinstance(val, (int, float)):
        res["status"] = "NONNUMERIC"
        res["detail"] = f"producer value {val!r} ({type(val).__name__}) is not a JSON number"
        res["hard_fail"] = True
        return res
    elif isinstance(val, float) and not math.isfinite(val):
        tag = "nan" if math.isnan(val) else ("+inf" if val > 0 else "-inf")
        allowed = [str(a).strip().lower() for a in (claim.get("allowed_nonfinite") or [])]
        if tag not in allowed:
            res["status"] = "NONFINITE"
            res["detail"] = (f"producer value is {tag}; the claim requires a finite number "
                             f"(allowed_nonfinite={allowed})")
            res["hard_fail"] = True
            return res

    # ci_method (sensitivity claims): read the produced ci_method next to 'fraction'. A claim that
    # names the method it expects (`ci_expected`) fails when the record carries another one, so a
    # record relabelled "none" cannot switch off the interval check.
    if claim.get("block") == "sensitivity" or claim.get("ci_expected") is not None:
        res["ci_method"] = _sibling_scalar(root, claim["producer_field"], claim.get("producer_selector", {}), "ci_method")
        if claim.get("ci_expected") is not None and res["ci_method"] != claim["ci_expected"]:
            res["status"] = "CI_METHOD_MISMATCH"
            res["detail"] = f"produced ci_method {res['ci_method']!r}, the claim expects {claim['ci_expected']!r}"
            res["hard_fail"] = True
            return res

    # Taxonomy kind cross-check (skip when the producer emits no _taxonomy).
    tax = root.get("_taxonomy") if isinstance(root, dict) else None
    if tax and not claim.get("producer_emits_no_taxonomy"):
        tk = taxonomy_kind_for(claim["producer_field"], claim.get("producer_selector", {}), tax)
        res["taxonomy_kind"] = tk
        if tk is not None and tk != kind:
            res["status"] = "TAXONOMY_CONFLICT"
            res["detail"] = f"claim kind {kind!r} != producer _taxonomy {tk!r} for {claim['producer_field']!r}"
            res["hard_fail"] = True
            return res

    # --- Check 1: producer-kind check ---
    if mode == "smoke":
        res["status"] = "SKIPPED_SMOKE"
    elif kind in BLOCKING_KINDS:
        exp, tol = claim.get("anchor_expected"), claim.get("anchor_tol")
        res["expected"], res["tol"] = exp, tol
        if exp is None or tol is None:
            res["status"] = "ANCHOR_MISSING"
            res["detail"] = "blocking kind lacks anchor_expected/anchor_tol"
            res["hard_fail"] = True
            return res
        ok = abs(float(val) - float(exp)) <= float(tol)
        res["status"] = "PASS" if ok else "FAIL"
        res["blocking_fail"] = not ok
    else:  # empirical_discovery -- non-blocking
        res["status"] = "OK_NONBLOCKING"
        if claim.get("anchor_expected") is not None:
            res["expected"] = claim["anchor_expected"]

    # --- Check 2: manuscript-consistency check (non-blocking) ---
    if claim.get("manuscript_tracked") and claim.get("quoted_value") is not None:
        if mode == "smoke":
            res["manuscript_change_required"] = None  # SKIPPED_SMOKE
            res["manuscript_status"] = "SKIPPED_SMOKE"
        else:
            q = float(claim["quoted_value"])
            dtol = float(claim.get("display_tol") or 0.0)
            mism = abs(float(val) - q) > dtol
            within_ci = False
            # CI-method branch (Wilson vs none). The cell's ci_low/ci_high/n are the interval of its
            # terminal `fraction` only; a claim on another field of the same cell (first passage,
            # censored-inclusive first passage) must not be excused by an interval of a different
            # quantity.
            if res.get("ci_method") == "wilson" and str(claim["producer_field"]).rstrip("/").endswith("/fraction"):
                sel = claim.get("producer_selector", {})
                lo = _sibling_scalar(root, claim["producer_field"], sel, "ci_low")
                hi = _sibling_scalar(root, claim["producer_field"], sel, "ci_high")
                n_den = _sibling_scalar(root, claim["producer_field"], sel, "n")
                n_tot = _sibling_scalar(root, claim["producer_field"], sel, "n_total")
                n_cen = _sibling_scalar(root, claim["producer_field"], sel, "n_censored")
                # The interval is used to excuse a quotation only after it
                # is shown to be the Wilson interval of the produced count over the recorded sample
                # (n = n_total - n_censored); an inconsistent record is a hard failure, never
                # "consistent".
                bad = wilson_record_error(val, lo, hi, n_den, n_tot, n_cen)
                if bad:
                    res["status"] = "CI_INCONSISTENT"
                    res["detail"] = bad
                    res["hard_fail"] = True
                    return res
                if mism and lo <= q <= hi:
                    mism = False
                    within_ci = True
            res["manuscript_change_required"] = bool(mism)
            # CONSISTENT_WITHIN_CI: statistically compatible with the produced fraction, not equal
            # to it at display precision (transcription agreement is plain CONSISTENT).
            res["manuscript_status"] = ("CHANGE_REQUIRED" if mism else
                                        ("CONSISTENT_WITHIN_CI" if within_ci else "CONSISTENT"))
    # --- Check 3: does the quoted number actually appear in the manuscript text? ---
    res["manuscript_quote"] = manuscript_quote_check(claim, _TEXDIR)
    return res


def wilson_record_error(frac: Any, lo: Any, hi: Any, n: Any, n_total: Any, n_censored: Any,
                        z: float = 1.96) -> Optional[str]:
    """None when a produced Wilson record is internally consistent, else the reason: finite
    ordered bounds in [0, 1] containing the produced fraction (up to 1e-12 of roundoff: the
    producer's wilson(n, n) has an upper bound of 1 - 2e-16), a positive integer denominator n
    equal to n_total - n_censored -- both counts are REQUIRED, a missing or null count is an
    error, never a skipped check -- an integer extinct count fraction*n, and bounds equal to the
    Wilson interval of that count (the producer's own formula, z = 1.96)."""
    def real(v: Any) -> bool:
        return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))
    def count(v: Any) -> bool:
        return isinstance(v, int) and not isinstance(v, bool) and v >= 0
    if not (real(frac) and real(lo) and real(hi)):
        return f"Wilson record is not finite numbers (fraction {frac!r}, ci_low {lo!r}, ci_high {hi!r})"
    tol = 1e-12
    if not (-tol <= lo <= hi <= 1.0 + tol and lo - tol <= frac <= hi + tol):
        return f"Wilson bounds [{lo}, {hi}] do not contain the produced fraction {frac} within [0, 1]"
    if not (isinstance(n, int) and not isinstance(n, bool) and n > 0):
        return f"Wilson denominator n={n!r} is not a positive integer"
    if not (count(n_total) and count(n_censored) and n_censored <= n_total and n == n_total - n_censored):
        return (f"Wilson denominator n={n!r} is not n_total - n_censored "
                f"(n_total {n_total!r}, n_censored {n_censored!r}; both must be recorded)")
    k = float(frac) * n
    if abs(k - round(k)) > 1e-6:
        return f"fraction {frac} x n {n} = {k} is not an integer count"
    k = round(k)
    p_ = k / n
    d = 1 + z * z / n
    c = p_ + z * z / (2 * n)
    h = z * math.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n))
    wlo, whi = (c - h) / d, (c + h) / d
    if abs(wlo - float(lo)) > 1e-9 or abs(whi - float(hi)) > 1e-9:
        return f"Wilson bounds [{lo}, {hi}] differ from the interval of {k}/{n}, [{wlo:.12g}, {whi:.12g}]"
    return None


_MISSING = object()   # a member reached through `*` that lacks the requested field


def _pointer_values(root: Any, pointer: str) -> List[Tuple[str, Any]]:
    """(path, value) pairs a '/a/*/b' pointer reaches; `*` spans every non-underscore dict key and
    every list item. A segment missing BEFORE any `*` yields no pair (the record is absent); one
    missing inside a member reached through `*` yields (path, _MISSING), so a member that lacks the
    flag fails instead of being skipped."""
    out: List[Tuple[str, Any]] = []

    def rec(node: Any, toks: List[str], path: str, wild: bool) -> None:
        if not toks:
            out.append((path, node))
            return
        tok, rest = toks[0], toks[1:]
        if tok == "*":
            if isinstance(node, dict):
                for k, v in node.items():
                    if not str(k).startswith("_"):
                        rec(v, rest, f"{path}/{k}", True)
            elif isinstance(node, list):
                for i, v in enumerate(node):
                    rec(v, rest, f"{path}/{i}", True)
        elif isinstance(node, dict) and tok in node:
            rec(node[tok], rest, f"{path}/{tok}", wild)
        elif isinstance(node, list) and tok.isdigit() and int(tok) < len(node):
            rec(node[int(tok)], rest, f"{path}/{tok}", wild)
        elif wild:
            out.append((f"{path}/{tok}", _MISSING))

    rec(root, [t for t in pointer.split("/") if t], "", False)
    return out


def essential_flags_gate(producers: Dict[str, Tuple[Dict[str, Any], str]], blocks: List[str], mode: str) -> List[str]:
    """HARD-fail messages for result records whose own success flags report a failure
    (ESSENTIAL_FLAGS). A flag must be the JSON bool true and a mismatch count the integer 0; in
    full mode an absent producer or flag fails as well. Both modes."""
    fails: List[str] = []
    for stem, owners, pointer, need in ESSENTIAL_FLAGS:
        if not any(b in blocks for b in owners):
            continue
        prod = producers.get(stem)
        if prod is None:
            if mode == "full":
                fails.append(f"[essential] {stem}: producer JSON absent ({pointer} must be {need})")
            continue
        hits = _pointer_values(prod[0], pointer)
        if not hits:
            if mode == "full":
                fails.append(f"[essential] {stem}{pointer}: absent (must be present and {need})")
            continue
        for path, v in hits:
            ok = (v is True) if need == "True" else (isinstance(v, int) and not isinstance(v, bool) and v == 0)
            if not ok:
                shown = "<missing>" if v is _MISSING else repr(v)
                fails.append(f"[essential] {stem}{path} = {shown}: required {need}")
    # The members themselves (both modes, whenever the producer is present).
    cert = producers.get(_CERT)
    if "two_cycle" in blocks and cert is not None:
        checks = (cert[0].get("certificate") or {}).get("checks") if isinstance(cert[0], dict) else None
        names = [c.get("check") if isinstance(c, dict) else None for c in checks] if isinstance(checks, list) else None
        if names is None or len(names) != len(EXPECTED_CERT_CHECKS) or sorted(map(str, names)) != sorted(EXPECTED_CERT_CHECKS):
            fails.append(f"[essential] {_CERT}/certificate/checks: members {names!r} are not the nine expected checks")
    dyn = producers.get(_DYN)
    if "exact_skew" in blocks and dyn is not None:
        per = (dyn[0].get("exact_skew") or {}).get("per_nu1") if isinstance(dyn[0], dict) else None
        keys = sorted(k for k in per if not str(k).startswith("_")) if isinstance(per, dict) else None
        if keys != sorted(EXPECTED_SKEW_NU1):
            fails.append(f"[essential] {_DYN}/exact_skew/per_nu1: members {keys!r}, expected {list(EXPECTED_SKEW_NU1)}")
    return fails


# What the coverage audit must NOT see: source that is never printed -- the command words and keys of
# labels, references and citations, file names, environment names and float placements, lengths and
# graphics options, the spacing argument of a line break, the invisible argument of \phantom, and the
# command word \caption. Everything printed stays visible: caption text, a citation's printed pinpoint
# (its optional argument), theorem titles, math such as \in[0.05, 0.95], and prose that follows any of
# these commands on the same line. (The first version masked the rest of the line after \ref or
# \caption and every bracketed argument, which hid printed numbers from the gate.)
_AUDIT_KEY_CMDS = r"label|ref|eqref|pageref|autoref|cref|Cref|cite|citep|citet"
_A_ARG = r"[ \t]*\{[^{}\n]*\}"      # one mandatory argument on the same line, without braces inside
_A_OPT = r"(?:\[[^\]\n]*\])?"        # at most one optional argument on the same line
# One left-to-right scan over TeX tokens. A backslash starts a control sequence only where TeX reads
# one: "\\" (a line break, with its optional spacing argument) and "\%", "\,", "\{" ... are consumed
# as units first, so "Line.\\hspace{...}" is a line break followed by printed text and "\\%" starts a
# comment. Blanked, because never printed: comments; the keys of labels, references and citations (a
# citation's optional argument, its printed pinpoint, stays); and the exact arguments of the layout
# commands -- \hspace/\vspace (starred too) one, \phantom/\hphantom/\vphantom/\input/\include one,
# \includegraphics (starred too) an optional and one, \setlength/\addtolength two, \rule an optional
# and two --, environment names with a float placement, \addlinespace[..] and the word \caption.
# Anything else, including an argument the scan does not recognise, stays visible: it can add numbers
# to check, never hide one. (Earlier versions blanked the rest of a line after \ref or \caption, every
# following brace group, every starred form and 10^{b} unless "\times" stood right before it.)
_AUDIT_TOKEN = re.compile(
    r"(?P<verb>\\verb\*?(?P<vd>[^A-Za-z\s*])(?:(?!(?P=vd))[^\n])*(?P=vd))"
    r"|(?P<venv>\\begin\{(?P<vname>verbatim\*?|Verbatim|lstlisting)\}[\s\S]*?\\end\{(?P=vname)\})"
    r"|(?P<bs>\\\\(?:\[[^\]\n]*\])?)"
    r"|(?P<esc>\\[^A-Za-z])"
    r"|(?P<com>%[^\n]*)"
    r"|(?P<key>\\(?:" + _AUDIT_KEY_CMDS + r")(?![A-Za-z])(?P<kopt>\[[^\]\n]*\])?" + _A_ARG + r")"
    r"|(?P<lay>\\(?:hspace|vspace)(?![A-Za-z])\*?" + _A_ARG
    + r"|\\(?:phantom|hphantom|vphantom|input|include)(?![A-Za-z])" + _A_ARG
    + r"|\\includegraphics(?![A-Za-z])\*?" + _A_OPT + _A_ARG
    + r"|\\(?:setlength|addtolength)(?![A-Za-z])" + _A_ARG + _A_ARG
    + r"|\\rule(?![A-Za-z])" + _A_OPT + _A_ARG + _A_ARG
    + r"|\\(?:begin|end)\{[^{}\n]*\}(?:\[[!htbpH]+\])?"
    + r"|\\addlinespace(?![A-Za-z])\[[^\]\n]*\]"
    + r"|\\caption(?![A-Za-z]))"
    r"|(?P<cs>\\[A-Za-z]+)")


def _audit_scan(text: str) -> Tuple[str, List[int]]:
    """(the text with every non-printed construct blanked, newlines kept so line numbers are
    unchanged; the lines where a number continues across a comment). TeX drops a comment together
    with its line end, so digits on the two sides join into one printed number that the blanked text
    would split: such a continuation is reported as unsupported instead of being read as two numbers."""
    joins: List[int] = []

    def blank(s: str) -> str:
        return re.sub(r"[^\n]", " ", s)

    def repl(m: "re.Match") -> str:
        kind, s = m.lastgroup, m.group(0)
        if kind in ("esc", "cs", "verb", "venv"):   # verbatim text is printed as it stands
            return s
        if kind == "key":
            if m.group("kopt") is None:
                return blank(s)
            a, b = m.start("kopt") - m.start(), m.end("kopt") - m.start()
            return blank(s[:a]) + s[a:b] + blank(s[b:])
        if kind == "com":
            st, en = m.start(), m.end()
            k = en
            while True:   # TeX also drops the line ends of whole-line comments that follow
                nxt = re.match(r"\n[ \t]*%[^\n]*", text[k:])
                if not nxt:
                    break
                k += nxt.end()
            if st > 0 and (text[st - 1].isdigit() or text[st - 1] in ".-") and re.match(r"\n[ \t]*[\d.]", text[k:k + 64]):
                joins.append(text.count("\n", 0, st) + 1)
        return blank(s)

    return _AUDIT_TOKEN.sub(repl, text), joins


def _audit_mask(text: str) -> str:
    """The section-file text with every non-printed construct blanked (see _audit_scan)."""
    return _audit_scan(text)[0]


@functools.lru_cache(maxsize=16)
def _printed_text(raw: str) -> str:
    """The masked text of a section file for the quote gate. Only printed text can quote a number: the gate
    used to search the raw source, so a number left in a comment, in a label key or in a graphics option
    satisfied a pin. Runs of blanks are collapsed (the number grammar allows any spacing, and the blanked
    constructs would slow its scan); line breaks, hence line numbers, are kept. Kept per text, because the
    gate reads a file once per claim."""
    return re.sub(r"[ \t]{2,}", " ", _audit_mask(raw))


# An exponent: a braced integer with TeX spacing anywhere (sign and digits may be separated), or a
# single unbraced digit. Its digits are read with the spacing removed.
_A_EXP = r"(?:\{" + _SP + r"*([+-]?)" + _SP + r"*(\d(?:" + _SP + r"*\d)*)" + _SP + r"*\}|(\d))"
_BARE_POW_TAIL = re.compile(_SP + r"*\^" + _SP + r"*" + _A_EXP)
_TIMES_POW_TAIL = re.compile(_SP + r"*\\(?:times|cdot)" + _SP + r"*10" + _SP + r"*\^" + _SP + r"*" + _A_EXP)
# A power whose exponent carries a digit but is not an integer we can read (10^{-1.5}, ...) is reported,
# never split into a mantissa and fragments.
_POW_WITH_DIGIT = re.compile(_SP + r"*\^" + _SP + r"*(?:\{[^{}\n]*\d|\d)")
_TIMES_WITH_DIGIT = re.compile(_SP + r"*\\(?:times|cdot)" + _SP + r"*10" + _SP + r"*\^" + _SP + r"*(?:\{[^{}\n]*\d|\d)")
# Comparisons run in a dedicated decimal context, exact for the supported domain: at most 60
# significant digits and exponents up to 9999 in magnitude. A number outside it is reported as
# unsupported instead of being rounded or underflowed into a match.
_AUDIT_DEC = decimal.Context(prec=200, Emax=10 ** 6, Emin=-(10 ** 6))
_AUDIT_MAX_DIGITS, _AUDIT_MAX_EXP = 60, 9999


def _exp_text(m: "re.Match", first: int) -> str:
    sign, digits, single = m.group(first), m.group(first + 1), m.group(first + 2)
    return (sign + re.sub(r"\D", "", digits)) if digits is not None else single


def _printed_numbers(raw: str) -> Tuple[List[Tuple[int, "Decimal", "Decimal", bool, str]], List[Tuple[int, str]]]:
    r"""Every printed number of a section file's source as (line, value, half unit, in scope, as
    printed), and the unsupported constructs as (line, reason): a number continued across a comment
    (chain), a power of ten whose exponent is not a readable integer, a number outside the supported
    precision or range. A number is its COMPLETE signed value with its exponent -- \times 10^{b} or
    \cdot 10^{b}, or e-notation -- and a bare power of ten 10^{b} is 1 x 10^{b}, with TeX spacing
    anywhere in the power; a consumed span is never read again, so exponent digits are not numbers of
    their own. A number's line is the line of its sign or first digit (not of spacing before it).
    In scope: an exponent, two decimals, or three integer digits (an unsigned year excepted). This one
    iterator serves the coverage scan and the audit-line check."""
    masked, join_lines = _audit_scan(raw.replace("{,}", ""))
    starts = [0] + [i + 1 for i, ch in enumerate(masked) if ch == "\n"]
    problems: List[Tuple[int, str]] = [(ln, "a number continued across a comment (TeX joins the digits; write it on one line)")
                                       for ln in join_lines]
    out: List[Tuple[int, Decimal, Decimal, bool, str]] = []
    skip_to = -1
    for m in _NUM_RE.finditer(masked):
        if m.start() < skip_to:
            continue
        sign, core, e = m.group(1), m.group(2), (m.group(3) or m.group(4))
        pos = m.start(1) if sign else m.start(2)
        line = bisect.bisect_right(starts, pos)
        end = m.end()
        if not e:
            tail = _TIMES_POW_TAIL.match(masked, end)
            if tail:
                e, end = _exp_text(tail, 1), tail.end()
            elif _TIMES_WITH_DIGIT.match(masked, end):
                problems.append((line, f"a power of ten whose exponent is not a readable integer: {masked[pos:end + 24].strip()[:40]}"))
            elif core == "10":
                tail = _BARE_POW_TAIL.match(masked, end)
                if tail:
                    core, e, end = "1", _exp_text(tail, 1), tail.end()
                elif _POW_WITH_DIGIT.match(masked, end):
                    problems.append((line, f"a power of ten whose exponent is not a readable integer: {masked[pos:end + 24].strip()[:40]}"))
        skip_to = end
        dec = core.split(".")[1] if "." in core else ""
        exp = int(e) if e else 0
        if len(core.replace(".", "").lstrip("0")) > _AUDIT_MAX_DIGITS or abs(exp) > _AUDIT_MAX_EXP:
            problems.append((line, f"a number outside the supported precision or range: {masked[pos:end].strip()[:40]}"))
            continue
        in_scope = bool(e) or len(dec) >= 2 or ("." not in core and len(core) >= 3)
        if not e and not sign and re.fullmatch(r"(19|20)\d\d", core):
            in_scope = False
        val = Decimal(core).scaleb(exp, context=_AUDIT_DEC)
        if sign:
            val = val.copy_negate()
        out.append((line, val, Decimal(5).scaleb(exp - len(dec) - 1, context=_AUDIT_DEC), in_scope,
                    masked[pos:end].strip()))
    return out, problems


def _audit_number_lines(raw: str) -> Dict[int, List["Decimal"]]:
    """Printed values of a section file by the line on which each starts."""
    by_line: Dict[int, List[Decimal]] = {}
    for ln, val, _half, _scope, _shown in _printed_numbers(raw)[0]:
        by_line.setdefault(ln, []).append(val)
    return by_line


def _token_half_unit(tok: str) -> "Decimal":
    """Half a unit of an audit token's last printed digit, e-notation included (1.40e-3 -> 5e-6)."""
    c = tok.lstrip("+-").lower()
    mant, ex = (c.split("e") + ["0"])[:2] if "e" in c else (c, "0")
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return Decimal(5).scaleb(int(ex) - dec - 1, context=_AUDIT_DEC)


def _token_matches(val: "Decimal", tok: str) -> bool:
    """A printed value equals an audit token at the token's printed precision, with the same sign."""
    tv = Decimal(tok)
    return (tv < 0) == (val < 0) and _AUDIT_DEC.compare(_AUDIT_DEC.abs(_AUDIT_DEC.subtract(val, tv)), _token_half_unit(tok)) <= 0


def audit_coverage_gaps(audit: List[Dict[str, Any]], texdir: str, sig_tokens: Any) -> List[str]:
    r"""Numbers printed in the section files (see _printed_numbers) that no audit row of the same file
    names: a row's numeric token covers a printed number only with the same sign and a value equal to
    within half a unit of the printed last digit -- a mantissa alone, an unsigned value or another
    exponent covers nothing. A number continued across a comment is reported as unsupported.
    """
    cover: Dict[str, List[Decimal]] = {}
    for a in audit:
        vals = cover.setdefault(os.path.basename(str(a.get("tex") or "")), [])
        for tok in sig_tokens(a.get("literal", "")):
            try:
                vals.append(Decimal(tok))
            except Exception:  # noqa: BLE001 -- a malformed token covers nothing
                pass
    gaps: List[str] = []
    for base in AUDIT_FILES:
        path = os.path.join(texdir, base)
        if not os.path.isfile(path):
            gaps.append(f"{base}: missing")
            continue
        nums, problems = _printed_numbers(open(path).read())
        for ln, why in problems:
            gaps.append(f"{base}:{ln}: unsupported -- {why}")
        cov = cover.get(base, [])
        for ln, val, half, in_scope, shown in nums:
            if in_scope and not any((tv < 0) == (val < 0)
                                    and _AUDIT_DEC.compare(_AUDIT_DEC.abs(_AUDIT_DEC.subtract(tv, val)), half) <= 0
                                    for tv in cov):
                gaps.append(f"{base}:{ln}: {shown}")
    return gaps


def inventory_errors(claims_doc: Any, mode: str, scope: List[str]) -> List[str]:
    """Schema errors of the claims inventory, checked BEFORE evaluation: a
    JSON object with a non-empty 'claims' list of objects carrying unique string ids and string
    producer / producer_field / producer_kind / block, a non-empty 'coverage_audit' list, and -- in a
    full unscoped run -- at least one claim for every block of REQUIRED_CLAIM_BLOCKS."""
    if not isinstance(claims_doc, dict):
        return [f"[inventory] the claims file holds a {type(claims_doc).__name__}, not a JSON object"]
    errs: List[str] = []
    claims = claims_doc.get("claims")
    if not isinstance(claims, list) or not claims:
        errs.append("[inventory] 'claims' is missing, not a list, or empty")
        claims = claims if isinstance(claims, list) else []
    audit = claims_doc.get("coverage_audit")
    if not isinstance(audit, list) or not audit:
        errs.append("[inventory] 'coverage_audit' is missing, not a list, or empty")
    def finite(v: Any) -> bool:
        return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))
    for i, a in enumerate(audit if isinstance(audit, list) else []):
        if not isinstance(a, dict):
            errs.append(f"[inventory] coverage_audit row #{i} is a {type(a).__name__}, not an object")
            continue
        if not (isinstance(a.get("tex"), str) and os.path.basename(a["tex"]) in SECTION_FILES):
            errs.append(f"[inventory] coverage_audit row #{i}: 'tex' {a.get('tex')!r} is not a section file")
        lines = a.get("lines")
        if not (isinstance(lines, list) and lines and all(isinstance(x, int) and not isinstance(x, bool) and x >= 1 for x in lines)):
            errs.append(f"[inventory] coverage_audit row #{i}: 'lines' must be a non-empty list of line numbers")
        if not (isinstance(a.get("literal"), str) and a["literal"].strip()):
            errs.append(f"[inventory] coverage_audit row #{i}: 'literal' is missing or empty")
        disp = a.get("disposition")
        if not (isinstance(disp, str) and disp.startswith(("claim:", "exclusion:"))):
            errs.append(f"[inventory] coverage_audit row #{i}: 'disposition' must start with claim: or exclusion:")
        else:
            cats = re.findall(r"exclusion:([A-Za-z_]+)", disp)
            if not re.findall(r"claim:([A-Za-z0-9_]+)", disp) and not cats:
                errs.append(f"[inventory] coverage_audit row #{i}: disposition {disp[:40]!r} names no claim and no exclusion category")
            unknown = sorted(set(cats) - AUDIT_EXCLUSION_CATEGORIES)
            if unknown:
                errs.append(f"[inventory] coverage_audit row #{i}: unknown exclusion category {unknown}")
    seen: set = set()
    for i, c in enumerate(claims):
        if not isinstance(c, dict):
            errs.append(f"[inventory] claim #{i} is a {type(c).__name__}, not an object")
            continue
        cid = c.get("id")
        if not isinstance(cid, str) or not cid:
            errs.append(f"[inventory] claim #{i} has no string id")
            continue
        if cid in seen:
            errs.append(f"[inventory] duplicate claim id {cid!r}")
        seen.add(cid)
        for k in ("producer", "producer_field", "producer_kind", "block"):
            if not isinstance(c.get(k), str) or not c.get(k):
                errs.append(f"[inventory] claim {cid!r}: {k!r} is missing or not a string")
        # A claim says explicitly whether the manuscript quotes it; a tracked claim carries what the
        # manuscript checks need (a missing field used to switch those checks off silently).
        tracked = c.get("manuscript_tracked")
        if not isinstance(tracked, bool):
            errs.append(f"[inventory] claim {cid!r}: 'manuscript_tracked' must be true or false")
        elif tracked:
            if not finite(c.get("quoted_value")):
                errs.append(f"[inventory] claim {cid!r}: a tracked claim needs a finite 'quoted_value'")
            if not (finite(c.get("display_tol")) and float(c["display_tol"]) >= 0.0):
                errs.append(f"[inventory] claim {cid!r}: a tracked claim needs a finite 'display_tol' >= 0")
            if not (isinstance(c.get("tex"), str) and os.path.basename(c["tex"]) in SECTION_FILES):
                errs.append(f"[inventory] claim {cid!r}: a tracked claim needs 'tex' naming a section file")
            qi = c.get("quoted_in")
            if not (isinstance(qi, dict) and qi and all(isinstance(k, str) and os.path.basename(k) in SECTION_FILES
                                                        and isinstance(v, int) and not isinstance(v, bool) and v >= 1
                                                        for k, v in qi.items())):
                errs.append(f"[inventory] claim {cid!r}: a tracked claim needs 'quoted_in' pins (section file -> count >= 1)")
    if mode == "full" and not scope:
        present = {c.get("block") for c in claims if isinstance(c, dict)}
        missing = [b for b in REQUIRED_CLAIM_BLOCKS if b not in present]
        if missing:
            errs.append(f"[inventory] no claim for required block(s) {missing}")
    return errs


# ---------------------------------------------------------------------------
# Full-mode completeness + fallback gates
# ---------------------------------------------------------------------------

def completeness_gate(producers: Dict[str, Tuple[Dict[str, Any], str]], figdir: str,
                      blocks: List[str], claim_stems: List[str]) -> List[str]:
    """HARD-fail messages for any missing block / field / _taxonomy entry / figure /
    C0-C10 config. ``blocks`` is the scoped (or full) block set."""
    fails: List[str] = []
    for block in blocks:
        spec = BLOCK_SPEC[block]
        prod = producers.get(spec["stem"])
        if prod is None:
            fails.append(f"[completeness] block {block!r}: producer {spec['stem']!r} JSON absent")
            continue
        root = prod[0]
        for field, sel in spec["require"]:
            try:
                vals = resolve_scalars(root, field, sel)  # presence
                if not vals:
                    fails.append(f"[completeness] block {block!r}: field {field} (sel {sel}) resolves to nothing")
                # Presence is not a value. A required field must be a finite JSON
                # number (a bool is accepted only for *_ok flags); a string, a NaN or an
                # infinity in a required field fails the gate.
                for v in vals:
                    if isinstance(v, bool):
                        if not str(field).endswith("_ok"):
                            fails.append(f"[completeness] block {block!r}: field {field} is a bool, a number is required")
                    elif not isinstance(v, (int, float)):
                        fails.append(f"[completeness] block {block!r}: field {field} holds {v!r} ({type(v).__name__}), a finite number is required")
                    elif isinstance(v, float) and not math.isfinite(v):
                        fails.append(f"[completeness] block {block!r}: field {field} is {v!r}, a finite number is required")
            except SelectorError as exc:
                fails.append(f"[completeness] block {block!r}: field {field}: {exc}")
        tax = root.get("_taxonomy", {}) if isinstance(root, dict) else {}
        for tk in spec.get("taxonomy_keys", []):
            if tk not in tax:
                fails.append(f"[completeness] block {block!r}: missing _taxonomy entry {tk!r}")
        for fig in spec.get("figures", []):
            if not os.path.exists(os.path.join(figdir, fig)):
                fails.append(f"[completeness] block {block!r}: figure {fig} missing in FIGDIR")
        if block == "sensitivity":
            sens = root.get("sensitivity", {}) if isinstance(root, dict) else {}
            n_cfg = sum(1 for k in sens if not str(k).startswith("_"))
            need = spec.get("sens_configs", 11)
            if n_cfg != need:
                fails.append(f"[completeness] sensitivity: found {n_cfg} configs, need C0-C{need - 1} = {need}")
    return fails


def fallback_gate(producers: Dict[str, Tuple[Dict[str, Any], str]], blocks: List[str]) -> List[str]:
    """HARD-block messages for any sensitivity config with _fallback_used=true (full
    submission runs must use nominal C0-C10 params)."""
    hard: List[str] = []
    if "sensitivity" not in blocks:
        return hard
    prod = producers.get(_DYN)
    if prod is None:
        return hard
    sens = prod[0].get("sensitivity", {}) if isinstance(prod[0], dict) else {}
    for cfg, body in sens.items():
        if str(cfg).startswith("_"):
            continue
        if isinstance(body, dict) and body.get("_fallback_used") is True:
            hard.append(f"[fallback] sensitivity config {cfg!r} used a wall-clock fallback "
                        f"(_fallback_used=true) -- a submission rerun must use nominal params")
    return hard


def qp_source_check(qp_root: Dict[str, Any]) -> List[str]:
    """Verify the QP diagnostics' (mu0, nu1) match the env-resolved
    source regime. Reuses ``_resolve_qp_source_regime()`` from the QP module rather
    than a naive ``qp_root['nu1'] == source['nu1']`` (the current source JSON has NO
    'nu1' key). It runs in full mode only: a record whose provenance says smoke fails, and nu1 is
    always compared -- for status default_nu1 / missing with the resolver's documented default
    nu1 (4.90) (the default and the smoke label used to skip the comparison). Both
    parameters must be finite JSON numbers (a string or NaN used to bypass the comparison);
    ``mu0`` is compared against the resolver's u5_mu0, not any hardcoded echo. Best-effort:
    an import failure is a note, not a crash."""
    issues: List[str] = []
    here = os.path.dirname(os.path.abspath(__file__))
    qp_py = os.path.join(here, _QP + ".py")
    if not os.path.exists(qp_py):
        return [f"[qp-source] QP module not found for source-regime check: {qp_py}"]
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(_QP + "_forvalidate", qp_py)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    except Exception as exc:  # noqa: BLE001 -- best-effort import
        return [f"[qp-source] could not import QP module ({exc}); skipping source-regime check"]
    if not hasattr(mod, "_resolve_qp_source_regime"):
        return ["[qp-source] QP module lacks importable _resolve_qp_source_regime()"]
    try:
        u5_mu0, u5_nu1, _src, status = mod._resolve_qp_source_regime()
    except Exception as exc:  # noqa: BLE001
        return [f"[qp-source] _resolve_qp_source_regime() raised {exc}"]
    # Prefer the status the QP run itself recorded (smoke/default_nu1/missing/read).
    rec_status = (qp_root.get("_provenance") or {}).get("qp_source_status", status)
    def _finite_real(v: Any) -> bool:
        return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))
    out_mu0 = qp_root.get("mu0")
    if not _finite_real(out_mu0):
        issues.append(f"[qp-source] QP output mu0 is {out_mu0!r}; a finite number is required")
    elif abs(float(out_mu0) - float(u5_mu0)) > 1e-9:
        issues.append(f"[qp-source] mu0 mismatch: QP output {out_mu0} vs resolver u5_mu0 {u5_mu0}")
    # This check runs in full mode only: a record that says it came from a smoke run is not a full
    # record, and the exemption is never taken from producer metadata.
    if rec_status == "smoke":
        issues.append("[qp-source] the QP record says it came from a smoke run (qp_source_status='smoke'); "
                      "a full validation needs the full run's record")
    out_nu1 = qp_root.get("nu1")
    if not _finite_real(out_nu1):
        issues.append(f"[qp-source] QP output nu1 is {out_nu1!r}; a finite number is required")
    elif abs(float(out_nu1) - float(u5_nu1)) > 1e-9:
        issues.append(f"[qp-source] nu1 mismatch: QP output {out_nu1} vs resolver u5_nu1 "
                      f"{u5_nu1} (status={rec_status})")
    return issues


# ---------------------------------------------------------------------------
# Report / VALJSON writers
# ---------------------------------------------------------------------------

def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    return "" if v is None else str(v)


def write_report(path: str, header: Dict[str, Any], rows: List[Dict[str, Any]],
                 audit: List[Dict[str, Any]], gate: Dict[str, Any]) -> None:
    L: List[str] = []
    L.append("# Design-G validation report\n")
    L.append(f"- generated: {header['timestamp']}")
    L.append(f"- mode: **{header['mode']}**   scope: {header['scope'] or 'ALL'}")
    L.append(f"- claims file: {header['claims_in']}")
    L.append(f"- producer JSON dir: {header['data_in']}   figure dir: {header['figdir']}")
    L.append(f"- producers found: {', '.join(header['producers_found']) or '(none)'}")
    L.append(f"- overall (structural): **{header['overall']}** (exit {header['exit_code']})")
    L.append(f"- manuscript readiness: **{header.get('readiness', 'n/a')}**"
             f"{' [strict mode: a required change fails the run]' if header.get('strict') else ''}\n")

    L.append("## Per-claim checks\n")
    L.append("| kind | name | actual | expected | tol | status | manuscript_change | ci | tex_ref |")
    L.append("|------|------|--------|----------|-----|--------|-------------------|----|---------|")
    for r in rows:
        mc = r.get("manuscript_change_required")
        mc_s = "-" if mc is None else ("YES" if mc else "no")
        L.append("| {kind} | {name} | {actual} | {expected} | {tol} | {status} | {mc} | {ci} | {tex} |".format(
            kind=r.get("kind"), name=r.get("name"), actual=_fmt(r.get("actual")),
            expected=_fmt(r.get("expected")), tol=_fmt(r.get("tol")), status=r.get("status"),
            mc=mc_s, ci=_fmt(r.get("ci_method")), tex=r.get("tex_ref")))
    L.append("")

    if gate.get("manuscript_changes"):
        L.append("## Manuscript-change triggers (non-blocking)\n")
        for m in gate["manuscript_changes"]:
            L.append(f"- `{m['name']}` at {m['tex_ref']}: quoted {m['quoted']} vs produced {m['actual']}")
        L.append("")

    if gate.get("blocking"):
        L.append("## BLOCKING failures (full mode)\n")
        for m in gate["blocking"]:
            L.append(f"- {m}")
        L.append("")

    # Coverage audit.
    L.append("## Claim-coverage audit (every Design-G number -> claim id | exclusion)\n")
    n_claim = sum(1 for a in audit if str(a.get("disposition", "")).startswith("claim:"))
    n_excl = sum(1 for a in audit if str(a.get("disposition", "")).startswith("exclusion:"))
    L.append(f"- audit rows: {len(audit)}  ({n_claim} -> claim id, {n_excl} -> documented exclusion)")
    L.append(f"- audit integrity: {gate.get('audit_integrity', 'n/a')}\n")
    L.append("| tex | line(s) | literal | disposition |")
    L.append("|-----|---------|---------|-------------|")
    for a in audit:
        lines = a.get("lines") or ([a["line"]] if a.get("line") else [])
        L.append("| {tex} | {ln} | {lit} | {disp} |".format(
            tex=os.path.basename(str(a.get("tex", ""))), ln=",".join(str(x) for x in lines),
            lit=str(a.get("literal", "")).replace("|", "\\|"),
            disp=str(a.get("disposition", "")).replace("|", "\\|")))
    L.append("")

    L.append("## Gates (full mode)\n")
    L.append(f"- claims inventory: {gate.get('inventory_status', 'n/a')}")
    for f in gate.get("inventory", []):
        L.append(f"  - {f}")
    L.append(f"- essential success flags: {gate.get('essential_flags_status', 'n/a')}")
    for f in gate.get("essential_flags", []):
        L.append(f"  - {f}")
    L.append(f"- completeness: {gate.get('completeness_status', 'n/a')}")
    for f in gate.get("completeness", []):
        L.append(f"  - {f}")
    L.append(f"- fallback: {gate.get('fallback_status', 'n/a')}")
    for f in gate.get("fallback", []):
        L.append(f"  - {f}")
    L.append(f"- qp_source (mu0,nu1 regime): {gate.get('qp_source_status', 'n/a')}")
    for f in gate.get("qp_source", []):
        L.append(f"  - {f}")
    L.append(f"- manifest: {gate.get('manifest_status', 'n/a')}")
    for f in gate.get("manifest", []):
        L.append(f"  - {f}")
    L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def _resolve_scope(raw: Optional[str]) -> List[str]:
    if not raw:
        return []
    out: List[str] = []
    for tok in raw.replace(",", " ").split():
        block = _U_ALIAS.get(tok.upper(), tok)
        if block not in BLOCK_SPEC:
            sys.exit(f"[scope] unknown DESIGNG_VALIDATE_SCOPE token {tok!r} "
                     f"(known blocks: {sorted(BLOCK_SPEC)} or U1-U5)")
        out.append(block)
    return out


def _write_failed_record(report_out: str, valjson_out: str, exc: BaseException, exit_code: int = 1) -> None:
    """Replace the previous report and VALJSON by FAILED ones when the run aborts before its normal
    writers run (a truncated producer JSON, a null `_taxonomy` or `sensitivity` block,
    or a list where a dict is expected raised before any report was written, so the previous
    SUCCESSFUL artifacts stayed on disk advertising SUBMISSION_READY)."""
    import traceback
    msg = f"{type(exc).__name__}: {exc}"
    body = {"_meta": {"task": "Design-G validator", "generated": datetime.now(timezone.utc).isoformat(),
                      "aborted": True, "error": msg, "traceback": traceback.format_exc()[-4000:]},
            "overall": "FAIL", "exit_code": exit_code, "readiness": "FAIL (run aborted before evaluation completed)",
            "strict_mode": os.environ.get("DESIGNG_VALIDATE_STRICT", "").strip().lower() in ("1", "true", "yes"),
            "summary": {"n_claims_total": None, "n_evaluated": 0, "n_pass": 0, "n_fail": 0, "n_manuscript_change": None},
            "checks": [], "gates": {}, "manuscript_changes": []}
    for path, text in ((valjson_out, json.dumps(body, indent=2)),
                       (report_out, f"# Design-G validation report\n\n- overall (structural): **FAIL** (exit {exit_code})\n"
                                    f"- manuscript readiness: **FAIL** -- the run ABORTED before evaluation completed: {msg}\n\n"
                                    "This report replaces the previous one on purpose: a validator that crashes must not leave a "
                                    "stale PASS behind.\n")):
        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            with open(path, "w") as fh:
                fh.write(text)
        except OSError:
            pass


_STATE = {"reports_written": False}


def main(argv=None) -> int:
    """Guarded entry point. A run that ends with a non-zero status
    before its normal writers have produced the report and the VALJSON -- an exception, an
    invalid DESIGNG_VALIDATE_MODE, an unknown DESIGNG_VALIDATE_SCOPE token, a missing required
    variable -- replaces the previous report and VALJSON by FAILED records whenever the two
    output paths are known, so a stale PASS never outlives the run that superseded it. A run
    that evaluated the claims and wrote its own FAIL report keeps that report."""
    _STATE["reports_written"] = False
    report_out, valjson_out = os.environ.get("DESIGNG_REPORT_OUT"), os.environ.get("DESIGNG_VALJSON_OUT")

    def _record(reason: BaseException, code: int) -> None:
        if report_out and valjson_out and not _STATE["reports_written"]:
            _write_failed_record(report_out, valjson_out, reason, exit_code=code)

    try:
        rc = _main_unguarded(argv)
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
        if code != 0:
            _record(RuntimeError(f"configuration abort before evaluation: {exc.code}"), code)
        raise
    except Exception as exc:  # noqa: BLE001 -- the whole point is to record it
        sys.stderr.write(f"[validate] ABORT: {type(exc).__name__}: {exc}\n")
        if report_out and valjson_out:
            _write_failed_record(report_out, valjson_out, exc, exit_code=1)
        return 1
    if rc != 0:
        _record(RuntimeError(f"configuration abort before evaluation (exit {rc}); "
                             "see the [validate] ABORT line on stderr"), rc)
    return rc


def _main_unguarded(argv=None) -> int:
    # --- 6-env-var abort (every invocation must pass the complete set) ---
    missing = [v for v in REQUIRED_ENV if not os.environ.get(v)]
    if missing:
        sys.stderr.write(
            "[validate] ABORT: required env var(s) unset: " + ", ".join(missing) + "\n"
            "  every validator invocation must pass all six of: " + ", ".join(REQUIRED_ENV) + "\n")
        return 2

    mode = os.environ.get("DESIGNG_VALIDATE_MODE", "full").strip().lower()
    if mode not in ("smoke", "full"):
        sys.stderr.write(f"[validate] ABORT: DESIGNG_VALIDATE_MODE must be smoke|full, got {mode!r}\n")
        return 2
    scope = _resolve_scope(os.environ.get("DESIGNG_VALIDATE_SCOPE"))

    data_in = os.environ["DESIGNG_VALIDATE_IN"]
    figdir = os.environ["DESIGNG_VALIDATE_FIGDIR"]
    claims_in = os.environ["DESIGNG_CLAIMS_IN"]
    report_out = os.environ["DESIGNG_REPORT_OUT"]
    valjson_out = os.environ["DESIGNG_VALJSON_OUT"]
    manifest = os.environ["DESIGNG_MANIFEST"]
    global _TEXDIR
    _TEXDIR = os.environ.get("DESIGNG_TEXDIR") or _BUNDLE_ROOT

    with open(claims_in) as fh:
        claims_doc = json.load(fh)
    # The inventory is schema-checked before anything is evaluated; the
    # well-formed part is still evaluated so the report shows what was found.
    inv_errs = inventory_errors(claims_doc, mode, scope)
    _raw_claims = claims_doc.get("claims", []) if isinstance(claims_doc, dict) else []
    claims = [c for c in (_raw_claims if isinstance(_raw_claims, list) else [])
              if isinstance(c, dict) and isinstance(c.get("id"), str) and c.get("id")]
    _raw_audit = claims_doc.get("coverage_audit", []) if isinstance(claims_doc, dict) else []
    audit = [a for a in (_raw_audit if isinstance(_raw_audit, list) else []) if isinstance(a, dict)]

    # Audit integrity: every 'claim:<id>' reference resolves to a known claim id.
    known_ids = {c["id"] for c in claims}
    import re
    dangling = sorted({m for a in audit for m in re.findall(r"claim:([A-Za-z0-9_]+)", str(a.get("disposition", "")))
                       if m not in known_ids})
    # A row that says a number is "claim-covered" must cite a claim that actually
    # checks that section file (its tex or one of its quoted_in pins); otherwise the audit
    # claims coverage that no executable check provides.
    claim_by_id = {c["id"]: c for c in claims}
    unpinned: List[str] = []
    untracked: List[str] = []   # coverage attributed to a claim whose manuscript check is switched off
    for a in audit:
        atex = os.path.basename(str(a.get("tex") or ""))
        for cid in re.findall(r"claim:([A-Za-z0-9_]+)", str(a.get("disposition", ""))):
            c = claim_by_id.get(cid)
            if c is None or not atex:
                continue
            if c.get("manuscript_tracked") is not True:
                untracked.append(f"{atex}:{a.get('lines')} -> {cid}")
                continue
            pinned = set(c.get("quoted_in") or {}) | ({os.path.basename(str(c["tex"]))} if c.get("tex") else set())
            if atex not in pinned:
                unpinned.append(f"{atex}:{a.get('lines')} -> {cid}")
    # The line numbers an audit row lists are re-verified against the
    # CURRENT section files -- every listed line must still carry one of the row's numeric tokens
    # (>= 3 significant characters; the same complete-number grammar as the quote gate). A stale
    # row therefore fails instead of passing silently.
    stale: List[str] = []
    # Placement rule (2026-11-10): a row keeps a line only where a printed number WITHIN THE COVERAGE SCOPE equals one
    # of the row's tokens at the finer of the two printed precisions. A line that carries a row's number only because
    # another number rounds to a coarser token (0.0040 against the token 0.00), or only through a number outside the
    # scope (the 0 of a subscript), is reported as a weak placement and fails like a stale one.
    weak: List[str] = []
    _tok_re = re.compile(r"-?\d+(?:\.\d+)?(?:[eE]-?\d+)?")
    def _sig_tokens(lit: str) -> List[str]:
        """Numeric tokens of a literal with >= 3 significant characters, SIGN KEPT: a token
        matches only a number of the same sign (an unsigned token means a positive number)."""
        out = []
        for t in _tok_re.findall(str(lit)):
            core = t.lstrip("-"); dec = core.split(".")[1] if "." in core else ""
            if "e" in core.lower() or len(dec) >= 2 or ("." not in core and len(core) >= 3):
                out.append(t)
        return out
    _tex_cache: Dict[str, List[str]] = {}
    _num_cache: Dict[str, Dict[int, List[Decimal]]] = {}
    _full_cache: Dict[str, Dict[int, List[Tuple[Decimal, Decimal, bool]]]] = {}
    uncheckable: List[str] = []
    for a in audit:
        toks = _sig_tokens(a.get("literal", ""))
        if not _TEXDIR:
            continue
        atex = os.path.basename(str(a.get("tex") or ""))
        apath = os.path.join(_TEXDIR, atex)
        if atex not in _tex_cache:   # only a real file is read (a row with no tex names the directory)
            raw = open(apath).read() if (atex and os.path.isfile(apath)) else ""
            _tex_cache[atex] = raw.replace("{,}", "").split("\n") if raw else []
            _num_cache[atex] = _audit_number_lines(raw) if raw else {}
            _full_cache[atex] = {}
            for _ln, _v, _half, _scope, _shown in (_printed_numbers(raw)[0] if raw else []):
                _full_cache[atex].setdefault(_ln, []).append((_v, _half, _scope))
        lines_txt = _tex_cache[atex]
        if not toks:
            # A literal with no checkable number may only document a symbol written out, never claim
            # coverage; its listed lines must still exist.
            disp = str(a.get("disposition", ""))
            cats = set(re.findall(r"exclusion:([A-Za-z_]+)", disp))
            if re.findall(r"claim:([A-Za-z0-9_]+)", disp) or not cats or not cats <= AUDIT_UNCHECKABLE:
                uncheckable.append(f"{atex}:{a.get('lines')} {str(a.get('literal'))[:40]!r}")
            for ln in a.get("lines") or []:
                if not isinstance(ln, int) or ln < 1 or ln > len(lines_txt):
                    stale.append(f"{atex}:{ln} (out of range)")
            continue
        found_tokens: set = set()
        for ln in a.get("lines") or []:
            if not isinstance(ln, int) or ln < 1 or ln > len(lines_txt):
                stale.append(f"{atex}:{ln} (out of range)"); continue
            # The printed values starting on this line, read by the same iterator as the coverage scan
            # (complete signed values, exponents consumed, comments and keys masked). A token matches at
            # its own printed precision (a literal may legitimately round what the line prints, e.g.
            # -0.6931 for a printed -0.69314); a line edited from 3.958 to 3.956 still matches a coarser
            # token 3.96 of the same row, which the direct claim on the finer value and the all-tokens
            # rule below catch. Every token found on this line is recorded.
            ok = False
            for v in _num_cache[atex].get(ln, []):
                for t in toks:
                    if _token_matches(v, t):
                        ok = True; found_tokens.add(t)
            if not ok:
                stale.append(f"{atex}:{ln} <- {str(a.get('literal'))[:40]!r}")
            elif not any(scope_ok and (Decimal(t) < 0) == (v < 0)
                         and _AUDIT_DEC.compare(_AUDIT_DEC.abs(_AUDIT_DEC.subtract(Decimal(t), v)),
                                                min(half, _token_half_unit(t))) <= 0
                         for v, half, scope_ok in _full_cache[atex].get(ln, []) for t in toks):
                weak.append(f"{atex}:{ln} <- {str(a.get('literal'))[:40]!r}")
        # Every numeric token of the row's literal must be found on at least one listed
        # line; a row whose literal names a number no listed line carries is stale.
        missing_tok = [t for t in toks if t not in found_tokens]
        if missing_tok and (a.get("lines") or []):
            stale.append(f"{atex}:{a.get('lines')} token(s) {missing_tok} of {str(a.get('literal'))[:40]!r} on no listed line")
    # Full mode: every number the section files print must be accounted for by an audit row.
    uncovered: List[str] = audit_coverage_gaps(audit, _TEXDIR, _sig_tokens) if (mode == "full" and _TEXDIR) else []
    msgs = []
    if dangling:
        msgs.append(f"DANGLING claim refs: {dangling}")
    if unpinned:
        msgs.append(f"UNPINNED coverage: {len(unpinned)} row(s) cite a claim that does not check that file: {unpinned[:6]}")
    if untracked:
        msgs.append(f"UNTRACKED coverage: {len(untracked)} row(s) cite a claim whose manuscript check is off: {untracked[:6]}")
    if stale:
        msgs.append(f"STALE audit lines: {len(stale)} listed line(s) no longer carry the row's number: {stale[:6]}")
    if weak:
        msgs.append(f"WEAK audit placements: {len(weak)} listed line(s) print no number of the coverage scope equal to a "
                    f"token of their row at the printed precision: {weak[:6]}")
    if uncheckable:
        msgs.append(f"UNCHECKABLE rows: {len(uncheckable)} row(s) carry no checkable number and are not symbolic exclusions: {uncheckable[:6]}")
    if uncovered:
        msgs.append(f"UNCOVERED numbers: {len(uncovered)} number(s) printed in the section files are in no audit row: {uncovered[:6]}")
    audit_integrity = "OK" if not msgs else "; ".join(msgs)
    dangling = dangling or unpinned or untracked or stale or weak or uncheckable or uncovered  # all block below

    stems = sorted({c.get("producer") for c in claims if c.get("producer")})
    producers = load_producers(data_in, list(stems), mode)

    # Which claims to evaluate: scope-filter in both modes when scope is set.
    def in_scope(c: Dict[str, Any]) -> bool:
        return (not scope) or (c.get("block") in scope)

    def _safe_eval(c):
        # A malformed value that raised inside evaluate_claim aborted the run before
        # any report was written, leaving the previous SUCCESSFUL report and VALJSON on disk.
        # Every exception is now a SCHEMA_ERROR hard-fail row, so the run completes and writes
        # a FAILED report that replaces the stale one.
        try:
            return evaluate_claim(c, producers, mode)
        except Exception as exc:  # noqa: BLE001 -- converted into a blocking row
            return {"kind": c.get("producer_kind"), "name": c.get("id", "?"), "producer": c.get("producer"),
                    "producer_field": c.get("producer_field"), "selector": c.get("producer_selector"),
                    "block": c.get("block"), "manuscript_tracked": bool(c.get("manuscript_tracked")),
                    "tex_ref": f"{c.get('tex')}:{c.get('line')}" if c.get("tex") else "(untracked)",
                    "actual": None, "expected": None, "tol": None, "status": "SCHEMA_ERROR",
                    "manuscript_change_required": None, "ci_method": None, "hard_fail": True,
                    "blocking_fail": False, "detail": f"{type(exc).__name__}: {exc}", "manuscript_quote": None}
    rows = [_safe_eval(c) for c in claims if in_scope(c)]

    # Collect gate outcomes.
    blocking: List[str] = []
    manuscript_changes: List[Dict[str, Any]] = []
    for r in rows:
        if r.get("hard_fail"):
            blocking.append(f"{r['name']}: {r['status']} -- {r.get('detail')}")
        if r.get("blocking_fail"):
            blocking.append(f"{r['name']}: anchor FAIL (actual {r.get('actual')} vs {r.get('expected')} +/- {r.get('tol')})")
        if r.get("manuscript_change_required"):
            manuscript_changes.append({"name": r["name"], "tex_ref": r["tex_ref"],
                                       "quoted": next((c.get("quoted_value") for c in claims if c["id"] == r["name"]), None),
                                       "actual": r.get("actual")})
    if dangling:
        blocking.append(f"coverage-audit integrity: {audit_integrity}")

    gate: Dict[str, Any] = {
        "blocking": blocking, "manuscript_changes": manuscript_changes,
        "audit_integrity": audit_integrity,
        "completeness": [], "fallback": [], "manifest": [],
        "inventory": inv_errs, "inventory_status": "OK" if not inv_errs else f"{len(inv_errs)} error(s)",
    }
    blocking.extend(inv_errs)
    # Essential success flags, both modes: the blocks in scope, or in a full
    # unscoped run every required block, never a set derived from the inventory under check.
    flag_blocks = scope if scope else (REQUIRED_CLAIM_BLOCKS if mode == "full"
                                       else sorted({c.get("block") for c in claims}))
    gate["essential_flags"] = essential_flags_gate(producers, flag_blocks, mode)
    gate["essential_flags_status"] = "OK" if not gate["essential_flags"] else f"{len(gate['essential_flags'])} FAILED"
    blocking.extend(gate["essential_flags"])

    # Manifest verification: report in both modes; only BLOCK in full.
    gate["manifest"] = verify_manifest(manifest, _TEXDIR, figdir)
    gate["manifest_status"] = "OK" if not gate["manifest"] else f"{len(gate['manifest'])} issue(s)"
    # Manuscript-quote check: a tracked number absent from the text blocks in full mode.
    quote_missing = [quote_shortfall_message(r["name"], next((c.get("quoted_value") for c in claims if c["id"] == r["name"]), None),
                                             r["manuscript_quote"])
                     for r in rows if r.get("manuscript_quote") and r["manuscript_quote"]["status"] != "FOUND"]
    gate["manuscript_quotes"] = quote_missing
    gate["manuscript_quotes_status"] = "OK" if not quote_missing else f"{len(quote_missing)} tracked number(s) not found in the manuscript"
    if mode == "full" and quote_missing:
        blocking.extend(quote_missing)

    if mode == "full":
        blocks = scope if scope else list(BLOCK_SPEC)   # never derived from the inventory
        gate["completeness"] = completeness_gate(producers, figdir, blocks, list(stems))
        gate["completeness_status"] = "OK" if not gate["completeness"] else f"{len(gate['completeness'])} missing"
        gate["fallback"] = fallback_gate(producers, blocks)
        gate["fallback_status"] = "OK" if not gate["fallback"] else f"{len(gate['fallback'])} HARD-BLOCK"
        # QP (mu0, nu1) source-regime consistency (only when in scope + present).
        gate["qp_source"] = qp_source_check(producers[_QP][0]) if ("quasiperiodic" in blocks and _QP in producers) else []
        gate["qp_source_status"] = "OK" if not gate["qp_source"] else f"{len(gate['qp_source'])} issue(s)"
    else:
        gate["completeness_status"] = "SKIPPED_SMOKE (coverage-of-present only)"
        gate["fallback_status"] = "SKIPPED_SMOKE"
        gate["qp_source"] = []
        gate["qp_source_status"] = "SKIPPED_SMOKE"

    # --- Exit-code decision ---
    exit_code = 0
    if any(r.get("hard_fail") for r in rows) or dangling or inv_errs or gate["essential_flags"]:
        exit_code = 1  # schema / selector / taxonomy / audit / inventory / success-flag failure, BOTH modes
    if mode == "full":
        if (any(r.get("blocking_fail") for r in rows) or gate["completeness"]
                or gate["fallback"] or gate["manifest"] or gate.get("qp_source")
                or gate["manuscript_quotes"]
                or any(r.get("status") == "ABSENT" for r in rows)):
            exit_code = 1
    # Separate the STRUCTURAL verdict from manuscript READINESS.
    # ``overall`` is the structural verdict (schema, anchors, gates); ``readiness`` says whether
    # any tracked claim still requires a manuscript change. DESIGNG_VALIDATE_STRICT=1 makes a
    # required change fail the run (exit 1), so a report with CHANGE_REQUIRED rows cannot be
    # presented as an unconditional submission-ready PASS.
    strict = os.environ.get("DESIGNG_VALIDATE_STRICT", "").strip().lower() in ("1", "true", "yes")
    n_change = len(manuscript_changes)
    structural_exit = exit_code   # every row, inventory, audit and gate failure, before strict mode adds its own
    if strict and mode == "full" and n_change:
        exit_code = 1
    overall = "PASS" if exit_code == 0 else "FAIL"
    if mode != "full":
        readiness = "SMOKE_ONLY (no readiness statement in smoke mode)"
    elif structural_exit != 0:
        # A failed gate is a failure whatever else is pending. (In strict mode, with a required manuscript
        # change and no failing claim row, a failed gate used to leave the label CHANGES_REQUIRED.)
        readiness = "FAIL"
    elif n_change:
        readiness = (f"CHANGES_REQUIRED ({n_change} tracked claim(s) disagree with the manuscript)"
                     + (f" [scoped to {','.join(scope)}]" if scope else ""))
    elif scope:
        # A run restricted by DESIGNG_VALIDATE_SCOPE checked only part of the
        # submission, so it never states whole-submission readiness.
        readiness = f"SCOPED_PASS ({','.join(scope)}): not a whole-submission readiness statement"
    else:
        readiness = "SUBMISSION_READY"

    header = {
        "timestamp": datetime.now(timezone.utc).isoformat(), "mode": mode,
        "scope": ",".join(scope) if scope else None, "claims_in": claims_in,
        "data_in": data_in, "figdir": figdir,
        "producers_found": sorted(producers.keys()), "overall": overall, "exit_code": exit_code,
        "readiness": readiness, "strict": strict,
    }

    # Emit VALJSON (machine-readable).
    valjson = {
        "_meta": {"task": "Design-G validator", "generated": header["timestamp"],
                  "mode": mode, "scope": header["scope"], "claims_in": claims_in,
                  "data_in": data_in, "figdir": figdir, "manifest": manifest,
                  "producers_found": header["producers_found"]},
        "overall": overall, "exit_code": exit_code,
        "readiness": readiness, "strict_mode": strict,
        "summary": {
            "n_claims_total": len(claims), "n_evaluated": len(rows),
            "n_blocking_kinds": sum(1 for r in rows if r.get("kind") in BLOCKING_KINDS),
            "n_pass": sum(1 for r in rows if r.get("status") == "PASS"),
            "n_fail": sum(1 for r in rows if r.get("status") == "FAIL"),
            "n_skipped_smoke": sum(1 for r in rows if r.get("status") == "SKIPPED_SMOKE"),
            "n_manuscript_change": len(manuscript_changes),
            "n_audit_rows": len(audit),
            "n_audit_to_claim": sum(1 for a in audit if str(a.get("disposition", "")).startswith("claim:")),
            "n_audit_exclusion": sum(1 for a in audit if str(a.get("disposition", "")).startswith("exclusion:")),
        },
        "ci_methods": {r["name"]: r.get("ci_method") for r in rows if r.get("ci_method")},
        "checks": rows,
        "gates": {k: gate.get(k) for k in ("inventory", "inventory_status", "essential_flags",
                                           "essential_flags_status", "completeness", "completeness_status",
                                           "fallback", "fallback_status", "qp_source", "qp_source_status",
                                           "manifest", "manifest_status", "manuscript_quotes",
                                           "manuscript_quotes_status", "audit_integrity")},
        "manuscript_changes": manuscript_changes,
    }
    os.makedirs(os.path.dirname(os.path.abspath(valjson_out)), exist_ok=True)
    with open(valjson_out, "w") as fh:
        json.dump(valjson, fh, indent=2)

    os.makedirs(os.path.dirname(os.path.abspath(report_out)), exist_ok=True)
    write_report(report_out, header, rows, audit, gate)
    _STATE["reports_written"] = True

    print(f"[validate] mode={mode} scope={header['scope'] or 'ALL'} -> structural={overall} "
          f"readiness={readiness} strict={'on' if strict else 'off'} (exit {exit_code}); "
          f"{len(rows)} checks, {len(manuscript_changes)} manuscript-change trigger(s); "
          f"report={report_out} valjson={valjson_out}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
