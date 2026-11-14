# Design-G reproducibility bundle (CPU)

This bundle reproduces the Design-G numerical results reported in the manuscript
(the transverse exponent lambda_perp = ln(nu1/4), the extinction threshold,
terminal extinction fractions, sensitivity intervals, the exact skew-product
exponent, the mu0 sweep and its onset ensemble, the invasion-exponent comparison)
and the ergodic-quadrature checks of the Supplementary Material (arcsine moments, equilibrium potential
-ln4, sub-action ln(3/4), periodic enumeration, I-family), and the
quasiperiodic illustration at mu0 = 3.5 (Figures 13 and 14 of Section 6: phase portraits and spectrum; the U5
invariant-circle diagnostics of Supplementary Section S10). Three validated blocks are retained from earlier
versions without a figure in the manuscript: the finite-amplitude decay diagnostic and the exact-skew block
(the validator's completeness gate still requires their two PNGs; the exact-skew numbers are those of Table S4),
and the family-unification block, which has no quotation in the manuscript. All of these are numerical
consistency checks, not proofs.

## Requirements

Python 3.10.x and the PINNED packages in `requirements.txt` (numpy 1.26.4,
numba 0.61.2, matplotlib 3.9.2). The preflights of `designG_dynamics_v4.py` and
`designG_quasiperiodic_v2.py` ENFORCE Python 3.10.x, NumPy 1.26.4 and Numba 0.61.2 (numba
is version-sensitive) and abort on a mismatch; Matplotlib 3.9.2 is pinned for
reproduction but is not checked by these preflights, and the other current producers do not
check versions, so the pinned environment is required for a faithful reproduction
whichever script is run.

From the bundle root (the directory that holds `codes/`):

    python3.10 -m venv .venv
    . .venv/bin/activate
    pip install -r codes/requirements.txt

Validated environment: Python 3.10.12, numpy 1.26.4, numba 0.61.2,
matplotlib 3.9.2.

## Contents

Producers (write `*_values.json` + figures):

    designG_dynamics_v4.py            # U1-U4: lambda_perp=ln(nu1/4), TERMINAL extinction fractions,
                                      #        coexistence/Jury band, sensitivity, finite-amp, exact-skew
    designG_ergodic_quadrature_v2.py  # arcsine moments, equilibrium potential -ln4,
                                      #        sub-action ln(3/4), periodic enumeration, I-family
    designG_bifurcation_v2.py         # mu0 sweep at nu1=2.1: Figure 10 (prey, predator, the orbit's exponents) + the companion M0 exponents (Figure S3)
    designG_onset_ensemble_v1.py      # 24-seed x 3-resolution onset ensemble (the onset numbers of Section 6 and Supplementary Section S12)
    designG_invasion_exponents_v1.py  # closed-form invasion exponents vs nu1 (arcsine, period-3,
                                      #        mu0=2 control, mean-variance bound): the controlled-comparison figure
    designG_stable_leaf_v1.py         # stable-leaf / matched-product ILLUSTRATION (2026-10-21): the local stable
                                      #        leaf through an axis point as the fixed point of the manuscript's
                                      #        Lyapunov--Perron operator, and cor:rare-predator's matched product;
                                      #        not an instance of the theorem (see the section below)
    designG_schematics_v2.py          # schematic figures (2026-10-17 redesign)
    designG_quasiperiodic_v2.py       # Figures 13-14 phase portraits + spectrum; S10 U5 quasiperiodic diagnostics
    certify_period_two_v1.py          # exact rational interval certificate of the attracting two-cycle at
                                      #        (mu0, nu1) = (4, 46/25) [prop:two-cycle] + numerical basin fraction
                                      #        and nu1 scan (standard library for the certificate; numpy for the scan)
    designG_dynamics_v3.py            # SUPERSEDED first-passage producer, kept as the record of what was computed
                                      #        for the July 2026 version of the manuscript; not run by run_designG_figs.sh and not to be
                                      #        used to reproduce this paper (its output: data/designG/superseded/)
    (the retired designG_bifurcation_v1.py and the former Neimark--Sacker producer are
     no longer shipped: the latter carried an externally supplied normal-form coefficient
     that the manuscript does not use or claim)

If Figure 10 leaves no room for its insets (they sit where the orbit rests at the coexistence point left of the onset line, so there is none when DESIGNG_ZOOM_LO, 3.80 by default, is set at or above the onset), or its layout check finds a fault in the drawn figure (an overlap of text, insets or data; an inset, a label of panel (c) or the legend outside its place; anything beyond the canvas; or lettering under 7 pt at the printed width, possible with another Matplotlib version or font), designG_bifurcation_v2.py prints each violation, leaves the existing figures/designG_bifurcation_mu0.png unchanged, writes its other outputs and exits 1, which stops run_designG_figs.sh.

Validator + reproduction:

    validate_designG_v1.py            # claim/JSON/figure/manifest validator (smoke|full modes)
    reproduce_designG_v1.py           # proves the SUBMITTED numbers reproduce (see below)

Claim inventories, report, manifest, data:

    current_claims_designG_v1.json    # validator-native claim inventory (98 claims, 94 of them
                                      #        quoted in the manuscript and quote-checked; refreshed 2026-11-12)
    submitted_claims_designG_v1.json  # FROZEN inventory of the 77 quoted numbers of the July 2026 version of the manuscript
                                      #        (regression reference for reproduce_designG_v1.py; it still lists claims no longer
                                      #        in the manuscript, the row numbers in its notes refer to the validation report
                                      #        of 2026-07, and the files and lines that its entries name are those of the
                                      #        manuscript of 2026-07);
                                      #        the LIVE inventory of the current manuscript is
                                      #        current_claims_designG_v1.json
    report_calculation_v2.md          # human-readable validation report (full-mode PASS)
    data/designG/validation_designG_v1.json  # the same validator run, machine-readable (VALJSON)
    data/designG/*.json               # shipped reference producer JSON (regenerated values)
    data/designG/ics/                 # ARCHIVED random initial conditions of the two Monte-Carlo sensitivity
                                      #        configs (ic.runif, ic.rlogy; N = 40,000, seed 20260709): exact
                                      #        float64 arrays (.npy) + sidecar JSON with dtype, shape, order,
                                      #        generator, seed, sha256 and the NumPy CPU-dispatch state

Regression tests: `test_designG_v4.py` (dynamics/quasiperiodic/validator; targets
the interpreter running pytest, `sys.executable`; set `DESIGNG_PYTHON` to override), `test_designG_bifurcation.py` (including, since 2026-10-21, the two-vector
largest-exponent contract on the invariant axis and the onset helper's exact-zero-Jacobian corner, and, since 2026-11-08, the exponents of panel (c) of Figure 10 and the figure's layout check, with, since 2026-11-09, its refusal when no room is left for the insets), `test_designG_invasion_exponents.py` (exact
rational checks of the controlled-comparison figure) and `test_designG_stable_leaf.py` (the
stable-leaf illustration: convergence of the Lyapunov--Perron iteration with observed ratios of
successive changes below one, tangency to the slope-series direction (Richardson-extrapolated
secant), the 1/107 tracking band, a refused non-converging height, the exact axis orbit against the
float64 orbit, a finite-N comparison of the exact long-run rate with the closed form, and, since 2026-10-21,
the frozen inputs: the base point and the 25 heights agree with their defining formulas and their
SHA-256, and an edited literal is refused).

`test_designG_v4.py` resolves the scripts it exercises through `codes/.designg_bindings.sh`
rather than hard-coding names, so that file must be **sourced with export** first;
without it the tests abort with `DYN not exported` (`validate_designG_v1.py` reads only
its own `DESIGNG_*` variables, listed below). From the bundle root:

    set -a; source codes/.designg_bindings.sh; set +a
    python3 -m pytest codes/test_designG_v4.py codes/test_designG_bifurcation.py codes/test_designG_invasion_exponents.py codes/test_designG_stable_leaf.py -q
                                                      # 75 passed (38 + 25 + 4 + 8), ~5 min, CPU only

The 2026-10-18 additions are targeted regression tests with independent oracles: the
initial-condition archive round-trip and a one-ULP tamper, the exceptional-seed
pre-filter against the exact rational orbit 1/4 -> 3/4 -> 3/4, the exact fixed
point of the logistic map at mu0 = 2 whose tangential exponent is -inf, the
periodic-window transverse exponents at mu0 = 3.83 and 3.84, the validator's
finiteness gate, readiness statuses and strict mode on corrupted copies of the shipped
data and section files -- a NaN, a wrong-sign empirical exponent and an edited duplicate
quotation each fail the run -- the exact-angle enumeration's measured coverage,
independent of the maximising value, and the two-cycle certificate, including a
displaced centre and an off-band parameter that must NOT certify. Later additions: the log-space
continuation of the predator to the common horizon against a synchronous float64 oracle; the
validator against JSON strings, a NaN in a completeness-required field and a crash that would have left a
stale report; a truncated or malformed producer JSON must leave FAILED report/VALJSON files, never the
previous PASS, and so must the configuration aborts (an invalid mode, an unknown scope token, a missing
required variable). Nine more followed: both assurance tools
reject a record whose own success flags report a failure (an uncertified two-cycle or a failed
certificate check, a failed collapse guard, an incomplete enumeration); the validator rejects a vacuous
or malformed claims inventory and never calls a scoped run submission-ready; the QP regime parameters
must be finite numbers equal to the resolved regime, its documented default included; a Wilson interval
excuses a quotation only after it is shown to be the interval of the produced count; a checkpointed
configuration is reused only when its configuration and input identity match the request; an exact
predator zero in the telescoping demo stays zero (its logarithm -inf); the period-one enumeration
passes its own checks; and an orbit that reaches the boundary inside the retained window of the M0
exponent routine gets the exact extended exponents (transverse -inf) instead of a shortened average.
The last ones: deleted success flags or deleted certificate checks and exact-skew parameters fail
(both tools); malformed audit rows, an untracked claim cited as coverage, an audit of meaningless rows
and an audit that leaves printed numbers uncovered fail; a Wilson record without its sample counts, or
with another interval method than its claim expects, fails; a full validation fails a QP record
labelled as a smoke run; and a checkpoint cell with a missing field, an inconsistent sample or an
unusable value is recomputed, not reused.

The bindings ship bundle-relative (`DYN=codes/...`); the author tree's copy is
project-root-relative, so the two files differ by design; `run_designG_figs.sh`
differs between the two trees only in its default data directory.

## I/O contract (env vars)

Outside `--smoke`, `designG_dynamics_v4.py` and `designG_quasiperiodic_v2.py` require
explicit output directories (they never write shipped paths implicitly); the other
producers honour the same variables but fall back to their own default directories when
the variables are unset (see each script's header). `run_designG_figs.sh` sets both
variables explicitly, so a README run never relies on those defaults:

    DESIGNG_FIGDIR    figure output dir
    DESIGNG_DATADIR   *_values.json output dir
    DESIGNG_MANIFEST  (optional) manifest file the producers APPEND artifact rows to
                      (the onset ensemble, which writes one JSON and no figure, appends none)
    DESIGNG_QP_SOURCE (quasiperiodic full run) the U5 source-values JSON

The validator requires six vars, `DESIGNG_VALIDATE_IN`, `DESIGNG_VALIDATE_FIGDIR`,
`DESIGNG_CLAIMS_IN`, `DESIGNG_REPORT_OUT`, `DESIGNG_VALJSON_OUT` and `DESIGNG_MANIFEST`
(it exits with status 2 if one is unset), and reads four optional ones:
`DESIGNG_VALIDATE_MODE` (smoke|full; default full), `DESIGNG_TEXDIR` (default: the bundle
root, where the `*_clean.tex` files live), and `DESIGNG_VALIDATE_SCOPE` and
`DESIGNG_VALIDATE_STRICT`, both described below.
Two gates read the manuscript itself. The manifest gate refuses an empty
manifest and requires, for every figure the manuscript includes
(`\includegraphics{figures/...}` in `*_clean.tex`), a row whose path is that
figure AT `DESIGNG_VALIDATE_FIGDIR` (real-path equality) with a matching sha256,
so a correctly hashed row for another copy of a same-named file does not count. The
gate checks artifact path, byte integrity and that the row's producer column names a
script shipped in `codes/`; it does not authenticate source provenance by digest.
The quote gate requires every tracked `quoted_value` in
`current_claims_designG_v1.json` to occur as a number in the named `*_clean.tex`
file (within the larger of its display tolerance and half a unit in the last digit of the
quoted value; a claim that pins `manuscript_token` counts only that exact printed token) at
least as many times as the inventory pins
for that file (`quoted_in`: one count per section file, main-text and Supplement
tables alike, so a quotation edited at any pinned site is caught). The gate reads the printed
text only, through the mask of the coverage scan described below, so that a number in a comment, in a
label, reference or citation key, or in an argument of a length or graphics command is no occurrence.
The matcher parses
complete signed numbers in the TeX forms it supports (an optional minus across
whitespace, line breaks, thin spaces or `~`; a mantissa, whose thousands may be grouped
by the TeX separator `{,}`, so that `40{,}000` is the one number 40000 and never an
occurrence of 40; an exponent written as `\times 10^{b}`, `\cdot 10^{b}` with an
optional sign on `b`, or e-notation): a sign before a number and a decimal,
thousands-group or exponent continuation after it do not count as a match, and
`test_designG_v4.py` keeps these mutations as negative controls (the grouped form since
2026-10-24). The coverage audit is
gated too: a row that calls a number claim-covered must cite a claim that checks that
section file, every line a row lists is re-read from the current section file and
must still carry that number (a stale line number fails the run instead of passing),
and every numeric token of a row's literal must be found on at least one listed line.
This is a presence/occurrence check of the quoted number in
the section file -- it does not bind the number to a particular sentence or
label -- and it is stated as such. In full mode both gates block.

Three further rules were added on 2026-10-18:

- **Non-finite values are hard failures.** A producer value that is NaN, or an infinity
  the claim does not declare in `allowed_nonfinite`, is reported as `NONFINITE` and fails
  the run in both modes (previously a NaN compared "consistent" with any quoted number).
- **Structural verdict and manuscript readiness are separate.** `overall` is the
  structural verdict (schema, anchors, gates); `readiness` is `SUBMISSION_READY` only
  for a full run of all blocks in which nothing failed and no tracked claim requires a
  manuscript change; otherwise it is `FAIL` (a structural failure), `CHANGES_REQUIRED (n ...)`,
  `SCOPED_PASS (...)` (a run restricted by `DESIGNG_VALIDATE_SCOPE`, below) or, in smoke
  mode, `SMOKE_ONLY (...)`. Both are printed on the final line, written to the VALJSON and to the
  report header. A PASS with `CHANGES_REQUIRED` is a structural pass, not readiness.
- **Strict mode.** `DESIGNG_VALIDATE_STRICT=1` makes any `CHANGE_REQUIRED` row fail the
  run (exit 1), so a report that still requires manuscript changes cannot be presented as
  an unconditional pass. The shipped record was produced in strict mode.

Four more rules were added later on 2026-10-18:

- **Success flags are results.** A certificate that is not certified or has a failed check,
  a failed exact-skew collapse guard, or an incomplete periodic-orbit enumeration (or a
  non-zero least-period count mismatch) fails the run in both modes, whatever the numbers
  say (`essential_flags` gate); the reproducer checks the same flags on both sides and
  compares every boolean leaf exactly.
- **The inventory is checked before it is used.** The claims file must hold a non-empty
  claim list with unique string ids and a non-empty coverage audit, and a full unscoped run
  needs at least one claim for every block the validator specifies (a list fixed in the
  validator, not derived from the inventory under check). A run restricted by
  `DESIGNG_VALIDATE_SCOPE` reports `SCOPED_PASS (...)`, never `SUBMISSION_READY`.
- **A Wilson interval must belong to the quoted number.** It excuses a quotation (status
  `CONSISTENT_WITHIN_CI`: statistically compatible, not equal at display precision) only
  for a claim on the cell's terminal `fraction`, and only after it is shown to be the
  Wilson interval of that fraction's integer count; an inconsistent record is
  `CI_INCONSISTENT`, a hard failure.
- **The QP regime parameters are checked as numbers.** `mu0` and `nu1` of the
  quasiperiodic output must be finite JSON numbers equal to the resolved source regime,
  including its documented default nu1 = 4.90.

The last additions of 2026-10-18 close the remaining shortcuts:

- **The success flags cover fixed members.** The certificate must carry its nine named checks and the
  exact-skew block its three parameters; a deleted flag, check or parameter fails (both tools).
- **The coverage audit is evidence, not a list.** Every row names a known claim or a known exclusion
  category and lines that exist; a row without a checkable number may only document a symbol written
  out; a claim cited as coverage must have its manuscript check switched on; and in full mode every
  number printed in the six section files must equal, at its printed precision, a number of some row
  of its file. A printed number is its complete signed value with its exponent (`\times 10^{b}` or
  `\cdot 10^{b}`, or e-notation; a bare `10^{b}` counts as 1 x 10^{b}) and is in scope when it has an
  exponent, at least two decimals, or is an integer of three or more digits (years excepted). Values are compared as exact
  decimals; a number with more than 60 significant digits or an exponent beyond 9999 in magnitude is
  reported as unsupported, so a tiny power keeps its magnitude and sign. TeX spacing around `\times`,
  `\cdot`, `10` and `^`, between an exponent's sign and its digits, and among the digits of a bare
  power's exponent changes nothing. The source is read as TeX tokens from left to right (`\\` and `\%`
  are units, so text after a line break is text; a star counts only for a command that takes one;
  verbatim text is kept as printed text), and only source that is
  never printed is masked: comments, label, reference and citation keys, file and environment names,
  and the exact arguments of length and graphics commands (`\hspace`, `\setlength`,
  `\includegraphics`, `\phantom`, ...); a brace group that merely follows such a command, caption
  text, citation pinpoints and prose after a reference are checked. Reported as unsupported rather than
  read as fragments: a number whose digits continue on the next printed line after a comment or a
  chain of whole-line comments (TeX joins the digits), and a bare power of ten whose braced exponent,
  on one line and without nested braces, is not a readable integer. A number belongs to the line of
  its sign or first digit. An argument the mask does not recognise (one spanning lines or containing
  braces) stays visible, so it can only add numbers to check. The coverage scan and the audit-line
  check read numbers through one iterator (`_printed_numbers`).
  Since 2026-11-10 an
  audit row keeps a line only where that line prints, within the scope of the scan, a number equal to
  one of the row's numbers at the finer of the two printed precisions; a line held only because
  another number rounds to a coarser number of the row, or only through a number outside the scope
  (the 0 of a subscript), fails the validation as a weak placement.
- **The coverage scan is a guard, not a TeX interpreter.** It reads the TeX described above, and it
  exists to catch a printed number left out of the audit. Source that writes a number through other
  constructs is outside its contract. It misreads the following forms in a way that can let an
  uncovered number pass:
  - a multiplication sign in braces (`5{\times}10^{-3}`, `5{\cdot}10^{-3}`), which is read as the
    factor (here 5, outside the scope) followed by the bare power `10^{-3}`;
  - after `\times` or `\cdot`, an exponent with a space or a decimal point among its digits
    (`\times10^{0 9}`, `\times10^{0.5}`), and any exponent that is nested or spans lines
    (`10^{\frac{1}{2}}`);
  - a number continued after a comment onto a line indented by 63 or more characters, or with a space
    before the `%` in math mode;
  - `\times`, `\cdot` or `^` inside verbatim text, which is read as notation
    (`\verb|0.00002727\times10^{4}|` is read as 0.2727);
  - a decimal whose digits reach about a million places after the point (it underflows to zero);
  - digits split by an empty group or a control sequence (`0.27{}99`, `0.27\relax 99`); a number that
    a macro defined in the document prints but that is not written where the macro is used, other than
    the automatic numbers of sections, theorem environments, equations, figures, tables, list items
    and pages (the scan expands no macro); and catcode changes, which it does not interpret.

  The delivered section files use the first of these forms on four lines of `supplementary_clean.tex`
  and none of the others (checked 2026-11-12): line 2517 prints the grid resolutions `5{\times}10^{-3}`
  and `2{\times}10^{-3}` of the onset ensemble, which the scan reads there as `10^{-3}` (audit row A253
  records this), and lines 2649-2651 (Table S8) print the horizon `7.5{\times}10^4`, which it reads
  there as `10^4`. On those lines these numbers are not read as printed; the two resolutions are read
  as printed on line 2522, and the horizon on lines 2200, 2415, 2453 and 2612. The only command defined
  in the document that they call is `\suppref` (`preamble_local.tex`), which prints its second argument
  as written; the theorem environments declared in `main.tex` and `supplementary.tex` add no number
  other than their automatic numbers.
- **A Wilson record carries its sample.** `n_total` and `n_censored` are required and
  n = n_total - n_censored; a claim that names its expected interval method fails on any other.
- **A smoke record is not a full record.** A full validation fails a QP record whose provenance says
  it came from a smoke run.

The onset-ensemble numbers quoted in Section 6 and the Supplementary Material (medians, ranges and the plotted seed's
onsets at the three resolutions) and the two-cycle certificate numbers are value-checked by
direct claims against `designG_onset_ensemble_v1.json` and
`designG_period_two_certificate_v1.json`, which are validator producers since 2026-10-18.

## Quick smoke run (fast, schema only)

From `codes/` (from the bundle root: `cd codes`). The dynamics smoke run writes its figures and data
and then exits with status 1: with its reduced counts several of its small-sample diagnostics do not pass,
which is expected and is not an acceptance test (the smoke validation of the third command passes):

    S=$(mktemp -d)
    DESIGNG_FIGDIR="$S" DESIGNG_DATADIR="$S" DESIGNG_MANIFEST="$S/manifest.md" \
      python3 designG_dynamics_v4.py --smoke
    DESIGNG_FIGDIR="$S" DESIGNG_DATADIR="$S" DESIGNG_MANIFEST="$S/manifest.md" \
      python3 designG_quasiperiodic_v2.py --smoke
    DESIGNG_VALIDATE_MODE=smoke DESIGNG_VALIDATE_IN="$S" DESIGNG_VALIDATE_FIGDIR="$S" \
      DESIGNG_CLAIMS_IN=current_claims_designG_v1.json \
      DESIGNG_REPORT_OUT="$S/r.md" DESIGNG_VALJSON_OUT="$S/v.json" \
      DESIGNG_MANIFEST="$S/manifest.md" python3 validate_designG_v1.py

## FULL reproduction of the submitted numbers (non-smoke)

Regenerate every producer JSON + figures, validate all blocks, and confirm the
submitted numbers reproduce. Run from `codes/` (from the bundle root: `cd codes`), the
directory the shipped record was produced in:

    S=$(mktemp -d); mkdir -p "$S/data" "$S/figs"
    # every manuscript figure + its values JSON + manifest rows (dynamics, ergodic,
    # bifurcation, onset ensemble, invasion exponents, the stable-leaf illustration, schematics, the
    # quasiperiodic producer (Figures 13-14) with its U5 block, and the two-cycle certificate; several
    # minutes). PYTHON names the interpreter EXPLICITLY (absolute path) so the runner
    # cannot pick up another environment ahead of the one you installed:
    DESIGNG_FIGDIR="$S/figs" DESIGNG_DATADIR="$S/data" DESIGNG_MANIFEST="$S/manifest.md" \
      PYTHON="$(command -v python)" bash run_designG_figs.sh

    # (a) validator, all blocks, against the tracked anchors/invariants (strict: a
    #     tracked claim that still requires a manuscript change fails the run):
    DESIGNG_VALIDATE_MODE=full DESIGNG_VALIDATE_STRICT=1 \
      DESIGNG_VALIDATE_IN="$S/data" DESIGNG_VALIDATE_FIGDIR="$S/figs" \
      DESIGNG_CLAIMS_IN=current_claims_designG_v1.json \
      DESIGNG_REPORT_OUT="$S/report.md" DESIGNG_VALJSON_OUT="$S/val.json" \
      DESIGNG_MANIFEST="$S/manifest.md" python3 validate_designG_v1.py
                                   # -> structural=PASS readiness=SUBMISSION_READY strict=on (exit 0)

    # (b) confirm the SUBMITTED manuscript numbers reproduce (submitted_claims):
    python3 reproduce_designG_v1.py --regen "$S/data"      # -> OVERALL: PASS (exit 0)

Step (a) reads the outputs of the producer step. Run instead on the shipped `data/designG/`
reference (`DESIGNG_VALIDATE_IN=data/designG DESIGNG_VALIDATE_FIGDIR=../figures`, with no producer
step before it), the validator prints `structural=FAIL readiness=FAIL strict=on` (exit 1). Two of
its eight gates fail because no producer run preceded it: `manifest`, since the manifest is written
by a producer run, and `completeness`, since `designG_lambda_convergence.png` and
`designG_exact_skew.png`, two figures the producer draws and the manuscript does not use, are not
shipped. The other six gates report OK, and none of the 98 checks fails (67 PASS and 31
OK_NONBLOCKING).

`reproduce_designG_v1.py` checks four things and reports a fifth: (A) a deterministic
regression -- every numeric leaf of each regenerated producer JSON (the nine result
files of the manuscript figures and data: dynamics, quasiperiodic, ergodic quadrature,
onset ensemble, invasion exponents, the mu0-sweep onset file and the nu1-sweep window file, the stable-leaf
illustration, and the two-cycle certificate) matches the shipped `data/designG/` reference within tolerance
(`rtol=1e-6`, `atol=1e-9`; the pass is not a claim of bitwise identity), with one stated
exception: `designG_stable_leaf_v1_values.json:/nu1_*/contraction_max`, the largest observed ratio
of successive changes of the Lyapunov--Perron iteration (an illustration diagnostic, quoted nowhere
in the paper), is a quotient of changes near 1e-13 that rounding in the inputs moves by up to
3.2e-5, so it is compared at an absolute 1e-3 with `0 <= value < 1` required on both sides, and the
gate prints a `[FIELD-TOL]` line for it (it was measured 1.4e-5 off on an AVX-512 CPU;
the height grid and base point of that producer are now frozen binary64 literals; the pattern's `*`
stands for one path component), and a mismatch, a field-rule violation or a missing file in any one
producer fails (A) and the run; (B) a direct
resolution -- each machine-parseable `json_source_field` in the frozen
`submitted_claims_designG_v1.json` (65 of its 77 entries) is looked up in the
regenerated JSON and compared with the recorded value; (B) fails unless the file holds these 77
entries and 65 of them resolve (`--expect-claims ENTRIES RESOLVED` states the counts of another
inventory). The remaining 12 entries
are derived or hand-annotated quantities that (B) does not value-check; they are
covered only by the regression in (A). (A') compares every boolean leaf of the same files
exactly (pass flags, stability and chaos labels, certificate checks), and (D) requires, on
both sides, that the result records report success: the two-cycle certified with every
check ok, every exact-skew collapse guard ok, the enumeration complete with zero
least-period count mismatches (numeric agreement alone had passed records stating that a
certificate or a guard failed). (C) is informational: for the two random
sensitivity configurations it reports whether the regenerated run used the SAME
initial conditions as the shipped run (`STRICT`: the archive was loaded, so (A) is a
strict dataset regression) or different ones (`PLATFORM`: a fresh draw or a resample
whose inputs differ at the floating-point level), so the reader knows which comparison
(A) performed. Neither (A) nor (B) reads the current manuscript:
agreement between the manuscript text and the producer JSON is the job of
`validate_designG_v1.py` (claim anchors + the quote gate above).

### Random initial conditions and CPU dispatch

The two Monte-Carlo sensitivity configurations (`ic.runif`, `ic.rlogy`) draw 40,000
seeds from `numpy.random.default_rng(20260709)`; the log-uniform predator seeds are
`10**u` via `np.power`, which NumPy dispatches to a CPU-feature-specific kernel. An
external measurement found that on an AVX-512 host 11,147 of the 40,000 predator
seeds differ by one unit in the last place from the values obtained with the scalar
libm path (which the authors' AVX2 host uses, and which `NPY_DISABLE_CPU_FEATURES=
AVX512F,AVX2,FMA3` reproduced on that AVX-512 host); three of the 19,022 valid
trajectories then change terminal classification and the extinct count at nu1 = 3
reads 5189 instead of 5188. Version pins alone therefore do not fix these inputs.
The producer now LOADS the archived arrays (`data/designG/ics/`, verified against
their sidecar sha256) by default and records, under
`sensitivity.<config>._params.ic_provenance`, what it used, the fresh draw's ULP
distance from the archive and the NumPy CPU-dispatch state of the run.
`DESIGNG_IC_MODE=resample` uses the fresh seeded draw instead (a statistical
comparison, platform-dependent at the ULP level); `DESIGNG_IC_ARCHIVE_DIR` points the
loader at another archive. A resumed run (a checkpoint
`designG_dynamics_v4_sensitivity_checkpoint.json` in `DESIGNG_DATADIR`) reuses a
checkpointed configuration only when its configuration (grid, horizon, threshold, sample
size, nu1 sweep, seed, protocol, terminal window, fallback flag) and, for a random
configuration, its input identity (`DESIGNG_IC_MODE` and the SHA-256 of the array) match
the request, and every requested cell is present with usable, consistent values: the sample its inputs
define (the total and its left-censored part) with a non-empty valid part; finite fractions in [0, 1];
terminal and first-passage fractions that are counts over the valid sample, the terminal count not above
the first-passage count; the censored-inclusive fraction equal to (first-passage count + censored) /
total; the Wilson interval of the terminal count for a random configuration and none for a grid;
otherwise that
configuration is recomputed, so `DESIGNG_IC_MODE=resample`
always resamples (the FULL recipe starts from a fresh directory and never resumes).
Three things are kept apart: the fixed-input regression (archive loaded), the resampled
statistical comparison (Wilson intervals), and the mathematics, which none of this
touches.

The shipped validation record -- `report_calculation_v2.md` and
`data/designG/validation_designG_v1.json` -- is the validator's `$S/report.md` and `$S/val.json`
from step (a) above. Two separate runs stand behind it. The producer step of the FULL sequence
was run from `codes/` on 2026-11-10 (UTC). The validator was run afterwards on those producer
outputs, with the claims inventory re-pinned to the text of this version (its main text and
Supplementary Material) and in its latest amended state, so that the record's manuscript line
references match the shipped text; the record states the time of that validation run,
2026-11-14T04:59:50.983119+00:00 (its `_meta.generated`, repeated on line 3 of
`report_calculation_v2.md`). Every value, status and gate in it equals that of the record of
2026-10-24, run on the outputs behind the shipped `data/designG/` reference (whose dynamics and
ergodic values files were regenerated later with every result unchanged; see the next two
paragraphs); only the manuscript
quotation sites (file, line, printed form and number of occurrences), the audit-row counts and the
tolerances of six count ratios (the five extinction fractions of Table S6 and the base fraction at
nu1 = 2 of Table S7, which this version checks at the resolution of their grids) differ. Its
header paths are printed verbatim. The runner writes its log next to the
manifest (`$S/run_designG_figs.log`): the documented FULL sequence writes its scientific
outputs, manifest, validation reports and runner log into the scratch directory, while
interpreter and Numba compilation caches (`__pycache__/`) may still be created beside the
scripts. The
validator re-resolves the U5 source JSON the way the producer does (`DESIGNG_QP_SOURCE`, else
`LATEST_QP_VALUES`, else `data/designG/designG_quasiperiodic_values.json` of either tree
layout), so that variable need not be repeated for step (a). That source file was written on
2026-06-05 by `designG_quasiperiodic_v1.py`, a retired producer that is not shipped. The current
producer, `designG_quasiperiodic_v2.py`, reads only its `mu0` (3.5) and uses the nominal
nu1 = 4.90, recording the status `default_nu1` under `_provenance.qp_source_status`; the file's
`rho_at_4p9` and `theta_over_2pi` equal the values that `designG_quasiperiodic_v2.py` computes and
writes to `designG_quasiperiodic_v2_values.json` (0.2241 and 0.2301).

Of the shipped `data/designG/` reference, `designG_dynamics_v4_values.json` was regenerated on
2026-11-12 (UTC), after the closing sentence of the producer's `provenance.honesty` note was
reworded: "Closed forms are mu0=4-specific." now reads "All checks run at mu0=4; ln(nu1/4) rests on
the arcsine law, which holds at mu0=4 only." The file it replaced was produced on 2026-10-21. The
producer ran without a checkpoint, so every result was recomputed, the eleven sensitivity
configurations included; the checkpoint it wrote, which a bare run of `run_designG_figs.sh` reuses
(see "Figures only" below), replaced the shipped one of 2026-10-21 and differs from it only in the
run time recorded for each configuration (`wall_s`). Apart from the provenance block (its time, that
note and the two output directories), the list of figure paths and those run times, the new values
file equals the replaced one value for value. Of the six figures that this run drew, the four that
the manuscript uses are byte-identical to those in `../figures`. The run also wrote the archive of
initial conditions, `data/designG/ics/`, again: its two `.npy` arrays are byte-identical to those
they replaced, and its two `.json` descriptions differ from theirs only in `written_utc`, which
read 2026-10-17T16:44:08.873241 (`ic.runif`) and 2026-10-17T16:44:08.875720 (`ic.rlogy`) and now
read 2026-11-12T14:23:12.027586 and 2026-11-12T14:23:29.172145. The superseded
`designG_dynamics_v3.py` ended its note with the same sentence and now reads the same way; its
archived output in `data/designG/superseded/` was not regenerated and keeps the earlier sentence. On 2026-11-14 the provenance entry `round` was removed from the four values files that carried it (this file, the ergodic and the quasiperiodic values file and that archived output), one line each; no other byte changed, and no file was regenerated for it. On the same day two note strings of this file were changed by hand to read as `designG_dynamics_v4.py` now writes them (lines 191 and 2318: `legacy_first_passage_note`, which now gives the terminal count that Remark S12.1 of the Supplement prints at nu1 = 2, and the note of the check "overall max y (nu1<=4) <= 1"), and that archived output was changed by hand to what `designG_dynamics_v3.py` now writes (the note of the same check; the two entries `spec_figure_nu1_2` and `spec_figure_note` deleted, two lines: they held a figure that the manuscript does not print, not a result of the run); no computed value changed, and no file was regenerated for it.

`designG_ergodic_quadrature_v2_values.json` was regenerated on 2026-11-12 (UTC) as well, after the
producer's `provenance.honesty` note was made to say what the note it prints says. It read
"Closed forms (arcsine density, I(4)=-2ln2, nu1*(4)=4) are mu0=4-specific; nu1*(mu0) for mu0<4 is
a per-parameter orbit average"; it now says that these three hold at mu0=4 only, that nu1*(mu0) has
a closed form on an attracting cycle of the prey (mu0/(mu0-1) for 1<mu0<=3, mu0/sqrt(mu0+1) for
3<mu0<=1+sqrt6), and that the six values at mu0=3.60 to 3.99 are orbit averages. The file it
replaced was produced on 2026-10-21. Apart from the provenance block (its time, that note and the
two output directories) and the list of figure paths, the new file equals the replaced one value
for value, and the three figures this run drew are byte-identical to those in `../figures`.

### Stable-leaf illustration (`designG_stable_leaf_v1.py`, added 2026-10-21)

Writes `figures/designG_stable_leaf.png` (manuscript Figure 5: panels (a) and (b) below),
`figures/designG_stable_leaf_rates.png` (manuscript Figure 12: panel (c) below) and
`data/designG/designG_stable_leaf_v1_values.json` (computed with NumPy and the standard library
only, the figures drawn with Matplotlib; deterministic, no random numbers, a few seconds). Both figures are drawn on a 6.5-inch canvas;
Figure 12 is included at the 6.27-inch text width of the A4 page and Figure 5 at 0.86 of it, so
their 7.5--8.5 pt labels print at about 7.2--8.2 pt and 6.2--7.1 pt.
The base point and the 25 heights are frozen binary64 literals (an import-time SHA-256 check
refuses an edited literal): np.logspace differed in the last bit on an AVX-512 CPU. At mu0 = 4,
for nu1 = 2.1 and 3, and at the base point p = sin^2(pi theta0/2) with theta0 = sqrt(2) - 1
(p ~ 0.3669; its angle is known to every bit, so its exact orbit can be computed), it shows
(a) the local stable leaf through (p, 0) as the fixed point of the manuscript's Lyapunov--Perron
operator (eq:LP-stable, eq:LP-unstable) on a 60-step horizon, plotted as the deviation
|psi(t)| = |X(t) - p - s0 t| from the slope-series tangent (eq:slope-series), consistent at
small heights with |psi_2| t^2 (the C^{1,1} structure of eq:leaf); (b) along the orbit of the leaf point of height
t = 1e-6, the predator against the matched reference product t V_n(p) of cor:rare-predator and
their ratio, which the corollary bounds by 1/107 on the theorem's set (40 steps; float64
rounding, roughly doubled at every step, ends the comparison beyond that); and (c) the running
rate (1/N) ln V_N(p) of the reference product along the EXACT orbit of p up to N = 1e5, computed
through the conjugacy x = sin^2(pi theta/2) with the tent map on the integer
floor(theta * 2^K) (exact; float64 only for each x_N), against the closed form ln(nu1/4).
The JSON records, per nu1: `s0`, the converged heights, leaf points and `psi`, the quadratic
coefficient `psi_over_t2_small_t`, the largest converged height and the first refused one; under
`matched`: `max_abs_ratio_minus_1`, the 40-step rates and the itinerary agreement; and the exact
running rate (`reference_rate_exact`, `reference_rate_at_40`, `reference_rate_N_long`).
It is an ILLUSTRATION of the objects the proof constructs, not an instance of the theorem: the
base point and the heights are not certified to lie in the theorem's block B x [0, c0(nu1)], no
constant of the proof is computed, and nothing here is used in any proof. A height at which the
iteration does not converge is refused (RuntimeError) and never plotted; the recorded ratios of
successive changes are observations, not a certificate that the operator contracts. The
finite-time decay rates are finite Birkhoff averages at this base point (after 40 steps -0.839
and -0.482; -0.638 and -0.281 at N = 1e5), compared with the almost-everywhere limit ln(nu1/4)
(-0.644 and -0.288); the genericity of the base point, and the convergence of its averages to
that limit, are not established.

## Figures only

    bash run_designG_figs.sh    # every manuscript figure -> ../figures (the manuscript's
                                # figure directory), values JSON -> data/designG, and a
                                # fresh producer manifest run_designG_figs_manifest.md

Run bare like this, the sensitivity sweep is not recomputed: `data/designG/` holds the shipped
checkpoint, whose eleven configurations match the request and are reused, and the shipped reference
JSON and figures are overwritten in place. Use the FULL recipe above (a fresh directory) to recompute
them. A producer that fails stops the script; the dynamics and the ergodic producer, when they cannot
draw their figures, still write their values JSON, print `[ERROR] figure generation FAILED` and exit
with status 1.

## Companion repository

The same bundle is published in the companion repository
`github.com/weishanlee/LotkaVolterraChaos-DesignG-Extinction` (its `codes/` is this
directory, its `figures/` the producers' outputs); `code_availability.md` lists what
ships with the submission.
