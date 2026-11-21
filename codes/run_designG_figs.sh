#!/usr/bin/env bash
# Regenerate EVERY figure the manuscript includes (CPU; about twelve minutes on the authors'
# host) and the two-cycle certificate. Writes figures + values JSON + a producer
# manifest; see README.md for the validator step that checks the manifest
# against the manuscript's \includegraphics list.
#
#   bash run_designG_figs.sh
#
# Override the interpreter with PYTHON=/path/to/python (defaults to the pinned
# ~/.venv/nolds if present, else python3).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -n "${PYTHON:-}" ]; then PY="$PYTHON";
elif [ -x "$HOME/.venv/nolds/bin/python" ]; then PY="$HOME/.venv/nolds/bin/python";
else PY="python3"; fi

# output/ layout: the manuscript's \includegraphics{figures/...} resolves against
# output/, so figures belong in output/figures -- one level ABOVE this script.
# (This is also the Python scripts' own default, _OUT/figures.)  The
# bundle has the same codes/ -> ../figures relation, so both layouts write the
# manuscript figure directory; only the data directory differs (the bundle keeps
# its reference JSON under codes/data/designG).  Pre-set DESIGNG_* variables are
# honoured, so the FULL reproduction in README.md can route everything to one
# scratch directory and one manifest.
export DESIGNG_FIGDIR="${DESIGNG_FIGDIR:-$(cd "$HERE/.." && pwd)/figures}"
export DESIGNG_DATADIR="${DESIGNG_DATADIR:-$HERE/data/designG}"
export DESIGNG_QP_SOURCE="${DESIGNG_QP_SOURCE:-$DESIGNG_DATADIR/designG_quasiperiodic_values.json}"
if [ -z "${DESIGNG_MANIFEST:-}" ]; then
  export DESIGNG_MANIFEST="$HERE/run_designG_figs_manifest.md"
  : > "$DESIGNG_MANIFEST"   # a fresh manifest per run; the validator refuses an empty one
fi
mkdir -p "$DESIGNG_FIGDIR" "$DESIGNG_DATADIR"

# The run log sits next to the manifest: the scratch directory of a README FULL run,
# or this directory for a bare run -- the outputs, manifest and log of a README run go to
# the scratch directory (Numba compilation caches may still appear in __pycache__/).
LOG="${DESIGNG_RUNLOG:-$(dirname "$DESIGNG_MANIFEST")/run_designG_figs.log}"
{
  echo "=== Design-G manuscript-figure regeneration $(date -u) ==="
  echo "interpreter: $PY"
  echo "figures -> $DESIGNG_FIGDIR ; data -> $DESIGNG_DATADIR ; manifest -> $DESIGNG_MANIFEST"
  echo "--- designG_schematics_v2.py (architecture / mechanism / workflow schematics) ---"
  "$PY" "$HERE/designG_schematics_v2.py"
  echo "--- designG_dynamics_v4.py (lambda_perp vs nu1, terminal extinction fractions, Jury band, sensitivity) ---"
  "$PY" "$HERE/designG_dynamics_v4.py"
  echo "--- designG_ergodic_quadrature_v2.py (equilibrium potential, sub-action, I(mu0) family) ---"
  "$PY" "$HERE/designG_ergodic_quadrature_v2.py"
  echo "--- designG_bifurcation_v2.py (mu0 sweep at nu1=2.1 + companion M0 Lyapunov figure) ---"
  "$PY" "$HERE/designG_bifurcation_v2.py"
  echo "--- designG_onset_ensemble_v1.py (24-seed x 3-resolution onset ensemble; the onset numbers of Section 6 and Supplementary Section S12) ---"
  "$PY" "$HERE/designG_onset_ensemble_v1.py"
  echo "--- designG_invasion_exponents_v1.py (closed-form invasion exponents vs nu1: arcsine, period-3,"
  echo "    mu0=2 control, mean-variance bound; prop:controls / lem:meanvar) ---"
  "$PY" "$HERE/designG_invasion_exponents_v1.py"
  echo "--- designG_stable_leaf_v1.py (stable-leaf illustration, Figures 5 and 12: LP fixed point + cor:rare-predator ratio; the two decay rates) ---"
  "$PY" "$HERE/designG_stable_leaf_v1.py"
  echo "--- designG_quasiperiodic_v2.py (Figures 13-14 phase portraits + spectrum; U5 quasiperiodic diagnostics JSON;"
  echo "    the U5 source regime is resolved by the producer from DESIGNG_QP_SOURCE or the shipped data/designG copy) ---"
  "$PY" "$HERE/designG_quasiperiodic_v2.py"
  echo "--- certify_period_two_v1.py (exact rational interval certificate of the attracting two-cycle at nu1=46/25,"
  echo "    prop:two-cycle, plus the numerical basin fraction and nu1 scan; writes a result JSON, no figure) ---"
  "$PY" "$HERE/certify_period_two_v1.py" --out "$DESIGNG_DATADIR"
  echo "=== done -> $DESIGNG_FIGDIR (manifest rows in $DESIGNG_MANIFEST) ==="
  echo "Every current manuscript producer has been run (the retained historical designG_dynamics_v3.py is not one)."
} 2>&1 | tee "$LOG"
