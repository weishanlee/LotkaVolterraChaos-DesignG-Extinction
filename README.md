# LotkaVolterraChaos-DesignG-Extinction -- companion code and data

Companion repository for *"Predator Extinction from Prey Chaos via Transverse Stability: A
Constructive Positive-Area Extinction Basin with Exact Rate in an Engineered Discrete
Predator--Prey Map"* (Wei Shan Lee, Kei Chon Sio and On Kio Sio). It holds the complete code and data behind
the paper's numerical results:
the CPU producers that regenerate every manuscript figure, the regenerated producer data,
the claim inventories, the validator with its validation record, the reproduction check
and the regression tests. The same bundle accompanies the submission. These are numerical
consistency checks, not proofs; the paper's results rest on its theorems.

## Version

This is the bundle of manuscript revision v31 (November 2026): a 36-page main text with 14
figures and a 56-page Supplementary Material with 3 figures. `codes/README.md` documents the
bundle.

## Layout

    codes/                  producers, validator, reproduce_designG_v1.py, regression tests,
                            claim inventories, validation report (see codes/README.md)
    codes/data/designG/     regenerated producer JSON (the reference the reproduction check
                            compares against) and the validator's machine-readable record;
                            ics/ holds the archived random initial conditions (with SHA-256)
                            that the two Monte-Carlo sensitivity runs load
    figures/                the figures the producers write (the manuscript's figures and the
                            retained diagnostics)
    requirements.txt        pinned CPU dependencies: Python 3.10.x, numpy 1.26.4, numba 0.61.2,
                            matplotlib 3.9.2 (see codes/README.md for what the preflights check)
    payload_manifest.json   every other file of this repository and every file of the
                            submission bundle, each with its SHA-256 and size
    LICENSE                 MIT

`codes/README.md` is the full documentation: the I/O contract (`DESIGNG_*` variables), the
validator gates, the quick smoke run and the FULL reproduction sequence.

## Quick start

    python3.10 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
    set -a; source codes/.designg_bindings.sh; set +a
    python -m pytest codes/test_designG_v4.py codes/test_designG_bifurcation.py \
                     codes/test_designG_invasion_exponents.py \
                     codes/test_designG_stable_leaf.py -q                 # 75 passed (see below)

    cd codes && S=$(mktemp -d) && mkdir -p "$S/data" "$S/figs"
    DESIGNG_FIGDIR="$S/figs" DESIGNG_DATADIR="$S/data" DESIGNG_MANIFEST="$S/manifest.md" \
      PYTHON="$(command -v python)" bash run_designG_figs.sh             # every manuscript figure
    python reproduce_designG_v1.py --regen "$S/data"                      # -> OVERALL: PASS

The validator's manuscript-quotation gate reads the manuscript's `*_clean.tex` section files,
which belong to the submission rather than to this repository. Eleven of the 75 tests run the
validator or its quotation check and look for these files at the repository root: with the
submission's seven `*_clean.tex` files copied there, all 75 pass; without them, those eleven fail
and the other 64 pass. The figure run and `reproduce_designG_v1.py` do not read them. To run
`validate_designG_v1.py` in full mode, point `DESIGNG_TEXDIR` at a directory holding them
(the sequence is given in `codes/README.md`).
