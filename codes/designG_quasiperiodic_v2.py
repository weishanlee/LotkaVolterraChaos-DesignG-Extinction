#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design-G quasiperiodic illustration (Figures 13 and 14 of Section 6; Supplementary Section S10) and U5 diagnostics.

Design-G map: x' = mu0*x*(1-x)*(1-y), y' = nu1*x*y*(1-y). At the coexistence
equilibrium (x*, y*) = (mu0/(mu0+nu1), 1 - 1/mu0 - 1/nu1) one has
det J = 1 - mu0*x*(1-2y*), so a non-real unit multiplier pair occurs exactly on
the curve y* = 1/2, i.e. nu1 = 2*mu0/(mu0-2) (manuscript Proposition prop:ns-locus).
At mu0 = 3.5 that unit-circle locus is nu1_NS = 14/3 ~ 4.6667, with linear angle
cos(theta) = 1 - mu0/4 = 1/8.

What this script establishes NUMERICALLY (and only that): below the locus
(nu1 = 4.6) the post-transient orbit sits at the equilibrium -- the panel shows the
transient spiral of the first 400 iterates from (x*+0.15, y*+0.10); above it
(nu1 = 4.8, 5.1) the post-transient orbit fills a closed invariant curve; at nu1 = 5.5
that curve has broken up. Whether the bifurcation is supercritical, the sign of the
first Lyapunov (normal-form) coefficient, any Chenciner-type degeneracy, and the
growth law of the circle are NOT computed or claimed here; the manuscript states the
same limitation. The circle diagnostics (U5 block) are finite-time quantities on a
computed orbit, not proofs of irrational rotation or of an invariant-circle theorem.

Outputs (DESIGNG_FIGDIR / DESIGNG_DATADIR):
  designG_quasiperiodic_panel.png  -- four phase portraits at mu0 = 3.5 (above)
  designG_quasiperiodic_fft.png    -- Fourier amplitude |X(f)| of the Hann-windowed
                                      prey deviation x_n - x* at nu1 = 4.9, N = 8192
  designG_quasiperiodic_v2_values.json -- nu1_NS, theta/2pi, FFT-peak frequency, and
                                      the U5 block (phase-unwrapped rotation number,
                                      radial-spread residual, radius-series FFT
                                      peak-to-median ratio, Benettin spectrum)

Protocol: orbits start at (x*+1e-3, y*+1e-3); portraits discard 60000 iterates and
keep 4000; the spectrum keeps 8192 after the same burn-in; the U5 block uses a
dedicated orbit with 5000 burn-in and 20000 samples (see _provenance). Every saved
artifact appends a row to DESIGNG_MANIFEST (``_manifest_append``). ``--smoke`` runs
the full schema at reduced sizes (``_smoke_params``); the preflight enforces the
pinned Python 3.10 / numpy 1.26.4 / numba 0.61.2 environment.
"""
# === Revision v2 (2026-07-09, Task B0 scaffold) =================================
# Additive scaffold ONLY -- the v1 orbit()/designG_step()/coexistence_fixed_point()
# math and the mu0=3.5, nu1 in {4.60,4.80,5.10,5.50} (panel) / nu1=4.90 (FFT)
# regime are UNCHANGED; a DEFAULT (non-smoke, non-only) run reproduces the same
# figures/JSON values as v1 bit-for-bit (same n_burn/n_keep passed explicitly
# where v1 relied on orbit()'s own defaults). Added: _preflight() (interpreter/
# package guard), resolve_paths(smoke) (env-var I/O contract), _smoke_params()
# (KEEPS mu0=3.5/nu1=4.9; reduces ONLY n_burn 60000->500 and the retained
# sample count -- panel's 4000 and the FFT's previously-hardcoded 8192 both
# ->2000 -- never switches to the mu0=4/nu1=2 extinction axis), --smoke/--only/
# --public CLI dispatch, _manifest_append() (DESIGNG_MANIFEST producer rows),
# and top-level "_taxonomy"/"_provenance" JSON keys (B0 JSON-schema convention;
# "_taxonomy" is left EMPTY here -- rho/largest_lyap/... + their taxonomy
# entries are added by Task B5). The previously-hardcoded FFT window (8192) is
# now the variable ``nfft`` so --smoke (n_keep=2000) does not shape-mismatch
# against a hardcoded np.hanning(8192); the DEFAULT path still uses nfft=8192,
# so non-smoke behaviour is bit-for-bit identical to v1.
#
from __future__ import annotations
import argparse
import os, math, json, sys
from datetime import datetime, timezone
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__)); _OUT = os.path.dirname(_HERE)

def _rel(p):
    """Record paths RELATIVE to the bundle/authoring root (the parent of codes/), never an
    absolute author-local path; non-string or empty values pass through."""
    if not isinstance(p, str) or not p:
        return p
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        return os.path.relpath(os.path.abspath(p), root)
    except ValueError:
        return os.path.basename(p)

FIGDIR = os.environ.get("DESIGNG_FIGDIR", os.path.join(_OUT, "figures"))
DATADIR = os.environ.get("DESIGNG_DATADIR", os.path.join(_OUT, "data", "designG"))
MU0 = 3.5
NU1_NS = 2.0*MU0/(MU0-2.0)   # = 14/3 ~ 4.6667

def designG_step(x, y, mu0, nu1):
    return mu0*x*(1.0-x)*(1.0-y), nu1*x*y*(1.0-y)
def coexistence_fixed_point(mu0, nu1):
    return mu0/(mu0+nu1), 1.0 - 1.0/mu0 - 1.0/nu1

def orbit(nu1, n_burn=60000, n_keep=4000, mu0=MU0):
    x0, y0 = coexistence_fixed_point(mu0, nu1)
    x, y = x0+1e-3, y0+1e-3; xs, ys = [], []
    for n in range(n_burn+n_keep):
        x, y = designG_step(x, y, mu0, nu1)
        if not (0 <= x <= 1 and 0 <= y <= 1): return None, None, (x0, y0)
        if n >= n_burn: xs.append(x); ys.append(y)
    return np.array(xs), np.array(ys), (x0, y0)

def orbit_radius(nu1, mu0=MU0):
    xs, ys, (x0, y0) = orbit(nu1, mu0=mu0)
    if xs is None: return float("nan")
    return float(np.max(np.hypot(xs-x0, ys-y0)))

def angular_gap(nu1, mu0=MU0):
    xs, ys, (x0, y0) = orbit(nu1, mu0=mu0)
    if xs is None: return float("nan")
    ang = np.sort(np.arctan2(ys-y0, xs-x0))
    return float(np.diff(np.concatenate([ang, [ang[0]+2*np.pi]])).max())

def dominant_frequency(nu1, mu0=MU0):
    xs, ys, (x0, y0) = orbit(nu1, n_keep=8192, mu0=mu0)
    if xs is None: return float("nan")
    sig = (xs - x0)[:8192]*np.hanning(8192)
    sp = np.abs(np.fft.rfft(sig)); fr = np.fft.rfftfreq(8192)
    return float(fr[1 + np.argmax(sp[1:])])


# =============================================================================
# Task B5 (U5): quasiperiodic invariant-circle diagnostics (2026-07-09).
# Four
# empirical diagnostics + one regression-anchor rotation number on the
# nu1=4.9 circle (above the unit-circle locus 14/3), computed from a DEDICATED orbit (burn-in
# 5e3 / sample N=2e4 in the full run, 500/2000 under --smoke via
# _smoke_params() -- distinct from the panel's 60000/4000 and the FFT block's
# 60000/8192 above, both of which stay v1-bit-for-bit UNCHANGED):
#   rho                        -- phase-unwrap rotation number (regression_anchor)
#   invariant_circle_residual  -- relative radial spread across 360 phase bins
#   fft_peak_snr               -- FFT peak-to-median power ratio of the radius signal
#   largest_lyap / normal_lyap -- Benettin QR Lyapunov spectrum of the FULL 2x2
#     map. The LARGER exponent is the tangent/neutral direction (~0 on a
#     genuine quasiperiodic circle -- a finite-time numerical diagnostic, NOT
#     an exact invariant; tagged empirical_discovery and NON-blocking, no
#     PASS/FAIL tolerance gate); the SMALLER is the transverse/normal exponent
#     (<0 => the circle attracts transversally). The larger one is the
#     neutral/tangent direction, not the transverse one.
# =============================================================================


def _jacobian(x, y, mu0, nu1):
    """2x2 Jacobian of the full map x'=mu0*x*(1-x)*(1-y), y'=nu1*x*y*(1-y),
    written out explicitly."""
    return np.array([
        [mu0*(1.0-2.0*x)*(1.0-y), -mu0*x*(1.0-x)],
        [nu1*y*(1.0-y),            nu1*x*(1.0-2.0*y)],
    ])


def rotation_number(xs, ys, x0, y0):
    """rho = mean(diff(unwrap(atan2(y-y0, x-x0)))) / (2*pi): mean phase
    advance per iterate around the fixed point, in revolutions/step."""
    phi = np.unwrap(np.arctan2(ys - y0, xs - x0))
    return float(np.mean(np.diff(phi)) / (2.0*np.pi))


def invariant_circle_residual(xs, ys, x0, y0, n_bins=360):
    """Relative radial spread across n_bins phase bins: mean_over_bins(std(r
    in bin)) / mean(r). Small => the orbit hugs a clean circle. Bins with
    fewer than 2 samples carry no std estimate and are excluded from the
    mean (undefined std, not a zero-spread claim)."""
    r = np.hypot(xs - x0, ys - y0)
    theta = np.arctan2(ys - y0, xs - x0)
    bin_idx = np.clip(((theta + np.pi)/(2.0*np.pi)*n_bins).astype(int), 0, n_bins - 1)
    stds = [np.std(r[bin_idx == b]) for b in range(n_bins) if np.count_nonzero(bin_idx == b) >= 2]
    if not stds:
        return float("nan")
    return float(np.mean(stds) / np.mean(r))


def fft_peak_snr(r):
    """P = |rfft(r - mean(r))|^2; fft_peak_snr = max(P[1:]) / median(P[1:])
    (the DC bin P[0] is excluded from both, per the Task B5 brief)."""
    P = np.abs(np.fft.rfft(r - np.mean(r)))**2
    return float(np.max(P[1:]) / np.median(P[1:]))


def benettin_spectrum(mu0, nu1, n_burn=5000, n_iter=20000):
    """Benettin QR Lyapunov spectrum of the full 2x2 map, starting from the
    same small-offset IC as orbit(). Carries an orthonormal frame Q (2x2);
    at each step forms M = J(x_n,y_n) @ Q, QR-decomposes M = Q_new R
    (numpy.linalg.qr), and accumulates log|R[i,i]| (abs() guards QR's
    sign-ambiguous diagonal). Returns (largest_lyap, normal_lyap) = (max,
    min) of the two exponents -- see the section banner above for which is
    tangent/neutral vs transverse/normal."""
    x0, y0 = coexistence_fixed_point(mu0, nu1)
    x, y = x0 + 1e-3, y0 + 1e-3
    for _ in range(n_burn):
        x, y = designG_step(x, y, mu0, nu1)
    Q = np.eye(2)
    log_sum = np.zeros(2)
    for _ in range(n_iter):
        J = _jacobian(x, y, mu0, nu1)
        x, y = designG_step(x, y, mu0, nu1)
        Q, R = np.linalg.qr(J @ Q)
        log_sum += np.log(np.abs(np.diag(R)))
    lyaps = log_sum / n_iter
    return float(np.max(lyaps)), float(np.min(lyaps))


# =============================================================================
# Scaffold: preflight, path resolution, --smoke params, manifest (Task B0,
# 2026-07-09 -- shared contract, duplicated verbatim-in-spirit in the sibling
# script designG_dynamics_v3.py; this project duplicates rather than shares
# utilities across CPU scripts by design)
# =============================================================================


def _preflight() -> None:
    """Abort early on an interpreter/package mismatch (never a silent wrong-env run).

    Enforces the lockfile -- Python 3.10.x, numpy 1.26.4, numba 0.61.2 -- for EVERY
    run. (This script does not itself use numba, but the lockfile is shared across
    the Design-G suite.)

    DEFAULT INVERTED 2026-07-17, identically to
    ``designG_dynamics_v4.py::_preflight``; see that docstring for the reasoning.
    In short: the old default additionally demanded ``sys.prefix == ~/.venv/nolds``,
    an author-local path, so the bundle's own README commands aborted for every
    third party. The version pins already catch the wrong-environment case strictly
    earlier, so the path check was redundant for the author and a hard blocker for
    everyone else. ``LV_LOCAL=1`` restores it; ``--public``/``LV_PUBLIC=1`` are now
    accepted NO-OPS so existing documented commands keep working.
    """
    import numpy, numba

    if sys.version_info[:2] != (3, 10):
        sys.exit(f"[preflight] Python 3.10.x required (numba 0.61.2 is version-sensitive), got {sys.version.split()[0]}")
    want = {"numpy": "1.26.4", "numba": "0.61.2"}
    got = {"numpy": numpy.__version__, "numba": numba.__version__}
    bad = {k: (got[k], want[k]) for k in want if got[k] != want[k]}
    if bad:
        sys.exit(f"[preflight] package mismatch vs requirements.txt: {bad}")
    if os.environ.get("LV_LOCAL") == "1":
        if sys.prefix.rstrip("/") != os.path.expanduser("~/.venv/nolds"):
            sys.exit(f"[preflight] LV_LOCAL=1 asserts the author venv ~/.venv/nolds, got {sys.prefix}")


def resolve_paths(smoke: bool) -> "dict[str, str]":
    """Env-var I/O contract: DESIGNG_FIGDIR/DESIGNG_DATADIR are REQUIRED outside --smoke."""
    figd, datad = os.environ.get("DESIGNG_FIGDIR"), os.environ.get("DESIGNG_DATADIR")
    if not smoke and (not figd or not datad):
        sys.exit("[paths] DESIGNG_FIGDIR/DESIGNG_DATADIR required outside --smoke")
    return dict(figdir=figd or "./smoke", datadir=datad or "./smoke")


def _smoke_params() -> "dict[str, object]":
    """QP smoke: KEEPS the U5 regime mu0=3.5, nu1=4.9; reduces ONLY burn-in and sample size.

    n_burn 60000(panel default)/60000(fft default) -> 500 (the plan's "5e3->500"
    describes the eventual U5 (Task B5) full-run burn-in of 5e3, which differs
    from this v1/v2-inherited 60000 default -- B0 does not change that default,
    only the smoke reduction target); n_keep (panel 4000 default / fft's
    previously-hardcoded 8192) -> 2000. The regime is not switched to mu0=4/nu1=2.
    """
    return dict(mu0=MU0, nu1=4.90, n_burn=500, n_keep=2000)


def _manifest_append(path: str) -> None:
    """Append a producer-authored ``| artifact | version | sha256 |`` row to DESIGNG_MANIFEST (if set)."""
    manifest = os.environ.get("DESIGNG_MANIFEST")
    if not manifest:
        return
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    with open(manifest, "a") as fh:
        fh.write(f"| {path} | {os.path.basename(__file__)} | {h.hexdigest()} |\n")


def _resolve_qp_source_regime():
    """Task B5 (U5) full-run env contract: (mu0, nu1) for the U5 diagnostics
    come from the JSON at DESIGNG_QP_SOURCE (falls back to the LATEST_QP_VALUES
    binding -- exported by `.designg_bindings.sh` -- then to the default data path of
    either tree layout (output/data/designG in the authoring tree, codes/data/designG
    in the submitted bundle), so a bare run without bindings sourced still works).

    The only source JSON that currently exists (v1's
    designG_quasiperiodic_values.json) predates an explicit "nu1" field (only
    "nu1_NS"/"rho_at_4p9" -- the latter's name records evaluation AT nu1=4.9
    but is not itself a machine-readable "nu1" key), and a companion-repo
    clone may not bundle a source file at all -- both cases fall back to the
    U5 nominal circle regime nu1=4.90, documented via the returned status
    (never a silent mismatch).

    Returns (mu0, nu1, source_path, status) with status in
    {"read", "missing", "default_nu1"}.
    """
    # Candidates in priority order: the explicit variable, the bindings variable, then the
    # default data path of either tree layout (the authoring tree keeps its data one level
    # above codes/, output/data/designG; the submitted bundle keeps it under codes/data/designG).
    # The first EXISTING candidate wins, so a bindings path exported relative to the bundle
    # root still resolves when the README FULL sequence is run from codes/.
    cands = [os.environ.get("DESIGNG_QP_SOURCE"), os.environ.get("LATEST_QP_VALUES"),
             os.path.join(_OUT, "data", "designG", "designG_quasiperiodic_values.json"),
             os.path.join(_HERE, "data", "designG", "designG_quasiperiodic_values.json")]
    cands = [c for c in cands if c]
    src = next((c for c in cands if os.path.exists(c)), cands[0])
    try:
        with open(src) as fh:
            doc = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        print(f"[qp-source] WARNING: could not read DESIGNG_QP_SOURCE={src!r} ({exc}); "
              f"falling back to the hardcoded U5 regime mu0={MU0}, nu1=4.90", file=sys.stderr)
        return MU0, 4.90, src, "missing"
    mu0 = float(doc.get("mu0", MU0))
    if "nu1" in doc:
        return mu0, float(doc["nu1"]), src, "read"
    print(f"[qp-source] NOTE: {src!r} has no explicit 'nu1' field; using the U5 "
          f"nominal circle regime nu1=4.90 (rho_at_4p9's evaluation point)", file=sys.stderr)
    return mu0, 4.90, src, "default_nu1"


def main(argv=None):
    _preflight()

    parser = argparse.ArgumentParser(
        description="Design-G quasiperiodic illustration (mu0=3.5; unit-circle locus nu1=14/3; attracting circle observed numerically above it -- criticality not computed)."
    )
    parser.add_argument("--smoke", action="store_true", help="tiny/reduced run, full schema -> <name>_smoke.json")
    parser.add_argument("--only", choices=["quasiperiodic"], default=None,
                         help="run a single named block (only one exists pre-Task-B5/U5)")
    parser.add_argument("--public", action="store_true",
                         help="relax preflight to the public/lockfile check (companion-repo clone)")
    args = parser.parse_args(argv)

    paths = resolve_paths(args.smoke)
    figdir, datadir = paths["figdir"], paths["datadir"]
    os.makedirs(figdir, exist_ok=True); os.makedirs(datadir, exist_ok=True)

    # --smoke (B0): KEEPS mu0=3.5, nu1=4.9 (asserted below) -- only n_burn/n_keep
    # shrink, and only for THIS run; the default path below is untouched so a
    # non-smoke run stays bit-for-bit identical to v1 (n_burn=60000 both cases,
    # panel n_keep=4000, fft nfft=8192).
    if args.smoke:
        smoke = _smoke_params()
        assert smoke["mu0"] == MU0 and smoke["nu1"] == 4.90, "QP smoke must keep the U5 regime mu0=3.5, nu1=4.9"
        n_burn_panel, n_keep_panel = smoke["n_burn"], smoke["n_keep"]
        n_burn_fft, nfft = smoke["n_burn"], smoke["n_keep"]
    else:
        n_burn_panel, n_keep_panel = 60000, 4000   # orbit()'s own defaults (v1-identical)
        n_burn_fft, nfft = 60000, 8192              # v1's explicit FFT call (n_keep=8192)

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    # Type sizes (2026-10-17; widths as of 2026-11-12): the manuscript sets the panel at 0.72 and
    # the spectrum at 0.78 of its 6.27 in text width (the panel is 9 in wide at design size, the
    # spectrum 8 in), so 12 pt here prints at about 6 pt in the panel and 7.3 pt in the spectrum.
    plt.rcParams.update({"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12,
                         "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 11.5, "figure.titlesize": 13})
    cases = [(4.60, "stable focus (below $\\nu_1^{\\rm NS}$)"), (4.80, "circle (+3%)"),
             (5.10, "circle (+9%)"), (5.50, "invariant-circle breakup")]
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for ax, (nu1, label) in zip(axes.ravel(), cases):
        if nu1 < NU1_NS:
            # Stable focus (2026-10-17): the post-transient orbit sits at
            # the equilibrium to rounding, so plotting it gave a roundoff-scale cloud
            # with meaningless offset tick labels. Show the TRANSIENT instead: the
            # first 400 iterates from a visible offset spiral into the equilibrium
            # (checked to stay inside [0,1]^2; radius 0.18 -> 0.013). This changes
            # only this picture -- the data file (rho_at_4p9, theta_over_2pi, U5
            # block) and the other three panels are untouched.
            x0, y0 = coexistence_fixed_point(MU0, nu1)
            x, y = x0 + 0.15, y0 + 0.10; pts = []
            for _ in range(400):
                pts.append((x, y)); x, y = designG_step(x, y, MU0, nu1)
            xs, ys = np.array(pts).T
            ax.plot(xs, ys, "-", color="0.6", lw=0.5); ax.plot(xs, ys, "k.", ms=1.6)
            ax.plot([x0], [y0], "r+", ms=10)
            ax.set_title(f"$\\nu_1={nu1}$: {label}\n(transient: 400 iterates from $(x^*{{+}}0.15,\\ y^*{{+}}0.10)$)")
            ax.set_xlim(x0 - 0.20, x0 + 0.20); ax.set_ylim(y0 - 0.20, y0 + 0.20)
        else:
            xs, ys, (x0, y0) = orbit(nu1, n_burn=n_burn_panel, n_keep=n_keep_panel)
            if xs is None:
                ax.text(0.5, 0.5, f"escaped\n$\\nu_1$={nu1}", ha="center", transform=ax.transAxes); continue
            ax.plot(xs, ys, "k.", ms=0.4, alpha=0.5); ax.plot([x0], [y0], "r+", ms=10)
            ax.set_title(f"$\\nu_1={nu1}$: {label}")
            # Task B5 axis fix: tight (padded) data limits.
            x_span = xs.max() - xs.min(); y_span = ys.max() - ys.min()
            pad_x = 0.05*x_span if x_span > 0 else 1e-3
            pad_y = 0.05*y_span if y_span > 0 else 1e-3
            ax.set_xlim(xs.min()-pad_x, xs.max()+pad_x); ax.set_ylim(ys.min()-pad_y, ys.max()+pad_y)
        ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.grid(True, alpha=0.3)
        # equal aspect so a genuine circle LOOKS circular; plain tick labels (no offset).
        ax.set_aspect("equal", adjustable="box")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.tick_params(labelbottom=True, labelleft=True)
    fig.suptitle(r"At $\mu_0=3.5$: focus below $\nu_1^{\rm NS}=14/3$, invariant circle above, breakup at $\nu_1=5.5$")
    fig.tight_layout()
    panel_path = os.path.join(figdir, "designG_quasiperiodic_panel.png")
    fig.savefig(panel_path, dpi=200, bbox_inches="tight", pad_inches=0.02); plt.close(fig)  # 2026-10-17: trim the white margin
    _manifest_append(panel_path)

    xs, ys, (x0, y0) = orbit(4.90, n_burn=n_burn_fft, n_keep=nfft)
    sig = (xs - x0)[:nfft]; sp = np.abs(np.fft.rfft(sig*np.hanning(nfft))); fr = np.fft.rfftfreq(nfft)
    rho = fr[1 + np.argmax(sp[1:])]; theta2pi = math.acos((4.0-MU0)/4.0)/(2*math.pi)
    fig, ax = plt.subplots(figsize=(8, 4)); ax.semilogy(fr, sp+1e-12, "b-", lw=0.8)
    ax.axvline(rho, color="r", ls="--", lw=1, label=fr"measured $\rho\approx{rho:.4f}$")
    ax.axvline(theta2pi, color="green", ls=":", lw=1.2, label=fr"linear $\theta/2\pi={theta2pi:.4f}$")
    ax.set_xlabel("frequency (cycles per iterate)"); ax.set_ylabel(r"Fourier amplitude $|X(f)|$")
    ax.set_ylim(top=200.0 * float(sp.max())); ax.legend(loc="upper right")   # head room: the legend lies above every spectral line
    ax.set_title(r"Spectrum of $x_n-x^*$ at $\mu_0=3.5,\ \nu_1=4.9$: lines at $n\rho$")
    fig.tight_layout()
    fft_path = os.path.join(figdir, "designG_quasiperiodic_fft.png")
    fig.savefig(fft_path, dpi=200); plt.close(fig)
    _manifest_append(fft_path)

    # --- Task B5 (U5) diagnostics: a DEDICATED orbit (burn-in 5e3 / N=2e4 in
    # the full run; reuses the `smoke` dict from the --smoke branch above,
    # 500/2000, under --smoke) -- distinct from the panel's 60000/4000 and
    # this FFT block's 60000/8192, both of which stay UNCHANGED above
    # (v1-bit-for-bit for rho_at_4p9/theta_over_2pi and the 2 PNGs).
    if args.smoke:
        u5_mu0, u5_nu1 = float(smoke["mu0"]), float(smoke["nu1"])
        u5_n_burn, u5_n_keep = int(smoke["n_burn"]), int(smoke["n_keep"])
        u5_source_path, u5_source_status = None, "smoke"
    else:
        u5_mu0, u5_nu1, u5_source_path, u5_source_status = _resolve_qp_source_regime()
        u5_n_burn, u5_n_keep = 5000, 20000

    xs5, ys5, (x0_5, y0_5) = orbit(u5_nu1, n_burn=u5_n_burn, n_keep=u5_n_keep, mu0=u5_mu0)
    assert xs5 is not None, (
        f"[U5] orbit escaped [0,1]^2 at mu0={u5_mu0}, nu1={u5_nu1} "
        f"(n_burn={u5_n_burn}, n_keep={u5_n_keep}) -- cannot compute the QP diagnostics"
    )
    r5 = np.hypot(xs5 - x0_5, ys5 - y0_5)
    u5_rho = rotation_number(xs5, ys5, x0_5, y0_5)
    u5_residual = invariant_circle_residual(xs5, ys5, x0_5, y0_5, n_bins=360)
    u5_snr = fft_peak_snr(r5)
    u5_largest_lyap, u5_normal_lyap = benettin_spectrum(u5_mu0, u5_nu1, n_burn=u5_n_burn, n_iter=u5_n_keep)

    values = {
        "utc": datetime.now(timezone.utc).isoformat(),
        # Task B6: record the U5 regime mu0 (env-resolved by
        # _resolve_qp_source_regime()), not the hardcoded MU0 echo, so "mu0"
        # is authoritative alongside "nu1"=u5_nu1. No-op on current data
        # (u5_mu0 == MU0 == 3.5 in every path), but semantically correct and
        # what the U6 validator's qp_source_check compares against.
        "mu0": u5_mu0,
        "nu1_NS": NU1_NS,
        "rho_at_4p9": float(rho),
        "theta_over_2pi": float(theta2pi),
        # Task B5 (U5) diagnostics on the dedicated (mu0, nu1) circle -- see
        # the "Task B5 (U5)" comment block above `_jacobian` for the full
        # derivation. "nu1" records the regime actually used for these 5
        # fields (env-contract-resolved in full-run mode; hardcoded in
        # --smoke); "mu0" above equals it in every current path.
        "nu1": u5_nu1,
        "rho": u5_rho,
        "invariant_circle_residual": u5_residual,
        "fft_peak_snr": u5_snr,
        "largest_lyap": u5_largest_lyap,
        "normal_lyap": u5_normal_lyap,
        # B0 JSON-schema convention: top-level "_taxonomy" map (field name / wildcard
        # JSON-Pointer -> regression_anchor|theory_invariant|empirical_discovery).
        # Populated by Task B5 (U5): "rho" is a regression anchor (measured
        # rotation number); the other four are finite-sample/finite-time
        # empirical diagnostics -- NEVER theory_invariant, and NON-blocking
        # (no PASS/FAIL tolerance gate). largest_lyap in particular is the
        # tangent/neutral exponent (~0 on a genuine circle), not "transverse".
        "_taxonomy": {
            "rho": "regression_anchor",
            "largest_lyap": "empirical_discovery",
            "normal_lyap": "empirical_discovery",
            "invariant_circle_residual": "empirical_discovery",
            "fft_peak_snr": "empirical_discovery",
        },
        "_provenance": {
            "generated_by": "output/codes/" + os.path.basename(__file__),

            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "smoke": bool(args.smoke),
            "only_block": args.only,
            "n_burn_panel": n_burn_panel,
            "n_keep_panel": n_keep_panel,
            "n_burn_fft": n_burn_fft,
            "nfft": nfft,
            "n_burn_u5": u5_n_burn,
            "n_keep_u5": u5_n_keep,
            "qp_source_path": _rel(u5_source_path),
            "qp_source_status": u5_source_status,
            "output_fig_dir": _rel(figdir),
            "output_data_dir": _rel(datadir),
        },
    }
    stem = os.path.splitext(os.path.basename(__file__))[0]
    out_name = f"{stem}_smoke.json" if args.smoke else f"{stem}_values.json"
    out_json = os.path.join(datadir, out_name)
    json.dump(values, open(out_json, "w"), indent=2)
    _manifest_append(out_json)
    print(f"OK quasiperiodic figures written; rho(nu1=4.9)={rho:.4f}, theta/2pi={theta2pi:.4f}")
    print(f"OK U5 diagnostics (mu0={u5_mu0}, nu1={u5_nu1}, source={u5_source_status}): "
          f"rho={u5_rho:.6f}, largest_lyap={u5_largest_lyap:.6e}, normal_lyap={u5_normal_lyap:.6e}, "
          f"invariant_circle_residual={u5_residual:.6e}, fft_peak_snr={u5_snr:.4f}")
    print("wrote", out_json)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
