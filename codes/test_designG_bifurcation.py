"""Regression + sanity tests for the SHIPPED Design-G bifurcation generator.

As of 2026-07-02 every Design-G bifurcation figure is produced by the single script
designG_bifurcation_v2.py (the nu1-sweep diagnostics were folded in from the now-retired
designG_bifurcation_v1.py), so this suite imports only v2. The two test families below keep their
historical v1/v2 names but both now exercise functions that live in v2:
  * test_v2_* -- the paper's mu0-sweep extinction route (Figure 10 in v30, three panels since 2026-11-08; Figure 3 in v29, Figure 15 in v19, Figure 14 in v18; the number has moved with revisions)
    (sweep_mu0_xy / lyap_max / m0_lyapunov / detect_onsets). This is the
    load-bearing "prey chaos onset = predator extinction" generator.
  * test_v1_* -- the nu1-sweep generator across the unit-circle locus (lyapunov_spectrum / fp_seed),
    which produces designG_bifurcation_nu1.png and designG_lyapunov_nu1.png.

Run in the bundle's codes directory: python -m pytest test_designG_bifurcation.py -q (README.md gives the full command).
"""
import json, os
import numpy as np
import designG_bifurcation_v2 as bif2

_HERE = os.path.dirname(os.path.abspath(__file__))
# The shipped fixture lives one level above codes/ in the authoring tree (output/data/designG)
# and under codes/data/designG in the bundle; the first existing path wins and the
# regression assertions are REQUIRED (a missing fixture must fail, not skip).
_CANDIDATES = [os.path.join(os.path.dirname(_HERE), "data", "designG", "designG_extinction_bif_values.json"),
               os.path.join(_HERE, "data", "designG", "designG_extinction_bif_values.json")]
_STORED = next((p for p in _CANDIDATES if os.path.exists(p)), _CANDIDATES[0])
_GRID_STEP = 0.01   # detect_onsets scans mu0 on arange(1,4,0.01)

# ---------------------------------------------------------------------------
# v2: the shipped mu0-sweep extinction-route generator (previously UNTESTED)
# ---------------------------------------------------------------------------

def test_v2_lambda_perp_closed_form_at_mu4():
    # On M0={y=0} the transverse predator exponent is lambda_perp=<ln(nu1 x)>;
    # at mu0=4 the arcsine law gives the closed form ln(nu1/4). This anchors the
    # numeric m0_lyapunov against the manuscript's exact value.
    _, lperp = bif2.m0_lyapunov(4.0, 2.1, 0.2, 20_000, 40_000)
    assert abs(lperp - np.log(2.1 / 4.0)) < 0.02, \
        f"lambda_perp={lperp} != ln(2.1/4)={np.log(2.1/4.0)}"

def test_v2_onsets_coincide_and_match_stored():
    # Paper-facing headline: prey chaos onset == predator extinction onset, both
    # ~3.96 at nu1=2.1. Regression-lock detect_onsets against the shipped JSON
    # within the grid step -- the check the old (v1-only) suite never performed.
    chaos, ext = bif2.detect_onsets(2.1, 0.2, 1e-3)
    assert chaos is not None and ext is not None, \
        f"both onsets must be detected (got chaos={chaos}, ext={ext})"
    assert abs(chaos - ext) <= _GRID_STEP, \
        f"chaos onset {chaos} and extinction onset {ext} must coincide to the grid step"
    # Regression-lock against the shipped JSON (both tree layouts ship it; a missing fixture
    # is a packaging defect and FAILS here rather than silently skipping).
    assert os.path.exists(_STORED), f"shipped fixture missing: {_STORED}"
    stored = json.load(open(_STORED))
    assert abs(chaos - stored["chaos_onset_mu0"]) <= _GRID_STEP, \
        f"chaos onset {chaos} regressed from stored {stored['chaos_onset_mu0']}"
    assert abs(ext - stored["extinction_onset_mu0"]) <= _GRID_STEP, \
        f"extinction onset {ext} regressed from stored {stored['extinction_onset_mu0']}"

def test_v2_extinction_is_real_not_escape():
    # Anti-masquerade sanity (a guard against the seed defect):
    # in the extinction band the predator is genuinely extinguished (kept y all
    # below Y_EXT) while the prey is ALIVE and chaotic (x samples present, with
    # spread). Design-G keeps [0,1]^2 invariant, so a real orbit IS recorded --
    # extinction must be a true y-decline, never an empty/escaped-orbit artifact.
    _, X, Y = bif2.sweep_mu0_xy(np.array([3.98]), 2.1, 0.2, 1e-3, bif2.N_BURN, bif2.N_KEEP)
    assert X.size > 0, "prey samples must exist (orbit did not escape the box)"
    assert Y.size == X.size and Y.size > 0, "predator samples recorded alongside prey"
    assert Y.max() < bif2.Y_EXT, f"predator must be extinct in the band (y_max={Y.max()})"
    assert X.std() > 0.05, f"prey must be chaotic/alive in the band (x_std={X.std()})"

def test_v2_coexistence_below_onset():
    # Well below onset the rare predator invades (lambda_perp>0 on M0) and
    # REGULATES to a coexistence equilibrium: predator present above Y_EXT, i.e.
    # NOT extinct. Complements the extinction test so the onset is a real crossing.
    _, _, Y = bif2.sweep_mu0_xy(np.array([2.5]), 2.1, 0.2, 1e-3, bif2.N_BURN, bif2.N_KEEP)
    assert Y.size > 0 and Y.max() > bif2.Y_EXT, \
        f"predator must survive below onset (y_max={Y.max() if Y.size else 0.0})"

# ---------------------------------------------------------------------------
# v1: the shipped nu1-sweep generator (designG_bifurcation_nu1.png); observational checks only
# ---------------------------------------------------------------------------

def test_v1_logistic_lyapunov_at_mu4():
    # On M0 (y=0) the map reduces to the logistic x->mu0 x(1-x); at mu0=4, lambda_max=ln2.
    lam = bif2.lyapunov_spectrum(0.2, 0.0, 4.0, 3.9, n_burn=2000, n_iter=20000)
    assert abs(max(lam) - np.log(2.0)) < 0.02, f"logistic mu0=4 lambda_max={max(lam)} != ln2"

def test_v1_band_not_chaotic():
    # Just above the unit-circle locus at mu0=3 (nu1=6 -> 6.2) the attractor observed is the large
    # invariant ring: lambda_max ~ 0 (quasiperiodic) or < 0 (phase-locked), never > 0 (chaos).
    lam = bif2.lyapunov_spectrum(*bif2.fp_seed(3.0, 6.2), 3.0, 6.2, n_burn=75000, n_iter=20000)
    assert not np.isnan(lam).any(), "ring orbit must stay bounded at nu1=6.2"
    assert max(lam) < 5e-3, f"ring must be non-chaotic, got lambda_max={max(lam)}"

def test_v1_focus_below_ns_lambda_negative():
    lam = bif2.lyapunov_spectrum(*bif2.fp_seed(3.0, 5.8), 3.0, 5.8, n_burn=75000, n_iter=20000)
    assert max(lam) < -1e-3, f"stable focus should have lambda_max<0, got {max(lam)}"

def test_v1_mu0_3p5_band_quasiperiodic():
    # At mu0=3.5, just above the unit-circle locus 14/3 (nu1=4.76), the attractor observed is a
    # quasiperiodic circle => lambda_max ~ 0 (not chaotic, not a stable focus). This is an
    # observation; no criticality classification is asserted.
    lam = bif2.lyapunov_spectrum(*bif2.fp_seed(3.5, 4.76), 3.5, 4.76, n_burn=75000, n_iter=20000)
    assert not np.isnan(lam).any(), "ring must stay bounded at nu1=4.76"
    assert abs(max(lam)) < 5e-3, f"quasiperiodic lambda_max should be ~0, got {max(lam)}"



# ---------------------------------------------------------------------------
# 2026-10-18: zero-derivative propagation; periodic windows
# ---------------------------------------------------------------------------

def test_v2_m0_lyapunov_exact_critical_hit_is_minus_infinity():
    # Independent oracle: at mu0 = 2 the point x = 1/2 is the superstable fixed point of the
    # logistic map (f(1/2) = 1/2 exactly in float64 and f'(1/2) = 0), so the tangential
    # exponent is -inf, NOT 0 (the skipped-term bug); the transverse exponent is the
    # closed form ln(nu1 x*) = ln(2.1/2).
    lpar, lperp = bif2.m0_lyapunov(2.0, 2.1, 0.5, 20_000, 40_000)
    assert lpar == -np.inf, f"tangential exponent at an exact critical hit must be -inf, got {lpar}"
    assert abs(lperp - np.log(2.1 * 0.5)) < 1e-10, f"lambda_perp={lperp} != ln(1.05)={np.log(1.05)}"
    # a generic seed at the same mu0 converges to the same fixed point but never lands on it
    # exactly: its tangential exponent is finite (the orbit's derivative products stay > 0)
    lpar2, _ = bif2.m0_lyapunov(2.0, 2.1, 0.2, 20_000, 40_000)
    assert np.isfinite(lpar2) and lpar2 < 0.0


def test_v2_m0_lyapunov_negative_transverse_inside_periodic_windows():
    # Under the figure's axis protocol (nu1 = 2.1, x0 = 0.2, 20k burn-in, 40k retained)
    # the transverse exponent is negative at mu0 = 3.83 and 3.84, where the prey orbit is
    # PERIODIC (tangential exponent < 0) -- so "negative again only in the high-mu0 chaotic
    # band" was false; the manuscript now says "including periodic windows".
    for mu0, lpar_ref, lperp_ref in ((3.83, -0.3696402747005968, -0.1195264537578609),
                                     (3.84, -0.0443724829926532, -0.1446783469485086)):
        lpar, lperp = bif2.m0_lyapunov(mu0, 2.1, 0.2, 20_000, 40_000)
        assert lpar < 0.0 and lperp < 0.0, (mu0, lpar, lperp)
        assert abs(lpar - lpar_ref) < 1e-6 and abs(lperp - lperp_ref) < 1e-6, (mu0, lpar, lperp)


def test_v2_m0_lyapunov_boundary_hit_in_the_retained_window():
    # An orbit that reaches the boundary inside the retained window used to return
    # the average of the shortened prefix (+ln 1.05 here). Exact oracle at mu0 = 4, nu1 = 2.1 from
    # x0 = 1/2: the float orbit is 1/2 -> 1 -> 0 -> 0 ..., with transverse multipliers nu1*x =
    # 1.05, 2.1, 0, 0, ... and tangential |4(1 - 2x)| = 0, 4, 4, ...
    lpar, lperp = bif2.m0_lyapunov(4.0, 2.1, 0.5, 0, 1)
    assert lpar == -np.inf and abs(lperp - np.log(1.05)) < 1e-15
    lpar, lperp = bif2.m0_lyapunov(4.0, 2.1, 0.5, 0, 2)
    assert lpar == -np.inf and abs(lperp - 0.5 * (np.log(1.05) + np.log(2.1))) < 1e-15
    lpar, lperp = bif2.m0_lyapunov(4.0, 2.1, 0.5, 0, 10)
    assert lpar == -np.inf and lperp == -np.inf      # x = 0 from the third step on: a vanishing multiplier
    # from x0 = 1 every tangential term is ln 4 (x = 1, then 0) and the transverse exponent is -inf
    lpar, lperp = bif2.m0_lyapunov(4.0, 2.1, 1.0, 0, 5)
    assert abs(lpar - np.log(4.0)) < 1e-15 and lperp == -np.inf
    # an orbit that hits the boundary during burn-in has no retained window: nan, as before
    lpar, lperp = bif2.m0_lyapunov(4.0, 2.1, 0.5, 5, 10)
    assert np.isnan(lpar) and np.isnan(lperp)


# ---------------------------------------------------------------------------
# 2026-10-21: the "largest exponent" helpers must see BOTH
# directions. A single tangent started at (1,0) stays tangential on the invariant axis
# y = 0 and misses a larger transverse rate. Oracle: mu0 = 1.8, nu1 = 2.1, x0 = 4/9, y0 = 0
# (fixed point; multipliers 0.2 and 14/15, largest exponent ln(14/15)).
# ---------------------------------------------------------------------------
import math
import designG_onset_ensemble_v1 as ons


def test_v2_lyap_max_on_axis_fixed_point_is_the_larger_rate():
    assert abs(bif2.lyap_max(1.8, 2.1, 4/9, 0.0, 20_000, 20_000) - math.log(14/15)) < 1e-9
    l1, l2 = bif2.lyap_spectrum_2d(1.8, 2.1, 4/9, 0.0, 20_000, 20_000)
    assert abs(l1 - math.log(0.2)) < 1e-9 and abs(l2 - math.log(14/15)) < 1e-9


def test_v2_lyap_spectrum_on_axis_matches_m0_lyapunov():
    for mu0 in (3.83, 4.0):
        lpar, lperp = bif2.m0_lyapunov(mu0, 2.1, 0.2, 20_000, 40_000)
        l1, l2 = bif2.lyap_spectrum_2d(mu0, 2.1, 0.2, 0.0, 20_000, 40_000)
        assert abs(l1 - lpar) < 1e-9 and abs(l2 - lperp) < 1e-9


def test_v2_lyap_max_interior_first_exponent_unchanged():
    # measured with the earlier single-vector helper (2026-10-21, numpy 1.26.4 / numba 0.61.2)
    FROZEN_3P95 = -0.55868788642131
    l1, l2 = bif2.lyap_spectrum_2d(3.95, 2.1, 0.2, 1e-3, 20_000, 20_000)
    assert abs(l1 - FROZEN_3P95) < 1e-12 and l2 < l1
    assert bif2.lyap_max(3.95, 2.1, 0.2, 1e-3, 20_000, 20_000) == l1


def test_onset_lyap_max_2d_axis_oracle():
    assert abs(ons._lyap_max_2d(1.8, 2.1, 4/9, 0.0, ons.LYAP_BURN, ons.LYAP_ITER) - math.log(14/15)) < 1e-9


def test_v2_lyap_spectrum_exact_zero_tangential_derivative():
    # mu0 = 2: (1/2, 0) is a fixed point with tangential multiplier 0 and transverse multiplier 2.1/2
    # the second vector starts at (0, 1), not on the eigenline, so its first step carries a
    # transient of ln(sqrt(b^2 + d^2)/d) that the 5000-step average dilutes below 1e-4
    l1, l2 = bif2.lyap_spectrum_2d(2.0, 2.1, 0.5, 0.0, 0, 5000)
    assert l1 == -math.inf and abs(l2 - math.log(1.05)) < 1e-4
    assert bif2.lyap_max(2.0, 2.1, 0.5, 0.0, 0, 5000) == l2


def test_onset_lyap_spectrum_zero_jacobian_corner():
    # 2026-10-21: at (x0, y0) = (0, 1) the Jacobian is the zero matrix, so both
    # propagated directions are annihilated at the first step: both exponents are -inf (the convention of the
    # bifurcation helper), not "no data". A non-finite norm stays invalid (NaN); zero iterations stay no-data.
    l1, l2 = ons._lyap_spectrum_2d(4.0, 2.1, 0.0, 1.0, 0, 10)
    assert l1 == -math.inf and l2 == -math.inf
    assert ons._lyap_max_2d(4.0, 2.1, 0.0, 1.0, 0, 10) == -math.inf
    assert all(np.isnan(v) for v in ons._lyap_spectrum_2d(4.0, 2.1, float("nan"), 0.5, 0, 5))
    assert all(np.isnan(v) for v in ons._lyap_spectrum_2d(4.0, 2.1, 0.0, 1.0, 0, 0))


# ---- Figure 10, three panels (2026-11-08): panel (c) is the plotted orbit's own spectrum -------------
import matplotlib
import pytest

GZ = bif2._grid(3.80, 4.0, 0.0004)


def _at(m):
    return int(np.argmin(np.abs(GZ - m)))


@pytest.fixture(scope="module")
def zoom_spectrum():
    return bif2.orbit_spectrum(GZ)


@pytest.fixture(scope="module")
def zoom_rows():
    P, X, Y = bif2.sweep_mu0_xy(GZ, bif2.NU1, bif2.X0, bif2.Y0, bif2.N_BURN, bif2.N_KEEP)
    assert len(P) == len(GZ) * bif2.N_KEEP
    return X.reshape(len(GZ), bif2.N_KEEP), Y.reshape(len(GZ), bif2.N_KEEP)


@pytest.fixture(scope="module")
def fig10_d():
    g = bif2._grid(1.0, 4.0, 0.003)
    P, X, Y = bif2.sweep_mu0_xy(g, bif2.NU1, bif2.X0, bif2.Y0, bif2.N_BURN, bif2.N_KEEP)
    return bif2.fig10_data(P, X, Y, *bif2.detect_onsets(bif2.NU1, bif2.X0, bif2.Y0))


def test_fig10_panel_c_signs(zoom_spectrum):
    l1, _ = zoom_spectrum
    assert len(GZ) == 501
    assert l1[_at(3.9604)] < 0 and l1[_at(3.9608)] < 0 and l1[_at(3.9612)] > 0
    assert round(float(GZ[np.argmax(l1 > 0)]), 4) == 3.8496
    # On the 0.01 detector grid the largest exponent of this seed is positive at one point below the sustained
    # onset 3.96: at 3.85, in the chaotic part of the period-three window of the prey (caption of Figure 10 and
    # Section S12 of the Supplement).
    gl = bif2._grid(1.0, 3.95, 0.01)
    assert [round(float(m), 2) for m in gl
            if bif2.lyap_max(m, bif2.NU1, bif2.X0, bif2.Y0, 20_000, 20_000) > 0] == [3.85]


def test_fig10_protocols_agree_except_3p9088(zoom_spectrum, zoom_rows):
    l1, _ = zoom_spectrum
    lm = np.array([bif2.lyap_max(m, bif2.NU1, bif2.X0, bif2.Y0, 20_000, 20_000) for m in GZ])
    assert np.round(GZ[(l1 > 0) != (lm > 0)], 4).tolist() == [3.9088]
    i = _at(3.9088)
    assert l1[i] < 0 < lm[i]
    xs, ys = zoom_rows[0][i], zoom_rows[1][i]
    x, y, m = xs[-1], ys[-1], GZ[i]
    assert np.ptp(xs) < 1e-9 and np.ptp(ys) < 1e-9 and y > 1e-3
    assert abs(m * x * (1 - x) * (1 - y) - x) < 1e-9 and abs(bif2.NU1 * x * y * (1 - y) - y) < 1e-9


def test_fig10_panel_c_on_M0(zoom_spectrum, zoom_rows):
    l1, l2 = zoom_spectrum
    idx = np.flatnonzero(np.all(zoom_rows[1] < 1e-3, axis=1) & (GZ >= 3.96))
    assert len(idx) == 101
    # Left of the line, about half of the 50 display-grid points of [3.94, 3.96) are already extinct
    # (caption of Figure 10; 27 on the reference machine).
    w = (GZ >= 3.94 - 1e-9) & (GZ < 3.96 - 1e-9)
    assert int(w.sum()) == 50 and 20 <= int(np.all(zoom_rows[1] < 1e-3, axis=1)[w].sum()) <= 30
    for i in idx:
        lp, lt = bif2.m0_lyapunov(GZ[i], bif2.NU1, 0.2, 20_000, 40_000)
        assert abs(l1[i] - max(lp, lt)) <= 0.02 and abs(l2[i] - min(lp, lt)) <= 0.03
    i4 = _at(4.0)
    assert abs(l1[i4] - np.log(2)) <= 0.01 and abs(l2[i4] - np.log(2.1 / 4)) <= 0.01


def test_fig10_layout(fig10_d):
    fig, axes, ins = bif2.draw_fig10(fig10_d)
    assert len(axes) == 3 and sorted(ins) == ["a", "b", "c"]
    assert all(axes[0].get_shared_x_axes().joined(axes[0], a) for a in axes[1:])
    assert bif2.fig10_violations(fig, axes, ins) == []
    fig2, axes2, ins2 = bif2.draw_fig10(fig10_d, inset_bounds={"a": [0.07, 0.10, 0.40, 0.40]})   # the old place
    assert any(v.startswith("inset (a): covers") for v in bif2.fig10_violations(fig2, axes2, ins2))


def test_fig10_not_written_on_violation(tmp_path, monkeypatch, capsys):
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    monkeypatch.delenv("DESIGNG_MANIFEST", raising=False)
    out = tmp_path / "designG_bifurcation_mu0.png"
    out.write_bytes(b"the figure already in the manuscript")
    fig = plt.figure()
    assert bif2.save_fig10(fig, str(out), ["OVERLAP x <-> y"]) is False
    assert out.read_bytes() == b"the figure already in the manuscript"
    text = capsys.readouterr().out
    assert "Figure 10 layout check: 1 violation(s)" in text and "  - OVERLAP x <-> y" in text and "NOT written" in text
    assert bif2.save_fig10(fig, str(out), []) is True
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


# ---- Figure 10, room for the insets (2026-11-09): no room is a layout violation, not a crash ----------------------
def _fig10_made_up(gz, calm, onset):
    """Figure 10 input without a sweep: on the columns marked calm the orbit rests at one point with the predator
    alive, elsewhere it is spread; the full-range arrays are filler."""
    n, k = len(gz), bif2.N_KEEP
    rng = np.random.default_rng(0)
    Xm = np.where(calm[:, None], 0.55, rng.uniform(0.2, 0.9, (n, k)))
    Ym = np.where(calm[:, None], 0.12, rng.uniform(0.0, 0.2, (n, k)))
    gf = bif2._grid(1.0, 4.0, 0.03)
    return dict(P=np.repeat(gf, 5), X=rng.uniform(0.0, 1.0, 5 * len(gf)), Y=rng.uniform(0.0, 0.3, 5 * len(gf)),
                Pz=np.repeat(gz, k), Xz=Xm.ravel(), Yz=Ym.ravel(), gf=gf, gz=gz,
                l1f=np.full(len(gf), -0.3), l2f=np.full(len(gf), -0.6),
                l1z=np.where(calm, -0.2, 0.3), l2z=np.full(n, -0.5), chaos_onset=onset, ext_onset=onset)


def test_fig10_insets_skip_a_run_at_the_window_end(monkeypatch):
    # a calm run touching the left end of the zoom window has no data column on that side; the insets go in the
    # narrower run that has one on each side (until 2026-11-09 the wider run was taken and an assert stopped the program)
    import matplotlib.pyplot as plt
    monkeypatch.setattr(bif2, "ZOOM_LO", 3.80)
    g = bif2._grid(3.80, 4.0, 0.002)
    calm = np.zeros(len(g), bool)
    calm[0:40] = True
    calm[45:75] = True
    fig, axes, ins = bif2.draw_fig10(_fig10_made_up(g, calm, 3.96))
    for a, k in zip(axes, "abc"):
        (x0, _), (x1, _) = a.transData.inverted().transform(ins[k].get_window_extent().get_points())
        assert g[44] < x0 < x1 < g[75], (k, x0, x1)
    plt.close(fig)


def test_fig10_no_room_at_the_onset_is_refused_by_main(tmp_path, monkeypatch, capsys):
    # DESIGNG_ZOOM_LO=3.96, the onset itself: no column lies left of the onset line. main() refuses Figure 10 like any
    # layout violation (the reason printed, the file left unchanged), writes its other outputs and returns 1.
    monkeypatch.setattr(bif2, "ZOOM_LO", 3.96)
    monkeypatch.setattr(bif2, "FIGDIR", str(tmp_path / "figures"))
    monkeypatch.setattr(bif2, "DATADIR", str(tmp_path / "data"))
    monkeypatch.delenv("DESIGNG_MANIFEST", raising=False)
    (tmp_path / "figures").mkdir()
    fig10 = tmp_path / "figures" / "designG_bifurcation_mu0.png"
    fig10.write_bytes(b"the figure already in the manuscript")
    assert bif2.main() == 1
    assert fig10.read_bytes() == b"the figure already in the manuscript"
    assert sorted(os.listdir(tmp_path / "figures")) == sorted(
        ["designG_bifurcation_mu0.png", "designG_lyapunov_mu0.png", "designG_bifurcation_mu0_M0logistic.png",
         "designG_lyapunov_mu0_M0logistic.png", "designG_bifurcation_nu1.png", "designG_lyapunov_nu1.png"])
    assert sorted(os.listdir(tmp_path / "data")) == ["designG_bifurcation_values.json", "designG_extinction_bif_values.json"]
    text = capsys.readouterr().out
    assert "Figure 10 layout check: 1 violation(s)" in text and "  - no room for the insets" in text
    assert "NOT written" in text


def test_fig10_too_narrow_a_span_is_no_room(monkeypatch):
    # one calm column between two columns of data cannot hold an inset (until 2026-11-09: a negative inset width,
    # ValueError); draw_fig10 reports no room and leaves no figure open
    import matplotlib.pyplot as plt
    monkeypatch.setattr(bif2, "ZOOM_LO", 3.80)
    g = bif2._grid(3.80, 4.0, 0.004)
    calm = np.zeros(len(g), bool)
    calm[10] = True
    open_before = plt.get_fignums()
    with pytest.raises(bif2.Fig10NoRoom, match="no room for the insets"):
        bif2.draw_fig10(_fig10_made_up(g, calm, 3.96))
    assert plt.get_fignums() == open_before


# Script mode: one entry point, after every definition. All 25 tests run through pytest (the eight Figure 10
# tests need its fixtures), which prints one verdict at the end and gives the script its exit status.
if __name__ == "__main__":
    import sys
    sys.exit(pytest.main(["-q", "-p", "no:cacheprovider", os.path.abspath(__file__)]))
