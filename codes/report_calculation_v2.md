# Design-G validation report

- generated: 2026-11-14T04:59:50.983119+00:00
- mode: **full**   scope: ALL
- claims file: current_claims_designG_v1.json
- producer JSON dir: ../../repro_run/data   figure dir: ../../repro_run/figs
- producers found: designG_dynamics_v4, designG_ergodic_quadrature_v2, designG_extinction_bif, designG_onset_ensemble_v1, designG_period_two_certificate_v1, designG_quasiperiodic_v2, designG_stable_leaf_v1
- overall (structural): **PASS** (exit 0)
- manuscript readiness: **SUBMISSION_READY** [strict mode: a required change fails the run]

## Per-claim checks

| kind | name | actual | expected | tol | status | manuscript_change | ci | tex_ref |
|------|------|--------|----------|-----|--------|-------------------|----|---------|
| regression_anchor | ext_frac_nu1_1 | 1 | 1 | 0.0001 | PASS | no |  | supplementary_clean.tex:2580 |
| regression_anchor | ext_frac_nu1_1p5 | 1 | 1 | 0.0001 | PASS | no |  | supplementary_clean.tex:2581 |
| regression_anchor | ext_frac_nu1_2 | 0.226944 | 0.2269 | 0.0001 | PASS | no |  | supplementary_clean.tex:2247 |
| regression_anchor | ext_frac_nu1_3 | 0.0111111 | 0.0111 | 0.0001 | PASS | no |  | supplementary_clean.tex:2384 |
| regression_anchor | ext_frac_nu1_4 | 0 | 0 | 0.0001 | PASS | no |  | supplementary_clean.tex:2584 |
| theory_invariant | lambda_perp_nu1_2 | -0.693147 | -0.693147 | 1e-06 | PASS | no |  | supplementary_clean.tex:2167 |
| theory_invariant | lambda_perp_nu1_3 | -0.287682 | -0.287682 | 1e-06 | PASS | no |  | supplementary_clean.tex:2168 |
| theory_invariant | lambda_perp_nu1_4 | 3.55271e-15 | 0 | 1e-06 | PASS | no |  | supplementary_clean.tex:2169 |
| regression_anchor | jury_lower_edge | 1.85641 | 1.85641 | 1e-06 | PASS | no |  | supplementary_clean.tex:1518 |
| regression_anchor | arcsine_E_lnx | -1.38629 | -1.38629 | 0.0001 | PASS | no |  | supplementary_clean.tex:1058 |
| regression_anchor | jensen_gap_mln2 | -0.693144 | -0.693147 | 0.0001 | PASS | no |  | supplementary_clean.tex:1081 |
| theory_invariant | subaction_sup_ln34 | -0.287682 | -0.287682 | 1e-06 | PASS | no |  | supplementary_clean.tex:1283 |
| theory_invariant | nu1star_mu0_4 | 3.99999 | 4 | 0.001 | PASS | no |  | supplementary_clean.tex:1597 |
| regression_anchor | family_spread_nu1_2 | 5.84959e-05 | 5.8496e-05 | 1e-05 | PASS | - |  | (untracked) |
| regression_anchor | family_spread_nu1_4 | 0.00184731 | 0.00184731 | 0.0001 | PASS | - |  | (untracked) |
| regression_anchor | rho_quasiperiodic | 0.224101 | 0.2241 | 0.0005 | PASS | no |  | resultsAndDiscussions_clean.tex:178 |
| regression_anchor | sens_base_nu1_2 | 0.3287 | 0.3287 | 5e-05 | PASS | no | none | supplementary_clean.tex:2486 |
| empirical_discovery | sens_runif_nu1_2 | 0.325675 |  |  | OK_NONBLOCKING | - | wilson | (untracked) |
| empirical_discovery | finite_amp_decay_nu1_2_pre_underflow | -0.693893 |  |  | OK_NONBLOCKING | - |  | (untracked) |
| regression_anchor | absorb_before_tail_nu1_2 | 817 | 817 | 0.5 | PASS | no |  | supplementary_clean.tex:2439 |
| regression_anchor | first_zero_step_min_nu1_2 | 591 | 591 | 0.5 | PASS | no |  | supplementary_clean.tex:2440 |
| regression_anchor | first_zero_step_max_nu1_2 | 1567 | 1567 | 0.5 | PASS | no |  | supplementary_clean.tex:2440 |
| regression_anchor | absorb_before_tail_nu1_3 | 40 | 40 | 0.5 | PASS | no |  | supplementary_clean.tex:2440 |
| regression_anchor | first_zero_step_min_nu1_3 | 1650 | 1650 | 0.5 | PASS | no |  | supplementary_clean.tex:2440 |
| regression_anchor | first_zero_step_max_nu1_3 | 4244 | 4244 | 0.5 | PASS | no |  | supplementary_clean.tex:2441 |
| regression_anchor | sens_base_nu1_3 | 0.0459 | 0.0459 | 5e-05 | PASS | no | none | supplementary_clean.tex:2468 |
| regression_anchor | sens_base_nu1_4 | 0.0005 | 0.0005 | 5e-05 | PASS | no | none | supplementary_clean.tex:2420 |
| empirical_discovery | sens_eps1e3_nu1_4 | 0.0006 | 0.0006 |  | OK_NONBLOCKING | no | none | supplementary_clean.tex:2426 |
| empirical_discovery | sens_rlogy_nu1_2 | 0.809641 | 0.8096 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2477 |
| empirical_discovery | sens_rlogy_nu1_3 | 0.272737 | 0.2727 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2478 |
| empirical_discovery | sens_rlogy_nu1_4 | 0.00404794 | 0.004 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2478 |
| empirical_discovery | sens_fp_eps1e3_nu1_3 | 0.1181 | 0.1181 |  | OK_NONBLOCKING | no | none | supplementary_clean.tex:2466 |
| empirical_discovery | sens_fp_eps1e12_nu1_3 | 0.0474 | 0.0474 |  | OK_NONBLOCKING | no | none | supplementary_clean.tex:2466 |
| empirical_discovery | sens_rlogy_fp_censored_incl_nu1_2 | 0.9095 | 0.9095 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2481 |
| empirical_discovery | sens_rlogy_fp_censored_incl_nu1_3 | 0.695075 | 0.6951 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2481 |
| empirical_discovery | sens_rlogy_fp_censored_incl_nu1_4 | 0.607225 | 0.6072 |  | OK_NONBLOCKING | no | wilson | supplementary_clean.tex:2481 |
| regression_anchor | lambda_hat_nu1_1 | -1.38535 | -1.3853 | 5e-05 | PASS | no |  | supplementary_clean.tex:2356 |
| regression_anchor | lambda_hat_nu1_1p5 | -0.979883 | -0.9799 | 5e-05 | PASS | no |  | supplementary_clean.tex:2357 |
| regression_anchor | lambda_hat_nu1_2 | -0.692201 | -0.6922 | 5e-05 | PASS | no |  | supplementary_clean.tex:2358 |
| regression_anchor | lambda_hat_nu1_3 | -0.286736 | -0.2867 | 5e-05 | PASS | no |  | supplementary_clean.tex:2359 |
| regression_anchor | lambda_hat_nu1_3p5 | -0.132585 | -0.1326 | 5e-05 | PASS | no |  | supplementary_clean.tex:2360 |
| regression_anchor | lambda_hat_nu1_4 | 0.000945951 | 0.0009 | 5e-05 | PASS | no |  | supplementary_clean.tex:2361 |
| regression_anchor | lambda_hat_n_trajectories | 128 | 128 | 0.5 | PASS | no |  | methods_clean.tex:332 |
| regression_anchor | lambda_hat_n_noncollapsed | 125 | 125 | 0.5 | PASS | no |  | methods_clean.tex:332 |
| regression_anchor | lambda_hat_n_accepted | 118 | 118 | 0.5 | PASS | no |  | methods_clean.tex:332 |
| empirical_discovery | lambda_hat_sem | 0.00140198 | 0.0014 |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2340 |
| regression_anchor | sweep_chaos_onset_nu1_2p1 | 3.96 | 3.96 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:36 |
| regression_anchor | sweep_ext_onset_nu1_2p1 | 3.96 | 3.96 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:36 |
| empirical_discovery | lambda_hat_sd | 0.0152294 | 0.0152 |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2341 |
| regression_anchor | lambda_hat_nu1_1_6d | -1.38535 | -1.38535 | 5e-06 | PASS | no |  | supplementary_clean.tex:2691 |
| regression_anchor | lambda_hat_nu1_1p5_6d | -0.979883 | -0.979883 | 5e-06 | PASS | no |  | supplementary_clean.tex:2692 |
| regression_anchor | lambda_hat_nu1_2_6d | -0.692201 | -0.692201 | 5e-06 | PASS | no |  | supplementary_clean.tex:2693 |
| regression_anchor | lambda_hat_nu1_3_6d | -0.286736 | -0.286736 | 5e-06 | PASS | no |  | supplementary_clean.tex:2694 |
| regression_anchor | lambda_hat_nu1_3p5_6d | -0.132585 | -0.132585 | 5e-06 | PASS | no |  | supplementary_clean.tex:2695 |
| theory_invariant | nu1star_mu0_4_5d | 3.99999 | 3.99999 | 5e-06 | PASS | no |  | supplementary_clean.tex:2783 |
| theory_invariant | nu1_NS_mu0_3p5 | 4.66667 | 4.66667 | 1e-09 | PASS | no |  | supplementary_clean.tex:1957 |
| theory_invariant | theta_over_2pi_mu0_3p5 | 0.230053 | 0.230053 | 1e-09 | PASS | no |  | supplementary_clean.tex:1958 |
| empirical_discovery | qp_residual_nu1_4p9 | 0.00125677 | 0.00125677 |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:1998 |
| empirical_discovery | qp_largest_lyap_nu1_4p9 | -0.000112245 | -0.000112245 |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2000 |
| empirical_discovery | qp_normal_lyap_nu1_4p9 | -0.0495186 | -0.0495186 |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2001 |
| empirical_discovery | ons_median_sustained_0p01 | 3.96 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2519 |
| empirical_discovery | ons_median_sustained_0p002 | 3.962 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2519 |
| empirical_discovery | ons_min_sustained_0p01 | 3.94 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2519 |
| empirical_discovery | ons_max_sustained_0p01 | 4 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2519 |
| empirical_discovery | ons_min_sustained_0p002 | 3.944 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2520 |
| empirical_discovery | ons_max_sustained_0p002 | 3.996 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2520 |
| empirical_discovery | ons_seed0_sustained_0p01 | 3.96 |  |  | OK_NONBLOCKING | no |  | resultsAndDiscussions_clean.tex:36 |
| empirical_discovery | ons_seed0_sustained_0p005 | 3.955 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2522 |
| empirical_discovery | ons_seed0_sustained_0p002 | 3.958 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:2523 |
| empirical_discovery | ons_seed0_first_positive_0p01 | 3.85 |  |  | OK_NONBLOCKING | no |  | resultsAndDiscussions_clean.tex:55 |
| regression_anchor | tc_p_x | 0.502328 | 0.502328 | 1e-09 | PASS | no |  | methods_clean.tex:812 |
| regression_anchor | tc_p_y | 0.175974 | 0.175974 | 1e-09 | PASS | no |  | methods_clean.tex:812 |
| regression_anchor | tc_phip_x | 0.824009 | 0.824009 | 1e-09 | PASS | no |  | methods_clean.tex:812 |
| regression_anchor | tc_phip_y | 0.134027 | 0.134027 | 1e-09 | PASS | no |  | methods_clean.tex:812 |
| regression_anchor | tc_mult_lo | -0.632808 | -0.632808 | 1e-09 | PASS | no |  | methods_clean.tex:813 |
| regression_anchor | tc_mult_hi | 0.963659 | 0.963659 | 1e-09 | PASS | no |  | methods_clean.tex:814 |
| regression_anchor | tc_trace_lo | 0.330852 | 0.330852 | 1e-09 | PASS | no |  | supplementary_clean.tex:1423 |
| regression_anchor | tc_trace_hi | 0.330852 | 0.330852 | 1e-09 | PASS | no |  | supplementary_clean.tex:1423 |
| regression_anchor | tc_det_lo | -0.609811 | -0.609811 | 1e-09 | PASS | no |  | supplementary_clean.tex:1423 |
| regression_anchor | tc_det_hi | -0.609811 | -0.609811 | 1e-09 | PASS | no |  | supplementary_clean.tex:1423 |
| empirical_discovery | tc_basin_two_cycle | 0.585275 |  |  | OK_NONBLOCKING | no |  | supplementary_clean.tex:1439 |
| regression_anchor | enum_closest_gap | 0.0632898 | 0.0632898 | 1e-09 | PASS | no |  | supplementary_clean.tex:1286 |
| regression_anchor | bif_window_lperp_3p83 | -0.119526 | -0.119526 | 1e-09 | PASS | no |  | resultsAndDiscussions_clean.tex:78 |
| regression_anchor | bif_window_lperp_3p84 | -0.144678 | -0.144678 | 1e-09 | PASS | no |  | resultsAndDiscussions_clean.tex:78 |
| regression_anchor | tele_max_stop_step_1p00 | 570 | 570 | 0.5 | PASS | no |  | supplementary_clean.tex:2704 |
| regression_anchor | tele_max_log10_y_horizon_1p00 | -2122.17 | -2122.17 | 1e-06 | PASS | no |  | supplementary_clean.tex:2705 |
| regression_anchor | tele_max_stop_step_1p20 | 671 | 671 | 0.5 | PASS | no |  | supplementary_clean.tex:2707 |
| regression_anchor | tele_max_log10_y_horizon_1p20 | -1830.28 | -1830.28 | 1e-06 | PASS | no |  | supplementary_clean.tex:2708 |
| regression_anchor | tele_max_stop_step_1p30 | 752 | 752 | 0.5 | PASS | no |  | supplementary_clean.tex:2710 |
| regression_anchor | tele_max_log10_y_horizon_1p30 | -1666.08 | -1666.08 | 1e-06 | PASS | no |  | supplementary_clean.tex:2711 |
| empirical_discovery | leaf_psi2_nu1_2p1 | -424.66 |  |  | OK_NONBLOCKING | no |  | methods_clean.tex:483 |
| empirical_discovery | leaf_psi2_nu1_3 | -30286.9 |  |  | OK_NONBLOCKING | no |  | methods_clean.tex:483 |
| empirical_discovery | leaf_track_nu1_2p1 | 5.74122e-05 |  |  | OK_NONBLOCKING | no |  | methods_clean.tex:490 |
| empirical_discovery | leaf_track_nu1_3 | 0.00132484 |  |  | OK_NONBLOCKING | no |  | methods_clean.tex:490 |
| regression_anchor | leaf_rate40_nu1_2p1 | -0.839114 | -0.839114 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:140 |
| regression_anchor | leaf_rate40_nu1_3 | -0.482439 | -0.482439 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:141 |
| regression_anchor | leaf_rateN_nu1_2p1 | -0.637615 | -0.637615 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:141 |
| regression_anchor | leaf_rateN_nu1_3 | -0.28094 | -0.28094 | 1e-06 | PASS | no |  | resultsAndDiscussions_clean.tex:141 |

## Claim-coverage audit (every Design-G number -> claim id | exclusion)

- audit rows: 254  (73 -> claim id, 181 -> documented exclusion)
- audit integrity: OK

| tex | line(s) | literal | disposition |
|-----|---------|---------|-------------|
| supplementary_clean.tex | 2584 | 0.0000 | claim:ext_frac_nu1_4 |
| resultsAndDiscussions_clean.tex | 13 | -0.6931 (predicted ln(1/2)) | claim:lambda_perp_nu1_2 |
| supplementary_clean.tex | 2358 | -0.6931 (predicted ln(1/2)) | claim:lambda_perp_nu1_2 |
| supplementary_clean.tex | 1354,2359 | -0.2877 (predicted ln(nu1/4), nu1=3) | claim:lambda_perp_nu1_3 |
| supplementary_clean.tex | 2361 | 0.0000 (predicted ln(nu1/4), nu1=4) | claim:lambda_perp_nu1_4 |
| methods_clean.tex | 794 | -0.2877 / -0.287682 / ln(3/4) (worst-orbit / sub-action supremum) | claim:subaction_sup_ln34 |
| supplementary_clean.tex | 1272,1283,1354,2359 | -0.2877 / -0.287682 / ln(3/4) (worst-orbit / sub-action supremum) | claim:subaction_sup_ln34 |
| supplementary_clean.tex | 1283,2242 | -0.287682 / -0.2877 (ln(3/4) sub-action / all-orbit boundary) | claim:subaction_sup_ln34 |
| supplementary_clean.tex | 2168,2694,2697,2773,2774,2778 | -0.287682 (worst-orbit ln x / sub-action H) | claim:subaction_sup_ln34 |
| supplementary_clean.tex | 1046,1058,1063,1066 | -1.38629 / -1.386294 (-ln4 equilibrium-measure constant / arcsine E[ln x]) | claim:arcsine_E_lnx |
| supplementary_clean.tex | 1081,1082 | -0.69314 / -0.69315 (-ln2 Jensen gap) | claim:jensen_gap_mln2 |
| supplementary_clean.tex | 2691,2757,2760,2761,2762,2763,2764,2765,2766,2767,2768,2769,2770,2782 | -1.38629 (-ln4 arcsine / equilibrium potential on-support) | claim:arcsine_E_lnx |
| supplementary_clean.tex | 2167,2693,2759 | -0.693147/-0.693144 (-ln2 Jensen gap) | claim:jensen_gap_mln2 |
| abstract_clean.tex | 6 | ln2 (symbolic: the Jensen gap) | exclusion:symbolic_mention (jensen_gap_mln2; written symbolically here, the numeric value is checked where printed) |
| introduction_clean.tex | 39 | ln2 (symbolic: the Jensen gap) | exclusion:symbolic_mention (jensen_gap_mln2; written symbolically here, the numeric value is checked where printed) |
| methods_clean.tex | 308,677 | -ln4 (equilibrium potential) | exclusion:symbolic_mention (arcsine_E_lnx; written symbolically here, the numeric value is checked where printed) |
| supplementary_clean.tex | 950,954,1016,1163,1165,1167,1168,1174 | -ln4 (equilibrium potential) | exclusion:symbolic_mention (arcsine_E_lnx; written symbolically here, the numeric value is checked where printed) |
| conclusions_clean.tex | 10,33 | -ln2 (symbolic: the Jensen gap) | exclusion:symbolic_mention (jensen_gap_mln2; written symbolically here, the numeric value is checked where printed) |
| methods_clean.tex | 144,951 | 1.856 | claim:jury_lower_edge |
| supplementary_clean.tex | 1387 | 1.856 | claim:jury_lower_edge |
| methods_clean.tex | 154,155 | 1.856 / 1.85641 / 1.8564 (Jury lower edge) | claim:jury_lower_edge |
| supplementary_clean.tex | 1518 | 1.856 / 1.85641 / 1.8564 (Jury lower edge) | claim:jury_lower_edge |
| supplementary_clean.tex | 2701 | 1.85641 (Jury band lower edge predicted/computed) | claim:jury_lower_edge |
| methods_clean.tex | 154 | 3.999995 (Jury band upper edge measured) | exclusion:closed_form_band_edge (upper edge nu1->4; measured_upper_edge reproduced by the coexistence block; anchored via jury lower edge + block completeness) |
| supplementary_clean.tex | 1519 | 3.999995 (Jury band upper edge measured) | exclusion:closed_form_band_edge (upper edge nu1->4; measured_upper_edge reproduced by the coexistence block; anchored via jury lower edge + block completeness) |
| supplementary_clean.tex | 1516,2702 | 4 / 5e-6 (Jury upper edge) | exclusion:closed_form_band_edge |
| supplementary_clean.tex | 1597 | 4.000 (nu1*(4)) | claim:nu1star_mu0_4 |
| supplementary_clean.tex | 2783 | 4 / 3.99999 (nu1*(4)=exp(-I(4))) | claim:nu1star_mu0_4_5d (six-significant-digit form of the same I_family.nu1star_4 field; Section S8 quotes 4.000) |
| supplementary_clean.tex | 1613 | .../.../.../... (nu1* SRB thresholds at mu0=.../.../.../...), 3.74 non-chaotic | exclusion:frozen_verification_record (I_family SRB thresholds nu1star_ae; empirical many-IC SRB average, reproduced by the sibling I_family block; block-covered) |
| supplementary_clean.tex | 1599,1600,1601,1611 | 1.65/1.59/1.76/2.04/2.17/2.84 (nu1* many-initial-condition averages at mu0=3.60/3.70/3.74/3.90/3.95/3.99, the values that Figure 8 marks; 3.74 lies in a periodic window) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json I_family.family.<mu0>.nu1star_ae: 1.653286, 1.591643, 1.764682, 2.038466, 2.168660, 2.842489; reproduce (A)) |
| supplementary_clean.tex | 2787,2788,2789,2790 | 1.653/1.592/1.765/2.038 (nu1*(mu0) exp(-I) rows) | exclusion:frozen_verification_record (I_family SRB thresholds; sibling I_family.family.<mu0>.nu1star_ae; block-covered) |
| supplementary_clean.tex | 1050,1068,1069 | -0.51928/0.37645/0.376 (off-support potential quadrature vs closed form) | exclusion:frozen_verification_record (equilibrium_measure off-support potential; sibling block; block-covered) |
| supplementary_clean.tex | 2771,2772 | -0.51928/0.376453 (off-support potential rows) | exclusion:frozen_verification_record (equilibrium_measure off-support; block-covered) |
| supplementary_clean.tex | 1354,1360,1361 | -1.0002/-0.8541/-0.7642 (worst realised telescoping rate) & -0.1054/-0.0253 (bounds ln(3nu1/4)) | exclusion:frozen_verification_record (global_telescoping worst-realised-rate diagnostic; DYN global_telescoping block; bounds are closed-form ln(3 nu1/4); block-covered) |
| supplementary_clean.tex | 1217 | 0.9498/1.0 (max y / max x boundedness) | exclusion:frozen_verification_record (boundedness; DYN boundedness block; block-covered) |
| resultsAndDiscussions_clean.tex | 13 | measured Benettin lambda_hat col (.../.../-0.6922/.../.../0.0009) | claim:lambda_hat_nu1_2 + claim:lambda_hat_nu1_4 (4-decimal forms) |
| supplementary_clean.tex | 2356,2357,2358,2359,2360,2361 | measured Benettin lambda_hat col (-1.3853/-0.9799/-0.6922/-0.2867/-0.1326/0.0009) | claim:lambda_hat_nu1_1 + claim:lambda_hat_nu1_1p5 + claim:lambda_hat_nu1_2 + claim:lambda_hat_nu1_3 + claim:lambda_hat_nu1_3p5 + claim:lambda_hat_nu1_4 (4-decimal forms) |
| supplementary_clean.tex | 2691,2692,2693,2694,2695 | hat-lambda rows (-1.38535/-0.979883/-0.692201/-0.286736/-0.132585 computed vs predicted) | claim:lambda_hat_nu1_1_6d + claim:lambda_hat_nu1_1p5_6d + claim:lambda_hat_nu1_2_6d + claim:lambda_hat_nu1_3_6d + claim:lambda_hat_nu1_3p5_6d (six-digit forms of the same /lambda_hat/per_nu1/*/lambda_hat fields; Section 6 and the Supplement quote the four-decimal forms) |
| methods_clean.tex | 335 | 9.5e-4 (maximum absolute error of the Benettin estimates, in the caption of the on-axis exponent figure) | exclusion:derived_quantity (9.5e-4 = \|lambda_hat(4) - 0\|, the largest discrepancy, derived from the fields of the claim lambda_hat_nu1_4, which is quote-checked where Section 6 and the Supplement print it) |
| resultsAndDiscussions_clean.tex | 12 | 9.5e-4 / ... / ... (lambda abs-err & sem) | claim:lambda_hat_nu1_4 (0.0009); 9.5e-4 = \|lambda_hat(4) - 0\|, the largest discrepancy, derived from the same fields |
| supplementary_clean.tex | 2336,2340,2351,2356,2357,2358,2359,2360,2361 | 9.5e-4 / 1.40e-3 / 0.0009 (lambda abs-err & sem) | claim:lambda_hat_sem (1.40e-3 dispersion) + claim:lambda_hat_nu1_4 (0.0009); 9.5e-4 = \|lambda_hat(4) - 0\|, the largest discrepancy, derived from the same fields |
| supplementary_clean.tex | 1047,1062 | 5.6e-5 (largest on-support potential departure) | exclusion:diagnostic_error_bar (equilibrium_measure on-support departure; block-covered) |
| methods_clean.tex | 157 | 1.6309 / 0.667 / 0.250 (coexistence fixed point at mu0 = 4: modulus of its dominant multiplier at nu1 = 1.5, where it is unstable, and its coordinates (x*, y*) at nu1 = 2; at nu1 = 4 it is (0.5, 0.5)) | exclusion:frozen_verification_record (designG_dynamics_v4_values.json coexistence.per_nu1: eig_abs_max at 1.50; x_star and y_star at 2.00 and at 4.00, with x* = mu0/(mu0+nu1); closed-form-verified; block-covered) |
| supplementary_clean.tex | 1522,1524 | 1.6309 / 0.667 / 0.250 / 0.7676 (coexistence fixed point at mu0 = 4: modulus of its dominant multiplier at nu1 = 1.5, where it is unstable; its coordinates (x*, y*) and its dominant modulus at nu1 = 2) | exclusion:frozen_verification_record (designG_dynamics_v4_values.json coexistence.per_nu1: eig_abs_max at 1.50 and at 2.00; x_star and y_star at 2.00, with x* = mu0/(mu0+nu1); closed-form-verified; block-covered) |
| supplementary_clean.tex | 2698,2699,2700 | 1.249e-16/1.6653e-16 (FP residual, tr/det J err) | exclusion:diagnostic_error_bar (coexistence FP residual & trace/det closed-form error; DYN coexistence block; block-covered) |
| supplementary_clean.tex | 2704,2707,2710 | 9.9998e-251/9.9955e-251/9.9965e-251 (predator value at the stop step -> numerical floor) | exclusion:frozen_verification_record (global_telescoping predator underflow diagnostic; DYN global_telescoping block; block-covered) |
| supplementary_clean.tex | 2507 | -0.64 (ln(nu1/4) at nu1=2.1, lyapunov_mu0 figure) | exclusion:closed_form_figure_value (closed form ln(2.1/4) on the M0 Lyapunov figure; = lambda_perp theory at nu1=2.1) |
| supplementary_clean.tex | 47,175 | 3.5699456 (Feigenbaum logistic accumulation point) | exclusion:literature_constant (Feigenbaum period-doubling accumulation; not a Design-G measured output) |
| supplementary_clean.tex | 1596 | -1.386291 (I(4) arcsine quadrature) | claim:arcsine_E_lnx (I(4)=E_mu4[ln x]=-ln4; same quadrature value) |
| resultsAndDiscussions_clean.tex | 18,28,90 | 1e-3 (Design-G bifurcation & QP FIGURE parameters: nu1=2.1 coincident collapse, mu0=3.5 circle, y0=1e-3 seed, mu0 ranges; the row's other values moved to the appendix with the restructuring) | exclusion:figure_input_parameter (Design-G extinction-bifurcation / quasiperiodic figure INPUT parameters mu0,nu1,y0 axis ranges -- not measured claims) |
| supplementary_clean.tex | 2347 | 1e-3 (ceiling quoted for the difference between the Benettin estimates and ln(nu1/4) over the sampled conversion rates) | exclusion:rounded_ceiling (the recorded difference is 9.4595e-4 at every sampled nu1, the \|lambda_hat - ln(nu1/4)\| column of Table S9, derived from the lambda_hat fields; printed as 'better than 1e-3') |
| supplementary_clean.tex | 2691,2692,2693,2694,2695,2696,2701,2702,2703,2706,2709,2757,2758,2759,2760,2761,2762,2763,2764,2765,2766,2767,2768,2769,2770,2771,2772,2773,2775,2782,2787,2788,2789,2790 | 0.02 / 0.005 / 0.001 / 0.002 / 0.25 (tolerance column of Tables S9 and S10) | exclusion:input_parameter_or_tolerance (PASS/FAIL tolerances of the validation_table rows of designG_dynamics_v4_values.json and designG_ergodic_quadrature_v2_values.json; inputs, not measured values) |
| supplementary_clean.tex | 2760,2761,2762,2768,2769 | on-support potential -1.38629 / -1.38630 & abs-err 1.8373e-6/1.9595e-6 (equilibrium potential rows) | exclusion:diagnostic_error_bar (equilibrium_measure on-support potential = -ln4 (claim:arcsine_E_lnx); the ~1e-5 departures are integrable-singularity error bars; block-covered) |
| resultsAndDiscussions_clean.tex | 60 | 3.96 (plotted mu0-sweep detector: chaos onset = extinction onset, nu1=2.1; 0.01 grid, rounded to the nearest double) | claim:sweep_chaos_onset_nu1_2p1 + claim:sweep_ext_onset_nu1_2p1 |
| resultsAndDiscussions_clean.tex | 36,60 | 3.96 (plotted mu0-sweep detector: chaos onset = extinction onset, nu1=2.1; detector grid of step one hundredth, rounded to the nearest double) | claim:sweep_chaos_onset_nu1_2p1 + claim:sweep_ext_onset_nu1_2p1 |
| resultsAndDiscussions_clean.tex | 36,60 | 3.96 (sustained-chaos onset that the onset ensemble finds for the plotted seed at the coarsest grid resolution) | claim:ons_seed0_sustained_0p01 (direct value check, 2026-10-18; the ensemble file is a validator producer) |
| supplementary_clean.tex | 2519 | 3.96 / 4.00 (3.96 ([moved] at the finest resolution) median and [moved] to 4.00 ([moved] to [moved]) range of the 22 collapsing seeds' sustained-chaos onset; the row's other values moved to the appendix with the restructuring) | claim:ons_median_sustained_0p01 (direct value checks, 2026-10-18); the seed counts 22/24 remain regression-covered (reproduce (A)) |
| methods_clean.tex | 793 | -4.556 (second-derivative bound -41/9 of H) | exclusion:closed_form_constant (-41/9 = -4.5556) |
| supplementary_clean.tex | 1274 | -4.556 (second-derivative bound -41/9 of H) | exclusion:closed_form_constant (-41/9 = -4.5556) |
| supplementary_clean.tex | 2662,2777 | -4.5556 (second-derivative bound -41/9 of H) | exclusion:closed_form_constant (-41/9) |
| resultsAndDiscussions_clean.tex | 10 | 50000 / 5000 / 75000 / 1e-14 (Benettin protocol: retained iterates, acceptance window, burn-in, and the distance from x = 0 or x = 1 below which a trajectory is discarded) | exclusion:protocol_parameter (lambda_hat n_keep / window / n_burn and the endpoint guard 1e-14 of designG_dynamics_v4.py) |
| supplementary_clean.tex | 1047,1062,1064 | -1.38635 (on-support equilibrium-potential quadrature at r=0.5; departure 5.6e-5 from -ln 4) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json equilibrium_measure.on_support; reproduce (A)) |
| supplementary_clean.tex | 1050 | -0.519 (off-support potential at r=1.2, 3-decimal form) | exclusion:regression_covered (equilibrium_measure.off_support.1.20; reproduce (A)) |
| supplementary_clean.tex | 2356,2357,2360 | -1.3863 / -0.9808 / -0.1335 (predicted ln(nu1/4) column of the Benettin table) | exclusion:closed_form_constant (ln(nu1/4) at nu1=1, 1.5, 3.5) |
| resultsAndDiscussions_clean.tex | 51 | 1.909 (21/11: the mu0 at which the coexistence stock turns positive and ln(nu1(1-1/mu0)) changes sign on the fixed-point branch, at nu1=2.1) | exclusion:closed_form_constant (21/11) |
| supplementary_clean.tex | 2505 | 1.909 (21/11: sign change of ln(nu1(1-1/mu0)) on the fixed-point branch at nu1=2.1) | exclusion:closed_form_constant (21/11) |
| supplementary_clean.tex | 2299 | 0.062 / 0.001 (first-passage fraction at nu1=4 across the eps ladder) | exclusion:regression_covered (extinction_fraction.per_nu1.4.00.eps_sweep.*.fraction_first_passage; reproduce (A)) |
| methods_clean.tex | 439 | ... (1/107 < ... in the proof of cor:rare-predator) | exclusion:closed_form_constant |
| resultsAndDiscussions_clean.tex | 103,104 | ... / 0.2269 / 0.0111 / 0.0000 (extinction-fraction table copies) | claim:ext_frac_nu1_2 + claim:ext_frac_nu1_3 + claim:ext_frac_nu1_4 (the caption of Figure 11) |
| supplementary_clean.tex | 2383,2384,2393,2484,2580,2581,2582,2583,2584,2592 | 1.0000 / 0.2269 / 0.0111 / 0.0000 (extinction-fraction table copies) | claim:ext_frac_nu1_1 + claim:ext_frac_nu1_1p5 + claim:ext_frac_nu1_2 + claim:ext_frac_nu1_3 + claim:ext_frac_nu1_4 (Supplement copies of the values in the caption of Figure 11) |
| supplementary_clean.tex | 2420,2421,2425,2426,2430,2460,2461,2462,2468,2478,2486,2592,2607,2608,2609,2610,2611,2612,2613,2614,2615,2616,2617,2640,2641,2642,2644,2645,2646,2647,2648,2649,2650,2651,2652,2653,2654 | sensitivity table cells 0.3287 / 0.0459 / 0.0005 / 0.3336 / 0.0412 / 0.0000 / 0.3284 / 0.0488 / 0.0013 / 0.3278 / 0.0483 / 0.0009 / 0.0020 / 0.0006 / 0.3257 / 0.0482 / 0.0012 / 0.2727 / 0.0040 and Wilson bounds 0.3211 / 0.3303 / 0.0462 / 0.0504 / 0.0016 / 0.2665 / 0.2791 / 0.0032 / 0.0051 | claim:sens_base_nu1_2 + claim:sens_base_nu1_3 + claim:sens_base_nu1_4 + claim:sens_eps1e3_nu1_4 + claim:sens_rlogy_nu1_2 + claim:sens_rlogy_nu1_3 + claim:sens_rlogy_nu1_4 (table copies); every other cell and bound exclusion:regression_covered (sensitivity.*.*.fraction / ci_low / ci_high and convergence_band; reproduce (A)) |
| supplementary_clean.tex | 2144 | 3.55e-15 (\|lambda_perp^quad - ln(nu1/4)\|, the exact-skew quadrature error, 3-digit form) | exclusion:regression_covered (exact_skew.per_nu1.*; reproduce (A)) |
| supplementary_clean.tex | 2152,2167,2168,2169 | 3.5527e-15 / 7.0989e-4 / 2.5437e-4 / -0.287428 (exact-skew table: quadrature error, orbit s.e.m., orbit cross-check difference, orbit average at nu1=3) | exclusion:regression_covered (exact_skew.per_nu1.*; reproduce (A)) |
| supplementary_clean.tex | 2691,2692,2693,2694,2695,2696 | 9.4595e-4 (\|lambda_hat - ln(nu1/4)\| column, identical across nu1 because the prey orbits are shared) | exclusion:derived_quantity (from the lambda_hat claims) |
| supplementary_clean.tex | 2695 | -0.133531 (predicted ln(3.5/4)) | exclusion:closed_form_constant |
| supplementary_clean.tex | 2760,2761,2762,2763,2764,2765,2766,2767,2768,2769,2770 | 0.00 / 0.10 / 0.20 / 0.30 / 0.40 / 0.50 / 0.60 / 0.70 / 0.80 / 0.90 / 1.00 (the eleven abscissae r on the support [0, 1] of the equilibrium-potential rows of Table S10) | exclusion:protocol_parameter (table abscissae: the r values of the validation_table rows 'int ln\|x-r\| dmu4' of designG_ergodic_quadrature_v2_values.json) |
| supplementary_clean.tex | 2771 | 1.1102e-16 (recorded absolute error of the off-support potential row at r = 1.2 of Table S10: machine precision) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of the off-support potential row at r = 1.2; reproduce (A)) |
| supplementary_clean.tex | 2329,2559 | 1.856 (Jury lower edge, 3-decimal form) | claim:jury_lower_edge |
| supplementary_clean.tex | 2677,2678,2679,2680,2681,2682,2683,2684,2685,2686,2687,2688,2689,2690,2697,2704,2705,2707,2708,2710,2711,2774,2776,2777,2778,2779,2784 | 1e-9 / 1e-6 / 1e-12 (tolerance column of Tables S9 and S10) | exclusion:protocol_parameter (PASS/FAIL tolerances of the validation_table rows of designG_dynamics_v4_values.json and designG_ergodic_quadrature_v2_values.json) |
| methods_clean.tex | 332 | 128 / 125 / 118 (Benettin launches, non-collapsed, accepted) | claim:lambda_hat_n_trajectories + claim:lambda_hat_n_noncollapsed + claim:lambda_hat_n_accepted |
| supplementary_clean.tex | 1271 | 0.7500 (x* = 3/4, the sub-action maximiser) | exclusion:closed_form_constant (3/4) |
| methods_clean.tex | 797 | 8799 / ... (interior periodic orbits and periodic points of period <= 16, exact enumeration) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json periodic_enumeration.interior_orbits_exact / interior_periodic_points_exact; reproduce (A)) |
| supplementary_clean.tex | 1279 | 8799 / 130485 (interior periodic orbits and periodic points of period <= 16, exact enumeration) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json periodic_enumeration.interior_orbits_exact / interior_periodic_points_exact; reproduce (A)) |
| supplementary_clean.tex | 1279 | 8799 / 130485 (interior periodic orbits and periodic points of period <= 16) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json periodic_enumeration; reproduce (A)) |
| supplementary_clean.tex | 2780,2781 | 8799 / 130485 (interior periodic orbits and periodic points of period <= 16) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json periodic_enumeration; reproduce (A)) |
| supplementary_clean.tex | 1351 | 3600 (60 x 60 deterministic initial-condition grid) | exclusion:protocol_parameter (n_total of the extinction-fraction grid) |
| methods_clean.tex | 654,659 | 2.773 / 2.520 (4*3^(-1/3) and 4^(2/3)) | exclusion:closed_form_constant (invasion-exponent thresholds of prop:controls / lem:meanvar) |
| resultsAndDiscussions_clean.tex | 116 | 2.773 / ... (4*3^(-1/3) and 4^(2/3)) | exclusion:closed_form_constant (invasion-exponent thresholds of prop:controls / lem:meanvar) |
| methods_clean.tex | 344,357,416,437,439,442,444,446 | 3240 / 384 / 1536 / 108 / 107 / ... (integer constants of thm:basin-positive and cor:rare-predator; the row's other values moved to the appendix with the restructuring) | exclusion:closed_form_constant (c0 < 1/3240, D < 1/384, c0 <= ln(4/nu1)/1536, 1/768, 1/135 + 1/540 = 1/108, 1/107, eta < min/1000) |
| supplementary_clean.tex | 2202,2246,2260,2307,2326,2379,2565,2625 | 1e3 (terminal window of the terminal-persistence count: W(H) = min(10^3, H/2) iterates, equal to 10^3 at the nominal horizon) | exclusion:protocol_parameter (designG_dynamics_v4.py: the terminal window min(_EXT_TAIL, H/2) with _EXT_TAIL = 1000) |
| methods_clean.tex | 605,615,683,692 | 1.789 / 2.7734 / 2.5198 / -0.182 (4/sqrt5, 4*3^(-1/3), 4^(2/3), ln(2.1*4^(-2/3))) | exclusion:closed_form_constant (invasion-exponent laws, prop:controls / lem:meanvar) |
| supplementary_clean.tex | 1597 | -1.386294 (-ln 4, equilibrium potential) | claim:arcsine_E_lnx (6-digit form) |
| supplementary_clean.tex | 2247,2249,2293 | 0.2269 / 0.2294 (terminal vs first-passage fraction at nu1=2, rmk:first-passage) | claim:ext_frac_nu1_2 (0.2269); 0.2294 exclusion:regression_covered (per_nu1.2.00.fraction_first_passage; reproduce (A)) |
| supplementary_clean.tex | 1440,1445,2245,2261,2276 | 2e4 / 3600 / 400 (Iteration protocols: horizons, grid sizes; the row's other values moved to the appendix with the restructuring) | exclusion:protocol_parameter |
| supplementary_clean.tex | 2240,2241 | 1.000 / 0.9498 / -0.2867 (Iteration protocols: largest prey and predator values, Benettin exponent at nu1=3) | exclusion:regression_covered (boundedness.overall_max_x / overall_max_y) + claim:lambda_hat_nu1_3 |
| resultsAndDiscussions_clean.tex | 104 | ... / 3600 (n_extinct / n_total behind the extinction fractions) | claim:ext_frac_nu1_2 + claim:ext_frac_nu1_3 + claim:ext_frac_nu1_4 + exclusion:protocol_parameter (3600 = 60 x 60 grid) |
| supplementary_clean.tex | 2330,2385,2439,2565,2580,2581,2582,2583,2584 | 817 / 3600 (n_extinct / n_total behind the extinction fractions) | claim:ext_frac_nu1_2 + claim:ext_frac_nu1_3 + claim:ext_frac_nu1_4 + exclusion:protocol_parameter (3600 = 60 x 60 grid) |
| supplementary_clean.tex | 2200,2261,2414,2415,2422,2426,2430,2452,2453,2469,2472,2476,2592,2607,2608,2609,2610,2611,2612,2613,2614,2615,2616,2617,2625,2646,2647,2648,2649,2650,2651 | 100 / 200 / 400 / 2500 / 10000 / 40000 / 160000 / 1000 / 75000 (sensitivity-table protocol columns: grid side, cell count, horizon) | exclusion:protocol_parameter (sensitivity _params in designG_dynamics_v4_values.json) |
| supplementary_clean.tex | 2430,2477,2592,2617,2625,2643,2822 | 19022 / 0.8096 (valid log-uniform seeds and fraction, nu1=2) | claim:sens_rlogy_nu1_2 + exclusion:regression_covered (ic.rlogy n) |
| supplementary_clean.tex | 2430,2477,2617,2643 | 0.8096 / 0.8040 / 0.8152 (log-uniform fraction with Wilson interval, nu1=2) | claim:sens_rlogy_nu1_2 + exclusion:regression_covered (ic.rlogy ci_low/ci_high) |
| supplementary_clean.tex | 2167 | -0.692893 (orbit Birkhoff average of the exact-skew exponent, nu1=2) | exclusion:regression_covered (exact_skew.per_nu1.2.lambda_perp_orbit; reproduce (A)) |
| supplementary_clean.tex | 2692 | -0.980829 (predicted ln(1.5/4)) | exclusion:closed_form_constant (ln(nu1/4) at nu1=1.5) |
| supplementary_clean.tex | 2706,2707,2708,2709,2710,2711 | 1.20 / 1.30 (nu1 sample values of the all-orbit table) | exclusion:protocol_parameter (parameter values) |
| supplementary_clean.tex | 2727,2728,2732,2787,2788,2789,2790 | 3.60 / 3.70 / 3.74 / 3.90 (mu0 sample values of the I-family table) | exclusion:protocol_parameter (parameter values) |
| resultsAndDiscussions_clean.tex | 90 | 0.05 (smallest predator density of the deterministic grid, y0 >= 0.05) | exclusion:protocol_parameter (headline linspace grid on [0.05, 0.95], 60 x 60, _make_grid in designG_dynamics_v4.py; the midpoint lattice belongs to the sensitivity protocol) |
| supplementary_clean.tex | 1530 | 1.7e-16 (largest trace and determinant error of the closed-form Jacobian across the band) | exclusion:regression_covered (validation_table rows 'tr J vs closed form (max err)' / 'det J vs closed form (max err)'; reproduce (A)) |
| introduction_clean.tex | 34 | 3240 (c0 < 1/3240) | exclusion:closed_form_constant (thm:basin-positive constant) |
| resultsAndDiscussions_clean.tex | 91 | 3240 (c0 < 1/3240) | exclusion:closed_form_constant (thm:basin-positive constant) |
| supplementary_clean.tex | 2341,2351 | 1.52e-2 (across-seed standard deviation of the accepted Benettin means) | claim:lambda_hat_sd |
| supplementary_clean.tex | 2481,2592 | 0.6072 (first-passage rule over a censoring-inclusive denominator, nu1=4) | claim:sens_rlogy_fp_censored_incl_nu1_4 |
| resultsAndDiscussions_clean.tex | 100 | 1e3 (terminal window of the extinction-fraction protocol of Figure 11: an orbit counts as extinct only if it stays below the threshold throughout its final 1000 iterates) | exclusion:protocol_parameter (designG_dynamics_v4.py: the terminal window min(_EXT_TAIL, H/2) = 1000 at the horizon 2e4) |
| supplementary_clean.tex | 2423,2592,2625 | 500 (terminal window W of the sensitivity row with the shortest horizon: half of that horizon) | exclusion:protocol_parameter (designG_dynamics_v4.py: W(H) = min(_EXT_TAIL, H/2) = 500 at the horizon of the row horizon.1e3) |
| supplementary_clean.tex | 2208,2209 | 20000 / 40000 (companion on-axis exponents of Figure S3: burn-in and retained iterations) | exclusion:protocol_parameter (designG_bifurcation_v2.py: m0_lyapunov(m, NU1, 0.2, 20_000, 40_000)) |
| resultsAndDiscussions_clean.tex | 36 | 3.96 (sampled classification of the plotted sweep seed near the onset: from this grid point on, the retained iterates lie on the prey manifold, that is, the seed is extinct) | exclusion:protocol_observation (per-mu0 classification of seed (0.2, 1e-3) under the shipped sweep protocol, reproducible with designG_bifurcation_v2.sweep_mu0_xy / lyap_max; not a stored field) |
| resultsAndDiscussions_clean.tex | 28 | 400 (retained iterates per mu0 of the sweep of Figure 10) | exclusion:protocol_constant (N_KEEP = 400 in designG_bifurcation_v2.py; not a measured quantity) |
| resultsAndDiscussions_clean.tex | 99 | 2e4 (horizon of the extinction-fraction protocol of Figure 11) | exclusion:protocol_constant (_EXT_ANCHOR_PROTOCOL horizon 20000 in designG_dynamics_v4.py; not a measured quantity) |
| resultsAndDiscussions_clean.tex | 162,175 | ... / ... / 8192 / ... / ... / ... / 400 (protocol of Figures 13-14 and Section S10: panel burn-in and plotted iterates, FFT length, U5 burn-in and samples, phase bins, transient iterates) | exclusion:protocol_constant (iteration counts of designG_quasiperiodic_v2.py, recorded in its _provenance block; not measured quantities) |
| supplementary_clean.tex | 1986,1987,1988,1991,1995,1998 | 6e4 / 4e3 / 8192 / 5e3 / 2e4 / 360 / 400 (protocol of Figures 13-14 and Section S10: panel burn-in and plotted iterates, FFT length, U5 burn-in and samples, phase bins, transient iterates) | exclusion:protocol_constant (iteration counts of designG_quasiperiodic_v2.py, recorded in its _provenance block; not measured quantities) |
| resultsAndDiscussions_clean.tex | 163 | 0.15 / 0.10 (initial offset of the stable-focus transient panel) | exclusion:protocol_constant (plotting initial condition of designG_quasiperiodic_v2.py's focus panel) |
| supplementary_clean.tex | 1992 | 0.15 / 0.10 (initial offset of the stable-focus transient panel) | exclusion:protocol_constant (plotting initial condition of designG_quasiperiodic_v2.py's focus panel) |
| resultsAndDiscussions_clean.tex | 152,178,179 | ... / 0.2301 / 0.2241 / ... / ... / ... (quasiperiodic numbers of Figures 13-14 and Section S10) | claim:theta_over_2pi_mu0_3p5 + claim:rho_quasiperiodic |
| supplementary_clean.tex | 1983,1984,1996,1998,2000,2001 | 4.667 / 0.2301 / 0.2241 / 1.26e-3 / -1.12e-4 / -4.95e-2 (quasiperiodic numbers of Figures 13-14 and Section S10) | claim:nu1_NS_mu0_3p5 + claim:theta_over_2pi_mu0_3p5 + claim:rho_quasiperiodic + claim:qp_residual_nu1_4p9 + claim:qp_largest_lyap_nu1_4p9 + claim:qp_normal_lyap_nu1_4p9 |
| supplementary_clean.tex | 1957,1958,1960 | 4.667 / 0.2301 / 0.2241 (Neimark-Sacker locus at mu0=3.5, linear angle, measured rotation number, quoted in Section S10) | claim:nu1_NS_mu0_3p5 + claim:theta_over_2pi_mu0_3p5 + claim:rho_quasiperiodic |
| methods_clean.tex | 812,813,814 | 0.502328 / 0.175974 / 0.824009 / 0.134027 / -0.6328 / 0.9637 / ... / ... / ... / ... (certified two-cycle: points, multipliers, trace/determinant enclosures) | claim:tc_p_x + claim:tc_p_y + claim:tc_phip_x + claim:tc_phip_y + claim:tc_mult_lo + claim:tc_mult_hi |
| supplementary_clean.tex | 1409 | 0.502328 / 0.175974 / ... / ... / ... / ... / ... / ... / ... / ... (certified two-cycle: points, multipliers, trace/determinant enclosures) | claim:tc_p_x + claim:tc_p_y |
| supplementary_clean.tex | 1423 | ... / ... / ... / ... / ... / ... / 0.330851 / 0.330852 / -0.609812 / -0.609811 (certified two-cycle: points, multipliers, trace/determinant enclosures) | claim:tc_trace_lo + claim:tc_trace_hi + claim:tc_det_lo + claim:tc_det_hi |
| supplementary_clean.tex | 1416 | 1.3e-7 (certificate: computed contraction bound of E, printed as a ceiling) | exclusion:rounded_ceiling (designG_period_two_certificate_v1.json certificate.contraction_bound_q = 1.2695e-07, a bound computed by certify_period_two_v1.py for the row-sum norm of E and printed rounded up; the field is required by the validator's certificate block and test_designG_v4.py asserts it below 1e-6; the box radius is the protocol input and has its own rows) |
| supplementary_clean.tex | 1439 | 0.585 (fraction of the 200x200 lattice attracted to the two-cycle at nu1=46/25) | claim:tc_basin_two_cycle |
| methods_clean.tex | 830 | 1.79 / 1.77 / 1.716 (nu1 scan of the observed interior attractor: periodic attractors down to about 1.79; non-extinct interior trajectories with no period at most 64 down to about 1.77; one seed of the lattice at 1.716; numerical observation) | exclusion:regression_covered (designG_period_two_certificate_v1.json nu1_scan rows: periods_observed non-empty from 1.792 up, fraction_interior > 0 at 1.772 to 1.788 and, with 0.0025 = 1/400, at 1.716; reproduce (A); observation on a 20x20 lattice, not a classification) |
| supplementary_clean.tex | 1442,1443,1444 | 1.836 / 1.796 / 1.832 / 1.792 / 1.772 / 1.788 / 1.77 / 1.79 (nu1 scan of the observed interior attractor: period 2, 4, 8, aperiodic, extinct-below; numerical observation) | exclusion:regression_covered (designG_period_two_certificate_v1.json nu1_scan rows; reproduce (A); observation on a 20x20 lattice, not a classification) |
| supplementary_clean.tex | 1439 | 200 (200x200 interior midpoint lattice of the basin-fraction observation) | exclusion:protocol_constant (certify_period_two_v1.py basin_fraction grid) |
| supplementary_clean.tex | 1286 | 0.063 (gap between ln(3/4) and the largest orbit average over least periods 2..16) | claim:enum_closest_gap |
| supplementary_clean.tex | 944,1589 | 1.11 (citation locators: Example 1.11 of Saff's potential-theory survey; Sect. 1.11.8 of the arXiv version of Dobbs and Todd; not a quantity) | exclusion:not_a_quantity (citation locators) |
| supplementary_clean.tex | 1286 | 0.063 (closest periodic-orbit competitor gap below ln(3/4), period 16) | claim:enum_closest_gap |
| supplementary_clean.tex | 2704,2705,2707,2708,2710,2711 | 570 / 671 / 752 / -2122 / -1830 / -1666 (telescoping demo: latest floor-crossing step and largest log10 y_N at the common horizon, nu1 = 1, 1.2, 1.3) | claim:tele_max_stop_step_1p00 + claim:tele_max_stop_step_1p20 + claim:tele_max_stop_step_1p30 + claim:tele_max_log10_y_horizon_1p00 + claim:tele_max_log10_y_horizon_1p20 + claim:tele_max_log10_y_horizon_1p30 |
| supplementary_clean.tex | 2704,2705,2707,2708,2710,2711 | 1e-250 / 4000 / 1e-9 (telescoping demo protocol: floor, horizon N, tolerance of the clipped excess row) | exclusion:protocol_constant (designG_dynamics_v4.py global_telescoping safe_floor, --tele-iter, row tolerance) |
| supplementary_clean.tex | 1356,1359 | 570 / -2122 (telescoping demo at nu1=1.0 quoted in Section S7: latest stopping step and largest log10 y_N at the common horizon) | claim:tele_max_stop_step_1p00 + claim:tele_max_log10_y_horizon_1p00 |
| supplementary_clean.tex | 1356,1358 | 4000 (horizon N of the telescoping demo, quoted in Section S7) | exclusion:protocol_constant (designG_dynamics_v4.py --tele-iter) |
| supplementary_clean.tex | 1439 | 0.585 (fraction of the certificate script's 200x200 lattice attracted to the two-cycle at nu1=46/25, quoted in Remark S7.4) | claim:tc_basin_two_cycle |
| supplementary_clean.tex | 1439,2609,2646,2647,2648 | 200 (side of two lattices: the 200x200 interior lattice of the two-cycle basin observation of Remark S7.4, and the 200x200 deterministic grid of the sensitivity study in Tables S7 and S8) | exclusion:protocol_constant (line 1439: certify_period_two_v1.py basin_fraction grid, n_side = 200; lines 2609 and 2646 to 2648: the grid.200 configuration of the sensitivity study, sensitivity _params in designG_dynamics_v4_values.json) |
| supplementary_clean.tex | 1444,1445,1446 | 1.716 (nu1 scan: the one value below 1.77 with a non-extinct seed, 1 of the 400, no period <= 64 detected) | exclusion:regression_covered (designG_period_two_certificate_v1.json nu1_scan row nu1=1.716: fraction_interior 0.0025 = 1/400, n_interior_no_period_le_64 = 1; reproduce (A); observation on a 20x20 lattice, not a classification) |
| supplementary_clean.tex | 1445 | 400 (seeds of the 20x20 scan lattice) | exclusion:protocol_constant (certify_period_two_v1.py nu1_scan lattice 20x20) |
| supplementary_clean.tex | 1351 | 3600 (interior points of the telescoping grid quoted in Section S7) | exclusion:protocol_constant (designG_dynamics_v4.py global_telescoping n_ic = 60x60, recorded in the JSON) |
| supplementary_clean.tex | 2783 | 3.99999 (arcsine-quadrature value of nu1*(4), Table S10) | claim:nu1star_mu0_4_5d |
| supplementary_clean.tex | 2342 | 4e-15 (ceiling of the exact-skew quadrature error, quoted in Section S12; Table S4 prints the value) | exclusion:rounded_ceiling (Table S4 \|err\| = 3.5527e-15, designG_dynamics_v4_values.json exact_skew abs_err; regression-covered by reproduce (A)) |
| supplementary_clean.tex | 2758,2778,2779 | 1.6653e-16 / 1.1102e-16 (recorded absolute errors of three rows of Table S10: E_mu4[x], max orbit-avg ln x, argmax orbit point x0) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of the rows 'E_mu4[x]', 'max orbit-avg ln x (all interior orbits, period<= 16)' and 'argmax orbit point x0'; reproduce (A)) |
| resultsAndDiscussions_clean.tex | 102 | 0.05 / 0.95 (y0 range of the interior extinction-fraction grid, Figure 11 caption) | exclusion:protocol_parameter (designG_dynamics_v4.py extinction-fraction grid: y0 in [0.05, 0.95], 60 x 60 interior points) |
| supplementary_clean.tex | 1350,2486,2565 | 0.05 / 0.95 (range, in each coordinate, of the interior grid of the extinction-fraction table and of the telescoping diagnostic) | exclusion:protocol_parameter (designG_dynamics_v4.py _make_grid: linspace on [0.05, 0.95], 60 x 60 interior points; the grid of the extinction fractions and of the telescoping diagnostic) |
| supplementary_clean.tex | 2559 | 0.2294 / 0.0289 / 0.0172 (first-passage counts on the same orbits at nu1 = 2, 3, 4) | exclusion:regression_covered (designG_dynamics_v4_values.json extinction_fraction.per_nu1.2.00/3.00/4.00.fraction_first_passage; reproduce (A)) |
| supplementary_clean.tex | 2427,2466,2592,2625 | 0.0474 / 0.1181 (first-passage fraction at nu1 = 3 for the tightest and the loosest threshold) | claim:sens_fp_eps1e12_nu1_3 claim:sens_fp_eps1e3_nu1_3 |
| supplementary_clean.tex | 2475,2592,2625 | 20978 (left-censored initial conditions of ic.rlogy) | exclusion:regression_covered (ic.rlogy cells n_censored in designG_dynamics_v4_values.json; reproduce (A)) |
| supplementary_clean.tex | 2152 | 256 (orbits M of the Monte-Carlo orbit average, Table S4 caption) | exclusion:protocol_parameter (designG_dynamics_v4_values.json exact_skew.M = 256) |
| supplementary_clean.tex | 2698,2699,2700 | 1e-10 (tolerance of the residual rows of the Table S9) | exclusion:input_parameter_or_tolerance (designG_dynamics_v4_values.json validation_table: tol = 1e-10 of the residual rows) |
| supplementary_clean.tex | 2701 | 3.5394e-6 (absolute error of the Jury lower-edge row) | exclusion:regression_covered (designG_dynamics_v4_values.json validation_table row 'Jury band lower edge' abs_err; reproduce (A)) |
| supplementary_clean.tex | 2757,2759,2760,2770,2782 | 3.4657e-6 (absolute error of the arcsine-moment, Jensen-gap, potential r = 0 and r = 1, and I(4) rows) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of those rows; reproduce (A)) |
| supplementary_clean.tex | 2763,2767 | 6.8765e-7 (absolute error of the potential rows at r = 3/10 and 7/10) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err; reproduce (A)) |
| supplementary_clean.tex | 2764,2766 | 2.6131e-8 (absolute error of the potential rows at r = 4/10 and 6/10) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err; reproduce (A)) |
| supplementary_clean.tex | 2765 | -1.38635 / 5.6128e-5 (potential at r = 1/2: quadrature value and absolute error) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json equilibrium_measure.on_support.0.50 and its validation_table abs_err; reproduce (A)) |
| supplementary_clean.tex | 2773 | 1.6681e-13 (absolute error of the sub-action row) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of sup_x H(x); reproduce (A)) |
| supplementary_clean.tex | 2783 | 1.3863e-5 (absolute error of the nu1*(4) row) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of nu1*(4) = exp(-I(4)); reproduce (A)) |
| resultsAndDiscussions_clean.tex | 36 | 1e-2 (step of the detector grid of the mu0 sweep) | exclusion:protocol_parameter (designG_bifurcation_v2.py detect_onsets: the grid _grid(1.0, 4.0, 0.01), which is the coarsest resolution of designG_onset_ensemble_v1.json) |
| supplementary_clean.tex | 2783 | 0.01 (tolerance of the nu1*(4) row) | exclusion:input_parameter_or_tolerance (designG_ergodic_quadrature_v2_values.json validation_table: tol = 0.01) |
| supplementary_clean.tex | 2787,2788,2789,2790 | 1.65329 / 1.59164 / 1.76468 / 2.03847 (many-initial-condition orbit-average estimates of nu1* at the four sampled mu0 of the Table S10 rows) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json I_family.family.3.60/3.70/3.74/3.90.nu1star_ae; reproduce (A)) |
| supplementary_clean.tex | 2787,2788,2789,2790 | 2.8642e-4 / 3.5654e-4 / 3.1828e-4 / 4.6591e-4 (absolute errors of those rows) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of the I_family rows; reproduce (A)) |
| methods_clean.tex | 821 | 1e-10 (radius of the certificate box B) | exclusion:protocol_parameter (certify_period_two_v1.py: box radius 1e-10 about the rational centre) |
| supplementary_clean.tex | 1411 | 1e-10 (radius of the certificate box B) | exclusion:protocol_parameter (certify_period_two_v1.py: box radius 1e-10 about the rational centre) |
| supplementary_clean.tex | 2296,2395 | 1e-2 / 1e-8 (largest and smallest threshold of the ladder quoted for the terminal fraction) | exclusion:protocol_parameter (designG_dynamics_v4_values.json extinction_fraction eps_sweep: thresholds 1e-2 ... 1e-8) |
| supplementary_clean.tex | 1355 | 1e-250 (stopping threshold of the telescoping demonstration) | exclusion:protocol_parameter (designG_dynamics_v4_values.json global_telescoping: safe floor 1e-250) |
| methods_clean.tex | 160 | 1e-16 (order of magnitude of the trace and determinant errors of the simulated Jacobian) | exclusion:diagnostic_error_bar (order-of-magnitude statement: designG_dynamics_v4_values.json validation_table rows 'tr J vs closed form (max err)' and 'det J vs closed form (max err)', abs_err 1.6653e-16 each) |
| supplementary_clean.tex | 1705,1741,1843 | 1.14 (Theorem 1.14 of Vidiella et al., a citation locator in the comparison subsection; not a quantity) | exclusion:not_a_quantity (citation locator: Theorem 1.14 of the predecessor chapter) |
| methods_clean.tex | 475 | 0.3669 (base point p = sin^2(pi(sqrt2-1)/2) of the stable-leaf illustration) | exclusion:figure_input_parameter (base point of Figure 5) |
| methods_clean.tex | 487 | 1e-6 (height t = 10^{-6} of the matched-product orbit in Figure 5) | exclusion:figure_input_parameter (height of the illustrated orbit in Figure 5) |
| resultsAndDiscussions_clean.tex | 139 | 1e-6 (height t = 10^{-6} of the matched-product orbit in Figure 5) | exclusion:figure_input_parameter (height of the illustrated orbit in Figure 5) |
| resultsAndDiscussions_clean.tex | 138,141 | 1e5 (exact-orbit horizon N = 10^5 of Figure 12) | exclusion:protocol_parameter (horizon of the exact-orbit running rate in Figure 12) |
| supplementary_clean.tex | 2543 | 1e5 (exact-orbit horizon N = 10^5 of Figure 12) | exclusion:protocol_parameter (horizon of the exact-orbit running rate in Figure 12) |
| resultsAndDiscussions_clean.tex | 122 | 1e-4 / 1e-5 (heights at which the Lyapunov--Perron iteration of Figure 5 stops converging) | exclusion:protocol_observation (convergence edge of the floating-point leaf computation; not an estimate of c_0) |
| resultsAndDiscussions_clean.tex | 142 | -0.644 (closed form ln(2.1/4) drawn in Figure 12) | exclusion:closed_form_figure_value (closed form ln(nu1/4) at nu1=2.1, the limit shown in Figure 12) |
| resultsAndDiscussions_clean.tex | 142 | -0.288 (closed form ln(3/4), the almost-everywhere value printed in the running-rate figure's caption) | exclusion:closed_form_figure_value (closed form ln(nu1/4) at nu1=3, the limit shown in the running-rate figure) |
| methods_clean.tex | 489 | 107 (the bound 1/107 of cor:rare-predator, drawn as the band of Figure 5(b)) | exclusion:closed_form_constant (the corollary's tracking bound 1/107) |
| supplementary_clean.tex | 6,1731,2661,2741 | 0.09 / 0.37 / 0.46 / 0.19 / 0.36 / 0.36 / 0.35 / 0.39 (column widths of the Supplement's tables as fractions of the line width; layout, not a quantity) | exclusion:not_a_quantity (LaTeX column widths of tab:correspondence, tab:prior-comparison and the two closed-form consistency tables) |
| methods_clean.tex | 483 | -425 (psi_2, nu1=2.1, Figure 5(a)) | claim:leaf_psi2_nu1_2p1 |
| methods_clean.tex | 483 | -3.03e4 (psi_2, nu1=3, Figure 5(a)) | claim:leaf_psi2_nu1_3 |
| methods_clean.tex | 490 | 5.7e-5 (max \|ratio - 1\|, nu1=2.1, Figure 5(b)) | claim:leaf_track_nu1_2p1 |
| methods_clean.tex | 490 | 1.3e-3 (max \|ratio - 1\|, nu1=3, Figure 5(b)) | claim:leaf_track_nu1_3 |
| resultsAndDiscussions_clean.tex | 140 | -0.839 (running rate at N=40, nu1=2.1, Figure 12) | claim:leaf_rate40_nu1_2p1 |
| resultsAndDiscussions_clean.tex | 141 | -0.482 (running rate at N=40, nu1=3, Figure 12) | claim:leaf_rate40_nu1_3 |
| resultsAndDiscussions_clean.tex | 141 | -0.638 (running rate at N = 10^5, nu1=2.1, Figure 12) | claim:leaf_rateN_nu1_2p1 |
| resultsAndDiscussions_clean.tex | 141 | -0.281 (running rate at N = 10^5, nu1=3, Figure 12) | claim:leaf_rateN_nu1_3 |
| resultsAndDiscussions_clean.tex | 108 | 1e-2 / 1e-8 (largest and smallest threshold of the ladder quoted for the terminal fraction) | exclusion:protocol_parameter (designG_dynamics_v4_values.json extinction_fraction eps_sweep: thresholds 1e-2 ... 1e-8) |
| supplementary_clean.tex | 2395 | 1e-8 (smallest threshold of the ladder quoted for the terminal fraction) | exclusion:protocol_parameter (designG_dynamics_v4_values.json extinction_fraction eps_sweep: thresholds 1e-2 ... 1e-8) |
| supplementary_clean.tex | 2439 | 19000 (the iterate after which the terminal window begins) | exclusion:protocol_parameter (tail_start = H - W recorded in designG_dynamics_v4_values.json; the loop advances the orbit and then tests the zero-based index k >= tail_start, so the window holds y_19001, ..., y_20000; a protocol constant, not a measurement) |
| supplementary_clean.tex | 2440 | 591 / 1567 (first-exact-zero iterate range, nu1=2) | claim:first_zero_step_min_nu1_2 + claim:first_zero_step_max_nu1_2 |
| supplementary_clean.tex | 2440,2441 | 1650 / 4244 (first-exact-zero iterate range, nu1=3; the count 40 shares the sentence) | claim:first_zero_step_min_nu1_3 + claim:first_zero_step_max_nu1_3 + claim:absorb_before_tail_nu1_3 |
| supplementary_clean.tex | 2479 | 0.007 (bound on the Wilson half-widths of the log-uniform fractions) | exclusion:derived_quantity ((ci_high-ci_low)/2 of the ic.rlogy cells: 0.0056/0.0063/0.0009) |
| supplementary_clean.tex | 2481,2592 | 0.9095 / 0.6951 / 0.6072 (first-passage rule over a censoring-inclusive denominator) | claim:sens_rlogy_fp_censored_incl_nu1_2 + claim:sens_rlogy_fp_censored_incl_nu1_3 + claim:sens_rlogy_fp_censored_incl_nu1_4 (each censored seed counted as a crossing at n=0; replaces the superseded v3 values 0.9095/0.6926/0.5955) |
| resultsAndDiscussions_clean.tex | 28,41,90 | 1e-3 / 3.80 (inputs of Figure 10: the predator seed y0 of the plotted mu0 sweep, and the left end of the mu0 range of the zoomed panels) | exclusion:figure_input_parameter (designG_bifurcation_v2.py: the plotted seed X0 = 0.2, Y0 = 1e-3 and ZOOM_LO = 3.80; inputs, not measured values) |
| supplementary_clean.tex | 1985 | 1e-3 (offset from the coexistence equilibrium, in each coordinate, at which the post-transient orbits of Figures 13 and 14 start) | exclusion:figure_input_parameter (designG_quasiperiodic_v2.py orbit(): the start (x* + 1e-3, y* + 1e-3); an input, not a measured value) |
| supplementary_clean.tex | 518,645,646,647,873,876,922,927,928,929 | 3240 / 384 / 1536 / 768 / 108 / 107 / 135 / 540 / 1000 (integer constants of thm:basin-positive and cor:rare-predator) | exclusion:closed_form_constant (c0 < 1/3240, D < 1/384, c0 <= ln(4/nu1)/1536, 1/768, 1/135 + 1/540 = 1/108, 1/107, eta < min/1000) |
| resultsAndDiscussions_clean.tex | 35,36,55 | 3.85 / 3.86 / 3.95 (sampled classification of the plotted sweep seed near the onset: extinct at 3.85, coexisting again at 3.86 to 3.95, extinct from 3.96) | exclusion:protocol_observation (per-mu0 classification of seed (0.2, 1e-3) under the shipped sweep protocol, reproducible with designG_bifurcation_v2.sweep_mu0_xy / lyap_max; not a stored field) |
| resultsAndDiscussions_clean.tex | 63,64 | 3.9604 / 3.9608 / 3.9612 (display-grid points right of the onset line: extinct with lambda_max < 0 at the first two, positive from the third) | exclusion:protocol_observation (seed (0.2, 1e-3) on the 4e-4 display grid under the shipped protocol; reproducible with designG_bifurcation_v2.sweep_mu0_xy / lyap_max; not a stored field) |
| supplementary_clean.tex | 2505 | 1.909 (21/11: sign change of ln(nu1(1-1/mu0)) on the fixed-point branch at nu1=2.1) | exclusion:closed_form_constant (21/11) |
| resultsAndDiscussions_clean.tex | 78 | -0.12 / -0.14 / 3.83 / 3.84 (transverse exponent inside periodic windows of the mu0 sweep) | claim:bif_window_lperp_3p83 + claim:bif_window_lperp_3p84 |
| supplementary_clean.tex | 2519 | 3.96 (median, over the 22 collapsing seeds of the onset ensemble, of the sustained-chaos onset on the coarsest grid) | claim:ons_median_sustained_0p01 (direct value check, 2026-10-18) |
| supplementary_clean.tex | 2519,2520 | 3.96 (3.962 at the finest resolution) median and 3.94 to 4.00 (3.944 to 3.996) range of the 22 collapsing seeds' sustained-chaos onset | claim:ons_median_sustained_0p01 claim:ons_median_sustained_0p002 claim:ons_min_sustained_0p01 claim:ons_max_sustained_0p01 claim:ons_min_sustained_0p002 claim:ons_max_sustained_0p002 (direct value checks, 2026-10-18); the seed counts 22/24 remain regression-covered (reproduce (A)) |
| resultsAndDiscussions_clean.tex | 36 | 3.96 / ... / ... (ensemble sustained-chaos onset of the plotted seed at the three grid resolutions) | claim:ons_seed0_sustained_0p01 (direct value checks, 2026-10-18; the ensemble file is now a validator producer) |
| supplementary_clean.tex | 2522,2523 | 3.955 / 3.958 (sustained-chaos onset that the ensemble finds for the plotted seed at the two finer grid resolutions; its onset at the coarsest resolution is printed in the Results) | claim:ons_seed0_sustained_0p005 claim:ons_seed0_sustained_0p002 (direct value checks, 2026-10-18; the ensemble file is a validator producer) |
| supplementary_clean.tex | 2526,2527 | 3.950000000000003 / 3.9499999999999797 (the two float64 values of the decimal grid point 3.95 built from the anchors 1.0 and 3.0) | exclusion:arithmetic_demonstration (double-precision representation; the scripts now round every sweep grid to the nearest double) |
| supplementary_clean.tex | 2526,2527 | 3.950000000000003 / 3.9499999999999797 (the two float64 values of the decimal grid point 3.95 built from the anchors 1.0 and 3.0) | exclusion:arithmetic_demonstration (double-precision representation; the scripts now round every sweep grid to the nearest double) |
| supplementary_clean.tex | 2337 | 5000 (Benettin protocol: the window of the stability test; its burn-in and retained length are printed in the Results) | exclusion:protocol_parameter (designG_dynamics_v4.py check_lambda_hat: its window, the --window option, default 5000) |
| supplementary_clean.tex | 2230 | 0.9498/1.0 (max y / max x boundedness) | exclusion:frozen_verification_record (boundedness; DYN boundedness block; block-covered) |
| supplementary_clean.tex | 2820,2824 | -512 (the AVX-512 instruction-set name in the reproducibility paragraph, parsed as a signed number; not a quantity) | exclusion:not_a_quantity (instruction-set name) |
| supplementary_clean.tex | 2430,2476,2477,2592,2616,2617,2625,2820,2821,2822,2823,2824 | 11147 / 40000 / 19022 / 5189 / 5188 (CPU-dispatch effect on the log-uniform predator seeds, external measurement) | exclusion:external_measurement (11,147 one-ULP changes and the 5189 count are external measurements on an AVX-512 host and are not reproducible on the authors' AVX2 host; 40,000 is the sample size, 19,022 the valid count and 5188 = fraction x n of the shipped ic.rlogy nu1=3 cell, all regression-covered by reproduce (A)) |
| methods_clean.tex | 485 | 5.25 / 4.25 (decades of height plotted in Figure 5(a): from the input grid's 10^{-9} to the last converged height) | exclusion:protocol_observation (span of the plotted heights of Figure 5(a); log10(1.778e-4/1e-9) = 5.25 at nu1=2.1, log10(1.778e-5/1e-9) = 4.25 at nu1=3) |
| supplementary_clean.tex | 1703,1741,1843 | 1.13 (Theorem 1.13 of Vidiella et al., a citation locator in the comparison subsection; not a quantity) | exclusion:not_a_quantity (citation locator: Theorem 1.13 of the predecessor chapter) |
| supplementary_clean.tex | 429 | 799 (page of the Proposition in Section 3 of Alexander, Yorke, You and Kan, a citation locator in Remark S2.2; not a quantity) | exclusion:not_a_quantity (citation locator: p. 799 of Alexander et al.) |
| supplementary_clean.tex | 432 | 811 (page of the Theorem in Section 7 of Alexander, Yorke, You and Kan, a citation locator in Remark S2.2; not a quantity) | exclusion:not_a_quantity (citation locator: p. 811 of Alexander et al.) |
| supplementary_clean.tex | 433 | 2.19 (Theorem 2.19 of Ashwin, Buescu and Stewart, a citation locator in Remark S2.2; not a quantity) | exclusion:not_a_quantity (citation locator: Theorem 2.19 of Ashwin, Buescu and Stewart) |
| supplementary_clean.tex | 1803 | 1.0851 / 1.0997 (bounds of the chaotic regime of the template antecedent published by Lopez-Ruiz and Fournier-Prunaret, Section 2, quoted in the paragraph 'Template antecedents'; not quantities of this paper) | exclusion:literature_constant (published regime bounds of the template antecedent, Lopez-Ruiz and Fournier-Prunaret, Section 2; quoted, not Design-G outputs) |
| methods_clean.tex | 855 | 832 / -833 (pages 832--833 of Avila and Moreira, a citation locator for the account of Benedicks and Carleson in the chaotic-family corollary; the en dash reads as -833; not quantities) | exclusion:not_a_quantity (citation locator: pp. 832--833 of Avila and Moreira) |
| resultsAndDiscussions_clean.tex | 40,62 | 4e-4 (mu0 step of the display grid of the zoomed main panels of the bifurcation figure) | exclusion:protocol_parameter (designG_bifurcation_v2.py: _grid(ZOOM_LO, 4.0, 0.0004); not a stored field) |
| supplementary_clean.tex | 2441 | 0.28 (lower bound of the geometric decay rate \|ln(nu1/4)\| at nu1 <= 3: ln(4/3) printed rounded down) | exclusion:closed_form_constant (\|ln(nu1/4)\| >= ln(4/3) = 0.2877 for nu1 <= 3; printed as the rounded-down bound 0.28) |
| resultsAndDiscussions_clean.tex | 58 | 1e-3 (extinction threshold of the mu0-sweep detector: the predator is deemed extinct when every retained iterate lies below it) | exclusion:protocol_parameter (designG_bifurcation_v2.py: Y_EXT = 1e-3; not a stored field) |
| resultsAndDiscussions_clean.tex | 30 | 20000 (Figure 10 panel (c): iterations averaged after the burn-in) | exclusion:protocol_parameter (designG_bifurcation_v2.py N_AVG, orbit_spectrum) |
| supplementary_clean.tex | 2207 | 20000 / 75000 / 501 / 3.9088 (Figure 10 panel (c) against the 20000-iteration burn-in) | exclusion:protocol_observation (orbit_spectrum and lyap_max on the 4e-4 display grid; recomputed by test_fig10_protocols_agree_except_3p9088) |
| methods_clean.tex | 490 | 5.0e-4 (max \|ratio - 1\| over the steps n <= 30, nu1=3, Figure 5(b)) | exclusion:derived_quantity (designG_stable_leaf_v1_values.json nu1_3.matched.ratio: the largest \|ratio - 1\| over n <= 30 is 4.9885e-4, at n = 28; recomputed by test_matched_product_stays_within_the_corollary_band; not a stored field) |
| resultsAndDiscussions_clean.tex | 28 | 75000 (burn-in of the mu0 sweep of Figure 10) | exclusion:protocol_constant (N_BURN = 75000 in designG_bifurcation_v2.py; not a measured quantity) |
| resultsAndDiscussions_clean.tex | 42 | 3.94 (left end of the stretch of the display grid of Figure 10 on which about half of the points left of the line are extinct) | exclusion:protocol_observation (sweep_mu0_xy of designG_bifurcation_v2.py on _grid(3.80, 4.0, 0.0004) for the seed (0.2, 1e-3): the kept window stays below Y_EXT at 27 of the 50 points of [3.94, 3.96); recomputed by test_fig10_panel_c_on_M0; not a stored field) |
| resultsAndDiscussions_clean.tex | 78 | 2.04 (the sampled threshold nu1* at mu0 = 3.90, the value Remark S8.2 lists) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json I_family.family.3.90.nu1star_ae = 2.038466; reproduce (A)) |
| resultsAndDiscussions_clean.tex | 99 | 1e-3 (extinction threshold of the protocol of Figure 11) | exclusion:protocol_constant (_EXT_ANCHOR_PROTOCOL eps = 1e-3 in designG_dynamics_v4.py; not a measured quantity) |
| supplementary_clean.tex | 928 | 0.01 (1/107 < 0.01 in the proof of cor:rare-predator) | exclusion:closed_form_constant |
| supplementary_clean.tex | 1602 | 75000 / 50000 / 5e-3 (protocol of the I(mu0) averages: burn-in, retained iterates, threshold of the chaotic flag) | exclusion:protocol_parameter (designG_ergodic_quadrature_v2.py I_family: --n-burn 75000, --n-keep 50000, chaotic = (the average of ln\|f'\| exceeds 5e-3)) |
| supplementary_clean.tex | 2209 | 1e-2 / 3e-3 (mu0 steps of the grids of Figure S3 and of the full-range insets of Figure 10) | exclusion:protocol_parameter (designG_bifurcation_v2.py: _grid(1.0, 4.0, 0.01) and _grid(1.0, 4.0, 0.003)) |
| supplementary_clean.tex | 2336 | 128 / 0.6180339887 / 1e-14 (Benettin launches of Figure 4: number of trajectories, the step of their seeds x0 = frac(1/2 + m phi), and the distance from x = 0 or x = 1 below which a trajectory is discarded) | exclusion:protocol_parameter (designG_dynamics_v4.py lambda_hat: 128 launches, x0 = (0.5 + 0.6180339887498949 (m + 1)) mod 1, endpoint guard 1e-14) |
| supplementary_clean.tex | 2337 | 5000 / 5e-3 (stability test of the Benettin protocol: window length and the floor of the allowed drift between consecutive windows) | exclusion:protocol_parameter (designG_dynamics_v4.py lambda_hat: window = 5000, stab_tol = max(2 standard errors, 5e-3)) |
| supplementary_clean.tex | 2421,2462 | 1e4 (order of magnitude: the nu1 = 4 entries of the sensitivity table are a few orbits in 10^4) | exclusion:derived_quantity (the fractions 0.0000 to 0.0013 of the nu1 = 4 cells, read as counts per 10^4 orbits; the cells themselves are regression-covered, reproduce (A)) |
| supplementary_clean.tex | 2516 | 0.05 / 0.95 / 1e-4 / 1e-2 (laws of the 23 drawn seeds of the onset ensemble: x0 uniform on [0.05, 0.95], y0 log-uniform on [1e-4, 1e-2]) | exclusion:protocol_parameter (designG_onset_ensemble_v1.py: rng.uniform(0.05, 0.95) and 10 ** rng.uniform(-4, -2)) |
| supplementary_clean.tex | 2529 | 3.85 (the least mu0 of the detector grid at which the largest exponent of the plotted sweep seed is positive) | exclusion:protocol_observation (lyap_max of designG_bifurcation_v2.py on _grid(1.0, 4.0, 0.01) for the seed (0.2, 1e-3) at nu1 = 2.1: positive at 3.85 and at no other point below 3.96; recomputed by test_fig10_panel_c_signs; not a stored field) |
| supplementary_clean.tex | 2545 | 1e-14 / 1e-15 / 200 / 1e-6 (leaf computation of Figure 5: relative and absolute stopping tolerance, largest number of sweeps, height of the matched orbit) | exclusion:protocol_parameter (designG_stable_leaf_v1.py lp_fixed_point: tol = 1e-14, floor = 1e-15, max_sweeps = 200; T_MATCH = 1e-6) |
| supplementary_clean.tex | 2565 | 1e-3 (extinction threshold eps of the extinction-fraction protocol of Table S6, as its caption states it) | exclusion:protocol_constant (_EXT_ANCHOR_PROTOCOL eps = 1e-3 in designG_dynamics_v4.py; not a measured quantity) |
| supplementary_clean.tex | 2466 | 1e-3 (loosest extinction threshold eps of the threshold rows of the sensitivity table, quoted with the first-passage count at nu1 = 3) | exclusion:protocol_parameter (designG_dynamics_v4.py sensitivity configurations: eps = 1e-3 in the row eps.1e-3) |
| supplementary_clean.tex | 2517,2522 | 1e-2 (coarsest of the three mu0 grid resolutions of the onset ensemble) | exclusion:protocol_parameter (designG_onset_ensemble_v1.py RESOLUTIONS = (0.01, 0.005, 0.002)) |
| resultsAndDiscussions_clean.tex | 20,21 | 1.05 / 0.53 (predator multipliers at nu1 = 2.1 under full prey chaos: the arithmetic mean nu1<x> = 2.1/2 and the geometric mean nu1/4) | exclusion:closed_form_constant (2.1 x 1/2 = 1.05 with the arcsine mean <x> = 1/2; 2.1/4 = 0.525, printed as 0.53) |
| supplementary_clean.tex | 2580,2581,2582,2583,2584,2677,2678,2679,2680,2681,2682,2683,2684,2685,2686,2687,2688,2691,2692,2693,2694,2695,2696,2703,2704,2705 | 1.00 / 1.50 / 2.00 / 3.00 / 3.50 / 4.00 (nu1 values that label the rows of Tables S6 and S9) | exclusion:protocol_parameter (parameter values: the nu1 grids of designG_dynamics_v4.py) |
| supplementary_clean.tex | 2775,2779 | 0.75 (x = 3/4 in Table S10: the maximiser of the sub-action H and the starting point x0 of the orbit with the largest average, predicted and computed columns) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table rows 'argmax H(x)' and 'argmax orbit point x0'; the predicted value is the closed form 3/4; reproduce (A)) |
| supplementary_clean.tex | 2784 | 1.33333 (all-orbit bound mu0/(mu0-1) = 4/3 at mu0 = 4, predicted and computed columns of Table S10) | exclusion:closed_form_constant (4/3; the computed column is the validation_table row 'all-orbit bound mu0/(mu0-1) at mu0=4' of designG_ergodic_quadrature_v2_values.json) |
| supplementary_clean.tex | 2246,2379 | 1e-3 (extinction threshold eps of the extinction-fraction protocol of Figure 11 and Table S6, where the running text states it) | exclusion:protocol_constant (_EXT_ANCHOR_PROTOCOL eps = 1e-3 in designG_dynamics_v4.py; not a measured quantity) |
| supplementary_clean.tex | 2775 | 2.5e-7 (recorded absolute error of the argmax H(x) row of Table S10) | exclusion:regression_covered (designG_ergodic_quadrature_v2_values.json validation_table abs_err of the row 'argmax H(x)'; reproduce (A)) |
| supplementary_clean.tex | 2415,2426,2427,2428,2453,2454,2467,2592,2607,2608,2609,2610,2611,2612,2613,2614,2615,2616,2617,2652,2653,2654 | 1e-3 / 1e-6 / 1e-9 / 1e-12 (extinction thresholds eps of the sensitivity study of Tables S7 and S8: the base value 1e-6 and the three swept values) | exclusion:protocol_parameter (designG_dynamics_v4.py sensitivity configurations: eps = 1e-6 in the base, grid, horizon and initial-condition rows; 1e-3, 1e-9 and 1e-12 in the eps rows; sensitivity _params in designG_dynamics_v4_values.json) |
| supplementary_clean.tex | 1440,2322,2323 | 2e4 / 1e-12 (two-cycle basin observation of the certificate script: number of steps after which the terminal state is classified, and the predator level below which an orbit counts as extinct) | exclusion:protocol_constant (certify_period_two_v1.py basin_fraction: n_iter = 20000, extinct = (y < 1e-12); recorded in the basin_fraction_nu1_46_25 block of designG_period_two_certificate_v1.json) |
| supplementary_clean.tex | 2771,2772 | 1.20 / 2.00 (the two abscissae r off the support [0, 1] of the equilibrium-potential rows of Table S10) | exclusion:protocol_parameter (table abscissae: the r values of the two off-support validation_table rows 'int ln\|x-r\| dmu4' of designG_ergodic_quadrature_v2_values.json) |
| supplementary_clean.tex | 2200,2245,2261,2380,2437,2565 | 2e4 (nominal horizon H of the extinction-fraction protocol of Figure 11 and Table S6) | exclusion:protocol_parameter (_EXT_ANCHOR_PROTOCOL horizon = 20000 in designG_dynamics_v4.py) |
| supplementary_clean.tex | 2203 | 20000 (largest exponent of the mu0 sweep and of the onset ensemble: burn-in and retained iterations) | exclusion:protocol_parameter (designG_bifurcation_v2.py: lyap_max(m, nu1, x0, y0, 20_000, 20_000); designG_onset_ensemble_v1.json lyap_burn = lyap_iter = 20000) |
| supplementary_clean.tex | 2517,2522 | 5e-3 / 2e-3 (the two finer mu0 grid resolutions of the onset ensemble; where the multiplication sign is braced the validator's number reader reads each as the power 1e-3 after a one-digit factor) | exclusion:protocol_parameter (designG_onset_ensemble_v1.py RESOLUTIONS = (0.01, 0.005, 0.002); on line 2517 the multiplication sign is braced, so the validator's number reader reads 5x10^-3 and 2x10^-3 there as the power 1e-3 after a one-digit factor) |
| supplementary_clean.tex | 2152 | 75000 / 1000 (Monte-Carlo orbit average of Table S4: horizon and the initial iterates discarded) | exclusion:protocol_parameter (designG_dynamics_v4_values.json exact_skew: largest horizon 75000, B = 1000 discarded iterates) |

## Gates (full mode)

- claims inventory: OK
- essential success flags: OK
- completeness: OK
- fallback: OK
- qp_source (mu0,nu1 regime): OK
- manifest: OK
