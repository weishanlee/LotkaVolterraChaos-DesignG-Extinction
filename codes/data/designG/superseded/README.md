# Superseded producer JSONs

`designG_dynamics_v3_values.json` is retained as the historical record of what
`designG_dynamics_v3.py` computed, and is moved OUT of the validator's producer
directory deliberately.

Why it had to move rather than merely be ignored: `validate_designG_v1.py`
resolves a claim's `producer` stem to `<stem>_values.json` by DISCOVERY in the
data dir. While this file sat there and the claim registry still named
`designG_dynamics_v3`, the validator resolved `ext_frac_nu1_2` to `0.229444` --
v3's FIRST-PASSAGE count -- and reported PASS. A green validator against a
superseded artifact is worse than no validator: it certifies the retracted
number. Repointing the registry alone would have left that file able to satisfy
a future claim that named v3 again.

The live producer is `designG_dynamics_v4_values.json` (terminal persistence +
left-censoring). See the SUPERSEDED banner at the head of
`designG_dynamics_v3.py` in the folder `codes/` of this bundle.

2026-10-17: the eight absolute path strings of the authors' machine in
`designG_dynamics_v3_values.json` were replaced by the placeholder `<project>/`; no numeric value
was touched (verified leaf-by-leaf), so the archived numbers are unchanged.

2026-11-14: the provenance entry `round` was removed from `designG_dynamics_v3_values.json` (one line);
no other byte changed.

2026-11-14: the entries `spec_figure_nu1_2` and `spec_figure_note` were deleted from
`designG_dynamics_v3_values.json` (two lines), and the note of its check "overall max y (nu1<=4) <= 1" was
shortened (one line), as `designG_dynamics_v3.py` now writes it. The two entries and the words removed from
the note held figures that the manuscript does not print, not results of the run; every computed value and
every other byte are unchanged.
