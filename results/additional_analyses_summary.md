# Additional observational analyses — flu-HI antigenic causal-discovery study

Five analyses requested after the Round-3 revision, all run on existing data (no new
measurements). #8 was run directly; #1, #4, #11/#13, #12 were run as parallel sub-agent tracks.
Numbering is common mature H3 (VHID: mature = alignment column; Bedford: mature = column − 10).

---

## #8 — Temporal (train-past / test-future) cross-validation  [Bedford H3N2 only]

VHID carries no collection-year metadata (0/2751), so this is H3N2-only (1968–2010).
Expanding-window forward prediction: train on all pairs up to year *t*, test on the next 5-year block.

| train ≤ | test window | n_test | Ridge R² | XGBoost R² |
|--------:|-------------|-------:|---------:|-----------:|
| 1990 | 1991–1995 | 927 | −4.70 | 0.43 |
| 1995 | 1996–2000 | 369 | 0.35 | **0.60** |
| 2000 | 2001–2005 | 3636 | −3.73 | **−0.37** |
| 2005 | 2006–2010 | 1297 | −0.59 | 0.25 |

**Finding.** Forward-in-time generalization is unstable and, in the worst window, negative:
XGBoost swings from R²=0.60 to **−0.37**, versus 0.613 under random-split CV. This is the
sharpest statement yet of the transportability limit the removed OOD test was meant to probe —
a model fit on past seasons does not reliably predict titers for future antigenic clusters,
consistent with drift crossing cluster boundaries the training data never saw.

## #1 — Cross-dataset structure replication with a permutation null

Both H3N2 panels run through identical discovery; sets mapped to mature numbering; cross-dataset
Jaccard compared against a 20,000-draw null (random same-size sets from each panel's node set).

| set | observed J | shared positions | null mean | p |
|-----|-----------:|-------------------|----------:|--:|
| PC | 0.167 | **158, 189** | 0.035 | 0.062 |
| GES | 0.182 | 133, 276 | 0.035 | 0.051 |
| Bootstrap (freq≥0.5) | 0.182 | **158, 189** | 0.034 | **0.047** |

**Finding.** Overlap is modest but ~5× the chance level for every set; **mature 158 and 189
replicate across both datasets and across constraint- and score-based discovery**, with the
bootstrap-set replication significant at p=0.047. This converts "the drivers coincide" into a
replication statistic with a null.

## #4 — Selection stability vs. adjusted effect (driver / hitchhiker separation)

Wilson 95% CI on each bootstrap selection frequency (B=200) paired with the shipped bootstrap CI
on the adjusted effect. (True BCa would need the stored per-resample draws, which are not shipped;
this is a binomial-proportion CI, stated as such.)

**Finding.** The quadrant plot cleanly separates defensible drivers from hitchhikers.
VHID drivers (stable AND large-effect): **158, 189, 289**; 156 is a stable *hitchhiker*
(freq=1.0 but below-median |β|); 144 is genuinely unstable (freq=0.505, Wilson 0.44–0.57).
Bedford drivers: **11, 133, 189**; 157/158/193 are stable-but-small; 190 (0.65) and 278 (0.82)
fall below the 0.9 stability line despite non-trivial effects.

## #12 — Formal test of the second-order KAN's epistatic pairs

For each KAN-nominated pair, an OLS interaction term x_a·x_b (Grantham) controlling for the causal
parents, HC3 robust SE, BH-corrected within dataset; plus a distance-correlation nonparametric check.

**Finding.** **10 of 16 nominated pairs** show a significant linear interaction (q<0.05):
4/8 in VHID, 6/8 in Bedford. Coefficients are small (|β| ~1e-4–1e-3 log2-titer per Grantham²)
but precisely estimated at these n. Critically, **KAN interaction norm does not track statistical
significance** — the norm ranks model curvature, not the presence of a formally testable epistatic
term — so the KAN surfaces should be read as an interpretability aid, not an inference. The
distance-correlation test flags nearly every pair, reflecting its power against any nonlinearity
rather than a specifically multiplicative effect.

## #11 + #13 — Encoding comparison & per-property decomposition  [exploratory]

PC titer-adjacency discovery (fisherz, α=0.01) under binary / Grantham / 12-property-L2 encodings.
(Full Bedford PC under all three encodings was computationally intractable — the dense skeleton did
not complete — so the encoding comparison is reported for VHID; the caveat is documented.)

**Finding (encoding).** Encoding choice materially changes selection. In VHID only **144, 156,
189, 289 survive under all three encodings** (pairwise Jaccard 0.50–0.63); five positions are
encoding-fragile — most notably **mature 158, selected only under the binary flag and dropping out
under both continuous encodings**, confirming the study's flagged fragility of 158.

**Finding (per-property) — strictly exploratory.** A single substitution sets all 12 property axes
at once, so per-property causal claims are **non-identifiable**: at every driver the property-change
vectors span only rank 1–4 of the active axes (289 is rank 1 — one substitution type), and partial
correlations conditioning one axis on the others collapse to ~0. The marginals only *rank which axis
co-varies most* (e.g. H-bond-acceptor/β-sheet at 158, β-sheet at 189). Reported as "which axis
co-varies," never as a property-level cause.

---

## What these add to the manuscript

- **158 and 189 are now the two most defensible drivers** — they replicate across datasets (p≈0.05),
  and 189 additionally sits in the driver quadrant of #4. 158's encoding-fragility (#11/#13) is a
  caveat to state alongside its replication.
- **The transportability limit is now quantified** (#8): honest forward-prediction R² is unstable,
  negative in the hardest era — the strongest available answer to the removed OOD test.
- **The KAN epistasis story is disciplined** (#12): real interactions exist (10/16 pairs), but the
  KAN norm is not an inference statistic — exactly the "capacity vs. signal" caveat the reviewers raised.
