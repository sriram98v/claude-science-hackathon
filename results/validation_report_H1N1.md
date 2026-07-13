# Causal DAG validation report

_Step 4 of the flu-HI antigenic workflow. Tests 1-3 establish *not inconsistent* + confidence; only Test 4 tests causality._

## 1. Goodness-of-fit (Shipley d-sep basis-set test)
- Fisher C = 230.308 on df = 12, overall p = 0 -> DAG **REJECTED**
- basis claims: 6; individual violations (p<0.05): 5
- _Tests the equivalence class, not arrow directions. Passing = not falsified, not proven._

## 2. Bootstrap structural stability
- B = 200 resamples
- top group inclusion frequencies:
  - `G199 -> titer` : 1.00
  - `G53 -> titer` : 1.00
  - `G151 -> titer` : 1.00
  - `G135 -> titer` : 0.97
- _Read stability at the linkage-GROUP level; raw linked columns dilute each other._

## 3. Direct-effect estimates (given fixed structure)
- estimation n = 913, data-split inference: False

| term | coef | 95% CI | p |
|---|---|---|---|
| (intercept) | 8.526 | [8.396, 8.656] | 0 |
| pos53 | -0.01396 | [-0.01714, -0.01079] | 0 |
| pos135 | -0.008341 | [-0.01162, -0.005061] | 7.2e-07 |
| pos151 | -0.02723 | [-0.03268, -0.02178] | 0 |
| pos199 | -0.009954 | [-0.0127, -0.007209] | 2.25e-12 |

- _WARNING: no data split -> CIs are post-selection anti-conservative (too narrow)._

## 4. Out-of-distribution / interventional validation
- reference parents: pos135, pos151, pos199, pos53

| held-out group | n | RMSE | R2 | parent recovery (Jaccard) |
|---|---|---|---|---|
| 1975 | 184 | 1.95 | -0.107 | 1.00 |
| 1995 | 157 | 1.53 | 0.445 | 1.00 |
| 1980 | 101 | 1.77 | -0.161 | 1.00 |
| 2000 | 53 | 1.62 | 0.134 | 1.00 |
| 1985 | 44 | 1.98 | 0.0913 | 0.75 |
| 1990 | 9 | 2 | 0.372 | 1.00 |
| 2005 | 365 | 1.25 | 0.247 | 0.50 |
- _OOD success is necessary, not sufficient; a failure may also mean a real cross-group mechanism shift. Gold standard is a mutagenesis intervention at a predicted direct cause._
