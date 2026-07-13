# Causal DAG validation report

_Step 4 of the flu-HI antigenic workflow. Tests 1-3 establish *not inconsistent* + confidence; only Test 4 tests causality._

## 1. Goodness-of-fit (Shipley d-sep basis-set test)
- Fisher C = 2000.328 on df = 56, overall p = 0 -> DAG **REJECTED**
- basis claims: 28; individual violations (p<0.05): 28
- _Tests the equivalence class, not arrow directions. Passing = not falsified, not proven._

## 2. Bootstrap structural stability
- B = 200 resamples
- top group inclusion frequencies:
  - `G168 -> titer` : 1.00
  - `G288 -> titer` : 1.00
  - `G203 -> titer` : 1.00
  - `G143 -> titer` : 1.00
  - `G199 -> titer` : 1.00
  - `G11 -> titer` : 1.00
  - `G200 -> titer` : 1.00
  - `G167 -> titer` : 0.96
- _Read stability at the linkage-GROUP level; raw linked columns dilute each other._

## 3. Direct-effect estimates (given fixed structure)
- estimation n = 7808, data-split inference: False

| term | coef | 95% CI | p |
|---|---|---|---|
| (intercept) | 8.781 | [8.724, 8.838] | 0 |
| pos11 | -0.009084 | [-0.01073, -0.007435] | 0 |
| pos143 | -0.01116 | [-0.01372, -0.008602] | 0 |
| pos167 | -0.002062 | [-0.003023, -0.0011] | 2.66e-05 |
| pos168 | -0.006831 | [-0.008397, -0.005264] | 0 |
| pos199 | -0.01422 | [-0.01565, -0.01279] | 0 |
| pos200 | -0.006241 | [-0.007905, -0.004577] | 2.16e-13 |
| pos203 | -0.002195 | [-0.002885, -0.001505] | 4.64e-10 |
| pos288 | -0.009597 | [-0.01086, -0.008334] | 0 |

- _WARNING: no data split -> CIs are post-selection anti-conservative (too narrow)._

## 4. Out-of-distribution / interventional validation
_not run (the only test of causality as such)_
