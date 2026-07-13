# Causal DAG validation report

_Step 4 of the flu-HI antigenic workflow. Tests 1-3 establish *not inconsistent* + confidence; only Test 4 tests causality._

## 1. Goodness-of-fit (Shipley d-sep basis-set test)
- Fisher C = 418.649 on df = 20, overall p = 0 -> DAG **REJECTED**
- basis claims: 10; individual violations (p<0.05): 10
- _Tests the equivalence class, not arrow directions. Passing = not falsified, not proven._

## 2. Bootstrap structural stability
- B = 200 resamples
- top group inclusion frequencies:
  - `G144 -> titer` : 1.00
  - `G189 -> titer` : 1.00
  - `G156 -> titer` : 1.00
  - `G158 -> titer` : 1.00
  - `G289 -> titer` : 0.99
- _Read stability at the linkage-GROUP level; raw linked columns dilute each other._

## 3. Direct-effect estimates (given fixed structure)
- estimation n = 2751, data-split inference: False

| term | coef | 95% CI | p |
|---|---|---|---|
| (intercept) | 9.336 | [9.241, 9.431] | 0 |
| pos144 | -0.01131 | [-0.01279, -0.009829] | 0 |
| pos156 | -0.01338 | [-0.01692, -0.009841] | 1.65e-13 |
| pos158 | -0.01661 | [-0.01911, -0.01411] | 0 |
| pos189 | -0.034 | [-0.03594, -0.03205] | 0 |
| pos289 | 0.02388 | [0.01562, 0.03215] | 1.59e-08 |

- _WARNING: no data split -> CIs are post-selection anti-conservative (too narrow)._

## 4. Out-of-distribution / interventional validation
_not run (the only test of causality as such)_
