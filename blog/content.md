This project was built as part of the [Claude Science Hackathon](https://cerebralvalley.ai/e/built-with-claude-life-sciences).

*A reproducible research notebook.* Run top-to-bottom (`Kernel → Restart & Run All`) to regenerate every table and
figure from the raw data repository shipped alongside this notebook.

## Abstract

Identifying specific hemagglutinin (HA) mutations that causally alter antibody recognition remains a significant challenge because dense viral phylogenies tightly link functional escape drivers with passenger mutations. While sequence-based models predict antigenic distance with high accuracy, they often conflate the evolutionary dynamics of population-level sweeps with the mechanistic physics of antibody-binding disruption. To address this ambiguity, we define the target estimand as a provenance-independent, type-level interventional contrast at the antigen–antibody interface, characterizing linkage-driven resolution loss as an intrinsic structural feature of observational serology. We analyze two H3N2 [hemagglutination-inhibition](https://en.wikipedia.org/wiki/Hemagglutination_assay) (HI) datasets by collapsing co-evolving positions and applying target-oriented [causal discovery](https://en.wikipedia.org/wiki/Causal_inference) algorithms ([PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html), [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html), [FCI](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/FCI.html)) prioritized by a 200-resample [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) stability framework. This is paired with an interpretable [B-spline](https://en.wikipedia.org/wiki/B-spline) [Kolmogorov–Arnold Network](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Arnold_representation_theorem) (KAN) to evaluate per-position response curves and capture second-order position-by-position [epistasis](https://en.wikipedia.org/wiki/Epistasis). Under matched-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29), the first-order KAN performs slightly below black-box [gradient boosting](https://en.wikipedia.org/wiki/Gradient_boosting) (lower $R^2$ by $\approx 0.025\text{--}0.028$; [Wilcoxon](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) $p < 10^{-5}$); a second-order KAN mitigates most of this performance gap under a 5-fold validation protocol. Cross-method convergence evaluates agreement across causal and association-based frameworks, identifying candidate drivers that overlap classical antigenic sites A and B. In the VHID dataset, convergent signals localize to mature H3 positions 156 and 189, where position 156 exhibits characteristics consistent with a stable hitchhiker (high frequency, small non-robust effect) and position 289 emerges as a doubly-robust candidate. In the Bedford dataset, convergent positions include mature 133, 158, and 189, with 158 demonstrating sensitivity to feature encoding. Under cluster resampling, mature position 189 remains the most robust signal across both datasets. Backdoor-adjusted effect sizes systematically shrink relative to marginal associations, consistent with the mitigation of phylogenetic confounding, though a Shipley [d-separation](https://en.wikipedia.org/wiki/Bayesian_network) test rejects a simplified sink-star structure. Finally, rigorous grouped (leave-serum-out) [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) establishes a baseline generalization accuracy (median $R^2 \approx 0.615$ for VHID, $\approx 0.498$ for Bedford). This work establishes a stability-ranked, confounding-audited, and linkage-aware feature-selection framework that systematically isolates candidate biophysical drivers of antigenic drift from observational data.

## Introduction

Influenza viruses continuously evolve under selective pressure from population immunity, accumulating substitutions in surface glycoproteins that facilitate immune escape. This process of antigenic drift causes circulating strains to diverge from those recognized by prior immunity, reducing vaccine effectiveness when selected vaccine strains mismatch circulating variants. This divergence is quantified as antigenic distance, and measuring it accurately is essential for evaluating vaccine efficacy and optimizing antigen selection. Traditionally, antigenic distance is determined via the [hemagglutination-inhibition](https://en.wikipedia.org/wiki/Hemagglutination_assay) (HI) assay, where lower cross-titers indicate greater immune escape. While frameworks like [antigenic cartography](https://en.wikipedia.org/wiki/Antigenic_cartography) have been foundational in mapping these titers, generating the required serum panels is logistically demanding, costly, and prone to inter-laboratory variation. This bottleneck has motivated sequence-based predictive models designed to forecast HI titers directly from hemagglutinin (HA) sequences. 

However, while machine learning approaches achieve high predictive accuracy, they frequently operate as black boxes, identifying predictive correlates rather than isolating the underlying causal drivers of immune escape. Because influenza strains share a dense phylogenetic history, HA positions exhibit strong [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium). Consequently, passenger mutations riding along with functional escape drivers appear as predictive as the drivers themselves. To resolve this ambiguity, this study introduces a precise conceptual reframing. Disentangling the drivers of antigenic drift requires separating the evolutionary question (which substitutions were favored by natural selection and swept the population) from the mechanistic question (which substitutions physically disrupt antibody recognition when introduced into a given strain background). This study focuses explicitly on the second, mechanistic question. We define our target estimand not as a historical claim about viral evolution, but as a provenance-independent, type-level interventional contrast at the antigen-antibody interface. Ideally, this contrast reflects a controlled biophysical experiment: introducing a single residue change at position $p$ in reference virus B to match virus A, while holding the rest of the protein sequence fixed, and measuring the resulting change in HI titer. The magnitude of this effect depends strictly on the structural footprint and local chemistry, independent of whether the mutation arose via positive selection or neutral drift. While evolutionary provenance does not dictate the biophysical effect itself, it heavily constrains our capacity to identify it from observational data. The selective history of the virus introduces systematic phylogenetic confounding, clustering distinct mutations into tightly linked blocks. Acknowledging this architecture allows us to treat linkage-driven resolution loss as an inherent structural characteristic of observational HI data, which must be formally accommodated within the analytical pipeline.

**Related works.** Antigenic cartography revealed the punctuated cluster structure of H3N2 drift by embedding HI tables into low-dimensional maps, and subsequent models unified antigenic and genetic evolution within joint phylogenetic frameworks. Modern sequence-based predictors forecast cross-immunity between drifted strains from sequence data with high accuracy. While these approaches are primarily predictive or descriptive, they do not explicitly learn causal structure over individual HA positions with quantified stability. The antigenic sites encompassing these positions were originally defined structurally and serologically (sites A through E on the H3 head) and subsequently refined via substitution resolution and deep mutational scanning escape maps. Our approach evaluates the extent to which these positions can be recovered directly from HI titers without structural priors. This framework complements influenza fitness models by isolating per-position candidate drivers and builds on intelligible-model literature for pairwise interactions.

**Our contributions.** We combine several ingredients that are not usually applied together for HI data to build a conservative, interpretable, and self-audited feature-selection pipeline. First, we implement linkage collapse with target-oriented [causal discovery](https://en.wikipedia.org/wiki/Causal_inference): we merge near-deterministic co-evolving positions into representative loci, then learn the direct-cause neighborhood of the HI target using [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html), [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html), and [FCI](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/FCI.html) algorithms, ranking every candidate by 200-resample [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) stability. Second, we deploy a genuine [B-spline](https://en.wikipedia.org/wiki/B-spline) [Kolmogorov–Arnold Network](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Arnold_representation_theorem) (KAN), an interpretable non-linear predictor whose learned per-position response curves are directly inspectable, and extend it to second order to capture and visualize position-by-position [epistasis](https://en.wikipedia.org/wiki/Epistasis). Third, we assess cross-method convergence, treating it honestly as agreement between one causal screen and three *correlated* association screens (the KAN, [gradient boosting](https://en.wikipedia.org/wiki/Gradient_boosting), and univariate association) rather than four independent lines of evidence; positions flagged across these screens form our strongest candidate-driver claims. We then subject those claims to a battery of audits: a [permutation](https://en.wikipedia.org/wiki/Permutation_test) calibration showing that the [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) conditional-independence test holds near-nominal size on our binary, [left-censored](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) data; a cluster (by-virus and by-serum) [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) that separates a robust convergent core from an over-optimistic stability tier; a [left-censoring](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) sensitivity analysis of the adjusted effect sizes; and a token-level identifiability audit of the discovered structure. Throughout, we estimate backdoor-adjusted per-position effect sizes but report them as partial-regression coefficients, because the baseline adjustment assumptions are rejected in-sample.

The remainder of the notebook is organized as an executable paper. [Section 2 (Methods)](#Methods) introduces the two H3N2 HI datasets and explains how their feature matrices are derived from the raw data. [Section 3 (Results)](#Results) then carries the analysis in full, each subsection stating how a step is performed, presenting its output, and interpreting it: the predictive benchmark and its leakage-free [cross-validated](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) comparison (§ [3.1](#Predictive-Benchmark)–[3.2](#Cross-Validated-$R^2$-and-Rigorous-Error-Bounds)); target-oriented [causal discovery](https://en.wikipedia.org/wiki/Causal_inference) together with the pre-collapse linkage block sizes, a [permutation](https://en.wikipedia.org/wiki/Permutation_test) calibration of the [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) independence test, and a cluster [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) of candidate stability (§ [3.3](#Causal-discovery)), followed by the discovered dependency structure (§ [3.4](#Discovered-Dependency-Structure)); the interpretable [B-spline](https://en.wikipedia.org/wiki/B-spline) KAN and its second-order [epistasis](https://en.wikipedia.org/wiki/Epistasis) extension (§ [3.5](#Interpretable-B-Spline-KAN)–[3.6](#Capturing-Epistasis:-A-Second-Order-KAN)); cross-method convergence (§ [3.7](#Cross-Method-Convergence)); backdoor-adjusted effect sizes and their sensitivity to titer [left-censoring](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) (§ [3.8](#Backdoor-Adjusted-Effect-Sizes)); three DAG-validation tests (§ [3.9](#Validating-the-Causal-Structure)); and a continuous per-position encoding that re-examines the titer [Markov blanket](https://en.wikipedia.org/wiki/Markov_blanket) (§ [3.10](#Continuous-Per-Position-Encoding-and-the-Titer-Markov-Blanket)). The narrative moves from how much of the HI signal is learnable and how non-linear it is, through the discovered candidate drivers and their audited dependency structure, to per-position effect sizes and the validation — and in-sample rejection — of the discovered graph. [Section 4 (Conclusion)](#Conclusion) interprets the convergent positions biologically, states plainly where residue-level attribution is limited and what assumptions the causal framing rests on, and looks ahead. [Section 5](#References) lists references.

## Methods
This section outlines the data structure and feature extraction protocols. Analytical procedures are described alongside their respective outputs in the Results section to maintain context. Configuration parameters (random seeds, linkage thresholds, [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) counts, and tier cutoffs) are centralized in src/analysis.py. Computationally intensive steps, such as the causal [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) and repeated $k$-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29), are managed via environment flags and are defaulted to load precomputed results to ensure reproducibility.

### Datasets
The study evaluates two H3N2 virus $\times$ reference-strain HI panels. Each dataset is derived from previously published works: the VHID panel is derived from the DPCIPI dataset (Du et al., 2023), and the Bedford dataset is obtained from Bedford et al. (2014). While both panels represent H3N2 HI datasets, they differ in metadata completeness. The Bedford H3N2 panel is curated from Bedford et al. (2014), GenBank accessions, and isolate collection years, all of which are fully populated. We analyze the two panels as independent replications of the same underlying structure-learning task rather than assuming a calibrated joint assay protocol.
The reported positions correspond to mature H3 residue numbers. Feature matrices are indexed by HA1 alignment columns, which map to mature numbers via fixed per-panel offsets. The VHID reference is gapless from mature Q1, yielding a 0-residue offset. In contrast, the Bedford H3N2 alignment contains a 9-residue signal-peptide prefix and an internal gap at column 17; thus, mature Q1 maps to column 10, and for all reported positions ($\ge \text{column 143}$), $\text{mature residue} = \text{column} - 10$. This mapping was verified against column-wise consensus references and gapless VHID references to ensure structural interpretability.

#### Preprocessing
The sequence data is formatted as per-position HA1 feature matrices under two distinct encodings:
Binary Mismatch Encoding: Sets a position to 1 when the virus and reference residues differ at that HA1 alignment position and 0 otherwise; this serves as the primary input for [causal discovery](https://en.wikipedia.org/wiki/Causal_inference).
[Grantham](https://doi.org/10.1126/science.185.4154.862) Encoding: Quantifies the physicochemical distance between residue pairs based on [Grantham](https://doi.org/10.1126/science.185.4154.862) (1974) metrics (gap-aware); this is utilized for predictive modeling and KAN response curves.

The target variable is modeled as $\log_2(\text{HI\_titer})$. To ensure complete pipeline traceability, feature matrices are regenerated from raw virus $\times$ reference pair tables using the repository's build scripts (scripts/build_vhid_matrices.py, scripts/build_bedford_matrices.py). Matrix identities were verified against shipped derivatives using [SHA-256](https://en.wikipedia.org/wiki/SHA-2) checksums. Feature spaces are evaluated independently per panel to prevent artifacts from lineage-specific HA1 trimming.

#### Reproducing the feature matrices from raw data

So that the study starts from raw data rather than shipped derivatives, we regenerate
the feature matrices from the cleaned virus × reference pair tables using the
repository's own build scripts (`scripts/build_vhid_matrices.py`,
`scripts/build_bedford_matrices.py`), which depend only on numpy and pandas. Running
them here makes the matrices provably the shipped ones (verified by [SHA-256](https://en.wikipedia.org/wiki/SHA-2)).

{{ table(src="tables/table-01.html", n=1) }}

Each dataset is internally row-aligned across its binary matrix, [Grantham](https://doi.org/10.1126/science.185.4154.862) matrix, and
cleaned pair table. Feature spaces are comparable *within* a lineage but not across
lineages (each lineage was HA1-trimmed against its own reference), so we analyze the
two datasets independently and compare only which HA positions emerge.

#### The target distribution

Before any modeling we look at what is being predicted. The panel below shows the
distribution of log2 HI titer in each dataset; its spread sets the scale against which
every R² reported in the Results should be read.

{{ figure(src="figures/target_distribution.png", n=2, alt="target distribution") }}

## Results
After establishing baseline performance metrics in [Section 3.1](#Predictive-Benchmark), we describe the generalization capacity and error bounds of our sequence-to-antigenic maps under strict grouped and temporal [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) protocols in [Section 3.2](#Cross-Validated-$R^2$-and-Rigorous-Error-Bounds). In [Section 3.3](#Causal-discovery), we outline the structural [causal discovery](https://en.wikipedia.org/wiki/Causal_inference) pipeline, linkage-collapse dynamics, and test calibrations. In [Section 3.4](#Discovered-Dependency-Structure), we analyze the resulting parent dependency structures and check for intermediate mediation. We then present our interpretable modeling frameworks, detailing the 1-D response curves of the first-order [B-spline](https://en.wikipedia.org/wiki/B-spline) KAN in [Section 3.5](#Interpretable-B-Spline-KAN) and the bivariate tensor-product surfaces for capturing [epistasis](https://en.wikipedia.org/wiki/Epistasis) in [Section 3.6](#Capturing-Epistasis:-A-Second-Order-KAN). In [Section 3.7](#Cross-Method-Convergence), we evaluate the cross-method convergence of our feature screens and verify non-linear omissions. Finally, we describe the estimation of backdoor-adjusted effect sizes and driver-hitchhiker differentiation in [Section 3.8](#Backdoor-Adjusted-Effect-Sizes), the global [d-separation](https://en.wikipedia.org/wiki/Bayesian_network) validation tests in [Section 3.9](#Validating-the-Causal-Structure), and the continuous physicochemical encoding replication in [Section 3.10](#Continuous-Per-Position-Encoding-and-the-Titer-Markov-Blanket).

### Predictive Benchmark

To establish performance baselines and characterize the mathematical properties of the antigenic signal, we evaluated the held-out test $R^2$ (20% split) for the top-performing single-position (max univariate $R^2$), [LASSO](https://en.wikipedia.org/wiki/Lasso_%28statistics%29), [Ridge](https://en.wikipedia.org/wiki/Ridge_regression), and [XGBoost](https://en.wikipedia.org/wiki/XGBoost) models. The predictive performance across both datasets consistently follows the ordering:

$$
\begin{aligned}
&\text{[XGBoost](https://en.wikipedia.org/wiki/XGBoost)} \gtrsim \text{[LASSO](https://en.wikipedia.org/wiki/Lasso_%28statistics%29)} \approx \text{[Ridge](https://en.wikipedia.org/wiki/Ridge_regression)} \\\\
&\qquad \gg \text{Best Single Position}.
\end{aligned}
$$

{{ table(src="tables/table-03.html", n=3) }}

The substantial performance delta between the single-position baseline and multivariable models indicates that the antigenic signal is distributed across multiple positions. Furthermore, the performance margin achieved by [XGBoost](https://en.wikipedia.org/wiki/XGBoost) over linear models suggests the presence of underlying non-linear and interaction structures, motivating the deployment of the KAN framework detailed in [Section 3.5](#Interpretable-B-Spline-KAN).

### Cross-Validated $R^2$ and Rigorous Error Bounds
To provide robust uncertainty estimates and ensure rigorous model comparison, all models were evaluated under an identical $5 \times 4$ repeated $k$-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) protocol using matched folds. [XGBoost](https://en.wikipedia.org/wiki/XGBoost) parameters were optimized via [early stopping](https://en.wikipedia.org/wiki/Early_stopping), performed strictly on an inner validation split carved from the training fold, thereby protecting the test fold from data leakage. Models are compared using a matched-fold paired test ([Wilcoxon](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) signed-rank test on the 20 per-fold differences) rather than evaluating confidence-interval overlaps.

While a paired test on random splits effectively differentiates model architectures, random partitioning allows identical viruses and reference antisera to recur across training and testing folds. This pair-level recurrence can artificially inflate performance metrics because models can memorize strain-specific profiles rather than generalizing to unseen variants. We treat random-split metrics purely as a baseline and leverage leakage-free grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) as our primary generalization metric.

{{ figure(src="figures/cv_r2.png", n=4, alt="cv r2") }}

Under the random-split baseline, [XGBoost](https://en.wikipedia.org/wiki/XGBoost) yields the highest [cross-validated](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) $R^2$ across both datasets (VHID: 0.845, Bedford: 0.613), followed closely by the KAN (VHID: 0.820, Bedford: 0.585). The matched-fold paired test confirms that the performance gap is statistically robust: the KAN trails [XGBoost](https://en.wikipedia.org/wiki/XGBoost) by 0.025 $R^2$ on VHID ([Wilcoxon](https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) $p \approx 2 \times 10^{-6}$) and by 0.028 on Bedford ($p \approx 6 \times 10^{-6}$). The KAN's primary utility lies in its additive interpretability, recovering most of the tree-based model's performance while explicitly exposing per-position response curves.

To determine true generalization performance on unseen strains, we executed grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) via leave-virus-out and leave-serum-out protocols.

{{ table(src="tables/table-06.html", n=5) }}

The performance metrics degrade under grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29), confirming that random-split protocols are systematically influenced by strain/serum recurrence. Under the strict leave-serum-out protocol, [XGBoost](https://en.wikipedia.org/wiki/XGBoost) performance settles at a median $R^2$ of 0.615 on VHID and 0.498 on Bedford. These grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) medians represent our honest predictive headlines for sequence-to-antigenic maps operating outside the training distribution.

#### Temporal Transportability
Since the Bedford H3N2 dataset includes temporal metadata (1968–2010), we evaluated forward-in-time generalization using an expanding window strategy: training on all pairs up to year $t$ and testing on the subsequent 5-year block.

| train ≤ | test window | n_test | Ridge R² | XGBoost R² |
|--------:|-------------|-------:|---------:|-----------:|
| 1990 | 1991–1995 | 927 | −4.70 | 0.43 |
| 1995 | 1996–2000 | 369 | 0.35 | **0.60** |
| 2000 | 2001–2005 | 3636 | −3.73 | **−0.37** |
| 2005 | 2006–2010 | 1297 | −0.59 | 0.25 |

Forward-in-time generalization displays notable instability; in the 2001–2005 test window, [XGBoost](https://en.wikipedia.org/wiki/XGBoost) performance drops to $R^2 = -0.37$. This highlights a primary boundary of transportability: models trained exclusively on past seasons struggle to reliably predict titers for future antigenic clusters when drift crosses major structural boundaries that are absent from the training history.

{{ figure(src="figures/figure-03.png", n=6, alt="figure 03") }}

#### Cross-Cluster Transportability of Property Encodings
We tested whether encoding substitutions by their physicochemical property shifts (a 12-property $L_2$ scalar) rather than by raw amino acid identity enhances temporal transportability. Mapping unseen substitutions to their local shifts in charge, volume, or hydrophobicity could enable the model to generalize based on biophysical similarity.

Our empirical results do not support this hypothesis. Across the unbiased [XGBoost](https://en.wikipedia.org/wiki/XGBoost) models, the mean future $R^2$ across all test windows was 0.19 for the binary encoding, 0.22 for the [Grantham](https://doi.org/10.1126/science.185.4154.862) distance, and 0.18 for the 12-property $L_2$ vector. The single-scalar [Grantham](https://doi.org/10.1126/science.185.4154.862) distance yielded the most stable performance across windows, undermining the assumption that higher-dimensional property vectors improve generalization to distribution shifts. All three encodings systematically fail during the 2001–2005 window, confirming that substitution-based representations do not fully capture major shifts in antigenic distribution. 

Consequently, the utility of property encodings rests on their interpretability ([Section 3.10](#Continuous-Per-Position-Encoding-and-the-Titer-Markov-Blanket)) rather than cross-cluster predictive transport.

{{ figure(src="figures/figure-04.png", n=7, alt="figure 04") }}

### Causal discovery
We model the HI titer as a downstream causal sink: HA sequence variations cause variations in titer, orienting all feature-to-target edges into the target variable. The direct-cause candidates are defined as the immediate parents of this target node.

The structural pipeline proceeds as follows 

$$
\begin{aligned}
\text{Power Filter} &\rightarrow \text{Linkage Collapse } (\vert{}\phi\vert{} \ge 0.8) \\\\
&\rightarrow \text{Constraint/Score Discovery } (\text{[PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html), [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html), [FCI](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/FCI.html)}) \\\\
&\rightarrow \text{Bootstrap Stability Evaluation } (B=200).
\end{aligned}
$$

Linkage collapse is a critical prerequisite; near-deterministic co-evolution violates the faithfulness assumption and introduces structural singularities into constraint-based searches. Collapsing these blocks into single representative loci resolves these dependencies. Because causal claims apply to the entire co-evolving unit, we explicitly report block sizes throughout.

{{ table(src="tables/table-09.html", n=8) }}

{{ table(src="tables/table-10.html", n=9) }}

{{ table(src="tables/table-102.html", n=10) }}

Following linkage collapse, residual strong locus pairs ($\vert{}\phi\vert{} \ge 0.9$) drop to zero in both datasets. [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) and [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html) algorithms were executed on the collapsed feature space, with selection frequency across 200 [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) resamples used to categorize candidates into stability tiers: High Confidence ($\ge 0.9$) and Moderate Confidence ($0.5\text{--}0.9$).
Under standard i.i.d. resampling, the High Confidence parent sets encompass:
- VHID: {156, 158, 189, 289}
- Bedford (Mature): {2, 133, 157, 158, 189, 193}
Positions with a block_size > 1 carry structural claims for the entire linkage group rather than the isolated index alone (e.g., Bedford mature position 2 represents a 9-position co-evolving block; see [Section 3.3.1](#Pre-Collapse-Linkage-Block-Dynamics)).
We apply two critical caveats to these stability classifications:
1. **Cluster Dependencies:** The i.i.d. [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) represents an upper bound on stability. Under a rigorous cluster [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) that resamples at the level of whole viruses or whole sera to preserve data dependencies, the VHID high-confidence set contracts: only positions {156, 189} remain robust under virus-level resampling, and {189, 289} under serum-level resampling. Mature position 189 consistently retains its High Confidence classification across all three resampling strategies.
2. **Algorithmic Concordance:** In the Bedford dataset, [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html) independently corroborates the [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) parents (mature 133, 157, 158). However, in the VHID dataset, the two algorithms diverge significantly: [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) selects {144, 156, 158, 189, 289} while [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html) selects {50, 133, 144, 145, 262, 276}, intersecting exclusively at position 144. This divergence indicates structural sensitivities to modeling assumptions within the VHID panel. The [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) framework ranks candidates by selection stability rather than arbitrating algorithmic divergence.

#### Pre-Collapse Linkage Block Dynamics
Linkage collapse groups positions co-evolving at $\vert{}\phi\vert{} \ge 0.8$ into representative units. The distribution of these blocks is heavy-tailed: while the majority of positions remain singletons (61/71 in VHID; 114/123 in Bedford), a few large blocks absorb substantial portions of the feature space.

In the Bedford H3N2 dataset, the largest pre-collapse block spans 88 positions (representative alignment column 181), the second spans 55 positions (column 90), and the third spans 26 positions (column 50). These large, clade-linked units absorb numerous head positions, indicating that purely observational methods cannot structurally disentangle individual residue effects within these blocks. Conversely, the VHID dataset exhibits less linkage; its largest block encompasses only 14 positions (column 173), allowing its collapsed loci to map more directly to individual amino acid changes.

{{ figure(src="figures/block_size.png", n=11, alt="block size") }}

{{ table(src="tables/table-11.html", n=12) }}

{{ table(src="tables/table-112.html", n=13) }}

#### Calibration of the Fisher-Z Test Under Permutation Nulls
Conditional independence decisions within the causal pipeline rely on the [Pearson](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) [partial-correlation](https://en.wikipedia.org/wiki/Partial_correlation) [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) test (fisherz). Because the data consists of binary mismatch features and [left-censored](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) titers, the multivariate normality assumption is violated, making the operating threshold ($\alpha = 0.01$) nominal.

To evaluate true type-I error rates, we constructed an empirical null distribution by permuting the $\log_2$ titer column, thereby disrupting feature-to-target relationships while preserving feature-to-feature correlation structures. We evaluated 10,000 independent tests across conditioning sizes $\{0, 1, 2, 3\}$ using three distinct encodings: VHID collapsed binary loci, VHID continuous $L_2$ physicochemical distances, and a top-40 screened subset of Bedford binary loci.

At the target operating point ($\alpha = 0.01$), the empirical false-positive rates align closely with nominal expectations, landing within their respective 95% [Clopper–Pearson](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval) intervals:
- VHID Binary Loci: 0.0103 [0.0084, 0.0125]
- VHID $L_2$ Continuous Loci: 0.0102 [0.0083, 0.0124]
- Bedford Binary Loci: 0.0094 [0.0076, 0.0115]

At a looser threshold ($\alpha = 0.05$), marginal deviations occur (VHID binary: 0.054, Bedford binary: 0.045, continuous: 0.051). These results indicate that the [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) test maintains controlled size at the target operating threshold ($\alpha = 0.01$), preventing inflation of the false-positive edge rate on these data.

{{ figure(src="figures/fisherz_calibration.png", n=14, alt="fisherz calibration") }}

#### Cluster (block) bootstrap: is selection stability an artifact of i.i.d. resampling?

The 200× selection-stability [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) above resamples HI **pairs** independently. But VHID's 2751 pairs derive from only 246 viruses crossed with 45 reference sera, and the grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) in §[3.2](#Cross-Validated-$R^2$-and-Rigorous-Error-Bounds) shows this clustering is decisive (held-out R² drops ~0.23 when folds respect virus grouping). Resampling pairs i.i.d. treats correlated pairs as independent draws and can therefore *overstate* how reproducibly a position is selected.

To test this directly we re-ran the identical collapse → [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) parent-selection routine (`src/causal_helpers.py`), changing only the resampling **unit**: instead of drawing pairs, we draw whole **viruses** with replacement (leave-virus-out clusters), and separately whole **reference sera**. Everything else — linkage collapse, the screened 50-locus node set, α=0.01, terminal-target background knowledge — is held fixed, so the comparison isolates the effect of respecting clustering. All three schemes use B=200 on VHID.

**Finding.** The HIGH-stability set is **not** preserved under clustering. Under i.i.d. pairs it is {156, 158, 189, 289}; under **virus** clustering only {156, 189} remain HIGH (158 and 289 fall to moderate), and under **serum** clustering only {189, 289} remain HIGH (156 and 158 fall to moderate). Only **mature 189** stays HIGH in all three schemes. Frequencies of the i.i.d.-HIGH parents move systematically **downward** toward the moderate range (mean change −0.04 under virus resampling, −0.12 under serum resampling; pos 158 falls from 0.95 to 0.69 under serum clustering) — i.e. the pipeline is *less* confident once pair correlation is accounted for, never more. The convergent-driver headline is unaffected in substance — 156 and 189 are the VHID convergent pair and both survive at least one clustering scheme, 189 survives both — but the four-position HIGH tier reported in §[3.3](#Causal-discovery) rests partly on i.i.d. resampling and should be read as an upper bound on selection confidence. The Bedford H3N2 cluster [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) is heavier (≈73 s per [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) fit at 7808×50 vs ≈13 s for VHID, so a full 3×200 is ≈50 min CPU); it is left to the shipped i.i.d. run here and flagged as a recommended robustness check.

### Discovered Dependency Structure
We map the learned relationships into the titer sink as a directed graph, where edge widths represent [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) stability. To evaluate whether the system operates as a strict star graph (parallel, mutually independent features directed into a common child) or exhibits intermediate dependencies, we performed a [partial-correlation](https://en.wikipedia.org/wiki/Partial_correlation) skeleton check across the parent nodes:
1. **Intermediate Mediation Check:** Testing each parent against the titer conditional on all remaining parents. A pure intermediate node would become conditionally independent and drop out.
2. **Mutual Independence Check:** Testing all parent-parent pairs conditional on the remaining parent set (excluding the titer node to prevent collider-induced bias/explaining-away effects).

{{ figure(src="figures/dag.png", n=15, alt="dag") }}

Our structural analysis shows that every parent node maintains a statistically significant direct association with the titer when conditioning on all other parents, indicating that the feature set contains no pure intermediates. However, the parent nodes are densely interconnected: we identify 5/10 significant parent-parent adjacencies in VHID and 21/28 in Bedford.
The graph is therefore not a strict star. The observational data confirms strong structural dependencies among the causal parents. While the direction of these internal edges or the exact proportion of mediated vs. direct effects cannot be uniquely oriented without direct interventional data, this dense linear interdependence parallels the non-linear [epistatic](https://en.wikipedia.org/wiki/Epistasis) surfaces captured by our second-order networks.

### Interpretable B-Spline KAN
To maintain predictive transparency without reverting to strict linear assumptions, we deployed a genuine [B-spline](https://en.wikipedia.org/wiki/B-spline) [Kolmogorov–Arnold Network](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Arnold_representation_theorem) (order-3 splines + [SiLU](https://en.wikipedia.org/wiki/Swish_function) residual), which optimizes learnable univariate functions along each network edge. This allows each first-layer connection to be characterized as an inspectable 1-D response curve: $f(\text{Grantham distance})$. Models were optimized using [Adam](https://arxiv.org/abs/1412.6980), $L_1$ spline regularization, and [early stopping](https://en.wikipedia.org/wiki/Early_stopping) on pinned random seeds to ensure numerical reproducibility. Feature importance is calculated as the standard deviation of each position's partial contribution evaluated strictly over the observed data distribution to avoid extrapolation artifacts.

{{ figure(src="figures/benchmark.png", n=16, alt="benchmark") }}

The [B-spline](https://en.wikipedia.org/wiki/B-spline) KAN implementation was validated on a synthetic additive function, recovering known functional forms and yielding an $R^2 \approx 0.99$. When applied to the empirical datasets, the KAN achieved predictive performance that tracked closely to [gradient boosting](https://en.wikipedia.org/wiki/Gradient_boosting), trailing by a small but stable margin under matched-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) ([Section 3.2](#Cross-Validated-$R^2$-and-Rigorous-Error-Bounds)). The learned per-position response curves across both datasets are predominantly monotone-decreasing with respect to [Grantham](https://doi.org/10.1126/science.185.4154.862) distance, indicating that increasing physicochemical divergence maps directly to decreases in cross-titer (greater immune escape).

### Capturing Epistasis: A Second-Order KAN
A comparison between depth-1 (additive) and depth-6 tree models indicates that a notable portion of the predictive signal depends on interaction effects, or [epistasis](https://en.wikipedia.org/wiki/Epistasis). To capture these interactions within an interpretable framework, we extended the KAN to second order by incorporating bivariate tensor-product spline surfaces, $g(x_i, x_j)$, across an interaction pool composed of the causal parents and top predictive positions. Group-sparsity penalties were applied to prune inactive pairs systematically. Cross-validation was nested relative to pool selection, ensuring that interaction pools were selected strictly within training folds to prevent data leakage.

To verify whether pairwise terms sufficiently capture the non-linear signal, we computed an interaction-order ladder across tree models of increasing depth:

{{ table(src="tables/table-16.html", n=17) }}

{{ table(src="tables/table-17.html", n=18) }}

The structural dynamics diverge between the two panels:
- On VHID, the interaction signal is primarily pairwise. The transition from additive to pairwise terms yields a major performance gain (+0.073), whereas higher-order interactions (3-way through 6-way) contribute a smaller cumulative increase (+0.050).
- On Bedford, pairwise interactions account for roughly half of the non-linear signal. The pairwise increment (+0.080) is closely matched by the cumulative higher-order contribution (+0.079). This indicates that higher-order [epistatic](https://en.wikipedia.org/wiki/Epistasis) complexes are prominent in the Bedford dataset.

While a third-order KAN could theoretically map these higher-order relationships, the exponential expansion of the parameter space reduces model interpretability. We treat the second-order KAN as a pragmatic baseline that visualizes the pairwise component. Evaluating the optimized tensor-product surfaces reveals clear biophysical patterns: synergistic escape regions (where co-occurring substitutions reduce titer beyond their additive expectations) and compensatory interaction surfaces.

#### Formal Verification of Epistatic Pairs

To formally cross-examine the top KAN-nominated interactions, we fitted standard [OLS](https://en.wikipedia.org/wiki/Ordinary_least_squares) interaction terms ($x_a \cdot x_b$) with [HC3](https://en.wikipedia.org/wiki/Heteroskedasticity-consistent_standard_errors) robust standard errors, applying [Benjamini–Hochberg](https://en.wikipedia.org/wiki/False_discovery_rate) (BH) correction alongside a non-parametric [distance-correlation](https://en.wikipedia.org/wiki/Distance_correlation) test:
- **Statistical Significance:** 10 out of 16 top KAN-nominated interaction pairs (4/8 in VHID; 6/8 in Bedford) demonstrate statistically significant linear interactions ($q < 0.05$). The estimated coefficients are small but precisely bounded ($\vert{}\beta\vert{} \approx 10^{-4}\text{--}10^{-3} \log_2\text{-titer}$ per [Grantham](https://doi.org/10.1126/science.185.4154.862) unit squared).
- **Structural Norms vs. Classical Significance:** KAN interaction norms do not systematically correlate with standard linear $p$-values; the KAN optimizes overall surface curvature rather than isolated multiplicative parameters. Thus, tensor surfaces serve primarily as an interpretability asset rather than formal inferential tests.

{{ figure(src="figures/figure-09.png", n=19, alt="figure 09") }}

### Cross-Method Convergence
We evaluated the alignment of our causal feature selection ([bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) frequency $\ge 0.5$) against three distinct association screens: top-15 KAN features, top-15 [XGBoost](https://en.wikipedia.org/wiki/XGBoost) features (by gain), and top-15 univariate associations (by $R^2$).

Because the three association metrics are computed over identical feature spaces, they function as correlated screens rather than independent lines of evidence. Mutual alignment represents a multi-perspective validation of structural relevance rather than independent replication. In the Bedford panel, the [causal discovery](https://en.wikipedia.org/wiki/Causal_inference) input was pre-screened using target association (filtering 123 loci to 60), meaning the univariate screen is mechanically related to the causal search input; the VHID panel (71 loci, unscreened) does not share this dependency.

Positions that consistently converge across all four independent and correlated screens represent our strongest candidate drivers:

{{ figure(src="figures/convergence.png", n=20, alt="convergence") }}

#### Summary of Convergent Driver Positions
- **VHID:** Mature positions 156 and 189 (both represent singleton linkage blocks with [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) frequencies $\approx 1.0$).
- **Bedford:** Mature positions 133, 158, and 189 (each maps to a single residue).
These positions map directly to classical HA head antigenic sites flanking the receptor-binding domain (Site A for 133; Site B for 156, 158, and 189). Features flagged by predictive models that show low causal [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) stability indicate predictive correlates rather than direct causes.

#### Functional Glycosylation Profiling
Since our per-position features utilize a symmetric encoding (evaluating changes at an isolated column), the model cannot explicitly represent the gain or loss of [N-linked glycosylation](https://en.wikipedia.org/wiki/N-linked_glycosylation) sequons ($\text{N-X-S/T}$ motifs, where $\text{X} \neq \text{P}$), which require a three-residue window. To incorporate this context, we mapped the convergent positions back to consensus HA1 sequences:
- **Sequon Mapping:** Among the primary convergent drivers, only Bedford mature position 133 forms the root of an active sequon ($\text{N-G-T}$ motif). Positions 156, 158, and 189 do not disrupt or form sequons across the consensus profiles.
- **Lineage Differences:** Position 144 carries an active sequon in the Bedford consensus profile but is absent in the VHID consensus. This highlights an encoding boundary: a driver call at a sequon-associated position indicates that column-specific variations track with titer shifts. However, the true physical mechanism may operate via glycan shielding rather than direct epitope contact.

#### Non-linear Verification of Omitted Features
The primary causal selection loop relies on the [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) conditional independence test, which detects linear [partial correlations](https://en.wikipedia.org/wiki/Partial_correlation). To verify whether functional features with purely non-linear dependencies were discarded, we re-tested positions flagged by $\ge 3$ predictive screens that exhibited low causal frequencies (Pattern-A positions) using the non-parametric [Kernel Conditional Independence](https://arxiv.org/abs/1202.3775) ([KCI](https://arxiv.org/abs/1202.3775)) test. We evaluated the null hypothesis:
$$H_0: P \perp \text{HI\_titer} \mid \text{Discovered Parents}$$

The results are summarized below (reporting the median $p$-values across 3 independent seeded runs on 1,000-row subsamples):
- **Pattern-A Validation:** Both VHID position 145 and Bedford position 159 exhibit structural independence under linear assumptions but demonstrate strong conditional dependence under the non-parametric [KCI](https://arxiv.org/abs/1202.3775) test. This alignment with curved KAN response splines identifies them as non-linear driver candidates that the linear causal framework omits.
- **Full Pipeline Scalability Limits:** Running a complete constraint-based discovery search with [KCI](https://arxiv.org/abs/1202.3775) substituted for [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) proved computationally intractable on the available infrastructure, exceeding a hard 180-second per-dataset execution threshold. Consequently, the global structural graphs remain conditional on the assumptions of the linear [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) test, which are validated for size control ([Section 3.3.2](#Calibration-of-the-Fisher-Z-Test-Under-Permutation-Nulls)) but may omit purely non-linear structures.

{{ table(src="tables/table-23.html", n=21) }}

#### Cross-dataset replication of the discovered structure

We evaluated the structural replication between the VHID and Bedford datasets by mapping discovered loci to shared mature H3 positions and calculating [Jaccard](https://en.wikipedia.org/wiki/Jaccard_index) similarity metrics against a 20,000-draw [permutation](https://en.wikipedia.org/wiki/Permutation_test) null:

| set | observed J | shared positions | null mean | p |
|-----|-----------:|-------------------|----------:|--:|
| PC | 0.167 | **158, 189** | 0.035 | 0.062 |
| GES | 0.182 | 133, 276 | 0.035 | 0.051 |
| Bootstrap (freq≥0.5) | 0.182 | **158, 189** | 0.034 | **0.047** |

The structural overlap exceeds chance expectations across all sets, demonstrating stable structural replication. However, the specific positions driving this replication depend on the algorithm: constraint-based [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) and the [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) selection replicate positions 158 and 189 ([bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) $p=0.047$), whereas score-based [GES](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Score-based%20causal%20discovery%20methods/GES.html) localizes its cross-dataset intersection at positions 133 and 276.

{{ figure(src="figures/figure-11.png", n=22, alt="figure 11") }}

### Backdoor-Adjusted Effect Sizes
If the target variable operates as a pure causal sink, the remaining selected parents could serve as a valid backdoor adjustment set to isolate a position's specific interventional effect. However, the [d-separation](https://en.wikipedia.org/wiki/Bayesian_network) goodness-of-fit test rejects the simple sink-star structure due to dense parent-parent adjacencies. This indicates that conditioning on co-parents can introduce confounding via unmodeled mediators or colliders. We therefore interpret these estimates as partial regression coefficients rather than as identified causal effects.

{{ table(src="tables/table-25.html", n=23) }}

{{ figure(src="figures/effect_sizes.png", n=24, alt="effect sizes") }}

The partial-regression coefficients were estimated using [ordinary least squares](https://en.wikipedia.org/wiki/Ordinary_least_squares), controlling for the co-parent set, paired with [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) 95% confidence intervals:
- **Systematic Shrinkage:** All adjusted effect sizes exhibit systematic shrinkage toward zero relative to their unadjusted marginal associations. This pattern is consistent with mitigating phylogenetic confounding, though it remains structurally confounded by internal parent dependencies.
- **Sign Consistency:** No feature displayed a sign reversal upon adjustment; all 13 primary candidates (5 in VHID, 8 in Bedford) maintained stable directional signs across both marginal and adjusted models.
- **Directional Dynamics:** The majority of estimated effects are negative, meaning larger physicochemical divergence maps to a lower $\log_2$ titer (greater immune escape). Position 289 (VHID) is a notable exception, showing a positive adjusted coefficient ($+0.024$). This indicates that larger substitutions at position 289 track with increases in cross-titer, consistent with a structural stabilization role rather than direct antibody evasion.

#### Driver versus Hitchhiker Differentiation. 
To separate functional drivers from passenger mutations, we mapped [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) selection frequencies against backdoor-adjusted effect sizes:
- **True Drivers (High Selection Stability + Substantial Effect Size):** Positions 158, 189, and 289 in VHID; positions 2, 133, and 189 in Bedford.
- **Stable Hitchhikers (High Selection Stability + Marginal Adjusted Effect Size):** Position 156 in VHID maintains a selection frequency of 1.0 but displays a minimal, non-robust effect size, indicating it is a passenger mutation tightly linked to functional drivers. Position 144 in VHID exhibits high effect variability, rendering its selection unstable.

{{ figure(src="figures/figure-13.png", n=25, alt="figure 13") }}

### Left-Censoring Sensitivity Analysis
The target variable is bounded by a [left-censoring](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) floor representing the assay's lower detection limit (undetectable titers $<10$ are encoded at a floor value of 5.0). This affects 493 pairs ($\le 10$) in VHID and 616 pairs in Bedford. Because these values concentrate in high-antigenic-distance regimes (52% of top-quartile [Grantham](https://doi.org/10.1126/science.185.4154.862) pairs in VHID are censored vs. 0% in the bottom quartile; [Mann–Whitney](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test) $p < 10^{-100}$), we evaluated whether effect rankings are artifacts of censoring configurations.

We recomputed the adjusted effects under three sensitivity states: (a) as shipped, (b) dropping all censored rows, and (c) recoding the floor value to 10.0. The estimated signs and overall performance ranks remained stable across all states. The only adjustments were single-step rank swaps among the lowest-impact positions (e.g., VHID 144$\leftrightarrow$156). At the same time, the primary convergent drivers (133, 158, 189) maintained stable parameters, confirming that [left-censoring](https://en.wikipedia.org/wiki/Censoring_%28statistics%29) does not systematically bias our classifications of primary drivers.

{{ figure(src="figures/figure-14.png", n=26, alt="figure 14") }}

#### Token-Level Attribution and Structural Identifiability
To transition from overall population-level effects to token-level (per-pair) credit attribution, we executed an identifiability audit:
- **Structural Block Resolution:** Because linkage collapse groups collinear positions prior to [causal discovery](https://en.wikipedia.org/wiki/Causal_inference), each selected driver represents an independent linkage block. There is no block-level confounding between selected drivers.
- **Natural Experiments:** Among virus-reference pairs that vary at $\ge 1$ driver, 35–36% qualify as single-driver natural experiments (655 pairs in VHID; 1395 in Bedford), where the sequence differs at exactly one driver block. Bedford position 193 (629 single-mismatch pairs) and VHID position 158 (255 pairs) represent major clear contrasts for prospective experimental validation. The remaining 64% are multi-block pairs requiring block-level resolution.
- **Counterfactual Stability Mapping:** Evaluating a leave-one-difference-out ensemble across [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) refits confirms that primary drivers maintain stable token-level credit boundaries. Consistent with our population-level analysis, VHID position 156 displays highly unstable token credit (sign consistency 0.58), confirming its characterization as a passenger mutation.

{{ figure(src="figures/figure-15.png", n=27, alt="figure 15") }}

{{ figure(src="figures/figure-16.png", n=28, alt="figure 16") }}

### Validating the Causal Structure
To evaluate the structural validity of the discovered titer-sink graph, we performed three complementary macro-validation tests:
1. **Global Goodness-of-Fit (Shipley's [d-Separation](https://en.wikipedia.org/wiki/Bayesian_network) Test):** Evaluates whether the implied conditional independence constraints are valid across the empirical joint distribution.
2. **Linkage-Group Bootstrap Stability:** Measures the structural reproducibility of edges when the entire selection pipeline is executed over independent data resamples.
3. **Direct Effect Bounds:** Quantifies the variance of the adjusted partial-regression coefficients.

{{ table(src="tables/table-31.html", n=29) }}

{{ table(src="tables/table-312.html", n=30) }}

The global goodness-of-fit test **systematically rejects** the simplified sink-star structure across both datasets ($p < 0.05$). This rejection occurs because the graph does not account for residual parent-parent dependencies arising from shared phylogenetic history.
Crucially, this does not invalidate the identification of the target's primary parent set; rather, it indicates that the graph functions as a localized feature-selection model for direct target parents rather than a complete generative model of the sequence-population structure. The linkage-group [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) stability scores approach 1.0 for the primary epitope groups, demonstrating that the identities of these driver blocks are highly reproducible under resampling.

### Continuous Per-Position Encoding and the Titer Markov Blanket
To evaluate whether preserving substitution magnitude impacts structural discovery, we executed an independent replication of the structural workflow on the VHID dataset using a continuous encoding scheme. This modification aligns with three structural revisions:
1. **Continuous Feature Space:** Positions are encoded as the continuous scalar $L_2$ norm of the z-standardized 12-property mutation vector, preserving mutation scale while collapsing property dimensions.
2. **Target Reframing:** The primary deliverable is the target's [Markov blanket](https://en.wikipedia.org/wiki/Markov_blanket) (adjacency) rather than directed parent sets, since edges into a pure sink cannot be uniquely oriented by conditional independence alone.
3. **Algorithmic Path:** Linkage collapse was performed using continuous [Pearson](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) correlation ($\vert{}r\vert{} \ge 0.8$), followed by screened [Fisher-Z](https://en.wikipedia.org/wiki/Fisher_transformation) [FCI](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/FCI.html) discovery, cross-checked against the [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) and [BOSS](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Permutation-based%20causal%20discovery%20methods/BOSS.html) algorithms over a 200-resample [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29).

#### Physicochemical Scalar Mapping
We validated the continuous 12-property $L_2$ scalar against standard [Grantham](https://doi.org/10.1126/science.185.4154.862) distances across the active alignment columns. The continuous metrics display high correlation while maintaining enhanced precision regarding residue-specific volume and charge trajectories.

{{ figure(src="figures/figure-17.png", n=31, alt="figure 17") }}

#### Continuous Markov Blanket Analysis
The continuous pipeline over the VHID dataset ($n=2751$) yields a High Confidence [Markov blanket](https://en.wikipedia.org/wiki/Markov_blanket) encompassing positions {156, 189, 278, 289}.

Comparing this with our binary results reveals a stable core: positions 156, 189, and 289 maintain high confidence across both encodings (e.g., position 289 yields a binary stability of 0.955 and a continuous stability of 0.960). The encodings diverge at two points: position 158 (High Confidence under binary flags) drops to unstable ($0.490$) under the continuous metric. In contrast, position 278 (Site C) is promoted from low binary stability ($0.40$) to high continuous stability ($0.970$). This divergence indicates that certain positions operate as binary switches, while others depend on the specific physicochemical distance of the substitution.

{{ table(src="tables/table-33.html", n=32) }}

{{ figure(src="figures/figure-18.png", n=33, alt="figure 18") }}

{{ table(src="tables/table-35.html", n=34) }}

#### Encoding Sensitivity Analysis
Rerunning [PC](https://causal-learn.readthedocs.io/en/latest/search_methods_index/Constraint-based%20causal%20discovery%20methods/PC.html) target-adjacency searches across binary, [Grantham](https://doi.org/10.1126/science.185.4154.862), and continuous $L_2$ encodings on the VHID panel confirms that positions 144, 156, 189, and 289 are robust across all three frameworks (pairwise [Jaccard](https://en.wikipedia.org/wiki/Jaccard_index) similarities $0.50\text{--}0.63$).

Exploratory analysis mapping individual property dimensions shows that because single amino acid substitutions modify all 12 property axes simultaneously, individual properties are structurally non-identifiable ([partial correlations](https://en.wikipedia.org/wiki/Partial_correlation) collapse to zero when conditioning axes on one another). Marginal correlations can rank which axis covaries most strongly with titer shifts (e.g., hydrogen-bond-acceptor properties at position 158; $\beta$-sheet preferences at position 189), but these cannot be interpreted as isolated causal effects.

{{ figure(src="figures/figure-19.png", n=35, alt="figure 19") }}

## Conclusion

By integrating target-oriented [causal discovery](https://en.wikipedia.org/wiki/Causal_inference) with interpretable machine learning, this study establishes a structural framework to separate HA positions that drive antigenic escape from linked passenger mutations. Evaluated strictly on [hemagglutination-inhibition](https://en.wikipedia.org/wiki/Hemagglutination_assay) data without structural or structural epitope priors, the pipeline maps its highest-confidence selections to classical HA head antigenic sites, rediscovering residues implicated in immune evasion across both H3N2 datasets.
The primary convergent core—encompassing mature positions 156 and 189 in VHID, and 133, 158, and 189 in Bedford—localizes to antigenic site B (flanking the receptor-binding domain) and site A. Position 189 is documented as a primary determinant of H3N2 cluster transitions; recovering this signal directly from observational titers indicates that the feature selection maps to verified antigenic mechanisms rather than dataset-specific artifacts.

Our validation battery highlights the structural boundaries of observational serology. Rejection of the simplified sink-star graph in global [d-separation](https://en.wikipedia.org/wiki/Bayesian_network) tests indicates that, while the pipeline isolates immediate target parents, it does not capture the dense network of phylogenetic dependencies among them. Consequently, adjusted effect sizes represent partial-regression coefficients rather than fully identified causal parameters.

Furthermore, our framework demonstrates that selection stability does not inherently imply functional causality; passenger mutations tightly linked to functional loci can achieve high [bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29) frequencies. This is illustrated by position 156 in the VHID dataset, which displays high selection frequency alongside small, non-robust effect sizes, characterizing it as a stable hitchhiker.

Phylogenetic linkage imposes physical limits on the resolution of individual residues in observational datasets. This is pronounced in the Bedford H3N2 panel, where the largest co-evolving linkage blocks encompass 88 and 55 positions, binding multiple head residues into single covarying units that cannot be resolved without interventional data.
Additionally, raw HI titers integrate multiple biophysical phenotypes, conflating head-epitope antibody binding with variations in receptor-binding avidity and unmodeled glycosylation structures. Because avidity-associated residues (including 145, 189, and 193) overlap our parent sets, individual residue attributions remain mechanistically complex under a raw-titer target.

These constraints guide the evaluation of sequence-to-antigenic maps for prospective surveillance. Grouped [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29) establishes realistic generalization boundaries: when forecasting titers against entirely unseen reference antisera, median predictive performance settles at $R^2 \approx 0.615$ for VHID and $\approx 0.498$ for Bedford, down from random-split baselines near $0.85$. Furthermore, temporal transport analysis indicates that forward-in-time predictions can become unstable when viral evolution crosses major cluster boundaries that are absent from the training data.

In conclusion, this pipeline provides a transparent, self-auditing framework that maps its own structural limits. While observational data can isolate co-evolving blocks and prioritize candidate drivers, resolving individual-residue causality within dense lineages requires integration with prospective [reverse genetics](https://en.wikipedia.org/wiki/Reverse_genetics), deep mutational scanning escape maps, and structurally isolated serological assays.

## References

1. **Smith, D. J., Lapedes, A. S., de Jong, J. C., Bestebroer, T. M., Rimmelzwaan,
   G. F., Osterhaus, A. D. M. E., Fouchier, R. A. M.** (2004). [Mapping the antigenic
   and genetic evolution of influenza virus](https://doi.org/10.1126/science.1097211).
   *Science* **305**(5682), 371–376.
2. **Bedford, T., Suchard, M. A., Lemey, P., Dudas, G., Gregory, V., Hay, A. J.,
   McCauley, J. W., Russell, C. A., Smith, D. J., Rambaut, A.** (2014). [Integrating
   influenza antigenic dynamics with molecular evolution](https://doi.org/10.7554/eLife.01914).
   *eLife* **3**, e01914.
3. **Du, E., Zhong, Z., Wang, P., et al.** (2023). [DPCIPI: A pre-trained deep learning
   model for predicting cross-immunity between drifted strains of Influenza
   A/H3N2](https://arxiv.org/abs/2302.00926). *arXiv:2302.00926.*
4. **Grantham, R.** (1974). [Amino acid difference formula to help explain protein
   evolution](https://doi.org/10.1126/science.185.4154.862). *Science* **185**(4154),
   862–864.
5. **Spirtes, P., Glymour, C., Scheines, R.** (2000). [*Causation, Prediction, and
   Search* (2nd ed.)](https://doi.org/10.7551/mitpress/1754.001.0001). MIT Press.
   (PC algorithm.)
6. **Chickering, D. M.** (2002). [Optimal structure identification with greedy
   search](https://www.jmlr.org/papers/v3/chickering02b.html). *Journal of Machine
   Learning Research* **3**, 507–554. (GES.)
7. **Zhang, J.** (2008). [On the completeness of orientation rules for causal discovery
   in the presence of latent confounders and selection
   bias](https://doi.org/10.1016/j.artint.2008.08.001). *Artificial Intelligence*
   **172**(16–17), 1873–1896. (FCI.)
8. **Liu, Z., Wang, Y., Vaidya, S., et al.** (2024). [KAN: Kolmogorov–Arnold
   Networks](https://arxiv.org/abs/2404.19756). *arXiv:2404.19756.*
9. **Pearl, J.** (2009). [*Causality: Models, Reasoning, and Inference* (2nd
   ed.)](https://doi.org/10.1017/CBO9780511803161). Cambridge University Press.
   (Backdoor adjustment / do-calculus.)
10. **Chen, T., Guestrin, C.** (2016). [XGBoost: A scalable tree boosting
    system](https://doi.org/10.1145/2939672.2939785). *KDD 2016*, 785–794.
11. **Koel, B. F., Burke, D. F., Bestebroer, T. M., van der Vliet, S., Zondag, G. C. M.,
    Vervaet, G., Skepner, E., Lewis, N. S., Spronken, M. I. J., Russell, C. A., Eropkin,
    M. Y., Hurt, A. C., Barr, I. G., de Jong, J. C., Rimmelzwaan, G. F., Osterhaus, A. D.
    M. E., Fouchier, R. A. M., Smith, D. J.** (2013). [Substitutions near the receptor
    binding site determine major antigenic change during influenza virus
    evolution](https://doi.org/10.1126/science.1244730). *Science* **342**(6161),
    976–979.
12. **Neher, R. A., Bedford, T., Daniels, R. S., Russell, C. A., Shraiman, B. I.** (2016).
    [Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza
    viruses](https://doi.org/10.1073/pnas.1525578113). *Proceedings of the National
    Academy of Sciences* **113**(12), E1701–E1709.
13. **Łuksza, M., Lässig, M.** (2014). [A predictive fitness model for
    influenza](https://doi.org/10.1038/nature13087). *Nature* **507**(7490), 57–61.
14. **Harvey, W. T., Benton, D. J., Gregory, V., Hall, J. P. J., Daniels, R. S., Bedford,
    T., Haydon, D. T., Hay, A. J., McCauley, J. W., Reeve, R.** (2016). [Identification of
    low- and high-impact hemagglutinin amino acid substitutions that drive antigenic drift
    of influenza A(H3N2) viruses](https://doi.org/10.1371/journal.ppat.1005526). *PLoS
    Pathogens* **12**(4), e1005526.
15. **Wiley, D. C., Wilson, I. A., Skehel, J. J.** (1981). [Structural identification of the
    antibody-binding sites of Hong Kong influenza haemagglutinin and their involvement in
    antigenic variation](https://doi.org/10.1038/289373a0). *Nature* **289**(5796),
    373–378.
16. **Caton, A. J., Brownlee, G. G., Yewdell, J. W., Gerhard, W.** (1982). [The antigenic
    structure of the influenza virus A/PR/8/34 hemagglutinin (H1
    subtype)](https://doi.org/10.1016/0092-8674%2882%2990135-0). *Cell* **31**(2),
    417–427.
17. **Doud, M. B., Lee, J. M., Bloom, J. D.** (2018). [How single mutations affect viral
    escape from broad and narrow antibodies to H1 influenza
    hemagglutinin](https://doi.org/10.1038/s41467-018-03665-3). *Nature
    Communications* **9**, 1386.
18. **Lee, J. M., Eguia, R., Zost, S. J., Choudhary, S., Wilson, P. C., Bedford, T.,
    Stevens-Ayers, T., Boeckh, M., Hurt, A. C., Lakdawala, S. S., Hensley, S. E., Bloom,
    J. D.** (2019). [Mapping person-to-person variation in viral mutations that escape
    polyclonal serum targeting influenza hemagglutinin](https://doi.org/10.7554/eLife.49324).
    *eLife* **8**, e49324.
19. **Lou, Y., Caruana, R., Gehrke, J., Hooker, G.** (2013). [Accurate intelligible models
    with pairwise interactions](https://doi.org/10.1145/2487575.2487579). *KDD 2013*,
    623–631.
20. **Hensley, S. E., Das, S. R., Bailey, A. L., Schmidt, L. M., Hickman, H. D., Jayaraman, A.,
    Viswanathan, K., Raman, R., Sasisekharan, R., Bennink, J. R., Yewdell, J. W.** (2009).
    [Hemagglutinin receptor binding avidity drives influenza A virus antigenic
    drift](https://doi.org/10.1126/science.1178258). *Science* **326**(5930), 734–736.

---

*Data: [`influenza-hi-antigenic-distance`](https://zenodo.org/records/21339272)
repository (CC-BY-4.0). Code and this notebook are
released alongside it. Causal discovery uses the [`causal-learn`](https://causal-learn.readthedocs.io/en/latest/)
library; the KAN is a
custom PyTorch implementation in `src/bspline_kan.py`. The full pipeline — linkage
collapse, causal discovery, B-spline KAN, and cross-method convergence — is packaged
as the reusable `kan-causal-antigenic-workflow` skill.*
