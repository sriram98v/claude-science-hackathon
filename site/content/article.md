+++
title = "Article"
path = "article"
template = "article.html"
+++

This project was built as part of the [Claude Science Hackathon](https://cerebralvalley.ai/e/built-with-claude-life-sciences).

*A reproducible research notebook.* Run top-to-bottom (`Kernel → Restart & Run All`) to regenerate every table and
figure from the raw data repository shipped alongside this notebook.

## Abstract

Identifying the specific hemagglutinin (HA) mutations that causally alter antibody recognition remains a major challenge because dense viral phylogenies tightly link true escape drivers with passenger mutations. While sequence-based models can predict antigenic distance accurately, they typically conflate the evolutionary question of which mutations sweep a population with the mechanistic question of which mutations physically disrupt antibody binding. To probe this ambiguity, we frame the target as a provenance-independent, type-level interventional contrast at the antigen–antibody interface, accepting linkage-driven resolution loss as an inherent structural feature of observational serology rather than a pipeline defect. We analyze two H3N2 hemagglutination-inhibition (HI) datasets by collapsing co-evolving positions and applying target-oriented causal discovery (PC, GES, FCI) ranked by 200-resample bootstrap stability, paired with an interpretable B-spline Kolmogorov–Arnold Network (KAN) whose per-position response curves are directly inspectable and which extends to second order to visualize position-by-position epistasis. Under matched-fold cross-validation the first-order KAN trails black-box gradient boosting by a small but statistically robust margin (≈ 0.025–0.028 lower R²; Wilcoxon p < 10⁻⁵); a second-order KAN closes most of this gap, but only under a different 5-fold protocol, so the closure is suggestive rather than established. Cross-method convergence — agreement between one causal screen and three correlated association screens rather than four independent tests — flags candidate drivers overlapping classical antigenic sites A and B: in VHID the convergent pair sits at mature H3 positions 156 and 189, though 156 behaves as a stable hitchhiker (near-fixed frequency, small non-robust effect) rather than a driver while position 289 is the one doubly-robust candidate; in Bedford the convergent positions are mature 133, 158, and 189, with 158 fragile to encoding. Under cluster resampling, mature position 189 is the single most robust signal. Adjusted partial-regression effect sizes shrink systematically relative to marginal associations — a pattern consistent with, but not identifying, removal of phylogenetic confounding — and a Shipley d-separation test in fact rejects the sink structure from which the backdoor argument is read. Honest grouped (leave-serum-out) cross-validation places attainable accuracy well below random-split estimates (median R² ≈ 0.615 for VHID, ≈ 0.498 for Bedford). We therefore present this work as a stability-ranked, confounding-audited, linkage-aware feature-selection and interpretable-prediction pipeline that flags candidate biophysical drivers rather than a method that identifies causation, and we leave whether the recovered structure transports to prospective vaccine-strain selection to future evaluation.

## 1. Introduction

Influenza viruses evolve continuously under the selective pressure of population immunity, accumulating substitutions in surface glycoproteins that enable immune escape. This process of antigenic drift causes circulating strains to gradually diverge from those recognized by prior immunity, leading to a decline in vaccine effectiveness when the selected vaccine strain no longer matches circulating variants. This mismatch is a graded quantity termed antigenic distance, and measuring it accurately is a prerequisite for interpreting vaccine effectiveness and designing broadly protective vaccines. Traditionally, antigenic distance is measured using the hemagglutination-inhibition (HI) assay, in which a lower cross-titer indicates greater immune escape. While frameworks like antigenic cartography have been indispensable for mapping these titers, generating the required serum panels is slow, costly, and difficult to reproduce across laboratories. This logistical bottleneck has motivated a long line of sequence-based predictive models that forecast HI titers directly from hemagglutinin (HA) sequences. However, while modern machine learning approaches can achieve high predictive accuracy, they generally fail to explain why a strain escapes, as they identify which positions predict titer rather than which positions causally drive the change.

Because influenza strains share a dense phylogenetic history, HA positions are strongly linked, meaning a passenger mutation that merely rides along with a true escape driver looks just as predictive as the driver itself. To resolve this ambiguity, this study introduces a strict conceptual reframing. The common phrasing of identifying mutations that drive antigenic drift often conflates an evolutionary question (which substitutions were favored by natural selection and swept through the population) with a mechanistic one (which substitutions, if introduced into a strain, physically alter how strongly a specific antiserum recognizes it). This study focuses explicitly on the second, mechanistic question. We define our target estimand not as a historical claim about viral evolution, but as a provenance-independent, type-level interventional contrast at the antigen-antibody interface. Ideally, this contrast represents a precise biophysical experiment: take reference virus B, introduce the residue carried by virus A at a single HA position p while holding the rest of the protein fixed, and measure the resulting change in HI titer. The magnitude of this effect is dictated entirely by the chemistry of the structural footprint, regardless of whether that mutation arose via positive selection or as a neutral spontaneous variant. While evolutionary provenance is irrelevant to the biophysical effect itself, it is highly decisive for whether that effect can be identified from observational data. The virus's selective history acts as a major phylogenetic confounder, causing distinct mutations to co-occur within linkage blocks. Acknowledging this distinction allows us to treat linkage-driven resolution loss as an inherent structural feature of observational HI data rather than a flaw in the pipeline.

**Related works.** Antigenic cartography revealed the punctuated cluster structure of H3N2 drift by embedding HI tables into a low-dimensional antigenic map, and subsequent work unified antigenic and genetic evolution within a joint phylogenetic model. More recent sequence-based predictors forecast cross-immunity between drifted strains from sequence alone with strong accuracy. These approaches are primarily predictive or descriptive. They map, model, or forecast antigenic distance, but none explicitly learn causal structure over individual HA positions with quantified stability. The antigenic sites that those positions fall in were originally defined structurally and serologically, such as sites A to E on the H3 head, and were later refined at substitution resolution and by deep-mutational-scanning escape maps. We aim to see how far these same positions can be recovered structure-free, from HI titers alone. This approach is complementary to influenza fitness and phenotype-prediction models, which forecast strain frequencies or antigenic phenotypes but do not isolate per-position candidate drivers, and it builds directly on the intelligible-model literature for pairwise interactions.

**Our contributions.** We combine several ingredients that are not usually applied together for HI data to build a conservative, interpretable, and self-audited feature-selection pipeline. First, we implement linkage collapse with target-oriented causal discovery: we merge near-deterministic co-evolving positions into representative loci, then learn the direct-cause neighborhood of the HI target using PC, GES, and FCI algorithms, ranking every candidate by 200-resample bootstrap stability. Second, we deploy a genuine B-spline Kolmogorov–Arnold Network (KAN), an interpretable non-linear predictor whose learned per-position response curves are directly inspectable, and extend it to second order to capture and visualize position-by-position epistasis. Third, we assess cross-method convergence, treating it honestly as agreement between one causal screen and three *correlated* association screens (the KAN, gradient boosting, and univariate association) rather than four independent lines of evidence; positions flagged across these screens form our strongest candidate-driver claims. We then subject those claims to a battery of audits: a permutation calibration showing that the Fisher-z conditional-independence test holds near-nominal size on our binary, left-censored data; a cluster (by-virus and by-serum) bootstrap that separates a robust convergent core from an over-optimistic stability tier; a left-censoring sensitivity analysis of the adjusted effect sizes; and a token-level identifiability audit of the discovered structure. Throughout, we estimate backdoor-adjusted per-position effect sizes but report them as partial-regression coefficients, because the baseline adjustment assumptions are rejected in-sample. The result is a stability-ranked, confounding-audited, linkage-aware feature-selection and interpretable-prediction pipeline that flags candidate biophysical drivers, not a demonstration of identified causation.

The remainder of the notebook is organized as an executable paper. Section 2 (Methods) introduces the two H3N2 HI datasets and explains how their feature matrices are derived from the raw data. Section 3 (Results) then carries the analysis in full, each subsection stating how a step is performed, presenting its output, and interpreting it: the predictive benchmark and its leakage-free cross-validated comparison (§ <a href="#3.1-Predictive-benchmark">3.1</a>–3.2); target-oriented causal discovery together with the pre-collapse linkage block sizes, a permutation calibration of the Fisher-z independence test, and a cluster bootstrap of candidate stability (§ <a href="#3.3-Causal-discovery">3.3</a>), followed by the discovered dependency structure (§ <a href="#3.4-The-discovered-dependency-structure">3.4</a>); the interpretable B-spline KAN and its second-order epistasis extension (§ <a href="#3.5-Interpretable-prediction:-a-B-spline-KAN">3.5</a>–3.6); cross-method convergence (§ <a href="#3.7-Cross-method-convergence">3.7</a>); backdoor-adjusted effect sizes and their sensitivity to titer left-censoring (§ <a href="#3.8-Backdoor-adjusted-effect-sizes">3.8</a>); three DAG-validation tests (§ <a href="#3.9-Validating-the-causal-structure">3.9</a>); and a continuous per-position encoding that re-examines the titer Markov blanket (§ <a href="#3.10-Continuous-per-position-encoding-and-the-titer-Markov-blanket">3.10</a>). The narrative moves from how much of the HI signal is learnable and how non-linear it is, through the discovered candidate drivers and their audited dependency structure, to per-position effect sizes and the validation — and in-sample rejection — of the discovered graph. Section 4 (Conclusion) interprets the convergent positions biologically, states plainly where residue-level attribution is limited and what assumptions the causal framing rests on, and looks ahead. Section 5 lists references.

## 2. Methods

This section describes the data. We introduce the two HI panels and how their
per-position feature matrices are prepared from raw sequence data; the analytical
procedures that consume these matrices are described step-by-step in the Results
section, where each is presented alongside the output it produces. Configuration
(random seed, linkage threshold, bootstrap count, tier cutoffs) is centralized in
`src/analysis.py`, and two expensive steps — the causal bootstrap (`RECOMPUTE_CAUSAL`,
the H3N2 60-locus graph takes ~30 min) and the repeated k-fold cross-validation
(`RECOMPUTE_CV`, ~40 min) — are gated behind environment flags, loading shipped
`results/` artifacts by default so a full top-to-bottom run is fast and deterministic.

### 2.1 Datasets

The study uses two H3N2 virus × reference-strain HI panels, each shipped in the data
repository as per-position HA1 feature matrices with two encodings. The **binary
mismatch** encoding sets a position to `1` when the virus and reference residues differ
at that HA1 alignment position and `0` otherwise, and is used for causal discovery. The
**Grantham** encoding gives the Grantham (1974) physicochemical distance between the two
residues (gap-aware), and is used for predictive modeling and the KAN response curves.
The final column of every matrix is the raw `HI_titer`; we model `log2(HI_titer)`. The
two study datasets map to the repository as **VHID H3N2** → `VHID/` (n=2751) and **Bedford
H3N2** → `Bedford/H3/` (n=7808). **Provenance.** The VHID panel derives from the DPCIPI
release (arXiv:2302.00926); its HA sequences carry GenBank accessions but no
collection-year metadata (0/2751 pairs), which bounds the temporal analyses below to the
year-stamped Bedford panel. We begin by importing the analysis module and fixing
configuration.

**Numbering note.** Throughout, results prose reports **mature H3 residue numbers**; the
feature matrices are indexed by HA1 alignment column, and the two differ by a fixed
per-panel offset. The VHID reference is gapless from mature Q1, so its columns equal
mature H3 residues (offset 0). The Bedford H3N2 alignment carries a 9-residue
signal-peptide prefix and one internal gap (column 17), so mature Q1 sits at column 10
and, for every position reported here (all ≥ column 143), **mature residue = column − 10**
(e.g. columns 143/168/199 → mature 133/158/189). The offset is a uniform −10 only from
column 18 onward; columns 10–16 carry −9 because of the single alignment gap at column 17,
but no reported position falls below column 143, so the constant −10 applies throughout.
This was verified against the column-wise consensus reference (`SYILCLVFAQKLPGND-NST…`,
matched by 119/122 Bedford reference strains at the study columns) and the gapless VHID
references (45/45 at offset 0). Feature spaces remain strictly per-panel; the mapping is
provided only for biological interpretation.

#### Provenance and assay conditions

Both panels are H3N2 hemagglutination-inhibition (HI) datasets, but they differ in
provenance and in how much of the assay context can be audited. The **Bedford H3N2**
panel is the curated HI table from Bedford et al. (2014, *eLife*; github.com/trvrb/flux),
with GISAID/GenBank accessions and isolate collection years populated for every pair.
The **VHID** panel is drawn from the DPCIPI dataset (Du et al., *DPCIPI*, arXiv:2302.00926);
each strain carries a GenBank protein accession (e.g. AY661038), but the cleaned pair
table has **no per-pair isolate collection year** (0/2751 populated), and neither panel
records the **RBC species** used in the assay or whether the sera are **ferret or human**.
Those assay conditions therefore cannot be audited for either panel, and the absence of
collection-year metadata means the VHID panel cannot be split or blocked by era. We report
titers as provided and treat the two panels as independent replications of the same
structure-learning question rather than as calibrated to a common assay protocol.

#### 2.1.1 Reproducing the feature matrices from raw data

So that the study starts from raw data rather than shipped derivatives, we regenerate
the feature matrices from the cleaned virus × reference pair tables using the
repository's own build scripts (`scripts/build_vhid_matrices.py`,
`scripts/build_bedford_matrices.py`), which depend only on numpy and pandas. Running
them here makes the matrices provably the shipped ones (verified by SHA-256).

{{ table(src="tables/table-01.html", n=1) }}

Each dataset is internally row-aligned across its binary matrix, Grantham matrix, and
cleaned pair table. Feature spaces are comparable *within* a lineage but not across
lineages (each lineage was HA1-trimmed against its own reference), so we analyze the
two datasets independently and compare only which HA positions emerge.

#### 2.1.2 The target distribution

Before any modeling we look at what is being predicted. The panel below shows the
distribution of log2 HI titer in each dataset; its spread sets the scale against which
every R² reported in the Results should be read.

{{ figure(src="figures/target_distribution.png", n=2, alt="target distribution") }}

## 3. Results

We now walk through each analytical step and report what it produced. Because the
Methods section covers only the data, each subsection below opens by stating *how* the
step is done, then presents its output and interprets it — moving from *prediction* to
*cause*: first how much of the HI signal is learnable and how nonlinear it is, which
fixes what any mechanistic account must explain; then the positions that plausibly
*drive* the signal and their dependency structure; then the interpretable predictor and
the epistasis it exposes; and finally the convergent driver claims, their
confounding-adjusted effect sizes, and the validation of the discovered graph. No number
or figure is presented without interpreting how it supports the point.

### 3.1 Predictive benchmark

We first establish how much of the HI signal is learnable and whether it is
multivariable and nonlinear. For each dataset we compare, on a held-out 20% test split,
the best single position (max univariate R²), LASSO, Ridge, and XGBoost.

The benchmark below reports, for each dataset, the held-out test R² of the best single
position, LASSO, Ridge, and XGBoost.

{{ table(src="tables/table-03.html", n=3) }}

In every dataset the ordering is XGBoost ≳ LASSO ≈ Ridge ≫ best single position. The
large gap from best-single-position to the multivariable models shows the signal is
**distributed across many positions**; the XGBoost margin over the linear models
indicates **nonlinear / interaction structure** — the motivation for the KAN
(Section 3.5). We add the KAN test R² to this comparison in Section 3.5.

### 3.2 Cross-validated R² with confidence intervals

A single split gives a point estimate but no uncertainty, so we repeat k-fold
cross-validation to attach 95% confidence intervals to each method's R², and — because a
method comparison is only fair when every method sees the same folds — we score all four
under one protocol. **All four methods use the same 5×4 repeated k-fold protocol** (Round 1
cross-validated the KAN on only 3 folds against the others' 20, making any overlap a
degrees-of-freedom artifact; that is fixed here). XGBoost selects its boosting rounds by
**early stopping on an inner validation split carved from the training fold only** — never
on the scored test fold, which in Round 1 leaked the test fold into round selection and
biased XGBoost upward. The fold scores are loaded from `results/cv_r2_folds.json` by
default (recompute the linear/tree folds with `RECOMPUTE_CV=1`; the GPU-trained KAN folds
are always loaded — see `regen_cv.py` for the KAN CV driver). Because both methods are now
scored under an identical protocol, we compare them with a **matched-fold paired test**
(Wilcoxon signed-rank on the 20 per-fold differences) rather than eyeballing interval
overlap.

That paired test compares *methods*, and is valid whatever the splitting scheme; the R²
*level*, however, is not. The repeated k-fold above partitions **pairs** (virus ×
reference-serum) at random, so the same virus and the same reference antiserum recur across
train and test folds — the model can memorize a strain's antigenic profile rather than
generalize to an unseen strain. This is the pair-level leakage Round 1 flagged (prior M3),
and it inflates every random-split R². We therefore treat the random-split numbers only as
the inflated baseline being corrected, and take the **leakage-free grouped-CV estimate**
(below) as the honest predictive headline.

{{ figure(src="figures/cv_r2.png", n=4, alt="cv r2") }}

Under this random-split baseline XGBoost has the highest cross-validated R² on both
datasets (VHID 0.845, H3N2 0.613), now that its test-fold early-stopping leak is removed,
and the KAN sits just below it (VHID 0.820, H3N2 0.585). These absolute levels are the
inflated ones — corrected downward by grouped CV below — but the *method ranking* they
support is leakage-robust, because the matched-fold paired test compares the two on
identical folds. That test is **decisive rather than ambiguous**: the KAN trails XGBoost by
0.025 R² on VHID (Wilcoxon p ≈ 2×10⁻⁶) and 0.028 on H3N2 (p ≈ 6×10⁻⁶) — a small but
statistically robust gap, not the interval-overlap "tie" the Round-1 3-vs-20-fold
comparison implied. The KAN's value is its additive interpretability (§<a href="#3.5-Interpretable-prediction:-a-B-spline-KAN">3.5</a>), not a
predictive edge; it recovers most of XGBoost's accuracy while exposing per-position
response curves XGBoost cannot. With the ranking settled, we turn to what these methods
actually achieve on the task they are for — predicting antigenic distance for a strain the
model has never seen — by removing the pair-level recurrence.

To get a leakage-free estimate we re-run all four methods under **grouped**
cross-validation: `leave-virus-out` (no virus appears in both train and test) and
`leave-serum-out` (no reference antiserum shared). These are the honest analogues of the
prediction task the model is actually for — reading antigenic distance for a *new* strain —
and the R² they yield, not the random-split value, is the number a downstream user of the
sequence-to-antigenic-distance map should expect. Precomputed in `results/cv_grouped.json`
(driver: `regen_cv.py`).

{{ table(src="tables/table-06.html", n=5) }}

The leakage-free numbers are **substantially lower** than the random-split baseline,
confirming that the random-CV R² was inflated by strain/serum recurrence. Under
leave-serum-out — the strictest, most realistic held-out task (predict titers against a
*new* antiserum) — XGBoost falls from its inflated ~0.85 to a median **0.615** on VHID and
from ~0.61 to **0.498** on H3N2. XGBoost remains the strongest method under every grouping;
the linear models degrade sharply and Ridge is unstable under leave-virus-out on H3N2 (a
single fold with R² = −12.5, so medians rather than means are reported throughout this
table). **These grouped medians — VHID ≈ 0.615, H3N2 ≈ 0.498 — are the honest predictive
headline** of this study, and they set the expectation for the map's use on strains outside
the training panel.

#### 3.2.1 Temporal transport: predicting future antigenic clusters

VHID carries no collection-year metadata (0/2751 pairs), so this train-past / test-future
analysis is **Bedford H3N2 only** (1968–2010). Using an expanding window, we train on all
pairs up to year *t* and test on the next 5-year block.

| train ≤ | test window | n_test | Ridge R² | XGBoost R² |
|--------:|-------------|-------:|---------:|-----------:|
| 1990 | 1991–1995 | 927 | −4.70 | 0.43 |
| 1995 | 1996–2000 | 369 | 0.35 | **0.60** |
| 2000 | 2001–2005 | 3636 | −3.73 | **−0.37** |
| 2005 | 2006–2010 | 1297 | −0.59 | 0.25 |

**Finding.** Forward-in-time generalization is unstable and, in the worst window, negative:
XGBoost swings from R²=0.60 to **−0.37**, versus 0.613 under random-split CV. This is the
sharpest statement yet of the study's transportability limit — a model fit on past seasons
does not reliably predict titers for future antigenic clusters, consistent with drift
crossing cluster boundaries the training data never saw.

{{ figure(src="figures/figure-03.png", n=6, alt="figure 03") }}

#### 3.2.2 Does property encoding transport better across antigenic clusters?

A natural hypothesis is that encoding substitutions by their **physicochemical property shift** (the 12-property $L_2$ scalar) rather than by raw amino-acid identity should *transport* better to future strains: a novel substitution the model never saw as a token can still be represented by the magnitude of its charge/volume/H-bond change, and so inherit an effect estimate from physicochemically similar substitutions. We tested this directly on Bedford H3N2 (the only year-stamped panel) by repeating the train-past / test-future protocol of §<a href="#3.2.1-Temporal-transport:-predicting-future-antigenic-clusters">3.2.1</a> under three encodings — binary (which-residue-changed), Grantham (scalar chemical distance), and the 12-property $L_2$ vector — holding the position set and model fixed.

**The hypothesis is not supported.** For the unbiased-XGBoost model the mean future-$R^2$ is 0.19 (binary), 0.22 (Grantham) and 0.18 (12-property $L_2$): property encoding is essentially tied with, or marginally below, raw binary, and the single scalar Grantham distance is the most consistent of the three — undercutting the specific claim that richer multi-axis property structure buys graceful degradation. The $L_2$ encoding does lead in the earliest window (train ≤ 1990: $R^2$ = 0.50 vs 0.34 binary) but is worst in the large 2001–2005 window ($R^2$ = −0.65 vs −0.41), which dominates any sample-weighted average. Encoding differences (≤ 0.04 in mean $R^2$) are within the noise of four windows, and all three encodings fail together in the 2000-cut window — a genuine antigenic distribution shift that no substitution-based representation recovers. The value of the property encoding therefore rests on its interpretability and its role in causal-structure stability (§<a href="#3.10-Continuous-per-position-encoding-and-the-titer-Markov-blanket">3.10</a>), **not** on cross-cluster predictive transport, which remains bounded by distribution shift for every encoding tested.

{{ figure(src="figures/figure-04.png", n=7, alt="figure 04") }}

### 3.3 Causal discovery

We treat HI titer as a causal **sink**: HA sequence differences cause the titer, not
vice-versa. Encoding this as background knowledge orients every feature→target edge
into the target, so the target's *parents* are its direct-cause candidates.

The pipeline is: **power filter** (drop invariant / rare positions) → **linkage
collapse** (merge positions with |φ|≥0.8 into one representative locus per co-evolving
block) → **PC** (constraint-based) and **GES** (score-based) discovery → **200×
bootstrap** for selection-stability tiers. Linkage collapse is the load-bearing step:
near-deterministic co-evolution violates faithfulness and makes search intractable, and
collapsing restores both. Every causal claim then attaches to a *block*, so we always
report block size. FCI (latent-aware) is run only for VHID, because causal-learn's FCI
*possible-d-sep* step does not terminate in practical time on the denser H3N2
collapsed graph (>30 min at ≥24 nodes); for those we report PC + GES + bootstrap, and
very wide graphs are screened to the top loci by target association first.

Running this pipeline yields, per dataset, the collapsed loci, the PC/GES parent sets, and the bootstrap-stability tiers below.

{{ table(src="tables/table-09.html", n=8) }}

Residual strong locus pairs (|φ|≥0.9) drop to **0** in both datasets, confirming
that linkage collapse succeeded and restored the conditions discovery needs. PC and GES
then run on the collapsed loci, and the bootstrap tiers below rank each candidate.

{{ table(src="tables/table-10.html", n=9) }}

{{ table(src="tables/table-102.html", n=10) }}

The bootstrap frequency tiers each candidate into **high** (≥0.9), **moderate**
(0.5–0.9), or unstable. Under iid resampling the high-stability parents are VHID
{156, 158, 189, 289} and Bedford mature {2, 133, 157, 158, 189, 193}. Positions with
`block_size > 1` (e.g. Bedford mature 2 = a 9-position block; the largest pre-collapse
blocks on this panel span 88 and 55 positions, §<a href="#3.3.1-Pre-collapse-linkage-block-sizes">3.3.1</a>) carry claims for the *whole* block,
not the single position. "High-stability" scopes to reproducibility of *selection*, not
causal correctness — a hitchhiker linked to a true driver is selected in nearly every
resample for the wrong reason.

Two cautions sharpen how far this tier can be pushed. First, the iid tier is an **upper
bound**: under a **cluster bootstrap** (§<a href="#3.3.3-Cluster-%28block%29-bootstrap:-is-selection-stability-an-artifact-of-i.i.d.-resampling?">3.3.3</a>, VHID) that resamples whole viruses or whole
sera rather than pairs, the four-position VHID set does *not* survive intact — only
{156, 189} persist by-virus and {189, 289} by-serum, with mature 189 HIGH in all three
schemes. The convergent pair 156/189 is robust to how the panel is resampled; the
four-position tier is not, and should be read as an upper bound. Second, the PC/GES
agreement is asymmetric across datasets: in Bedford H3N2, GES independently confirms the PC
parents mature 133/157/158, but in VHID the two algorithms **disagree almost entirely** —
per §<a href="#3.3-Causal-discovery">3.3</a>'s discovery output, VHID PC ({144, 156, 158, 189, 289}) and GES ({50, 133, 144,
145, 262, 276}) share only pos144. Under faithfulness and causal sufficiency two correct
learners should return the same equivalence class, so a disagreement this large diagnoses a
tension in those assumptions on VHID rather than a resolvable ranking. The high-stability
tiers therefore rest on the *bootstrap*, which **ranks** candidates by selection frequency
rather than **arbitrating** between the two algorithms — a distinction that matters most
exactly where PC/GES concordance is weak, as it is for VHID.

#### 3.3.1 Pre-collapse linkage block sizes

Linkage collapse merges positions co-evolving at |φ|≥0.8 into one representative locus, and every downstream causal claim attaches to a *block* rather than a single residue. The table and figure below quantify how much resolution that collapse absorbs, per dataset, and are the evidentiary source for the block sizes cited in the Conclusion.

The distributions are heavy-tailed: most blocks are singletons (VHID 61/71; Bedford H3N2 114/123), but a few very large blocks dominate. In **Bedford H3N2** the largest pre-collapse block spans **88 positions** (representative alignment column 181) and the second spans **55** (column 90), followed by a 26-position block (column 50); these single clade-linked units absorb many head positions that no observational method can resolve internally, which is also why H3N2 carries the lower predictive ceiling of the two panels. **VHID** is far less linked — its largest block is only **14 positions** (column 173) — so its collapsed loci map much more nearly one-to-one onto individual residues. `block_size > 1` therefore carries a claim for the whole block, not the single position (H3N2 alignment-column 11 = 9-position block is the worked example in §<a href="#3.4-The-discovered-dependency-structure">3.4</a>).

{{ figure(src="figures/block_size.png", n=11, alt="block size") }}

{{ table(src="tables/table-11.html", n=12) }}

{{ table(src="tables/table-112.html", n=13) }}

#### 3.3.2 Calibration of the Fisher-z test under the permutation null

Every edge in the causal discovery above is decided by a **Pearson partial-correlation Fisher-z** conditional-independence test (causal-learn's `fisherz`) applied to **0/1 binary mismatch loci** and a left-censored HI titer. That test's reference null is exact only under joint multivariate normality, which these data violate, so the operating threshold $\alpha = 0.01$ is *nominal* and its true size is not guaranteed. This subsection measures the test's actual size empirically.

**Protocol.** For each encoding we build a null in which the target is independent of every feature by *permuting* the log2-titer column (destroying any feature→titer dependence while preserving the feature–feature correlation structure that the conditioning sets see). We then run the **same** `fisherz` CIT the pipeline uses, over 250 permutations × 40 random (feature, conditioning-set) draws = **10,000 tests** per encoding, with conditioning-set size drawn uniformly from {0,1,2,3} and rows subsampled to 1,000 (fixed design). Under a correctly sized test the resulting p-values are Uniform(0,1) and the fraction below $\alpha$ equals $\alpha$. We test three encodings: the **VHID collapsed binary loci** (the actual pipeline input), the **same loci re-encoded with the continuous L2 physicochemical distance** (closer to the Gaussian working model, as a contrast), and a **top-40 screened subset of the Bedford H3N2 binary loci**.

**Result.** The empirical type-I error tracks the nominal level closely, and critically at the operating point $\alpha = 0.01$ every encoding — including the binary ones — lands within its 95% Clopper–Pearson interval of 0.01 (VHID binary 0.0103 [0.0084, 0.0125]; VHID L2 continuous 0.0102 [0.0083, 0.0124]; H3N2 binary 0.0094 [0.0076, 0.0115]). The QQ plots lie on the diagonal. At the looser $\alpha = 0.05$ the VHID binary encoding is very slightly liberal (0.054), while the H3N2 binary encoding is mildly conservative (0.045) and the continuous encoding is on-target (0.051); all deviations are small and the discovery is run at 0.01, not 0.05. **We therefore find no evidence that the Fisher-z test is anti-conservative at the operating threshold on these data: the binary encoding does not inflate the false-positive edge rate relative to nominal.** All structural claims in Section 3.3 should nonetheless be read as conditional on Fisher-z adequacy — the test detects only *linear* partial correlation, so it can miss purely nonlinear feature→titer dependence (the KCI cross-check in the sensitivity analysis addresses that separately); this calibration establishes only that its *size* is controlled, not that its *power* against every alternative is complete.

{{ figure(src="figures/fisherz_calibration.png", n=14, alt="fisherz calibration") }}

#### 3.3.3 Cluster (block) bootstrap: is selection stability an artifact of i.i.d. resampling?

The 200× selection-stability bootstrap above resamples HI **pairs** independently. But VHID's 2751 pairs derive from only 246 viruses crossed with 45 reference sera, and the grouped cross-validation in §<a href="#3.2-Cross-validated-R²-with-confidence-intervals">3.2</a> shows this clustering is decisive (held-out R² drops ~0.23 when folds respect virus grouping). Resampling pairs i.i.d. treats correlated pairs as independent draws and can therefore *overstate* how reproducibly a position is selected.

To test this directly we re-ran the identical collapse → PC parent-selection routine (`src/causal_helpers.py`), changing only the resampling **unit**: instead of drawing pairs, we draw whole **viruses** with replacement (leave-virus-out clusters), and separately whole **reference sera**. Everything else — linkage collapse, the screened 50-locus node set, α=0.01, terminal-target background knowledge — is held fixed, so the comparison isolates the effect of respecting clustering. All three schemes use B=200 on VHID.

**Finding.** The HIGH-stability set is **not** preserved under clustering. Under i.i.d. pairs it is {156, 158, 189, 289}; under **virus** clustering only {156, 189} remain HIGH (158 and 289 fall to moderate), and under **serum** clustering only {189, 289} remain HIGH (156 and 158 fall to moderate). Only **mature 189** stays HIGH in all three schemes. Frequencies of the i.i.d.-HIGH parents move systematically **downward** toward the moderate range (mean change −0.04 under virus resampling, −0.12 under serum resampling; pos 158 falls from 0.95 to 0.69 under serum clustering) — i.e. the pipeline is *less* confident once pair correlation is accounted for, never more. The convergent-driver headline is unaffected in substance — 156 and 189 are the VHID convergent pair and both survive at least one clustering scheme, 189 survives both — but the four-position HIGH tier reported in §<a href="#3.3-Causal-discovery">3.3</a> rests partly on i.i.d. resampling and should be read as an upper bound on selection confidence. The Bedford H3N2 cluster bootstrap is heavier (≈73 s per PC fit at 7808×50 vs ≈13 s for VHID, so a full 3×200 is ≈50 min CPU); it is left to the shipped i.i.d. run here and flagged as a recommended robustness check.

### 3.4 The discovered dependency structure

We visualize each dataset's learned structure as a **titer-sink graph**: every retained
direct-cause candidate is a colored arrow into the HI-titer node, with width and opacity
proportional to its bootstrap stability and color by tier. But is this a strict *star* —
positions that each affect titer while being mutually independent? We test that directly
rather than assume it, running a partial-correlation skeleton over each parent set plus
titer, which asks two questions:

1. **Are any parents pure intermediates?** A position that reaches titer only through
   another (posₐ → pos_P → titer, no direct edge) would become conditionally
   independent of titer once the other parents are conditioned on, and drop out. We test
   each parent against titer given all *other* parents.
2. **Are the parents mutually independent?** We test every parent–parent pair for
   adjacency given all the *other parents* (excluding titer — since titer is a common
   child of the parents, conditioning on it would induce spurious dependence by the
   collider / explaining-away effect). Any surviving edge means the star is incomplete.

The grey arcs overlaid on the graph are the significant parent–parent adjacencies from
test 2. They are drawn **without arrowheads**: a partial correlation cannot distinguish
posₐ → pos_b from shared-ancestry confounding (posₐ ← L → pos_b), so their direction is
not identified.

The titer-sink graphs below show each dataset's retained direct-cause candidates and the parent–parent adjacencies from the skeleton test.

> **Encoding note (see §<a href="#3.10-Continuous-per-position-encoding-and-the-titer-Markov-blanket">3.10</a>).** This skeleton uses the **binary** Hamming encoding. The revised
> *continuous* L2 encoding (§<a href="#3.10.2-Full-revised-pipeline-on-VHID-%28n-=-2751%29">3.10.2</a>), re-run on VHID with a B = 200 bootstrap, gives a HIGH-stability
> titer Markov blanket of {156, 189, 278, 289}. The two encodings share a HIGH core of
> **{156, 189, 289}** — pos_289 is HIGH under *both* (binary bootstrap 0.955, continuous 0.960),
> **not** a continuous-only promotion. They differ on two positions: pos_158, a HIGH site-B driver
> under the binary flag (0.95), collapses to UNSTABLE (0.490) under the continuous L2 scalar, while
> pos_278 (site C) is the one position genuinely *promoted* by the continuous encoding (binary 0.40
> → continuous 0.970). The difference is expected —
> the binary flag counts *whether* a residue changed while the L2 scalar weights *how far* it moved —
> and the two are reconciled position-by-position in §<a href="#3.10-Continuous-per-position-encoding-and-the-titer-Markov-blanket">3.10</a>, not treated as competing blankets.

{{ figure(src="figures/dag.png", n=15, alt="dag") }}

The result is consistent across both datasets: **every parent stays directly
associated with titer** conditional on the others — so the parent set contains no *pure*
intermediates, and each carries its own direct effect (matching the nonzero backdoor
estimates in §<a href="#3.8-Backdoor-adjusted-effect-sizes">3.8</a>). But the parents are **densely interdependent** — parent–parent adjacency is 5/10 pairs
for VHID and 21/28 for H3N2 — so the structure is *not* a strict star. This is precisely why the
goodness-of-fit test in §<a href="#3.9-Validating-the-causal-structure">3.9</a> rejects the bare star, and the linear parent–parent
adjacencies here are the first-order shadow of the same interdependence the second-order
KAN captures nonlinearly as the epistasis surfaces in §<a href="#3.6-Capturing-epistasis:-a-second-order-KAN">3.6</a> (the 144×189 pair, for
instance, is both a significant parent–parent adjacency here and a top-ranked KAN
interaction surface; note the two views need not coincide — a linear partial-correlation
edge and a nonlinear interaction surface measure related but distinct things, so most
pairs appear in only one). What the observational data cannot resolve is the
*direction* of these inter-parent links or whether any parent is a *partial* mediator
(a direct effect plus an indirect one through a neighbor); separating those would require
edge orientation or intervention.

### 3.5 Interpretable prediction: a B-spline KAN

Kolmogorov–Arnold Networks (KANs) replace fixed node activations with *learnable
univariate functions on each edge*. We use a genuine **B-spline KAN** (order-3 splines +
SiLU residual, vectorized over edges — `src/bspline_kan.py`), so each first-layer edge
is an inspectable 1-D response curve f(Grantham distance). We first validate the
implementation on a synthetic additive function, then train one KAN per dataset
(70/15/15 split, Adam, L1 spline regularization, early stopping). Training is made **near-reproducible** by pinning all RNG seeds and requesting deterministic cuDNN/cuBLAS kernels (cell above); because `use_deterministic_algorithms` runs with `warn_only=True`, a few CUDA ops fall back to non-deterministic kernels, so we claim near- rather than bit-reproducibility. In practice the headline top-15 importance ranking — consumed by the §<a href="#3.7-Cross-method-convergence">3.7</a> convergence table and the §<a href="#3.7.1-Nonlinear-re-test-of-the-predictive-but-not-causal-positions">3.7.1</a> curvature reconciliation — is stable across reruns; only membership right at the ranking cutoff can shift by a position (see §<a href="#3.7-Cross-method-convergence">3.7</a>). Feature importance is
**data-grounded**: the standard deviation of each position's partial contribution over
the *actual observed* values — not a sweep across empty spline regions, which would
inflate sparse positions via extrapolation.

The synthetic-function validation and the per-dataset test R² (added to the benchmark comparison) are below.

{{ figure(src="figures/benchmark.png", n=16, alt="benchmark") }}

The KAN validates at R²≈0.99 on the synthetic additive function (recovering the true
curve shapes and flagging irrelevant inputs as flat), and on the real data it recovers most of
XGBoost's accuracy while trailing it by a small, statistically robust margin under matched-fold CV (§<a href="#3.2-Cross-validated-R²-with-confidence-intervals">3.2</a>: −0.025 R² VHID, −0.028 H3N2; both Wilcoxon p < 10⁻⁵). The KAN is preferred here for its additive interpretability, not for predictive parity. Its learned per-position
curves (shipped in `results/kan_curves.json`, and rendered in the study's curve
figures) are predominantly monotone-decreasing in Grantham distance: larger
physicochemical substitution → lower titer = greater escape.

### 3.6 Capturing epistasis: a second-order KAN

The first-order KAN is **additive** — a sum of univariate curves — and a diagnostic
(fitting depth-1 vs depth-6 gradient-boosted trees) shows that much of the accuracy gap
between additive models and full XGBoost on these datasets is **epistasis**:
position×position interactions (the per-dataset breakdown is in §<a href="#3.6-Capturing-epistasis:-a-second-order-KAN">3.6</a>'s ladder — predominantly
pairwise on VHID, with non-negligible higher-order structure on H3N2). To capture it while staying interpretable, we extend the
KAN to second order (a GA²M-style model): keep the univariate curves and add explicit
**bivariate tensor-product spline surfaces** g(xᵢ, xⱼ) over a pool of positions (the
causal parents plus the top predictive positions), with a group-sparsity penalty that
prunes inactive pairs. Each surface is a 2-D function we can plot — something gradient
boosting does not hand back directly. CV is **nested with respect to pool selection**: the interaction pool is re-selected on each
training fold only, so the reported R² carries no *pool-selection* leakage. (The KAN
architecture and regularization are fixed a priori per dataset, not tuned on the test folds.) The pool size needed per dataset tracks how localized the epistasis is within
the causal parents (measured separately): in both VHID and H3N2 the interactions
spread across more positions than the causal-parent set alone (a wider pool).

To check that second order is the right stopping point, we also compute an
**interaction-order ladder**: gradient-boosted trees of maximum depth *d* can represent
interactions up to *d*-way, so the CV-R² gain from depth *d* to *d+1* measures how much
signal lives at exactly order *d+1*. If the curve flattens after depth 2, pairwise terms
capture essentially all the available interaction structure; if it keeps rising, higher-order structure remains. As the two datasets show below, they differ on this.

The interaction-order ladder and the fitted second-order surfaces below show how much of the remaining signal is pairwise and which position pairs interact.

{{ table(src="tables/table-16.html", n=17) }}

{{ table(src="tables/table-17.html", n=18) }}

The ladder tells different stories on the two datasets, so we read it **per dataset** rather than pooling. On **VHID** the interaction signal is predominantly pairwise: the additive→pairwise jump (depth 1→2) is +0.073, whereas everything at third order and above adds only +0.050 cumulatively (2→6-way), with the 3→6-way increment just +0.020 — so pairwise terms capture most of the available interaction structure. On **H3N2**, by contrast, pairwise interactions capture only **roughly half** the interaction signal: the pairwise jump is +0.080, but the cumulative third-order-and-higher gain is +0.079 — essentially equal to the pairwise jump — and the 3→6-way increment (+0.042) is *larger* than the single 3-way step (+0.037). Higher-order structure is therefore **not negligible on H3N2**; the earlier ladder run happened to be monotone non-decreasing in depth on both datasets, so there is no evidence here that higher-order terms *hurt* held-out R².

We still stop the KAN at second order, but as a deliberate interpretability trade-off rather than because higher-order signal is absent. A third-order KAN (adding trivariate spline *volumes* gᵢⱼₖ) is buildable in principle, and on H3N2 it would be chasing a real +0.079 of signal above the pairwise level; the cost is far more terms and far less interpretable objects to plot, and the second-order surfaces retain the pairwise component we can visualize cleanly. We do **not** fit a direct higher-order KAN here, so the split of the H3N2 higher-order gain into genuine ≥3-way epistasis versus tree-model artifact is not resolved in this notebook — it is left as a limitation. The second-order KAN is thus a well-motivated interpretable model for the pairwise structure, exactly right for VHID and a deliberate lower bound on the interaction order for H3N2.

The interaction terms close the gap to gradient boosting, and unlike the trees the
model exposes *which* pairs interact and *how*. We fit on the full data and render the
strongest surfaces (ranked by norm). Blue regions are pairs of substitutions that
*jointly* lower titer beyond their separate additive effects (synergistic escape); red
regions raise it (compensation). Some top pairs are **parent×parent** (both partners
already have stable direct effects — the cleanest epistasis claims), others involve a
predictive-only partner, consistent with much of the interaction signal living outside
the high-stability causal-parent set.

#### Formal test of the nominated epistatic pairs

For each KAN-nominated pair we fit an OLS interaction term $x_a \cdot x_b$ (Grantham) controlling
for the causal parents, with HC3 robust standard errors and BH correction within dataset, plus a
distance-correlation nonparametric check.

**Finding.** **10 of 16 nominated pairs** show a significant linear interaction (q < 0.05):
4/8 in VHID, 6/8 in Bedford. Coefficients are small (|β| ~ 1e-4–1e-3 log2-titer per Grantham²)
but precisely estimated at these n. Critically, **KAN interaction norm does not track statistical
significance** — the norm ranks model curvature, not the presence of a formally testable epistatic
term — so the KAN surfaces should be read as an interpretability aid, not an inference. The
distance-correlation test flags nearly every pair, reflecting its power against any nonlinearity
rather than a specifically multiplicative effect.

{{ figure(src="figures/figure-09.png", n=19, alt="figure 09") }}

### 3.7 Cross-method convergence

To test whether the causal drivers are corroborated by other views of the data, we
intersect the causal moderate+high set (bootstrap ≥0.5) with three **association screens**:
the KAN top-15, the XGBoost top-15 (by gain), and the univariate top-15 (by R²). Only causal
discovery is methodologically distinct here — the KAN, XGBoost-gain, and univariate-R²
rankings are all association measures computed on the same feature matrix, so they are
**correlated screens, not independent method families**, and their mutual agreement is
partly mechanical. The right reading is therefore *one causal method agreeing with three
correlated association screens*, i.e. convergent corroboration of the causal claim rather
than four independent confirmations. Under a random-independence null, three top-15 screens
would be expected to share only ~0.3 positions on VHID (N=102) and ~0.03 on H3N2 (N=312), so
the observed overlap far exceeds chance — but that excess reflects the screens being
correlated, not independent. One caveat is specific to H3N2: its causal discovery was itself
pre-screened by target association (§<a href="#3.3-Causal-discovery">3.3</a>, `screened=True`, 123→60 loci), so on H3N2 the
univariate screen overlaps the causal step's own input and is **not** an independent
converging family; VHID (71 loci, unscreened) does not have this dependence. Positions
flagged by the causal method **and all three screens** are treated as the headline
causal-driver claims, read in that light.

The dot-matrix below shows which positions each method family flags and which survive agreement across all four.

> **Selection caveat (enrichment).** The clean concentration of discovered positions on the known head epitope sites is partly a property of *what got sampled*, not solely of what the pipeline can discriminate. Because the substitutions observed here are the ones antigenic drift selected and swept, the data is already enriched for antigenically consequential changes — neutral changes at the same positions largely failed to fix and are under-observed. This is a favourable bias for detecting real drivers, but it is a bias: the epitope localization should be read as the pipeline working on a sample pre-filtered toward genuine effects.

{{ figure(src="figures/convergence.png", n=20, alt="convergence") }}

The positions on which the causal method and all three association screens **agree**
(printed above) are the strongest driver claims: **VHID mature {156, 189}** and **Bedford
mature {133, 158, 189}**. The high-stability agreements — 156, 189 in VHID; 133, 158, 189
in Bedford — are each a **singleton** linkage block with bootstrap frequency ≈1.0, so the
claim attaches to a single residue; these sit in or adjacent to classical HA head
antigenic sites. Positions that enter only as large blocks (e.g. Bedford mature 2, a
9-position block) are block-level claims. Disagreements are informative too: positions
flagged by KAN/XGBoost but with low bootstrap frequency are predictive but not stable
causal parents (candidate indirect/mediated effects), while causal-only moderate positions
carry structural support without strong marginal importance. Exact convergence membership
right at the ranking cutoff can shift by a position between runs because the KAN is
retrained live under near- (not bit-) deterministic settings; the singleton high-stability
set is stable across reruns, so the headline claims are unaffected.

#### Glycosylation-sequon check on the convergent drivers

The per-position feature encoding is symmetric (a Grantham distance between the two
strains at a fixed alignment column) and therefore cannot represent **gain or loss of
an N-linked glycosylation sequon** (the N-X-S/T motif, X\u2260P), which is a three-residue,
directional property: a substitution at position *i*, *i*+1, or *i*+2 can create or
destroy a glycan that shields a whole epitope patch, an effect the single-column code
does not see. We therefore annotate whether each convergent driver participates in a
sequon in the panel consensus HA1 sequences. Of the convergent drivers, only **Bedford
H3N2 mature-133 is itself the root (N) of a sequon** (motif N-G-T); VHID 156/189 and
Bedford 158/189/157/193 are not in a sequon in either consensus. The recovered
consensus sequon set (mature N-sites 8, 22, 38, 63, 122, 126, 133, 144, 165, 246, 285)
matches the canonical H3 HA1 glycosylation sites. Notably, the head site **144 carries a
sequon in the Bedford consensus but not in the VHID consensus** \u2014 a concrete case
where the symmetric encoding cannot express a glycosylation difference between the two
panels. This is a limitation of the encoding, not of the causal search: a driver call at
a sequon-forming position should be read as "distance at this column co-varies with
titer", which may act through glycan shielding rather than direct epitope contact.

#### 3.7.1 Nonlinear re-test of the predictive-but-not-causal positions

The causal step (§<a href="#3.3-Causal-discovery">3.3</a>–§<a href="#3.4-The-discovered-dependency-structure">3.4</a>) decides every edge with the **Fisher-z** conditional-independence
(CI) test, which sees only *linear* dependence after linearly regressing out the conditioning
set. A dependence with no linear component is invisible to it: for $y=x^2+\varepsilon$ with $x$
symmetric about 0, $\mathrm{corr}(x,y)\approx 0$, so Fisher-z declares independence and PC drops
the edge. A genuinely **nonlinear** direct cause of titer can therefore land in the
`causal_freq ≈ 0` bucket of §<a href="#3.7-Cross-method-convergence">3.7</a> while still being flagged by the predictive families (KAN,
XGBoost, univariate). This is a limit on *linearity*, not on *calibration*: §<a href="#3.3.2-Calibration-of-the-Fisher-z-test-under-the-permutation-null">3.3.2</a> confirms
the Fisher-z test holds approximately nominal size on these binary, left-censored data, so
the concern here is a missed nonlinear edge, not an anti-conservative one.

Here we re-test exactly those **Pattern-A** positions — flagged by $\ge 3$ predictive families but
with causal bootstrap frequency below `MOD_CONF` — with a **nonparametric** CI test (**KCI**,
kernel-based) that detects arbitrary dependence. The null is
$H_0:\; P \perp \text{HI\_titer} \mid \text{parents}$, run twice (Fisher-z and KCI) on the *same*
columns, so the contrast is apples-to-apples. A position independent under Fisher-z but **dependent
under KCI** is a nonlinear direct-cause candidate the linear pipeline discarded.

*Caveats.* KCI uses causal-learn's default kernel bandwidths; its p-values come from a
gamma/permutation approximation and are stochastic, so we seed and report the **median of 3 runs**.
KCI is $O(n^2)$–$O(n^3)$, so for large datasets we subsample to 1000 rows (stated per row). A
`nonlinear_flag=True` is a *candidate*, not a settled orientation — direction still rests on the
sink prior. We first reproduce the synthetic $y=x^2$ discrepancy as a unit test, so the reader sees
the test *can* detect nonlinearity before trusting its null results on real positions.

> **Selection caveat (weak nulls).** Reading a `nonlinear_flag=False` as evidence of *no* nonlinear cause rests on a thinner sample than the converse. Because drift makes “position changed but titer did not move” a rare outcome, negative controls are under-represented, so this test — like the causal step it re-checks — estimates the *presence* of an effect more sharply than its *absence*. A null here tempers, rather than closes, the case for a position.

**Reading Task 1 + Task 3 together.** The null is $H_0: P \perp \text{HI\_titer} \mid \text{parents}$.
`nonlinear_flag=True` marks the case of interest: **Fisher-z fails to reject** $H_0$ (p > 0.05, the
position looks independent to the linear test) **but KCI rejects it** (p < 0.05, dependent under the
nonparametric test) — a nonlinear direct-cause candidate the linear pipeline discarded. `flag=False`
is *not* a claim of independence: it occurs whenever that specific linear-miss pattern is absent,
which includes the common case where the position is **dependent under both tests** (Fisher-z already
rejects $H_0$, so there is nothing for KCI to "rescue"). Read the two p-value columns, not just the
flag. The `kan_spline_curvature` column is the independent corroboration from §<a href="#3.5-Interpretable-prediction:-a-B-spline-KAN">3.5</a>–§<a href="#3.6-Capturing-epistasis:-a-second-order-KAN">3.6</a>: the
second-difference norm of the KAN's learned univariate response for that position (all other inputs at
baseline), normalized by response scale, so a straight-line response scores ≈0 and a curved one higher.
**Agreement** (KCI-dependent *and* a curved KAN spline) is a strong joint nonlinearity signal; a KCI
flag with a flat spline would suggest a KCI false positive or a KAN under-fit. Only a position that is
independent under **both** tests would be a clean predictive hitchhiker rather than a direct cause;
here every Pattern-A position remains titer-dependent under at least one test, and Bedford
mature 159 (site B) is the one whose dependence is **only** visible nonlinearly.

**Task 2 — CI-test swap sensitivity (bounded).** Task 1 re-tests a fixed candidate against a
fixed parent set. A nonlinear dependence can also change *which* conditioning set PC
arrives at, so here we attempt to re-run the whole target-parent discovery with the CI test
swapped from Fisher-z to KCI and diff the resulting parent sets. KCI-based PC is in the
same cost class that made FCI intractable in §<a href="#3.3-Causal-discovery">3.3</a>: even restricted to the **top-25 screened
loci** and a **1000-row subsample**, the kernel CI test is evaluated over too many
conditioning sets to finish in a practical budget. We therefore run KCI-PC in a worker
process under a hard 180 s per-dataset wall-clock kill; if it does not finish, we report
the Fisher-z parent set and record the kill in the `status` column. This is the documented
fall-back from the handoff spec.

Two consequences follow, and we state them plainly. Task 1 already answers the headline
question for the specific Pattern-A positions — the result this section turns on — by
re-testing each against a fixed parent set with the nonparametric KCI. But because the
nonparametric *discovery* is killed at scale, the robustness of the whole parent **set** to
a nonlinear CI test remains untested here: every structural claim in §<a href="#3.3-Causal-discovery">3.3</a>–§<a href="#3.9-Validating-the-causal-structure">3.9</a> is
therefore **conditional on the adequacy of the Fisher-z test**. That adequacy is not merely
assumed — §<a href="#3.3.2-Calibration-of-the-Fisher-z-test-under-the-permutation-null">3.3.2</a> calibrates it directly, showing the Fisher-z CI test holds approximately
nominal size on these binary, left-censored data — but the conditionality is real and we
carry it forward rather than let it lapse.

{{ table(src="tables/table-23.html", n=21) }}

#### 3.7.2 Cross-dataset replication of the discovered structure

Both H3N2 panels (VHID H3N2 n=2751 vs Bedford H3N2 n=7808) were run through identical
discovery; the discovered sets were mapped to common mature H3 numbering and the
cross-dataset Jaccard was compared against a 20,000-draw permutation null (random same-size
sets drawn from each panel's node set).

| set | observed J | shared positions | null mean | p |
|-----|-----------:|-------------------|----------:|--:|
| PC | 0.167 | **158, 189** | 0.035 | 0.062 |
| GES | 0.182 | 133, 276 | 0.035 | 0.051 |
| Bootstrap (freq≥0.5) | 0.182 | **158, 189** | 0.034 | **0.047** |

**Finding.** Overlap is modest but ≈ 5× the chance level for every set. The replicated
positions, however, differ by discovery method: **mature 158 and 189 replicate across the
constraint-based PC set (p=0.062) and the bootstrap set (p=0.047, significant)**, whereas
the score-based GES set replicates a *different* pair, **mature 133 and 276 (p=0.051)**. So
158/189 are the pair that recurs across datasets under PC and the bootstrap, not under
every algorithm; GES agrees that *something* replicates above chance but localizes it
elsewhere. Read together, the three rows convert "the drivers coincide" into a replication
statistic with a null, while showing that which residues carry the cross-dataset signal is
itself method-dependent.

{{ figure(src="figures/figure-11.png", n=22, alt="figure 11") }}

### 3.8 Backdoor-adjusted effect sizes

Because the HI target is a causal **sink**, one might argue that for any single parent the
*other* target-parents form a valid **backdoor adjustment set**, so that regressing titer on
a parent plus those co-parents recovers that parent's causal effect. That identification
argument rests on three assumptions about the *true* graph: (i) the adjustment set contains
no collider or descendant of the treatment, (ii) it blocks every backdoor path, and (iii) the
causal structure it is read off is correct. **In-sample, these assumptions are violated.** The
d-separation goodness-of-fit test of §<a href="#3.9-Validating-the-causal-structure">3.9</a> *rejects* the sink-star structure on both datasets,
and the retained parents are densely inter-adjacent (VHID 5/10, H3N2 21/28 parent–parent pairs,
§<a href="#3.4-The-discovered-dependency-structure">3.4</a>), so we cannot certify that conditioning on the co-parents blocks backdoor paths without
also conditioning on a mediator or a collider between them. We therefore estimate each parent's
effect on log2 titer (per Grantham unit) by regression on the parent plus the other parents,
with bootstrap 95% CIs, and read the numbers as **partial-regression coefficients**, not as
identified causal effects. Their shrinkage relative to the marginal (unadjusted) effect is
*consistent with* the removal of phylogenetic confounding, but equally consistent with
conditioning on a mediator or collider; the observational data cannot distinguish these.

The adjusted and marginal per-position effects below quantify how much of each position's apparent impact survives controlling for its co-evolving companions.

{{ table(src="tables/table-25.html", n=23) }}

{{ figure(src="figures/effect_sizes.png", n=24, alt="effect sizes") }}

The adjusted effects are all significant (bootstrap CI excludes zero) and
**systematically smaller in magnitude than the marginal effects** — the ×'s sit
farther from zero than the dots — a shrinkage *consistent with* the removal of
phylogenetic confounding shared among co-evolving sites — but, because the underlying
structure is GOF-rejected and the parents are densely inter-adjacent (§<a href="#3.8-Backdoor-adjusted-effect-sizes">3.8</a>), equally
consistent with conditioning on a mediator or collider, so we do not read it as confounding
identified and removed. Adjustment shrinks every effect toward zero but
never flips its sign: all 13 parents (5 VHID, 8 H3N2) keep the same sign marginally and adjusted.
Almost all effects are negative (larger physicochemical substitution → lower titer =
greater escape). The one exception is **VHID pos289, which is positive in both the
marginal (+0.042) and adjusted (+0.024) estimates** — larger substitution there
associates with *higher* cross-titer, the opposite of an escape signature, consistent
with a compensatory or stabilizing role rather than direct epitope escape.

#### Selection stability vs. adjusted effect (driver / hitchhiker separation)

We pair a Wilson 95% CI on each bootstrap selection frequency (B=200) with the shipped bootstrap
CI on the adjusted effect. (A true BCa interval would require the stored per-resample draws, which
are not shipped; this is a binomial-proportion CI, stated as such.)

**Finding.** The quadrant plot cleanly separates defensible drivers from hitchhikers.
VHID drivers (stable **and** large-effect): **158, 189, 289**; 156 is a stable *hitchhiker*
(freq=1.0 but below-median |β|); 144 is genuinely unstable (freq=0.505, Wilson 0.44–0.57).
Bedford drivers: **2, 133, 189**; 157/158/193 are stable-but-small; 190 (0.65) and 278 (0.82)
fall below the 0.9 stability line despite non-trivial effects.

> **Selection caveat (weak nulls).** The *hitchhiker* label is a claim of small/absent adjusted effect at a stably-selected position, and that negative rests on the same asymmetry: drift under-samples substitutions that change without antigenic consequence, so the data constrains the *absence* of an effect less sharply than its presence. The driver-vs-hitchhiker split should therefore be read as well supported for the large-effect drivers and more provisional for the small-effect calls.

{{ figure(src="figures/figure-13.png", n=25, alt="figure 13") }}

### 3.8.2 Left-censoring sensitivity of the adjusted effects

The HI titer has a **left-censoring floor**: pairs whose titer falls below the assay
detection limit (the "<10" undetectable code) are recorded at a floor value of 5.0
rather than a true measurement. In the two panels this affects the low tail:
VHID has 83 pairs at the 5.0 floor and 410 at 10.0 (493 pairs \u2264 10);
Bedford H3N2 has 28 at 5.0 and 588 at 10.0 (616 pairs \u2264 10). Because the model
target is `log2(HI_titer)`, these censored values enter the regression as if they
were exact, and they do not fall at random: below-detection titers concentrate in
the **high-antigenic-distance / high-Grantham-distance** pairs \u2014 the escape
regime the model is trying to explain (panel a). This raises a fair question of
whether the \u00a73.8 backdoor-adjusted per-position effects are an artifact of the
censoring code.

We test this directly. For each panel we recompute the adjusted (backdoor)
per-position effect on `log2` titer three ways: **(a)** as shipped;
**(b)** with the censored floor rows dropped; and **(c)** with the floor recoded
from 5.0 to the detection limit 10.0. The signs and ranks of the driver effects are
stable across all three (panel b): every effect keeps its sign, and the headline
drivers keep their rank. The only rank changes are single-step swaps among the
**smallest-effect** positions (VHID 144\u2194156 at ranks 4/5; Bedford H3N2
mature-2\u2194mature-278 at ranks 3/4), none of which touch the convergent headline
drivers (VHID 189, 158; Bedford 189, 133). Below-detection pairs sit
overwhelmingly in the top Grantham quartile (VHID: 52% of top-quartile pairs are
\u2264 10 vs 0% in the bottom quartile; Bedford H3N2: 17% vs 0.3%;
Mann\u2013Whitney p < 1e-100 in both), so the censoring is informative about escape
rather than noise \u2014 but it does not drive the ranking of the epitope signal.

{{ figure(src="figures/figure-14.png", n=26, alt="figure 14") }}

#### 3.8.1 From type-level to token-level: per-pair attribution and its identifiability

The effects above are **type-level**: they describe, across the whole panel, whether changing a position moves titer. A different question is **token-level** — for one specific virus–reference escape pair that differs at several positions, which difference was responsible? We address it two ways, and the two agree.

**Structural identifiability audit.** Because the pipeline collapses near-deterministic co-evolving positions into representative loci *before* discovery, every surviving driver is its own linkage-block representative: the drivers occupy **distinct** blocks (5 of 5 in VHID, 8 of 8 in H3N2), so there is no driver-versus-driver confounding — **no escape pair is block-confounded**. Among pairs that differ at ≥ 1 driver, **35–36 % are *clean single-driver* natural experiments** (655 of 1849 VHID; 1395 of 3867 H3N2): the pair differs at exactly one driver, so the titer change is attributable to a single residue from observational data alone. Position 193 (H3N2; 629 clean pairs) and 158 (VHID; 255) are the richest reservoirs of such controlled contrasts — the natural place to point a validation effort. The remaining ~64 % are *multi-block*: attributable at block resolution, with several candidate residues.

**Model-based per-pair attribution.** For each pair we computed a leave-one-difference-out counterfactual (revert one differing driver to the reference residue, re-predict, record the swing), and ran it across an ensemble of bootstrap refits and the PC/GES equivalence class to obtain a *stability band* on each position's credit. Consistent with the audit, **most drivers carry stable per-pair credit** — sign-consistent across the ensemble with a band excluding zero (H3N2 2, 133, 157, 158, 189, 190, 193, 278; VHID 144, 158, 189, 289). The informative exceptions are position-specific: **VHID 156 is a credit-swapper** (sign consistency 0.58; band −0.32 to +0.55) — the per-pair confirmation of its §<a href="#3.8-Backdoor-adjusted-effect-sizes">3.8</a> status as a *stable hitchhiker* (high selection frequency, small and non-robust effect) — as are VHID 262 and H3N2 276. Two positions (H3N2 2, VHID 276) are flagged `credit-not-decomposable`: block representatives whose credit cannot be split from their block-mates, the correct narrow statement of within-block non-identifiability. Token-level attribution is therefore **identified at residue resolution for a third of escape pairs and for most individual drivers**, at block resolution otherwise, with hitchhikers surfacing exactly where credit is unstable.

{{ figure(src="figures/figure-15.png", n=27, alt="figure 15") }}

{{ figure(src="figures/figure-16.png", n=28, alt="figure 16") }}

### 3.9 Validating the causal structure

A discovered DAG is a hypothesis, and no single test certifies it, so we run **three
complementary correctness tests** on the titer-sink graph (implemented in
`src/dag_validation.py`), each answering a different question:

1. **Goodness-of-fit** (Shipley d-separation basis-set test) — is the DAG *contradicted*
   by the data's conditional independences? Falsification only, at the equivalence-class
   level.
2. **Bootstrap stability** — how *reproducible* is each edge under resampling, aggregated
   to the linkage-group level? This is a **separate, stronger computation** from the edge
   frequencies labeling the DAG figure in §<a href="#3.4-The-discovered-dependency-structure">3.4</a>: there, the numbers are the original
   PC/GES discovery bootstrap (`bootstrap_freq`); here, we re-run an independent
   partial-correlation parent-selection pipeline on each resample and aggregate to
   linkage groups, so the two sets of frequencies are not expected to be identical.
3. **Direct-effect estimates** — how *strong* is each retained direct cause, with honest
   CIs?

The three tests are reported together below. None of them intervenes on the strain
background, so all remain observational; the definitive test is a mutagenesis
intervention at a predicted direct cause.

{{ table(src="tables/table-31.html", n=29) }}

{{ table(src="tables/table-312.html", n=30) }}

**Reading the three tests together.** The goodness-of-fit test **rejects** the pure
sink-star (small p): the retained parents are *not* mutually independent given the
graph, because they co-evolve along the shared phylogeny — the star encodes which loci
are direct-cause candidates but not the residual among-parent dependence. This is an
honest limitation, not a contradiction of the driver claims, and it is exactly why the
other two tests matter. Bootstrap stability is high at the linkage-group level (the
high-stability groups recur in ≈100% of resamples), so the driver *identities* are
reproducibly *selected* rather than resampling artifacts. This is reproducibility of selection,
not evidence of causal correctness: a hitchhiker linked to a true driver recurs in ~100% of
resamples for the wrong reason, and no test here speaks to transportable causal content.
Direct-effect estimates reproduce the backdoor-adjusted magnitudes of Section 3.8, with
honest confidence intervals.

Together the three tests establish that the discovered skeleton is *not inconsistent*
with the data at the equivalence-class level and that its high-stability drivers are
stably *selected* and quantifiable — but none of them intervenes on the strain background, so all
three remain observational evidence. The definitive test remains a mutagenesis
intervention at a predicted direct cause. Full per-dataset reports are in
`results/validation_report_*.md`.

### 3.10 Continuous per-position encoding and the titer Markov blanket

Sections 3.3–3.9 decided every edge on the **binary** per-position encoding (Hamming change
flag) with PC/GES/FCI. This section asks what changes when substitution *magnitude* is
retained instead, and reports an independent replication of Step A on VHID under a
continuous encoding. Following the companion specification `revised_workflow.md`, three
changes are made together:

1. **Encoding** — one *continuous* scalar per HA1 position (the L2 norm of the z-standardized
   12-property change vector) replaces the binary flag, retaining substitution *magnitude* while
   collapsing the 12 property axes into one variable per position. The property axes re-enter only
   downstream (the A→B bridge and the additive KAN) as *functional* evidence, never as separately
   identifiable property-level causes.
2. **Discovery target** — the primary deliverable is the titer **Markov blanket (adjacency)**, not a
   directed parent set: HI titer is a pure sink, so edges *into* it cannot be oriented from
   conditional-independence structure. Directed "parents" are reported only with that caveat.
3. **Algorithm** — linkage collapse is a *mandatory* prerequisite (it removes the fisherz singularity
   that collinear linkage blocks create), after which top-*k*-screened fisherz-FCI runs in seconds;
   PC and BOSS cross-check. RFCI/GFCI live only in Tetrad and need a Java toolchain unavailable here,
   so collapse + screened FCI is the adopted path.

The heavy step (B = 200 bootstrap of the full collapse-and-discover pipeline) is precomputed and
shipped in `results/`; the cells below re-render from those tables, matching the load-by-default
policy used for the causal bootstrap in §<a href="#3.3-Causal-discovery">3.3</a>.

#### 3.10.1 The true 12-property L2 scalar vs. the Grantham stand-in

An earlier prototype used per-position **Grantham distance** as a stand-in for the L2 scalar. Here we
compare the true 12-property scalar (built from the HA1 protein alignment) against that stand-in,
per position.

{{ figure(src="figures/figure-17.png", n=31, alt="figure 17") }}

#### 3.10.2 Full revised pipeline on VHID (n = 2751)

The complete revised Step A on VHID — continuous L2 encoding, mandatory Pearson linkage collapse
(|r|≥0.8), top-40-screened fisherz-FCI with PC/BOSS cross-checks, and a B = 200 bootstrap of the
whole pipeline. Bootstrap *adjacency* frequency (not directed-parent status) tiers the blanket.

{{ table(src="tables/table-33.html", n=32) }}

{{ figure(src="figures/figure-18.png", n=33, alt="figure 18") }}

{{ table(src="tables/table-35.html", n=34) }}

#### 3.10.3 Encoding sensitivity and per-property decomposition (exploratory)

PC titer-adjacency discovery (Fisher-z, α=0.01) was rerun under binary / Grantham / 12-property-L2
encodings. (Full Bedford PC under all three encodings was computationally intractable — the dense
skeleton did not complete — so the encoding comparison is reported for VHID; the caveat is
documented.)

**Finding (encoding).** Encoding choice materially changes selection. In VHID only **144, 156,
189, 289 survive under all three encodings** (pairwise Jaccard 0.50–0.63); five positions are
encoding-fragile — most notably **mature 158, selected only under the binary flag and dropping out
under both continuous encodings**, confirming the study's flagged fragility of 158.

**Finding (per-property) — strictly exploratory.** A single substitution sets all 12 property axes
at once, so per-property causal claims are **non-identifiable**: at every driver the property-change
vectors span only rank 1–4 of the active axes (289 is rank 1 — a single substitution type), and
partial correlations conditioning one axis on the others collapse to ~0. The marginals only *rank
which axis co-varies most* (e.g. H-bond-acceptor / β-sheet at 158, β-sheet at 189). Reported as
"which axis co-varies," never as a property-level cause.

> **Caveat (read before interpreting panel b).** The per-property marginal correlations below are
> **not** identifiable causal effects. One amino-acid substitution moves all 12 property axes
> simultaneously; the partial correlations (conditioning each axis on the others) are ≈ 0
> everywhere, so no single property can be assigned a causal role. Panel (b) ranks *co-variation*
> only.

{{ figure(src="figures/figure-19.png", n=35, alt="figure 19") }}

## 4. Conclusion

We set out to separate the HA positions that cause antigenic drift from those that merely accompany it, and the honest reading of what we achieved reframes the contribution: not a more accurate predictor of antigenic distance, but a more honest one — a model that names which residues its predictions lean on, reports how encoding-sensitive and how stable each of those choices is, and separates the part of every association that survives controlling for co-evolving sites from the part that does not.

Run on hemagglutination-inhibition titers alone, with no structural prior beyond the sink orientation, the pipeline concentrates its selected positions on the classical head antigenic sites — the residues serology and structural biology already implicate in escape — and it does so reproducibly across two independent H3N2 datasets. The positions that recur as high-stability are mature 156 and 189 in VHID and mature 133, 158, and 189 in Bedford, falling in antigenic site B (around 156–160 and 187–198, flanking the receptor-binding site) and site A (133); mature 157 and 158 are both site B, which retires an earlier "site D" label that misread alignment columns as mature residue numbers. Position 189 is a well-known determinant of H3N2 cluster transitions (Koel et al. 2013), so recovering it and its neighbors from titers alone tells us the feature matrix is tracking real antigenic biology rather than a dataset idiosyncrasy. What it does not tell us — and what we are careful not to claim — is that the causal step added identification over association: HI titer is by construction a readout of antibody binding to these very sites, so any method that extracts antigenic signal lands on the same residues, and the coincidence with known epitopes is necessary evidence that the associations are real, but not evidence that the causal search resolved cause from co-evolution.

That reservation is grounded in our own validation. The Shipley d-separation goodness-of-fit test rejects the pure titer-sink star on both datasets, because the retained positions co-evolve along a shared phylogeny and the star does not capture their residual dependence; the per-position effects of the backdoor analysis are therefore reported as partial-regression coefficients, not identified causal effects. Read that way they remain informative — every adjusted effect shrinks relative to its marginal counterpart and none flips sign, a pattern consistent with, but not identifying, the removal of phylogenetic confounding. Throughout, "stability" is a property of the selection procedure, not of causal correctness: a hitchhiker tightly linked to a true driver is selected in nearly every bootstrap for exactly the wrong reason, which is why VHID 156 — present in every strain but small and non-robust in effect — is best read as a stable hitchhiker, while 189 and the doubly-robust mature 289 carry the signal. Cluster resampling sharpens the boundary: the convergent pair 156/189 is robust, whereas the four-position VHID high-stability tier is an upper bound rather than a settled set.

Where the data cannot support residue-level attribution at all, the limiting factor is linkage. In H3N2 the largest pre-collapse blocks of co-evolving positions span 88 and 55 residues (block-size distribution, §<a href="#3.3.1-Pre-collapse-linkage-block-sizes">3.3.1</a>); these absorb many head positions into single clade-linked units that no observational method can resolve internally, which is also why H3N2 has the lowest predictive ceiling. Beyond linkage, unmeasured drivers remain possible confounders — neuraminidase, glycosylation, non-HA1 regions, and receptor-binding avidity, which raw HI titer folds in alongside head-epitope binding (Hensley et al. 2009). Because the avidity-associated residues 145, 189, and 193 overlap our drivers, escape calls at 189 and 193 are mechanistically ambiguous under a raw-titer target: our readout is raw log2 HI titer, not the cartographic antigenic distance that separates serum potency, virus avidity, and true antigenic distance (Smith et al. 2004).

It follows that the promise most naturally attached to this work — anticipating escape for strains never seen serologically — is not something we have demonstrated, and we say so directly. No available test speaks to whether the discovered structure carries transportable causal content, so transportability of the causal structure, as opposed to the associations, is currently unevaluated. The closest present-day evidence of generalization is the leakage-free grouped cross-validation, and that is the number a downstream user should hold in mind: under the strictest grouping — predicting titers against a reference antiserum held out entirely — the strongest method reaches a median R² of roughly 0.615 on VHID and 0.50 on H3N2, well below the random-split figures near 0.85 and 0.61. That is the honest predictive headline, and it is modest.

Returning to the decision that opened this study — the twice-yearly bet on which strains next season's vaccine should target — the value of an honest model is that it marks its own edges. That the selected positions rediscover the classical head epitopes from titers alone is evidence the features are tracking real biology; that we cannot yet show the selection is causal, transportable, or robust to how substitutions are scored is a map of exactly where reverse-genetics single substitutions, denser linkage-breaking serum panels, and deep-mutational-scanning escape maps must take over. Extending the framework to other subtypes and lineages would test whether the convergence generalizes beyond H3N2. Reading antigenic distance from sequence remains within reach; reading its causes from sequence is the harder step that would let surveillance anticipate escape rather than only measure it, and this study marks honestly how far the current data can carry us.

## 5. References

1. **Smith, D. J., Lapedes, A. S., de Jong, J. C., Bestebroer, T. M., Rimmelzwaan,
   G. F., Osterhaus, A. D. M. E., Fouchier, R. A. M.** (2004). Mapping the antigenic
   and genetic evolution of influenza virus. *Science* **305**(5682), 371–376.
2. **Bedford, T., Suchard, M. A., Lemey, P., Dudas, G., Gregory, V., Hay, A. J.,
   McCauley, J. W., Russell, C. A., Smith, D. J., Rambaut, A.** (2014). Integrating
   influenza antigenic dynamics with molecular evolution. *eLife* **3**, e01914.
3. **Du, E., Zhong, Z., Wang, P., et al.** (2023). DPCIPI: A pre-trained deep learning
   model for predicting cross-immunity between drifted strains of Influenza A/H3N2.
   *arXiv:2302.00926.*
4. **Grantham, R.** (1974). Amino acid difference formula to help explain protein
   evolution. *Science* **185**(4154), 862–864.
5. **Spirtes, P., Glymour, C., Scheines, R.** (2000). *Causation, Prediction, and
   Search* (2nd ed.). MIT Press. (PC algorithm.)
6. **Chickering, D. M.** (2002). Optimal structure identification with greedy search.
   *Journal of Machine Learning Research* **3**, 507–554. (GES.)
7. **Zhang, J.** (2008). On the completeness of orientation rules for causal discovery
   in the presence of latent confounders and selection bias. *Artificial Intelligence*
   **172**(16–17), 1873–1896. (FCI.)
8. **Liu, Z., Wang, Y., Vaidya, S., et al.** (2024). KAN: Kolmogorov–Arnold Networks.
   *arXiv:2404.19756.*
9. **Pearl, J.** (2009). *Causality: Models, Reasoning, and Inference* (2nd ed.).
   Cambridge University Press. (Backdoor adjustment / do-calculus.)
10. **Chen, T., Guestrin, C.** (2016). XGBoost: A scalable tree boosting system.
    *KDD 2016*, 785–794.
11. **Koel, B. F., Burke, D. F., Bestebroer, T. M., van der Vliet, S., Zondag, G. C. M.,
    Vervaet, G., Skepner, E., Lewis, N. S., Spronken, M. I. J., Russell, C. A., Eropkin,
    M. Y., Hurt, A. C., Barr, I. G., de Jong, J. C., Rimmelzwaan, G. F., Osterhaus, A. D.
    M. E., Fouchier, R. A. M., Smith, D. J.** (2013). Substitutions near the receptor
    binding site determine major antigenic change during influenza virus evolution.
    *Science* **342**(6161), 976–979.
12. **Neher, R. A., Bedford, T., Daniels, R. S., Russell, C. A., Shraiman, B. I.** (2016).
    Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza
    viruses. *Proceedings of the National Academy of Sciences* **113**(12), E1701–E1709.
13. **Łuksza, M., Lässig, M.** (2014). A predictive fitness model for influenza. *Nature*
    **507**(7490), 57–61.
14. **Harvey, W. T., Benton, D. J., Gregory, V., Hall, J. P. J., Daniels, R. S., Bedford,
    T., Haydon, D. T., Hay, A. J., McCauley, J. W., Reeve, R.** (2016). Identification of
    low- and high-impact hemagglutinin amino acid substitutions that drive antigenic drift
    of influenza A(H3N2) viruses. *PLoS Pathogens* **12**(4), e1005526.
15. **Wiley, D. C., Wilson, I. A., Skehel, J. J.** (1981). Structural identification of the
    antibody-binding sites of Hong Kong influenza haemagglutinin and their involvement in
    antigenic variation. *Nature* **289**(5796), 373–378.
16. **Caton, A. J., Brownlee, G. G., Yewdell, J. W., Gerhard, W.** (1982). The antigenic
    structure of the influenza virus A/PR/8/34 hemagglutinin (H1 subtype). *Cell* **31**(2),
    417–427.
17. **Doud, M. B., Lee, J. M., Bloom, J. D.** (2018). How single mutations affect viral
    escape from broad and narrow antibodies to H1 influenza hemagglutinin. *Nature
    Communications* **9**, 1386.
18. **Lee, J. M., Eguia, R., Zost, S. J., Choudhary, S., Wilson, P. C., Bedford, T.,
    Stevens-Ayers, T., Boeckh, M., Hurt, A. C., Lakdawala, S. S., Hensley, S. E., Bloom,
    J. D.** (2019). Mapping person-to-person variation in viral mutations that escape
    polyclonal serum targeting influenza hemagglutinin. *eLife* **8**, e49324.
19. **Lou, Y., Caruana, R., Gehrke, J., Hooker, G.** (2013). Accurate intelligible models
    with pairwise interactions. *KDD 2013*, 623–631.
20. **Hensley, S. E., Das, S. R., Bailey, A. L., Schmidt, L. M., Hickman, H. D., Jayaraman, A.,
    Viswanathan, K., Raman, R., Sasisekharan, R., Bennink, J. R., Yewdell, J. W.** (2009).
    Hemagglutinin receptor binding avidity drives influenza A virus antigenic drift.
    *Science* **326**(5930), 734–736.

---

*Data: `influenza-hi-antigenic-distance` repository (CC0). Code and this notebook are
released alongside it. Causal discovery uses the `causal-learn` library; the KAN is a
custom PyTorch implementation in `src/bspline_kan.py`. The full pipeline — linkage
collapse, causal discovery, B-spline KAN, and cross-method convergence — is packaged
as the reusable `kan-causal-antigenic-workflow` skill.*
