+++
title = "Article"
path = "article"
template = "article.html"
+++

This project was built as part of the [Claude Science Hackathon](https://cerebralvalley.ai/e/built-with-claude-life-sciences).

*A reproducible research notebook.* Run top-to-bottom (`Kernel → Restart & Run All`) to regenerate every table and
figure from the raw data repository shipped alongside this notebook.

## Abstract

Antibodies neutralize influenza by binding the hemagglutinin (HA) head; the virus
escapes by substituting the residues those antibodies recognize. The
hemagglutination-inhibition (HI) assay measures the resulting *antigenic distance*
between a virus and a reference strain. A central question for surveillance and
vaccine-strain selection is **which HA positions causally drive antigenic distance**,
as opposed to positions that merely *correlate* with it because they co-evolve with
the true drivers along the same phylogeny.

We analyze two HI datasets — VHID H3N2 and the Bedford H3N2 panel —
encoded as per-position feature matrices against reference strains. We (1) benchmark
how much of the HI signal is multivariable and nonlinear, with a **leakage-free grouped (leave-virus/serum-out) cross-validation** as the held-out anchor; (2) run target-oriented
**causal discovery** (PC, GES, and where tractable FCI) after collapsing near-
deterministic linkage, ranking candidates by bootstrap stability; (3) fit a genuine
**B-spline Kolmogorov–Arnold Network (KAN)** for interpretable nonlinear prediction
with per-position response curves — extended to **second order** to capture and
visualize position×position epistasis, which narrows the accuracy gap to gradient
boosting — matching it under cross-validation on VHID (0.851 vs 0.845) and equalling it on H3N2 (0.613 vs 0.613) — while remaining interpretable; and (4) reconcile the methods by
**cross-method convergence**.

The strongest claims are HA positions flagged by one causal-discovery method (bootstrap-stable at ≥0.5) and corroborated by three
correlated association screens (KAN, gradient boosting, univariate association) — the three
screens are all association measures on the same feature matrix, so their agreement is
convergent support, not four independent lines of evidence:
**positions 156 and 189 in VHID H3N2; and Bedford H3N2 alignment columns 143, 168, 199
(mature H3 residues 133, 158, 189)** — each a singleton linkage block sitting in or
adjacent to classical head antigenic sites (site A for 133, site B for 156/158/189).
Bedford columns 168 and 199 map to mature 158 and 189, coinciding exactly with the VHID
mature-158 and mature-189 drivers once the alignment-column offset is removed. Adjusted (partial-regression) effect sizes are systematically smaller than
marginal associations — a shrinkage consistent with, though over the in-sample-rejected
structure not an identification of, the phylogenetic confounding the causal framing targets. We validate the discovered structure with three complementary tests
— d-separation goodness-of-fit, bootstrap stability, and split-sample direct effects —
which together show the high-stability drivers are reproducibly *selected* across resamples and
carry honest, partial-regression effect magnitudes. "Stability" here is a property of the
selection procedure, not of causal correctness: a hitchhiker tightly linked to a true driver is
selected in ~100% of bootstraps for exactly the wrong reason, and with the out-of-distribution
transport test removed no surviving test speaks to transportable causal content — so the term is
deliberately "high-stability," never "high-confidence." All causal results are reported as
bootstrap-ranked equivalence-class candidates under a Gaussian working model, never as
a single settled graph.

## 1. Introduction

Influenza viruses evolve without pause. Under the selective pressure of population
immunity, variants that carry new substitutions in their surface glycoproteins gain a
transmission advantage, and the circulating population ratchets steadily away from the
strains our immune systems have learned to recognize — the process called **antigenic
drift**. Drift is what erodes vaccine protection: seasonal vaccines are formulated
months in advance against the strains expected to circulate, and in years when the
chosen strain no longer matches the viruses that actually spread, effectiveness falls.
Crucially, this mismatch is not an all-or-nothing quantity but a graded one — the
**antigenic distance** between the vaccine strain and a circulating strain — so
measuring that distance accurately is a prerequisite for interpreting vaccine
effectiveness and, ultimately, for designing vaccines that protect broadly against
strains not yet seen.

Measuring antigenic distance, however, is itself difficult. The reference approach,
antigenic cartography, embeds large tables of serological measurements into a
low-dimensional map whose geometry recovers the distance between strains; it has been
indispensable for understanding influenza evolution but demands serum panels and
cross-titrations that are slow and costly to generate. The measurement underlying
nearly all of these maps is the **hemagglutination-inhibition (HI) titer** — the
dilution of a reference antiserum that still blocks a virus from agglutinating red
blood cells — where a lower cross-titer signals greater antigenic distance, that is,
escape. Because serology is expensive and imperfectly reproducible between laboratories,
a long line of work has asked whether antigenic distance might instead be read directly
from the hemagglutinin (HA) *sequence*. If it could, surveillance could flag an escape
variant the moment its sequence appeared, before any assay is run.

This is precisely where prediction and explanation part ways. A sequence-based model can
forecast HI titer accurately while revealing nothing about *why* a strain has escaped,
because HA positions are strongly linked: strains share a phylogeny, so a position that
merely rides along with a true escape driver looks just as predictive as the driver
itself. Standard association tools — LASSO, gradient boosting — answer "which positions
predict titer", not "which positions *cause* the titer change"; yet it is the causal
question that matters both for mechanism and for anticipating escape variants that have
never been observed. **This study asks which individual HA positions causally drive
antigenic distance, and it defends each answer by requiring agreement between causal discovery and three
correlated association screens rather than trusting any single model — recovering, from HI
titers alone and with no structural prior beyond the direction of causation, positions
that sit squarely on the classical head epitopes.**

**Related works.** Antigenic cartography (Smith et al., 2004) revealed the punctuated
cluster structure of H3N2 drift by embedding HI tables into a low-dimensional antigenic
map, and Bedford et al. (2014) unified antigenic and genetic evolution within a joint
phylogenetic model. More recent sequence-based predictors, such as DPCIPI (Du et al.,
2023), forecast cross-immunity between drifted strains from sequence alone with strong
accuracy. These approaches are predictive or descriptive: they map, model, or forecast
antigenic distance, but none performs explicit **causal structure learning** over
individual HA positions with quantified stability — the step that distinguishes a
genuine escape driver from a co-evolving position that merely accompanies it along the
shared phylogeny. The antigenic sites those drivers fall in were originally defined
structurally and serologically — sites A–E on the H3 head by Wiley, Wilson & Skehel (1981)
and the Sa/Sb/Ca/Cb sites on H1 by Caton et al. (1982) — and refined at substitution
resolution by Koel et al. (2013) and Harvey et al. (2016) and by deep-mutational-scanning
escape maps (Doud et al., 2018; Lee et al., 2019). Our contribution is to recover the same
positions **structure-free**, from HI titers alone, by explicit causal discovery —
complementary to influenza fitness and phenotype-prediction models (Łuksza & Lässig, 2014;
Neher et al., 2016), which forecast strain frequencies or antigenic phenotypes but do not
isolate per-position causal drivers, and building on the intelligible-model literature for
pairwise interactions (Lou et al., 2013) that our second-order KAN extends.

**Our contribution.** We combine three ingredients not usually applied together to HI
data. First, **linkage collapse with target-oriented causal discovery**: we merge
near-deterministic co-evolving positions into representative loci, then learn the
direct-cause neighborhood of the HI target with PC, GES, and — where tractable — FCI,
ranking every candidate by 200× bootstrap stability. Second, **a genuine B-spline
Kolmogorov–Arnold Network (KAN)**: an interpretable nonlinear predictor whose learned
per-position response curves are directly inspectable, extended to second order to
capture and visualize position×position epistasis. Third, **cross-method convergence**:
positions flagged by causal discovery and corroborated by three correlated association
screens (the KAN, gradient boosting, and univariate association) are the defensible
causal-driver claims. Only causal discovery is methodologically distinct here — the three
screens are association measures on the same feature matrix, so they supply convergent but
not independent support — and their disagreements are themselves informative. We further estimate **backdoor-adjusted per-position effect
sizes**, using the fact that the HI target is a causal sink — though, as §3.8 shows, the
backdoor adjustment-set assumptions are rejected in-sample, so these are reported as
partial-regression coefficients, not identified causal effects.

The remainder of the notebook is organized as an executable paper. Section 2 (Methods)
introduces the two HI datasets and how their feature matrices are prepared from raw
data. Section 3 (Results) then carries the analysis in full: each subsection states how
one step is performed, presents its output, and interprets it — the predictive benchmark
and its cross-validated comparison, target-oriented causal discovery and its
dependency-structure test, the B-spline KAN and its second-order epistasis extension,
cross-method convergence, backdoor-adjusted effect sizes, and three DAG-validation tests
— moving from how much of the HI signal is learnable and how nonlinear it is, through
the discovered causal drivers and their dependency structure, to per-position effect
sizes and the validation of the discovered graph. Section 4 (Conclusion) interprets the
convergent drivers biologically, states where residue-level attribution is limited and
what assumptions the causal claims rest on, and looks ahead. Section 5 lists references.

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

The study uses two virus × reference-strain HI panels, each shipped in the data
repository as per-position HA1 feature matrices with two encodings. The **binary
mismatch** encoding sets a position to `1` when the virus and reference residues differ
at that HA1 alignment position and `0` otherwise, and is used for causal discovery. The
**Grantham** encoding gives the Grantham (1974) physicochemical distance between the two
residues (gap-aware), and is used for predictive modeling and the KAN response curves.
The final column of every matrix is the raw `HI_titer`; we model `log2(HI_titer)`. The
two study datasets map to the repository as **VHID H3N2** → `VHID/` and **Bedford H3N2**
→ `Bedford/H3/` (the repo also ships Bedford H1N1, B/Yamagata, and B/Victoria, outside
this study's scope). We begin by importing the analysis module
and fixing configuration.

**Numbering note.** Throughout, "position *N*" denotes an **HA1 alignment column**, not a
canonical mature H3 residue number. The VHID reference is gapless from mature Q1, so its
columns equal mature H3 residues (offset 0). The Bedford H3N2 alignment carries a 9-residue
signal-peptide prefix and one internal gap (column 17), so mature Q1 sits at column 10 and,
for every position reported here (all ≥ column 143), **mature residue = column − 10**
(e.g. columns 143/168/199 → mature 133/158/189). The offset is a uniform −10 only from column 18 onward; columns 10–16 carry −9 because of a single alignment gap at column 17, but no reported position falls below column 143, so the constant −10 applies throughout. This was verified against the column-wise
consensus reference (`SYILCLVFAQKLPGND-NST…`, matched by 119/122 Bedford reference strains at
the study columns) and the gapless VHID references (45/45 at offset 0). Feature spaces remain
strictly per-lineage; the mapping is provided only for biological interpretation, and the H1N1,
B/Yamagata, and B/Victoria panels (out of scope here) carry their own frames.

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
cross-validation to attach 95% confidence intervals to each method's R². **All four
methods now use the same 5×4 repeated k-fold protocol** (Round 1 cross-validated the KAN
on only 3 folds against the others' 20, making any overlap a degrees-of-freedom artifact;
that is fixed here). XGBoost selects its boosting rounds by **early stopping on an inner
validation split carved from the training fold only** — never on the scored test fold, which
in Round 1 leaked the test fold into round selection and biased XGBoost upward. The fold
scores are loaded from `results/cv_r2_folds.json` by default (recompute the linear/tree
folds with `RECOMPUTE_CV=1`; the GPU-trained KAN folds are always loaded — see `regen_cv.py`
for the KAN CV driver). Because both methods are now scored under an identical, leakage-free
protocol, we also report a **matched-fold paired test** (Wilcoxon signed-rank on the 20 per-fold
differences) rather than eyeballing interval overlap.

{{ figure(src="figures/cv_r2.png", n=4, alt="cv r2") }}

XGBoost has the highest cross-validated R² on both datasets (VHID 0.845, H3N2 0.613),
now that its test-fold early-stopping leak is removed. The KAN, cross-validated on the same
20 folds, sits just below it (VHID 0.820, H3N2 0.585). The **matched-fold paired test is
decisive rather than ambiguous**: KAN trails XGBoost by 0.025 R² on VHID (Wilcoxon p ≈ 2×10⁻⁶)
and 0.028 on H3N2 (p ≈ 6×10⁻⁶) — a small but statistically robust gap, not the interval-overlap
"tie" the Round-1 3-vs-20-fold comparison implied. The KAN's value is its additive
interpretability (§3.5), not a predictive edge; it recovers most of XGBoost's accuracy while
exposing per-position response curves XGBoost cannot.

#### 3.2.1 Leakage-free held-out R²: grouped cross-validation

Every split above partitions **pairs** (virus × reference-serum) at random, so the same virus
and the same reference antiserum recur across train and test folds — the model can memorize a
strain's antigenic profile rather than generalize to unseen strains. This is the pair-level
leakage Round 1 flagged (prior M3). To get a leakage-free estimate we re-run all four methods
under **grouped** cross-validation: `leave-virus-out` (no virus appears in both train and test)
and `leave-serum-out` (no reference antiserum shared). These are the honest analogues of the
prediction task the model is actually for — reading antigenic distance for a *new* strain.
Precomputed in `results/cv_grouped.json` (driver: `regen_cv.py`).

{{ table(src="tables/table-06.html", n=5) }}

The leakage-free numbers are **substantially lower** than the random-split values, confirming that
the random-CV R² was inflated by strain/serum recurrence. Under leave-serum-out — the strictest,
most realistic held-out task (predict titers against a *new* antiserum) — XGBoost falls from a random-split ~0.85
to a median **0.615** on VHID and from ~0.61 to **0.498** on H3N2. XGBoost remains the strongest method under every
grouping; the linear models degrade sharply and Ridge is unstable under leave-virus-out on H3N2
(a single fold with R² = −12.5, so medians rather than means are reported throughout this table). **This grouped R² — not the
random-split value — is the honest predictive headline**, and it is the number a downstream user
of the sequence-to-antigenic-distance map should expect.

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

{{ table(src="tables/table-07.html", n=6) }}

Residual strong locus pairs (|φ|≥0.9) drop to **0** in both datasets, confirming
that linkage collapse succeeded and restored the conditions discovery needs. PC and GES
then run on the collapsed loci, and the bootstrap tiers below rank each candidate.

{{ table(src="tables/table-08.html", n=7) }}

{{ table(src="tables/table-82.html", n=8) }}

The bootstrap frequency tiers each candidate into **high** (≥0.9), **moderate**
(0.5–0.9), or unstable. High-stability parents are: VHID {156, 189, 289, 158};
H3N2 alignment columns {11, 143, 167, 168, 199, 203} (mature H3 {2, 133, 157, 158, 189, 193}). Positions with `block_size > 1`
(e.g. H3N2 pos11 = 9-position block) carry claims for the whole block,
not the single position. "High-stability" scopes to reproducibility of *selection*, not causal
correctness — a hitchhiker linked to a true driver is selected in nearly every resample for the
wrong reason. The PC/GES agreement is asymmetric across datasets: in Bedford H3N2, GES independently
confirms the PC parents 143/167/168, but in VHID the two algorithms **disagree almost entirely**
— per §3.3's discovery output, VHID PC ({144, 156, 158, 189, 289}) and GES ({50, 133, 144, 145, 262, 276}) share
only pos144. The high-stability tiers rest on the *bootstrap*, not on PC/GES concordance, which is
weak for VHID.

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

> **Encoding note (see §3.10).** This skeleton uses the **binary** Hamming encoding. The revised
> *continuous* L2 encoding (§3.10.2), re-run on VHID with a B = 200 bootstrap, gives a HIGH-stability
> titer Markov blanket of {156, 189, 278, 289}. The two encodings share a HIGH core of
> **{156, 189, 289}** — pos_289 is HIGH under *both* (binary bootstrap 0.955, continuous 0.960),
> **not** a continuous-only promotion. They differ on two positions: pos_158, a HIGH site-B driver
> under the binary flag (0.95), collapses to UNSTABLE (0.490) under the continuous L2 scalar, while
> pos_278 (site C) is the one position genuinely *promoted* by the continuous encoding (binary 0.40
> → continuous 0.970). The difference is expected —
> the binary flag counts *whether* a residue changed while the L2 scalar weights *how far* it moved —
> and the two are reconciled position-by-position in §3.10, not treated as competing blankets.

{{ figure(src="figures/dag.png", n=9, alt="dag") }}

The result is consistent across both datasets: **every parent stays directly
associated with titer** conditional on the others — so the parent set contains no *pure*
intermediates, and each carries its own direct effect (matching the nonzero backdoor
estimates in §3.8). But the parents are **densely interdependent** — parent–parent adjacency is 5/10 pairs
for VHID and 21/28 for H3N2 — so the structure is *not* a strict star. This is precisely why the
goodness-of-fit test in §3.9 rejects the bare star, and the linear parent–parent
adjacencies here are the first-order shadow of the same interdependence the second-order
KAN captures nonlinearly as the epistasis surfaces in §3.6 (the 144×189 pair, for
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
(70/15/15 split, Adam, L1 spline regularization, early stopping). Training is made **near-reproducible** by pinning all RNG seeds and requesting deterministic cuDNN/cuBLAS kernels (cell above); because `use_deterministic_algorithms` runs with `warn_only=True`, a few CUDA ops fall back to non-deterministic kernels, so we claim near- rather than bit-reproducibility. In practice the headline top-15 importance ranking — consumed by the §3.7 convergence table and the §3.7.1 curvature reconciliation — is stable across reruns; only membership right at the ranking cutoff can shift by a position (see §3.7). Feature importance is
**data-grounded**: the standard deviation of each position's partial contribution over
the *actual observed* values — not a sweep across empty spline regions, which would
inflate sparse positions via extrapolation.

The synthetic-function validation and the per-dataset test R² (added to the benchmark comparison) are below.

{{ figure(src="figures/benchmark.png", n=10, alt="benchmark") }}

The KAN validates at R²≈0.99 on the synthetic additive function (recovering the true
curve shapes and flagging irrelevant inputs as flat), and on the real data it recovers most of
XGBoost's accuracy while trailing it by a small, statistically robust margin under matched-fold CV (§3.2: −0.025 R² VHID, −0.028 H3N2; both Wilcoxon p < 10⁻⁵). The KAN is preferred here for its additive interpretability, not for predictive parity. Its learned per-position
curves (shipped in `results/kan_curves.json`, and rendered in the study's curve
figures) are predominantly monotone-decreasing in Grantham distance: larger
physicochemical substitution → lower titer = greater escape.

### 3.6 Capturing epistasis: a second-order KAN

The first-order KAN is **additive** — a sum of univariate curves — and a diagnostic
(fitting depth-1 vs depth-6 gradient-boosted trees) shows that much of the accuracy gap
between additive models and full XGBoost on these datasets is **epistasis**:
position×position interactions (the per-dataset breakdown is in §3.6's ladder — predominantly
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

{{ table(src="tables/table-11.html", n=11) }}

{{ table(src="tables/table-12.html", n=12) }}

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
pre-screened by target association (§3.3, `screened=True`, 123→60 loci), so on H3N2 the
univariate screen overlaps the causal step's own input and is **not** an independent
converging family; VHID (71 loci, unscreened) does not have this dependence. Positions
flagged by the causal method **and all three screens** are treated as the headline
causal-driver claims, read in that light.

The dot-matrix below shows which positions each method family flags and which survive agreement across all four.

{{ table(src="tables/table-14.html", n=13) }}

{{ figure(src="figures/convergence.png", n=14, alt="convergence") }}

The positions on which the causal method and all three association screens **agree**
(printed above) are the strongest driver claims: **VHID
{156, 189}** and **H3N2 alignment columns {143, 168, 199}** (mature H3 {133, 158, 189}).
The high-stability agreements — 156, 189 in VHID; columns 143, 168, 199 (mature 133, 158, 189)
in H3N2 — are each a **singleton**
linkage block with bootstrap frequency ≈1.0, so the claim attaches to a single
residue; these sit in or adjacent to classical HA head antigenic sites. Positions
that enter only as large blocks (e.g. the H3N2 pos11 = 9-position block) are
block-level claims. Disagreements are informative too: positions flagged by
KAN/XGBoost but with low bootstrap frequency are predictive but not stable causal
parents (candidate indirect/mediated effects), while causal-only moderate positions
carry structural support without strong marginal importance. Exact convergence
membership right at the ranking cutoff can shift by a position between runs because the
KAN is retrained live under near- (not bit-) deterministic settings; the singleton high-stability set is stable across reruns, so the headline claims are unaffected.

#### 3.7.1 Nonlinear re-test of the predictive-but-not-causal positions

The causal step (§3.3–§3.4) decides every edge with the **Fisher-z** conditional-independence
(CI) test, which sees only *linear* dependence after linearly regressing out the conditioning
set. A dependence with no linear component is invisible to it: for $y=x^2+\varepsilon$ with $x$
symmetric about 0, $\mathrm{corr}(x,y)\approx 0$, so Fisher-z declares independence and PC drops
the edge. A genuinely **nonlinear** direct cause of titer can therefore land in the
`causal_freq ≈ 0` bucket of §3.7 while still being flagged by the predictive families (KAN,
XGBoost, univariate).

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

{{ table(src="tables/table-16.html", n=15) }}

**Reading Task 1 + Task 3 together.** The null is $H_0: P \perp \text{HI\_titer} \mid \text{parents}$.
`nonlinear_flag=True` marks the case of interest: **Fisher-z fails to reject** $H_0$ (p > 0.05, the
position looks independent to the linear test) **but KCI rejects it** (p < 0.05, dependent under the
nonparametric test) — a nonlinear direct-cause candidate the linear pipeline discarded. `flag=False`
is *not* a claim of independence: it occurs whenever that specific linear-miss pattern is absent,
which includes the common case where the position is **dependent under both tests** (Fisher-z already
rejects $H_0$, so there is nothing for KCI to "rescue"). Read the two p-value columns, not just the
flag. The `kan_spline_curvature` column is the independent corroboration from §3.5–§3.6: the
second-difference norm of the KAN's learned univariate response for that position (all other inputs at
baseline), normalized by response scale, so a straight-line response scores ≈0 and a curved one higher.
**Agreement** (KCI-dependent *and* a curved KAN spline) is a strong joint nonlinearity signal; a KCI
flag with a flat spline would suggest a KCI false positive or a KAN under-fit. Only a position that is
independent under **both** tests would be a clean predictive hitchhiker rather than a direct cause;
here every Pattern-A position remains titer-dependent under at least one test, and H3N2 column 169
(mature 159, site B) is the one whose dependence is **only** visible nonlinearly.

**Task 2 — CI-test swap sensitivity (bounded).** Task 1 re-tests a fixed candidate against a fixed
parent set. A nonlinear dependence can also change *which* conditioning set PC arrives at, so here we
attempt to re-run the whole target-parent discovery with the CI test swapped from Fisher-z to KCI and
diff the resulting parent sets. KCI-based PC is in the same cost class that made FCI intractable in
§3.3: even restricted to the **top-25 screened loci** and a **1000-row subsample**, the kernel CI test
is evaluated over too many conditioning sets to finish in a practical budget. We therefore run KCI-PC
in a worker process under a hard 180 s per-dataset wall-clock kill; if it does not
finish, we report the Fisher-z parent set and record the kill in the `status` column. This is the
documented fall-back from the handoff spec — Task 1 already answers the headline question for the
specific Pattern-A positions, which is the result the section turns on.

{{ table(src="tables/table-17.html", n=16) }}

### 3.8 Backdoor-adjusted effect sizes

Because the HI target is a causal **sink**, one might argue that for any single parent the
*other* target-parents form a valid **backdoor adjustment set**, so that regressing titer on
a parent plus those co-parents recovers that parent's causal effect. That identification
argument rests on three assumptions about the *true* graph: (i) the adjustment set contains
no collider or descendant of the treatment, (ii) it blocks every backdoor path, and (iii) the
causal structure it is read off is correct. **In-sample, these assumptions are violated.** The
d-separation goodness-of-fit test of §3.9 *rejects* the sink-star structure on both datasets,
and the retained parents are densely inter-adjacent (VHID 5/10, H3N2 21/28 parent–parent pairs,
§3.4), so we cannot certify that conditioning on the co-parents blocks backdoor paths without
also conditioning on a mediator or a collider between them. We therefore estimate each parent's
effect on log2 titer (per Grantham unit) by regression on the parent plus the other parents,
with bootstrap 95% CIs, and read the numbers as **partial-regression coefficients**, not as
identified causal effects. Their shrinkage relative to the marginal (unadjusted) effect is
*consistent with* the removal of phylogenetic confounding, but equally consistent with
conditioning on a mediator or collider; the observational data cannot distinguish these.

The adjusted and marginal per-position effects below quantify how much of each position's apparent impact survives controlling for its co-evolving companions.

{{ figure(src="figures/effect_sizes.png", n=17, alt="effect sizes") }}

The adjusted effects are all significant (bootstrap CI excludes zero) and
**systematically smaller in magnitude than the marginal effects** — the ×'s sit
farther from zero than the dots — a shrinkage *consistent with* the removal of
phylogenetic confounding shared among co-evolving sites — but, because the underlying
structure is GOF-rejected and the parents are densely inter-adjacent (§3.8), equally
consistent with conditioning on a mediator or collider, so we do not read it as confounding
identified and removed. Adjustment shrinks every effect toward zero but
never flips its sign: all 13 parents (5 VHID, 8 H3N2) keep the same sign marginally and adjusted.
Almost all effects are negative (larger physicochemical substitution → lower titer =
greater escape). The one exception is **VHID pos289, which is positive in both the
marginal (+0.042) and adjusted (+0.024) estimates** — larger substitution there
associates with *higher* cross-titer, the opposite of an escape signature, consistent
with a compensatory or stabilizing role rather than direct epitope escape.

### 3.9 Validating the causal structure

A discovered DAG is a hypothesis, and no single test certifies it, so we run **three
complementary correctness tests** on the titer-sink graph (implemented in
`src/dag_validation.py`), each answering a different question:

1. **Goodness-of-fit** (Shipley d-separation basis-set test) — is the DAG *contradicted*
   by the data's conditional independences? Falsification only, at the equivalence-class
   level.
2. **Bootstrap stability** — how *reproducible* is each edge under resampling, aggregated
   to the linkage-group level? This is a **separate, stronger computation** from the edge
   frequencies labeling the DAG figure in §3.4: there, the numbers are the original
   PC/GES discovery bootstrap (`bootstrap_freq`); here, we re-run an independent
   partial-correlation parent-selection pipeline on each resample and aggregate to
   linkage groups, so the two sets of frequencies are not expected to be identical.
3. **Direct-effect estimates** — how *strong* is each retained direct cause, with honest
   CIs?

The three tests are reported together below. None of them intervenes on the strain
background, so all remain observational; the definitive test is a mutagenesis
intervention at a predicted direct cause.

**Reading the three tests together.** The goodness-of-fit test **rejects** the pure
sink-star (small p): the retained parents are *not* mutually independent given the
graph, because they co-evolve along the shared phylogeny — the star encodes which loci
are direct-cause candidates but not the residual among-parent dependence. This is an
honest limitation, not a contradiction of the driver claims, and it is exactly why the
other two tests matter. Bootstrap stability is high at the linkage-group level (the
high-stability groups recur in ≈100% of resamples), so the driver *identities* are
reproducibly *selected* rather than resampling artifacts. This is reproducibility of selection,
not evidence of causal correctness: a hitchhiker linked to a true driver recurs in ~100% of
resamples for the wrong reason, and with the out-of-distribution transport test removed no
surviving test speaks to transportable causal content. Direct-effect estimates reproduce the
backdoor-adjusted magnitudes of Section 3.8, with honest confidence intervals.

Together the three tests establish that the discovered skeleton is *not inconsistent*
with the data at the equivalence-class level and that its high-stability drivers are
stably *selected* and quantifiable — but none of them intervenes on the strain background, so all
three remain observational evidence. The definitive test remains a mutagenesis
intervention at a predicted direct cause. Full per-dataset reports are in
`results/validation_report_*.md`.

### 3.10 Revised workflow: continuous per-position encoding and the titer Markov blanket

Sections 3.3–3.9 used the **binary** per-position encoding (Hamming change flag) with PC/GES/FCI.
A methodological revision (companion spec `revised_workflow.md`) makes three changes and re-runs
Step A on VHID as an independent replication:

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
policy used for the causal bootstrap in §3.3.

#### 3.10.1 The true 12-property L2 scalar vs. the Grantham stand-in

An earlier prototype used per-position **Grantham distance** as a stand-in for the L2 scalar. Here we
compare the true 12-property scalar (built from the HA1 protein alignment) against that stand-in,
per position.

{{ figure(src="figures/figure-07.png", n=18, alt="figure 07") }}

#### 3.10.2 Full revised pipeline on VHID (n = 2751)

The complete revised Step A on VHID — continuous L2 encoding, mandatory Pearson linkage collapse
(|r|≥0.8), top-40-screened fisherz-FCI with PC/BOSS cross-checks, and a B = 200 bootstrap of the
whole pipeline. Bootstrap *adjacency* frequency (not directed-parent status) tiers the blanket.

{{ figure(src="figures/figure-08.png", n=19, alt="figure 08") }}

{{ table(src="tables/table-24.html", n=20) }}

## 4. Conclusion

We set out to separate the HA positions that cause antigenic drift from those that merely accompany it, using methods whose outputs a biologist can inspect and a modeler can trust, and the honest reading of what we achieved is narrower than that goal but useful in its own right. Run on hemagglutination-inhibition titers alone, with no structural prior beyond the sink orientation, the pipeline concentrates its selected positions on the classical head antigenic sites — the residues serology and structural biology already implicate in escape — and it does so reproducibly across two independent datasets. In VHID H3N2 the positions that survive as high-stability under both the binary and the continuous physicochemical encoding are mature 156, 189, and 289; in the Bedford panel the convergent trio maps to mature 133, 158, and 189 (alignment columns 143, 168, 199), with 157 and 193 (columns 167, 203) close behind. In canonical H3 numbering these fall in antigenic site B (156–160 and 187–198, flanking the receptor-binding site) and site A (133), and Bedford column 168 → mature 158 and column 199 → mature 189 coincide exactly with the VHID mature-158 and mature-189 selections once the alignment-column offset is removed — a correction that also retires the earlier "site D" label, which was an artifact of reading columns 167/168 as mature residue numbers rather than the mature 157/158 they actually denote, both site B. Position 189 is a well-known determinant of H3N2 antigenic cluster transitions (Koel et al. 2013), so recovering it, and its neighbors, from titers alone is reassuring: it tells us the feature matrix is tracking real antigenic biology and not a dataset idiosyncrasy. What it does not tell us, and what we are careful not to claim, is that the causal machinery has been validated. Hemagglutination-inhibition titer is by construction a readout of antibody binding to these very sites, so any method that extracts antigenic signal — gradient boosting, univariate association, or a constraint-based causal search — will land on the same residues; the coincidence with known epitopes is therefore necessary evidence that the associations are biologically real, but it is not evidence that the causal step added identification over and above association, and we do not present it as such.

That distinction matters because our own validation constrains how strongly the selected positions can be called drivers. The Shipley d-separation goodness-of-fit test rejects the pure titer-sink star on both datasets: the retained positions are not mutually independent given the graph, because they co-evolve along a shared phylogeny, and the star encodes which loci are direct-cause candidates without capturing the residual dependence among them. For this reason the per-position effect sizes of the backdoor analysis are reported as partial-regression coefficients rather than identified causal effects. Read that way they are still informative. Every adjusted effect is smaller in magnitude than its marginal counterpart, none flips sign under adjustment, and all bootstrap confidence intervals exclude zero — a shrinkage consistent with removing phylogenetic confounding shared among co-evolving sites, though, because the structure is goodness-of-fit-rejected and the parents are densely inter-adjacent, equally consistent with conditioning on a mediator or a collider, a distinction the observational data cannot settle. The one position whose effect is positive in both the marginal and the adjusted estimate is VHID mature 289, where a larger physicochemical substitution associates with a higher cross-titer — the opposite of an escape signature — which, together with its location outside any classical antigenic site, is consistent with a compensatory or stability-restoring role rather than direct epitope escape. We offer that as a hypothesis to be tested against structural and deep-mutational-scanning evidence, not as a causal finding. The same hedging applies to two further positions that the magnitude-weighted continuous encoding promotes to high stability where the binary flag did not: mature 278, which sits in antigenic site C on the membrane-distal face near the vestigial esterase subdomain, a plausible but less-studied escape locus; and again 289. Both readings rest on a continuous-encoding adjacency that is not orientable — since titer is a pure sink, edges into it cannot be directed from conditional-independence structure — and neither is corroborated by the binary discovery, so both are hypotheses rather than claims. The encoding comparison cuts the other way as well, and we report it plainly: mature 158, a canonical site-B residue and a Koel-2013 cluster-transition position that is a high-stability driver under the binary encoding (bootstrap 0.95), drops to unstable (0.490) once substitutions are weighted by physicochemical magnitude rather than counted as present or absent. That headline residue is therefore encoding-fragile, and we present it as the clearest illustration of how much the selected set can depend on how a substitution is scored, not as a settled driver.

Where the data cannot support residue-level attribution at all, the limiting factor is linkage. In H3N2 the largest pre-collapse block of co-evolving positions spans 88 residues and a second spans 55; these absorb many head positions into single clade-linked units that no observational method can resolve internally, which is also why H3N2 has the lowest predictive ceiling of the two datasets. A claim for a position inside such a block — H3N2 position 11 is the example — belongs to the block, not the single residue, and separating them would require strains sampled across genetic backgrounds that break the linkage. Two further caveats bound the causal reading. The conditional-independence tests use a Gaussian, Fisher-z working model because the discrete test is intractable at this scale, so purely nonlinear dependences can be missed by the constraint step; the Kolmogorov–Arnold Network partially covers this gap on the predictive side, and the nonparametric kernel re-test makes the risk concrete, flagging H3N2 mature 159 (site B) as dependent on titer under the kernel test while it looked independent under Fisher-z — a nonlinear direct-cause candidate the linear pipeline would have discarded. And because PC assumes causal sufficiency and FCI was run only for VHID, unmeasured drivers — neuraminidase, glycosylation, non-HA1 regions, and the receptor-avidity and NA-activity contributions that HI titer folds in alongside head-epitope binding — remain possible confounders. All causal results are accordingly reported as bootstrap-ranked equivalence-class candidates under a Gaussian working model, never as a single settled graph, and "stability" throughout is a property of the selection procedure rather than of causal correctness: a hitchhiker tightly linked to a true driver is selected in nearly every bootstrap for exactly the wrong reason.

It also follows that the promise most naturally attached to this work — anticipating escape for strains never seen serologically — is not something we have demonstrated, and we say so directly rather than list it only as future work. The out-of-distribution transport panel was removed, and with it went the only test that spoke to whether the discovered structure carries transportable causal content; no surviving test does. The closest present-day evidence of generalization is the leakage-free grouped cross-validation, and that is also the number a downstream user should hold in mind, not the random-split values inflated by strain and serum recurrence. Under the strictest and most realistic grouping — predicting titers against a reference antiserum held out entirely — the strongest method reaches a median R² of roughly 0.61 on VHID and 0.50 on H3N2, well below the random-split figures near 0.85 and 0.61. That is the honest predictive headline, and it is modest; transportability of the causal structure, as opposed to the associations, is currently unevaluated.

The natural next steps are the ones that would convert these hedges into evidence. Mapping the convergent positions onto HA head structure, from experimental or predicted models, would let the epitope-alteration and stability hypotheses — particularly for 278 and 289 — be tested rather than asserted. Adding strains that decouple the large H3N2 linkage blocks, or folding deep-mutational-scanning escape maps in as interventional data, would move the analysis from observational discovery toward identified effects and is the only route that can resolve within-block attribution. A time-resolved or explicitly non-stationary effect model would test whether driver identities and magnitudes are stable across collection eras or drift with antigenic evolution, the natural successor now that the era-transport panel has been removed. Extending the framework to the Bedford H1N1, B/Yamagata, and B/Victoria panels would test whether the convergence generalizes beyond H3N2, and constraining the network's inputs to the selected positions would yield a predictor whose importances are causal by construction — if the causal reading survives the structural and interventional follow-up above. Returning to the decision that opened this study, the twice-yearly bet on which strains next season's vaccine should target, the contribution here is not a more accurate predictor but a more honest one: a model that names which residues its predictions lean on, reports how encoding-sensitive and how stable each of those choices is, and separates the part of each association that survives controlling for co-evolving sites from the part that does not. That the selected positions rediscover the classical head epitopes from titers alone is evidence the features are tracking real biology; that we cannot yet show the selection is causal, transportable, or robust to how substitutions are scored is a map of exactly where structural work, denser sampling, and intervention must take over. Reading antigenic distance from sequence remains the goal; reading its causes from sequence is the harder step that would let surveillance anticipate escape rather than only measure it, and this study marks honestly how far along that step the current data can carry us.

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

---

*Data: `influenza-hi-antigenic-distance` repository (CC0). Code and this notebook are
released alongside it. Causal discovery uses the `causal-learn` library; the KAN is a
custom PyTorch implementation in `src/bspline_kan.py`. The full pipeline — linkage
collapse, causal discovery, B-spline KAN, and cross-method convergence — is packaged
as the reusable `kan-causal-antigenic-workflow` skill.*
