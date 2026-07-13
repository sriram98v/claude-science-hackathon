1. The problem

    Modern ML predicts influenza antigenic distance from HA sequence accurately, but it can't say why a strain escapes — it flags positions that predict titer, not positions that causally drive it.
    The obstacle is linkage: because H3N2 strains share a dense phylogeny, a passenger ("hitchhiker") mutation looks just as predictive as the true escape driver. (This is the core framing — the "driver vs. hitchhiker" separation problem.)

2. Why it's important

    Antigenic distance sets vaccine effectiveness and drives the twice-yearly vaccine-strain decision.
    The gold-standard HI assay is slow, costly, and hard to reproduce across labs — motivating sequence-based models.
    Separating causal drivers from hitchhikers is the step needed to move from measuring escape to anticipating it.

3. Data collection

    Two independent H3N2 HI datasets — this replication-across-two-panels design is itself a strength worth stating:
        VHID (from DPCIPI), n = 2,751 pairs
        Bedford et al. 2014, n = 7,808 pairs, with collection-year metadata (enables temporal tests)
    Both rebuilt from raw pair tables into per-position HA1 feature matrices, verified by SHA-256 (reproducibility).

4. Data processing

    Two encodings: binary mismatch (for causal discovery) and Grantham physicochemical distance (for prediction/KAN).
    Linkage collapse — the load-bearing step: merge co-evolving positions (|φ| ≥ 0.8) into representative blocks. Positive result: residual strong pairs (|φ| ≥ 0.9) drop to 0, restoring the conditions causal discovery needs.
    Target-oriented causal discovery (PC, GES, FCI) with HI titer as a causal sink, ranked by 200-resample bootstrap stability.
    Fisher-z test is well-calibrated: permutation test shows near-nominal type-I error at α = 0.01 on binary, left-censored data (all encodings within Clopper–Pearson intervals) — validates the whole discovery step.

5. Main results (the headliners)

Must-include:

    The pipeline rediscovers the known antigenic sites from HI titers alone, reproducibly across both datasets — high-stability positions land in classical head sites A and B (VHID: mature 156, 189; Bedford: 133, 158, 189).
    Position 189 is the single most robust signal — survives all three resampling schemes (i.i.d., by-virus, by-serum). It's a textbook H3N2 cluster-transition determinant (Koel et al. 2013), so recovering it structure-free confirms the features track real antigenic biology.
    Interpretable prediction nearly matches the black box: the B-spline KAN trails XGBoost by only 0.025 R² (VHID) / 0.028 (H3N2) — while exposing per-position response curves XGBoost can't. Validates at R² ≈ 0.99 on synthetic data; learned curves are monotone (larger substitution → more escape).
    Clean driver/hitchhiker separation: the quadrant analysis cleanly splits stable large-effect drivers (VHID 158/189/289; Bedford 2/133/189) from stable hitchhikers (VHID 156 — near-fixed, small non-robust effect).

Optional depth if time allows:

    Epistasis captured & tested: second-order KAN plots interacting position pairs; 10 of 16 nominated pairs show a significant interaction (4/8 VHID, 6/8 Bedford).
    Backdoor-adjusted effects all significant, shrink toward zero but never flip sign (all 13 parents) — consistent with removing phylogenetic confounding.
    Cross-dataset replication ≈ 5× chance; positions 158/189 replicate (bootstrap p = 0.047).
    Token-level identifiability: 35–36% of escape pairs are clean single-driver "natural experiments," giving residue-resolution attribution for a third of pairs and most drivers.
    Robustness: driver ranking stable to left-censoring treatment and to continuous L2 encoding (shared HIGH core {156, 189, 289}).

6. Future works

    Definitive causal test: reverse-genetics / mutagenesis single substitutions at predicted drivers — start with positions 193 and 158, the richest reservoirs of clean contrasts.
    Denser, linkage-breaking serum panels and deep-mutational-scanning escape maps to resolve within-block ambiguity.
    Target cartographic antigenic distance instead of raw HI titer (separates avidity from true antigenic distance).
    Extend to other subtypes/lineages to test whether convergence generalizes beyond H3N2.
    Evaluate whether the recovered structure transports to prospective vaccine-strain selection.

