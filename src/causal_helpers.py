"""Kernel helpers for target-oriented causal discovery on high-linkage
binary feature matrices (e.g. sequence-mismatch vs. a continuous readout).

Workflow (see SKILL.md): build hamming matrix -> power filter -> collapse
linkage (restores faithfulness + tractability) -> discover target neighborhood
(PC/FCI with terminal-target background knowledge) -> bootstrap edge stability
-> GES score-based cross-check.
"""

import json


def build_hamming_dataset(df, seq_a_col, seq_b_col, target_col, id_cols=None):
    """Binary mismatch (hamming) representation of two aligned sequence columns.

    Returns dict with:
      processed_df : DataFrame (id_cols + target + 'hamming_string')
      H            : uint8 ndarray (n_rows x L), 1 = mismatch, 0 = match
      length       : alignment length L
    Both sequence columns must be equal-length aligned strings.
    """
    import numpy as np, pandas as pd
    a = df[seq_a_col].astype(str).values
    b = df[seq_b_col].astype(str).values
    L = len(a[0])
    if not (all(len(s) == L for s in a) and all(len(s) == L for s in b)):
        raise ValueError("aligned sequences are not all the same length")
    n = len(df)
    A = np.frombuffer("".join(a).encode("ascii"), dtype="S1").reshape(n, L)
    B = np.frombuffer("".join(b).encode("ascii"), dtype="S1").reshape(n, L)
    H = (A != B).astype(np.uint8)
    hamming_str = ["".join(map(str, row)) for row in H]
    cols = {}
    if id_cols:
        for c in id_cols:
            cols[c] = df[c].values
    cols[target_col] = df[target_col].values
    cols["hamming_string"] = hamming_str
    return {"processed_df": pd.DataFrame(cols), "H": H, "length": L}


def variable_power_filter(H, min_minor=None):
    """Drop invariant positions and positions with a minor-class count below
    `min_minor` (default 30) so CI tests have adequate power.

    Returns dict: kept_idx (0-based cols into H), X (n x kept binary float),
    names (1-based position labels 'p<pos>'), n_dropped_invariant, n_dropped_lowpower.
    """
    import numpy as np
    if min_minor is None:
        min_minor = 30
    n = H.shape[0]
    col_sum = H.sum(0)
    variable = (col_sum > 0) & (col_sum < n)
    minor = np.minimum(col_sum, n - col_sum)
    keep = variable & (minor >= min_minor)
    kept_idx = np.where(keep)[0]
    names = ["p%d" % (i + 1) for i in kept_idx]
    return {"kept_idx": kept_idx, "X": H[:, kept_idx].astype(float),
            "names": names, "n_dropped_invariant": int((col_sum == 0).sum() + (col_sum == n).sum()),
            "n_dropped_lowpower": int((variable & (minor < min_minor)).sum())}


def linkage_abs_phi(X):
    """|phi| (= |Pearson| for binary) association matrix, diagonal zeroed."""
    import numpy as np
    n = X.shape[0]
    Xc = X - X.mean(0)
    cov = Xc.T @ Xc / n
    sd = np.sqrt(np.diag(cov))
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = cov / np.outer(sd, sd)
    absphi = np.abs(phi)
    np.fill_diagonal(absphi, 0.0)
    return np.nan_to_num(absphi)


def collapse_linkage(X, names, thr=None):
    """Collapse near-deterministic co-evolving positions into representative loci.

    Positions are linked if |phi| >= thr (default 0.8); each connected component
    becomes one locus represented by its highest-power member. This is the key
    step that restores approximate FAITHFULNESS and makes PC/FCI/GES tractable
    on high-linkage data. A locus stands in for its whole block: causal claims
    attach to the block, not the representative column.

    Returns dict: locus_names, rep_idx (cols into X), block_map {locus: [1-based pos,...]},
    Xr (collapsed matrix), n_loci, residual_strong_pairs (|phi|>=0.9 among loci).
    """
    import numpy as np, scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    if thr is None:
        thr = 0.8
    absphi = linkage_abs_phi(X)
    minor = np.minimum(X.sum(0), X.shape[0] - X.sum(0))
    ncomp, lab = connected_components(sp.csr_matrix(absphi >= thr), directed=False)
    reps, block_map = [], {}
    for c in range(ncomp):
        mem = np.where(lab == c)[0]
        rep = mem[np.argmax(minor[mem])]
        reps.append(int(rep))
        block_map[names[rep]] = [int(names[j][1:]) for j in mem]
    reps = sorted(reps)
    locus_names = [names[j] for j in reps]
    Xr = X[:, reps]
    ar = linkage_abs_phi(Xr)
    iu = np.triu_indices(len(reps), 1)
    return {"locus_names": locus_names, "rep_idx": np.array(reps), "block_map": block_map,
            "Xr": Xr, "n_loci": len(reps),
            "residual_strong_pairs": int((ar[iu] >= 0.9).sum())}


def terminal_background_knowledge(colnames, target_col):
    """BackgroundKnowledge placing every feature in tier 0 and the target in
    tier 1, so all feature<->target edges orient INTO the target (target is a
    causal sink). Its parents then equal its direct causes / Markov blanket."""
    from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
    from causallearn.graph.GraphNode import GraphNode
    nodes = {c: GraphNode(c) for c in colnames}
    bk = BackgroundKnowledge()
    for c in colnames:
        bk.add_node_to_tier(nodes[c], 0 if c != target_col else 1)
    return bk


def discover_target_parents(A_df, target_col, method=None, alpha=None, fci_depth=None,
                            ci_test=None):
    """PC or FCI with terminal-target background knowledge. Returns dict:
    parents (feature cols oriented into target), graph (adjacency), cols.

    method: 'pc' (default) or 'fci'. ci_test: conditional-independence test,
    'fisherz' (default, Gaussian working model) or 'kci' (kernel-based,
    nonparametric -- detects nonlinear dependence but O(n^2)-O(n^3), so use only
    on screened/small locus sets). The discrete G-square test is intractable
    beyond ~100 nodes here. Run this on the COLLAPSED locus matrix.
    """
    import numpy as np
    if method is None:
        method = "pc"
    if alpha is None:
        alpha = 0.01
    if fci_depth is None:
        fci_depth = 3
    if ci_test is None:
        ci_test = "fisherz"
    cols = list(A_df.columns)
    data = A_df.values.astype(np.float64)
    ti = cols.index(target_col)
    bk = terminal_background_knowledge(cols, target_col)
    if method == "pc":
        from causallearn.search.ConstraintBased.PC import pc
        cg = pc(data, alpha=alpha, indep_test=ci_test, stable=True, uc_rule=0,
                uc_priority=2, background_knowledge=bk, node_names=cols, show_progress=False)
        G = cg.G.graph
    elif method == "fci":
        from causallearn.search.ConstraintBased.FCI import fci
        g, _ = fci(data, independence_test_method=ci_test, alpha=alpha, depth=fci_depth,
                   max_path_length=fci_depth, background_knowledge=bk, node_names=cols,
                   show_progress=False)
        G = g.graph
    else:
        raise ValueError("method must be 'pc' or 'fci'")
    parents = [cols[j] for j in range(len(cols)) if j != ti and G[j, ti] == -1 and G[ti, j] == 1]
    return {"parents": parents, "graph": G, "cols": cols}


def patch_ges_bic():
    """Patch causal-learn's Gaussian BIC score for a numpy>=1.25 incompatibility
    (float() on a 1x1 array raises). Call ONCE before ges(). Idempotent."""
    import numpy as np
    import causallearn.score.LocalScoreFunction as L
    import causallearn.score.LocalScoreFunctionClass as LC
    import causallearn.search.ScoreBased.GES as G

    def bic(Data, i, PAi, parameters=None):
        cov, n = Data
        lam = 0.5 if parameters is None else parameters["lambda_value"]
        sigma = cov[i, i]
        if len(PAi) > 0:
            yX = cov[np.ix_([i], PAi)]; XX = cov[np.ix_(PAi, PAi)]
            try:
                XXi = np.linalg.inv(XX)
            except np.linalg.LinAlgError:
                XXi = np.linalg.pinv(XX)
            sigma = float(np.squeeze(np.asarray(cov[i, i] - yX @ XXi @ yX.T)))
        if sigma <= 0:
            sigma = np.finfo(float).eps
        return -0.5 * n * (1 + np.log(sigma)) - lam * (len(PAi) + 1) * np.log(n)
    bic.__name__ = "local_score_BIC_from_cov"
    for mod in (L, LC, G):
        if hasattr(mod, "local_score_BIC_from_cov"):
            setattr(mod, "local_score_BIC_from_cov", bic)
    return True


def screen_top_features(A_df, target_col, k=None):
    """Rank features by |point-biserial| with the (real-valued) target and return
    the top-k names. Method-independent screen for focusing slow searches (GES)
    or the bootstrap node set."""
    import numpy as np
    from scipy import stats
    if k is None:
        k = 50
    feats = [c for c in A_df.columns if c != target_col]
    y = A_df[target_col].values.astype(float)
    X = A_df[feats].values.astype(float)
    r = np.array([abs(stats.pointbiserialr(X[:, j], y)[0]) for j in range(X.shape[1])])
    order = np.argsort(r)[::-1][:k]
    return [feats[j] for j in order]


def ges_target_neighborhood(A_df, target_col, screen_k=None, must_include=None):
    """GES (Gaussian BIC) score-based cross-check on a SCREENED node set.
    GES in causal-learn is slow beyond ~40 nodes here, so restrict to the top
    `screen_k` associated features plus any `must_include` (e.g. PC parents).
    Returns dict: neighborhood {feature: '->'|'<-'|'--'|'<->'}, cols, graph.
    """
    import numpy as np
    if screen_k is None:
        screen_k = 20
    patch_ges_bic()
    from causallearn.search.ScoreBased.GES import ges
    top = set(screen_top_features(A_df, target_col, screen_k))
    if must_include:
        top |= set(must_include)
    sub_cols = sorted(top, key=lambda c: int("".join(ch for ch in c if ch.isdigit()) or 0))
    sub_cols = sub_cols + [target_col]
    sub = A_df[sub_cols].values.astype(np.float64)
    ti = len(sub_cols) - 1
    rec = ges(sub, score_func="local_score_BIC", node_names=sub_cols)
    Gg = rec["G"].graph
    code = {(-1, 1): "->", (1, -1): "<-", (-1, -1): "--", (1, 1): "<->"}
    nb = {}
    for j in range(len(sub_cols)):
        if j == ti:
            continue
        t = code.get((Gg[j, ti], Gg[ti, j]))
        if t:
            nb[sub_cols[j]] = t
    return {"neighborhood": nb, "cols": sub_cols, "graph": Gg}


def bootstrap_target_parents(A_df, target_col, B=None, alpha=None, n_jobs=None,
                             screen_k=None, must_include=None, seed0=0):
    """Bootstrap edge-selection stability for target parents under PC + terminal
    background knowledge. Resample rows with replacement B times; report the
    fraction of resamples in which each feature is a direct parent of the target.

    Runs on a screened node set (top `screen_k` associated features + must_include)
    for speed. Returns dict: freq {feature: selection_frequency}, B_effective, node_set.
    """
    import numpy as np
    from joblib import Parallel, delayed
    if B is None:
        B = 200
    if alpha is None:
        alpha = 0.01
    if n_jobs is None:
        n_jobs = 12
    if screen_k is None:
        screen_k = 50
    node_set = set(screen_top_features(A_df, target_col, screen_k))
    if must_include:
        node_set |= set(must_include)
    cols = sorted(node_set) + [target_col]
    sub = A_df[cols].values.astype(np.float64)
    ti = len(cols) - 1
    N = sub.shape[0]

    def one(seed):
        from causallearn.search.ConstraintBased.PC import pc
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, N, N)
        dat = sub[idx] + rng.normal(0, 1e-8, (N, len(cols)))
        bk = terminal_background_knowledge(cols, target_col)
        try:
            cg = pc(dat, alpha=alpha, indep_test="fisherz", stable=True, uc_rule=0,
                    uc_priority=2, background_knowledge=bk, node_names=cols, show_progress=False)
            G = cg.G.graph
            return [cols[j] for j in range(len(cols)) if j != ti and G[j, ti] == -1 and G[ti, j] == 1]
        except Exception:
            return None
    res = Parallel(n_jobs=n_jobs, backend="loky")(delayed(one)(seed0 + s) for s in range(B))
    ok = [r for r in res if r is not None]
    freq = {}
    for r in ok:
        for p in r:
            freq[p] = freq.get(p, 0) + 1
    freq = {p: c / len(ok) for p, c in freq.items()}
    return {"freq": freq, "B_effective": len(ok), "node_set": cols[:-1]}
