
# kernel.py for causal-dag-validation
# Step-4 validation of a titer-sink DAG for the flu-HI antigenic workflow.
# Helpers: dsep_basis_test, bootstrap_stability, effect_estimates,
# ood_validation, validation_report (+ small utilities).
# Third-party imports are deferred into function bodies (numpy/scipy/pandas
# are in the starter set; nothing else is required).

VALIDATION_VERSION = "1.0"


def parents_from_edges(edges):
    """Build a {node: [parents]} dict from an iterable of (parent, child) tuples.
    Every node that appears (as parent or child) gets a key."""
    parents = {}
    for p, c in edges:
        parents.setdefault(c, [])
        parents.setdefault(p, [])
        if p not in parents[c]:
            parents[c].append(p)
    return parents


def pcor_test(df, x, y, cond):
    """Fisher-z partial-correlation CI test of x _||_ y | cond.
    Returns (partial_r, p_value). cond is a list/iterable of column names."""
    import numpy as np
    from scipy import stats
    cond = [c for c in cond if c != x and c != y]
    n = len(df)
    xv = df[x].to_numpy(dtype=float)
    yv = df[y].to_numpy(dtype=float)
    k = len(cond)
    if k == 0:
        rx, ry = xv, yv
    else:
        Z = np.column_stack([np.ones(n)] + [df[c].to_numpy(dtype=float) for c in cond])
        bx, _, _, _ = np.linalg.lstsq(Z, xv, rcond=None)
        by, _, _, _ = np.linalg.lstsq(Z, yv, rcond=None)
        rx = xv - Z @ bx
        ry = yv - Z @ by
    sx, sy = rx.std(), ry.std()
    if sx == 0 or sy == 0:
        return 0.0, 1.0
    r = float(np.corrcoef(rx, ry)[0, 1])
    r = max(min(r, 0.9999999), -0.9999999)
    dof = n - k - 3
    if dof < 1:
        return r, float("nan")
    z = 0.5 * np.log((1 + r) / (1 - r))
    stat = np.sqrt(dof) * abs(z)
    p = 2.0 * (1.0 - stats.norm.cdf(stat))
    return r, float(p)


def dsep_basis_test(parents, data, ci_test=None, alpha=0.05):
    """Shipley d-separation basis-set goodness-of-fit test for a DAG.

    parents : {node: [parents]} dict (use parents_from_edges to build it).
    data    : pandas DataFrame with a column per node.
    ci_test : callable(df, x, y, cond) -> (stat, p); defaults to Fisher-z
              partial correlation (pcor_test).
    Returns dict with the per-claim table, Fisher's combined C, its df and
    overall p-value. Small overall p  ==>  the DAG's implied independences
    are contradicted (the DAG is falsified at the equivalence-class level)."""
    import numpy as np
    from scipy import stats
    if ci_test is None:
        ci_test = pcor_test
    nodes = list(parents.keys())

    def adjacent(a, b):
        return (a in parents.get(b, [])) or (b in parents.get(a, []))

    claims = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            a, b = nodes[i], nodes[j]
            if adjacent(a, b):
                continue
            cond = sorted((set(parents.get(a, [])) | set(parents.get(b, []))) - {a, b})
            r, p = ci_test(data, a, b, cond)
            claims.append({"x": a, "y": b, "cond": cond, "stat": r, "p": p})

    ps = [c["p"] for c in claims if c["p"] == c["p"]]  # drop nan
    k = len(ps)
    if k == 0:
        return {"claims": claims, "n_claims": 0, "C": float("nan"),
                "df": 0, "p_overall": float("nan"),
                "note": "no testable basis claims (fully connected DAG?)"}
    ps_clip = np.clip(np.array(ps), 1e-16, 1.0)
    C = float(-2.0 * np.sum(np.log(ps_clip)))
    dfree = 2 * k
    p_overall = float(1.0 - stats.chi2.cdf(C, dfree))
    return {"claims": claims, "n_claims": k, "C": C, "df": dfree,
            "p_overall": p_overall,
            "rejected_at_alpha": bool(p_overall < alpha),
            "n_individual_violations": int(sum(1 for p in ps if p < alpha))}


def bootstrap_stability(pipeline_fn, data, B=200, collapse_groups=None,
                        seed=0):
    """Bootstrap structural stability of a discovery pipeline.

    pipeline_fn : callable(df) -> iterable of (parent, child) edges.
    Re-runs the WHOLE pipeline on B row-resamples and returns per-edge
    inclusion frequencies. If collapse_groups ({feature: group_label}) is
    given, ALSO returns group-level frequencies computed per-bootstrap
    (so tightly-linked features do not dilute each other)."""
    import numpy as np
    from collections import Counter
    rng = np.random.default_rng(seed)
    n = len(data)
    edge_counts = Counter()
    group_counts = Counter()
    for _b in range(B):
        idx = rng.integers(0, n, n)
        sample = data.iloc[idx].reset_index(drop=True)
        edges = set(tuple(e) for e in pipeline_fn(sample))
        edge_counts.update(edges)
        if collapse_groups is not None:
            gedges = set()
            for (s, d) in edges:
                gedges.add((collapse_groups.get(s, s), collapse_groups.get(d, d)))
            group_counts.update(gedges)
    out = {"B": B,
           "edge_freq": {e: c / B for e, c in edge_counts.items()},
           "group_freq": None}
    if collapse_groups is not None:
        out["group_freq"] = {e: c / B for e, c in group_counts.items()}
    return out


def effect_estimates(data, target="titer", parents=None, discovery_fn=None,
                     split=True, frac=0.5, seed=0):
    """Direct-effect estimates for the parents of a SINK target, with honest
    (data-split) inference by default.

    Because target is a sink with no unobserved confounders, its parents are
    their own sufficient adjustment set, so a single OLS of target on all
    parents yields each direct effect. If discovery_fn is given and split is
    True, structure is discovered on a random `frac` of rows and effects are
    estimated on the disjoint remainder (avoids post-selection over-confidence).
    Returns per-term coef / se / t / p / 95% CI."""
    import numpy as np
    from scipy import stats
    df = data.reset_index(drop=True)
    used_split = bool(split and discovery_fn is not None)
    if used_split:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(df))
        cut = max(1, int(len(df) * frac))
        train = df.iloc[perm[:cut]]
        est = df.iloc[perm[cut:]].reset_index(drop=True)
        parents = [p for (p, c) in discovery_fn(train) if c == target]
    else:
        if parents is None:
            raise ValueError("pass parents=[...] or a discovery_fn")
        est = df
    parents = list(dict.fromkeys(parents))  # dedupe, keep order
    if len(parents) == 0:
        return {"parents": [], "estimates": [], "split_used": used_split,
                "note": "no parents of target"}
    X = np.column_stack([np.ones(len(est))] +
                        [est[p].to_numpy(dtype=float) for p in parents])
    y = est[target].to_numpy(dtype=float)
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n, kk = X.shape
    dof = n - kk
    if dof < 1:
        raise ValueError("not enough estimation rows for %d parents" % len(parents))
    sigma2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(sigma2 * XtX_inv))
    tvals = beta / se
    pvals = 2.0 * (1.0 - stats.t.cdf(np.abs(tvals), dof))
    tcrit = stats.t.ppf(0.975, dof)
    names = ["(intercept)"] + parents
    est_rows = []
    for i, nm in enumerate(names):
        est_rows.append({"term": nm, "coef": float(beta[i]), "se": float(se[i]),
                         "t": float(tvals[i]), "p": float(pvals[i]),
                         "ci_low": float(beta[i] - tcrit * se[i]),
                         "ci_high": float(beta[i] + tcrit * se[i])})
    return {"parents": parents, "n_estimation": int(n), "df": int(dof),
            "split_used": used_split, "estimates": est_rows}


def ood_validation(data, group_col, discovery_fn, predict_fn=None,
                   target="titer"):
    """Leave-one-group-out (e.g. clade/season) transport check.

    For each group: discover + fit on the other groups, predict target on the
    held-out group. Also reports parent-set recovery (Jaccard vs the parents
    discovered on the full data). discovery_fn(df)->edges; predict_fn defaults
    to OLS on the discovered parents.
    A causal parent set should transport across groups; clade-specific
    hitchhiker correlations should degrade. OOD success is necessary, not
    sufficient; a failure can also mean a genuine cross-group mechanism shift."""
    import numpy as np

    def _ols_predict(train, test, parents):
        Xtr = np.column_stack([np.ones(len(train))] +
                              [train[p].to_numpy(dtype=float) for p in parents])
        ytr = train[target].to_numpy(dtype=float)
        beta, _, _, _ = np.linalg.lstsq(Xtr, ytr, rcond=None)
        Xte = np.column_stack([np.ones(len(test))] +
                              [test[p].to_numpy(dtype=float) for p in parents])
        return Xte @ beta

    if predict_fn is None:
        predict_fn = _ols_predict
    ref_parents = set(p for (p, c) in discovery_fn(data) if c == target)
    rows = []
    for g in list(dict.fromkeys(data[group_col].tolist())):
        train = data[data[group_col] != g]
        test = data[data[group_col] == g]
        if len(test) == 0 or len(train) == 0:
            continue
        parents = [p for (p, c) in discovery_fn(train) if c == target]
        if len(parents) == 0:
            rows.append({"held_out_group": g, "n_test": int(len(test)),
                         "rmse": float("nan"), "r2": float("nan"),
                         "parent_recovery_jaccard": 0.0, "n_parents": 0})
            continue
        pred = np.asarray(predict_fn(train, test, parents), dtype=float)
        y = test[target].to_numpy(dtype=float)
        rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
        pj = set(parents)
        union = pj | ref_parents
        jac = len(pj & ref_parents) / len(union) if union else float("nan")
        rows.append({"held_out_group": g, "n_test": int(len(test)),
                     "rmse": rmse, "r2": r2,
                     "parent_recovery_jaccard": float(jac),
                     "n_parents": len(parents)})
    return {"reference_parents": sorted(str(p) for p in ref_parents),
            "per_group": rows}


def validation_report(gof=None, stability=None, effects=None, ood=None,
                      path="validation_report.md"):
    """Assemble the four test outputs into one markdown correctness report.
    Pass whichever of the four result dicts you produced; returns `path`."""
    L = ["# Causal DAG validation report", "",
         "_Step 4 of the flu-HI antigenic workflow. Tests 1-3 establish "
         "*not inconsistent* + confidence; only Test 4 tests causality._", ""]

    L.append("## 1. Goodness-of-fit (Shipley d-sep basis-set test)")
    if gof is None:
        L.append("_not run_")
    elif gof.get("n_claims", 0) == 0:
        L.append("- %s" % gof.get("note", "no basis claims"))
    else:
        verdict = "REJECTED" if gof.get("rejected_at_alpha") else "not rejected"
        L += ["- Fisher C = %.3f on df = %d, overall p = %.4g -> DAG **%s**"
              % (gof["C"], gof["df"], gof["p_overall"], verdict),
              "- basis claims: %d; individual violations (p<0.05): %d"
              % (gof["n_claims"], gof.get("n_individual_violations", 0)),
              "- _Tests the equivalence class, not arrow directions. "
              "Passing = not falsified, not proven._"]
    L.append("")

    L.append("## 2. Bootstrap structural stability")
    if stability is None:
        L.append("_not run_")
    else:
        L.append("- B = %d resamples" % stability["B"])
        gf = stability.get("group_freq") or stability.get("edge_freq", {})
        which = "group" if stability.get("group_freq") else "edge"
        top = sorted(gf.items(), key=lambda kv: -kv[1])[:20]
        L.append("- top %s inclusion frequencies:" % which)
        for e, f in top:
            L.append("  - `%s -> %s` : %.2f" % (e[0], e[1], f))
        L.append("- _Read stability at the linkage-GROUP level; raw linked "
                 "columns dilute each other._")
    L.append("")

    L.append("## 3. Direct-effect estimates (given fixed structure)")
    if effects is None:
        L.append("_not run_")
    elif not effects.get("estimates"):
        L.append("- %s" % effects.get("note", "no parents"))
    else:
        L.append("- estimation n = %d, data-split inference: %s"
                 % (effects.get("n_estimation", 0), effects.get("split_used")))
        L.append("")
        L.append("| term | coef | 95% CI | p |")
        L.append("|---|---|---|---|")
        for r in effects["estimates"]:
            L.append("| %s | %.4g | [%.4g, %.4g] | %.3g |"
                     % (r["term"], r["coef"], r["ci_low"], r["ci_high"], r["p"]))
        if not effects.get("split_used"):
            L.append("")
            L.append("- _WARNING: no data split -> CIs are post-selection "
                     "anti-conservative (too narrow)._")
    L.append("")

    L.append("## 4. Out-of-distribution / interventional validation")
    if ood is None:
        L.append("_not run (the only test of causality as such)_")
    else:
        L.append("- reference parents: %s" % ", ".join(ood["reference_parents"]))
        L.append("")
        L.append("| held-out group | n | RMSE | R2 | parent recovery (Jaccard) |")
        L.append("|---|---|---|---|---|")
        for r in ood["per_group"]:
            L.append("| %s | %d | %.3g | %.3g | %.2f |"
                     % (r["held_out_group"], r["n_test"], r["rmse"],
                        r["r2"], r["parent_recovery_jaccard"]))
        L.append("- _OOD success is necessary, not sufficient; a failure may "
                 "also mean a real cross-group mechanism shift. Gold standard "
                 "is a mutagenesis intervention at a predicted direct cause._")
    L.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(L))
    return path
