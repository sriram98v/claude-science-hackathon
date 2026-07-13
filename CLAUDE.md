# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A reproducible research paper delivered as a single executable notebook, **`antigenic_study.ipynb`**. It separates *causal* HA escape drivers from *linked hitchhikers* across three hemagglutination-inhibition (HI) datasets and predicts antigenic distance with an interpretable B-spline Kolmogorov–Arnold Network (KAN). The notebook is paper-structured (Abstract → Background → Datasets → Results → Conclusion → References) and regenerates every table and figure from raw data when run top-to-bottom. It ships *with* its output cells so the blog (a Zola site built from an exported markdown+media contract) can be published without a GPU.

## Environment & commands

Environment is pinned via `uv` (Python 3.11, `pyproject.toml` + `uv.lock`). `src/` is on the path (hatch wheel packages `src`).

```bash
uv sync                    # create pinned env
uv run jupyter lab         # open antigenic_study.ipynb, Run All
```

Headless execution:

```bash
# Standard Jupyter path
uv run jupyter nbconvert --to notebook --execute --inplace antigenic_study.ipynb

# In-process (no Jupyter kernel sockets) — use in sandboxed/restricted envs
# python run_nb.py <in.ipynb> <out.ipynb> <repo_root> [recompute]
uv run python run_nb.py antigenic_study.ipynb antigenic_study.ipynb . 1
```

`run_nb.py` runs cells through an IPython `InteractiveShell` with the inline backend, captures figures/DataFrames/stdout as real nbformat outputs, checkpoints after every cell, prints per-cell timing, and exits non-zero on any cell error. Its trailing `1` sets both recompute switches; `0` (or omitted) uses the fast cached path.

## Recompute switches (important)

Two env vars gate the slow steps. Default (both off) loads shipped artifacts from `results/` so a run finishes in minutes; everything else is recomputed live.

- `RECOMPUTE_CAUSAL=1` — recompute causal discovery + the ~30 min H3N2 60-locus bootstrap (200 resamples) instead of loading `results/causal_results.json`.
- `RECOMPUTE_CV=1` — recompute the heavy 5×4 repeated k-fold cross-validation.

Setting both reproduces the shipped executed notebook (~107 min on an RTX 3090). KAN training uses CUDA if available, else CPU.

## Data dependency

The raw data repo **`influenza-hi-antigenic-distance/`** must sit at the repository root (place or symlink it there before running). The notebook rebuilds feature matrices from the raw cleaned pair tables via that repo's own build scripts; reproduction is verified byte-for-byte by SHA-256. `src/analysis.py` hardcodes `DATA_DIR`, `DATASETS = ["vhid_HA1", "H3N2", "H1N1"]`, and per-dataset paths — this is the integration point with the data repo.

## Architecture

The notebook is a thin orchestration layer; all logic lives in `src/`. Pipeline stages (see README "Method summary"): linkage collapse → target-oriented causal discovery → B-spline KAN → cross-method convergence → backdoor-adjusted effect sizes.

- **`src/analysis.py`** — the hub. Global config/constants (`SEED`, `MIN_MINOR`, `COLLAPSE_THR=0.8`, `BOOTSTRAP_B=200`, confidence tiers), data loading (`load_all`, `rebuild_matrices_from_raw`), and per-stage routines (`benchmark`, `collapse_to_loci`, `tiered_parents`, `convergence_from_tops`, `sokan_cv`, `backdoor_effects`, `draw_titer_dag`, `build_validation_frame`). `load_result` reads cached `results/*.json`. Edit config here.
- **`src/causal_helpers.py`** — linkage collapse (`collapse_linkage`, |φ|≥`COLLAPSE_THR`) and target-oriented causal discovery: PC/GES/FCI wrappers with the HI target pinned as a causal *sink* (`terminal_background_knowledge`, `discover_target_parents`), plus `bootstrap_target_parents` for the 200× stability ranking.
- **`src/bspline_kan.py`** — first-order interpretable KAN (`KANLinear`, `BSplineKAN`) with per-position learned response curves.
- **`src/second_order_kan.py`** — pairwise-interaction (epistasis) KAN (`SecondOrderKAN`, `train_sokan`, `top_interactions`, `interaction_surface`).
- **`src/dag_validation.py`** — structural validation: d-separation basis tests (`dsep_basis_test`), `bootstrap_stability`, `effect_estimates`, `ood_validation`, and `validation_report` (writes `results/validation_report_*.md`).

Key modeling premise that shapes the code: because the HI target is a graph sink, the other target-parents form a valid backdoor adjustment set — this is why `backdoor_effects` can produce confounding-corrected per-position effects.

## Outputs

- `results/` — cached artifacts (`causal_results.json`, `kan_results.json`, `sokan_results.json`, `*_univariate_fdr.csv`, `kan_model_*.pt`, validation reports). These are the fast-path inputs; regenerated only under the recompute switches.
- `figures/` — written by the notebook on each run.
- **Blog pipeline (two stages, notebook-independent after export).** `export_blog.py` turns the executed notebook into a portable contract in `blog/`: `content.md` (all prose + `figure(...)`/`table(...)` shortcodes), `media/figures/*.png`, `media/tables/*.html`, and `full.html` (complete notebook render). `assemble_site.py` materializes that contract into the Zola project in `site/` (never reads the notebook); `zola build` produces the two-column reader (`/`) and single-column article (`/article/`). CI: `.github/workflows/export-blog.yml` (test gate → export → commit `blog/`) then `.github/workflows/deploy-web.yml` (zola build → publish `site/public` to the **`web`** branch, which is the Pages source). Local: `python export_blog.py && python assemble_site.py && (cd site && zola build)`. The blog's curated media set is controlled by `SKIP_NEVER_FOCUSED` in `export_blog.py`; the notebook-cell-output gate is `tests/test_notebook_outputs.py`.

## Conventions

Determinism matters: a fixed `SEED = 0` is threaded through benchmarks, CV, and bootstraps — preserve seeding when adding routines. Power filtering drops positions with minor-class count `< MIN_MINOR` (30). Code is licensed under MIT (see `LICENSE`); the raw HI dataset is CC-BY-4.0 (archived on Zenodo).
