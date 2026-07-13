# Causal & Interpretable Modeling of Influenza Antigenic Distance

> **Built with Claude — Life Sciences.** This project was created with Claude Science
> as an entry in the [Built with Claude: Life Sciences](https://cerebralvalley.ai/e/built-with-claude-life-sciences)
> hackathon.

A reproducible research notebook that separates **causal** HA escape drivers from
**linked hitchhikers** across three hemagglutination-inhibition (HI) datasets, and
predicts antigenic distance with an interpretable B-spline Kolmogorov–Arnold Network.

The single deliverable is **`antigenic_study.ipynb`** — a paper-structured notebook
(Abstract → Background → Datasets → Results → Conclusion → References) that regenerates
every table and figure from raw data when run top-to-bottom.

## Layout

```
antigenic_study.ipynb              the study, as an executable paper (ships executed)
run_nb.py                          headless in-process notebook executor (see below)
pyproject.toml / uv.lock           pinned environment
src/
  analysis.py                      config + data loading + analysis routines
  bspline_kan.py                   the B-spline KAN implementation
  causal_helpers.py                linkage collapse + PC/GES/FCI + bootstrap
  dag_validation.py                d-separation / stability / OOD validation tests
  second_order_kan.py              pairwise-interaction (epistasis) KAN
results/                           precomputed result artifacts (see note below)
figures/                          (written by the notebook on run)
export_blog.py                     notebook → portable blog contract (see below)
assemble_site.py                   blog/ → Zola site inputs (md + media only)
blog/                              the portable blog contract (committed)
  content.md                       all prose + figure/table shortcodes (single source)
  media/figures/*.png              every figure, extracted from the notebook
  media/tables/*.html              every kept table, as a styled-HTML fragment
  full.html                        full executed notebook render, code included
site/                              Zola project (builds from blog/ alone)
  config.toml / templates/ / static/
.github/workflows/export-blog.yml  CI A: test gate → export blog/ → commit
.github/workflows/deploy-web.yml   CI B: build Zola site → publish to `web` branch
influenza-hi-antigenic-distance/   the raw data repository (download via download_data.sh; gitignored)
```

**The dataset directory sits at the repository root, alongside the notebook.** The raw
HI dataset is not tracked in git — download it from Zenodo before running:

```bash
./download_data.sh          # fetches + verifies + extracts to influenza-hi-antigenic-distance/
```

The script pulls the archive from Zenodo (<https://zenodo.org/records/21339272>),
verifies its MD5 checksum, and extracts `influenza-hi-antigenic-distance/` into the repo
root. To place it manually instead, download and unpack the tarball there yourself (or
symlink an existing copy). The notebook then regenerates the feature matrices from the raw
cleaned pair tables via the data repo's own build scripts — the reproduction is
byte-for-byte identical to the shipped matrices (verified by SHA-256).

## Running

```bash
uv sync                    # create the pinned environment
uv run jupyter lab         # open antigenic_study.ipynb, Run All
```

or execute headless:

```bash
uv run jupyter nbconvert --to notebook --execute --inplace antigenic_study.ipynb
```

### The one slow step

The causal **bootstrap** on the H3N2 60-locus graph takes ~30 min (200 resamples). By
default the notebook loads the shipped `results/causal_results.json` for the
bootstrap-stability tiers, so a full run completes in a few minutes. To recompute
causal discovery from scratch:

```bash
RECOMPUTE_CAUSAL=1 uv run jupyter nbconvert --to notebook --execute --inplace antigenic_study.ipynb
```

Everything else — matrix reproduction, predictive benchmark, cross-validated R²,
KAN training, cross-method convergence, backdoor effect sizes — is recomputed live.
A CUDA GPU is used for the KAN if available (falls back to CPU).

To recompute **everything** from raw data — both the causal bootstrap and the heavy
5×4 repeated-k-fold cross-validation — set both switches (this is how the shipped
executed notebook was produced, on an RTX 3090, ~107 min total):

```bash
RECOMPUTE_CAUSAL=1 RECOMPUTE_CV=1 \
  uv run jupyter nbconvert --to notebook --execute --inplace antigenic_study.ipynb
```

### Headless execution without a Jupyter kernel

`run_nb.py` executes the notebook **in-process** via IPython — no out-of-process
Jupyter kernel, no ZeroMQ sockets. Use it in restricted/sandboxed environments where
`nbconvert --execute` cannot bind kernel sockets, or for scripted re-execution with
per-cell timing:

```bash
# python run_nb.py <in.ipynb> <out.ipynb> <repo_root> [recompute]
uv run python run_nb.py antigenic_study.ipynb antigenic_study.ipynb . 1
```

It captures stdout, inline figures, and rich DataFrame/HTML outputs as real nbformat
cells, checkpoints after every cell, prints per-cell START/OK timing, and exits
non-zero if any cell errors. Passing `1` as the last argument sets both recompute
switches; omit it (or pass `0`) for the fast cached path.

## Publishing to GitHub Pages (blog)

The blog is built in **two stages** from a small, portable contract rather than straight
from the notebook. The notebook ships **with its output cells**, so nothing here needs a
GPU, a recompute, or the project env.

**The contract** — everything the blog needs, committed and human-inspectable:

- `blog/content.md` — one markdown file with all prose, and a `figure(...)` / `table(...)`
  shortcode wherever a figure or table appears. The single source of the blog.
- `blog/media/figures/*.png`, `blog/media/tables/*.html` — every figure and kept table.
- `blog/full.html` — the complete executed notebook, code included, linked from the blog.

Both published pages — the two-column synchronized **reader** and the single-column
**article** (print/PDF friendly) — are generated by [Zola](https://www.getzola.org/) from
`blog/content.md` + `blog/media/` **alone**; `site/` never reads the notebook.

Rebuild locally:

```bash
uv run --group dev python -m pytest tests/   # gate: every notebook cell has output
uv run --group dev python export_blog.py     # notebook → blog/content.md + media/ + full.html
python assemble_site.py                       # blog/ → site/content + site/static
cd site && zola build && zola serve           # preview reader (/) and article (/article/)
```

**Stage A — `export-blog.yml`:** on a notebook change, runs the test gate, regenerates
`blog/` from the notebook, commits it to `main`, then dispatches Stage B.

**Stage B — `deploy-web.yml`:** on a `blog/` or `site/` change, runs `assemble_site.py`,
`zola build`, and publishes `site/public` to the **`web`** branch.

One-time GitHub setup:

1. Push this repo to `github.com/sriram98v/claude-science-hackathon` (the repo **name**
   determines the URL path).
2. **Settings → Pages → Build and deployment → Source → Deploy from a branch → `web` /
   `(root)`.** (The `web` branch is created by the first `deploy-web` run.)
3. Push to `main`. Editing the notebook triggers Stage A → Stage B; editing `blog/` or
   `site/` triggers Stage B directly.

The site is then live at **`https://sriram98v.github.io/claude-science-hackathon/`** (the
two-column reader), with the single-column article at **`/article/`**.

## Method summary

1. **Linkage collapse** merges positions with |φ|≥0.8 into representative loci,
   restoring approximate faithfulness and making causal search tractable.
2. **Target-oriented causal discovery** (PC, GES; FCI where tractable) with the HI
   target as a causal sink; candidates ranked by **200× bootstrap stability**.
3. **B-spline KAN** for interpretable nonlinear prediction with per-position learned
   response curves and data-grounded importance.
4. **Cross-method convergence** — positions flagged by all four method families
   (causal, KAN, gradient boosting, univariate) are the headline driver claims.
5. **Backdoor-adjusted effect sizes** — the other target-parents form a valid
   adjustment set (target is a sink), giving confounding-corrected per-position effects.

The full pipeline is also packaged as the reusable `kan-causal-antigenic-workflow`
skill.

## License

Code released under the MIT License (see [LICENSE](LICENSE)). The raw HI dataset is
archived on Zenodo (<https://zenodo.org/records/21339272>) under CC-BY-4.0.
