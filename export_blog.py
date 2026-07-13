"""Export antigenic_study.ipynb into the portable blog contract: one markdown file
plus a media directory, and a full-notebook reference render.

Emits three things under blog/:
  content.md            all markdown prose, in order, with a `figure(...)` /
                        `table(...)` Zola shortcode at each place a figure or table
                        appears. This is the single canonical source of the blog.
  media/figures/*.png   every figure, decoded from the notebook's image/png outputs
                        (named after the notebook's own savefig filenames).
  media/tables/*.html   every kept table, as a lossless styled-HTML fragment.
  full.html             the complete executed notebook (code + all outputs), rendered
                        with nbconvert's lab template, as a reference page.

The Zola site under site/ consumes ONLY content.md + media/ -- it never reads the
notebook. So both the single-column article and the two-column reader are, by
construction, reproducible from content.md + media/ alone.

The notebook ships already executed; this script only *reads* it -- no kernel, no
recompute. Usage:  python export_blog.py
"""
import base64
import os
import re

import nbformat
from nbconvert import HTMLExporter

HERE = os.path.dirname(os.path.abspath(__file__))
NOTEBOOK = os.path.join(HERE, "antigenic_study.ipynb")
BLOG = os.path.join(HERE, "blog")
MEDIA = os.path.join(BLOG, "media")
FIG_DIR = os.path.join(MEDIA, "figures")
TAB_DIR = os.path.join(MEDIA, "tables")

RICH = ("image/png", "image/jpeg", "image/svg+xml", "text/html")

# Media cells (counted by order among media-producing cells, matching the two-column
# reader's focus scan) that never become the focused block -- each is a table adjacent
# to a figure showing the same result, so it is dropped from the curated blog. They
# remain visible in full.html. Verified by scrolling the whole reader and recording
# which blocks ever gain the .active class. Keep in sync with tests/.
SKIP_NEVER_FOCUSED = {4, 13, 18, 20, 22}


def _cell_has_media(cell):
    """True if a code cell carries at least one figure/table (rich) output."""
    if cell.cell_type != "code":
        return False
    for out in cell.get("outputs", []):
        if out.get("output_type") in ("execute_result", "display_data"):
            if any(m in out.get("data", {}) for m in RICH):
                return True
    return False


def _png_names(cell):
    """Figure filenames this cell saves, in source order (basenames of savefig args)."""
    return re.findall(r'["\']([^"\'/]+\.png)["\']', cell.source or "")


def _as_text(value):
    """nbformat rich payloads are str or list[str]; normalise to one string."""
    return "".join(value) if isinstance(value, list) else value


def _dedup_title(source):
    """Drop the duplicated H1 + bold subtitle from the notebook's first markdown cell;
    the site masthead already shows them. Keep the rest (hackathon / repro note)."""
    idx = source.find("This project was built")
    return source[idx:] if idx != -1 else source


def export_content(nb):
    """Walk cells in order, writing prose to content.md and each figure/table to
    media/. Returns (markdown_text, n_figures, n_tables, n_skipped)."""
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(TAB_DIR, exist_ok=True)

    parts = []
    n = 0                 # global media index -> data-media / anchor ids in the reader
    media_ordinal = 0     # per-media-cell counter (matches SKIP_NEVER_FOCUSED)
    n_fig = n_tab = n_skip = 0
    fig_seen = 0
    title_done = False

    for cell in nb.cells:
        if cell.cell_type == "markdown":
            src = cell.source
            if not title_done:
                title_done = True
                src = _dedup_title(src)
            parts.append(src.strip())
            continue

        if not _cell_has_media(cell):
            continue

        media_ordinal += 1
        if media_ordinal in SKIP_NEVER_FOCUSED:
            n_skip += 1
            continue

        png_names = _png_names(cell)
        png_i = 0
        tab_i = 0
        for out in cell.get("outputs", []):
            if out.get("output_type") not in ("execute_result", "display_data"):
                continue
            data = out.get("data", {})

            if "image/png" in data:
                if png_i < len(png_names):
                    fname = png_names[png_i]
                else:
                    fname = "figure-%02d.png" % (fig_seen + 1)
                png_i += 1
                fig_seen += 1
                with open(os.path.join(FIG_DIR, fname), "wb") as f:
                    f.write(base64.b64decode(_as_text(data["image/png"])))
                n += 1
                n_fig += 1
                alt = os.path.splitext(fname)[0].replace("_", " ").replace("-", " ")
                parts.append(
                    '{{ figure(src="figures/%s", n=%d, alt="%s") }}' % (fname, n, alt))

            elif "text/html" in data:
                tab_i += 1
                fname = "table-%02d.html" % (media_ordinal if tab_i == 1
                                             else int("%d%d" % (media_ordinal, tab_i)))
                with open(os.path.join(TAB_DIR, fname), "w", encoding="utf-8") as f:
                    f.write(_as_text(data["text/html"]).strip())
                n += 1
                n_tab += 1
                parts.append('{{ table(src="tables/%s", n=%d) }}' % (fname, n))

    return "\n\n".join(p for p in parts if p) + "\n", n_fig, n_tab, n_skip


def export_full(nb):
    """Complete executed notebook (code + outputs) -> blog/full.html, self-contained."""
    exporter = HTMLExporter(template_name="lab")
    exporter.embed_images = True
    html, _ = exporter.from_notebook_node(nb)
    with open(os.path.join(BLOG, "full.html"), "w", encoding="utf-8") as f:
        f.write(html)


def main():
    nb = nbformat.read(NOTEBOOK, as_version=4)
    os.makedirs(BLOG, exist_ok=True)

    md, n_fig, n_tab, n_skip = export_content(nb)
    with open(os.path.join(BLOG, "content.md"), "w", encoding="utf-8") as f:
        f.write(md)

    export_full(nb)

    # Every figure/table shortcode must resolve to a file that now exists.
    missing = []
    for src in re.findall(r'src="((?:figures|tables)/[^"]+)"', md):
        if not os.path.exists(os.path.join(MEDIA, src)):
            missing.append(src)

    print(f"[export_blog] wrote blog/content.md ({len(md):,} bytes), "
          f"{n_fig} figures, {n_tab} tables, {n_skip} tables skipped; blog/full.html")
    print(f"[export_blog] dangling media references: {len(missing)} (expect 0)"
          + (" -> " + ", ".join(missing) if missing else ""))
    if missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
