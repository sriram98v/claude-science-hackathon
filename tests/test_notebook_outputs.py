"""Gate for the blog export pipeline: verify the shipped, executed notebook actually
produced output for every cell before it is turned into blog content.

This runs against the committed executed notebook -- no kernel, no recompute, no GPU.
It is Action A's precondition: if a cell ran clean but emitted nothing (or emitted an
error), the blog would silently lose a figure/table, so we fail loudly here instead.

Reuses export_blog's own RICH / _cell_has_media / SKIP_NEVER_FOCUSED so the checks stay
in lockstep with what the exporter will actually emit.
"""
import os

import nbformat
import pytest

import export_blog as E

NOTEBOOK = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "antigenic_study.ipynb")

# The executed notebook currently ships exactly this many figure (image/png) outputs
# and this many media-producing code cells; the exporter turns them into the blog.
EXPECTED_FIGURES = 9
EXPECTED_MEDIA_CELLS = 24   # cells with >=1 rich output (before SKIP_NEVER_FOCUSED)


@pytest.fixture(scope="module")
def nb():
    return nbformat.read(NOTEBOOK, as_version=4)


def _code_cells(nb):
    return [(i, c) for i, c in enumerate(nb.cells) if c.cell_type == "code"]


def test_every_code_cell_has_output(nb):
    empty = [i for i, c in _code_cells(nb) if len(c.get("outputs", [])) == 0]
    assert not empty, f"code cells produced no output at all: {empty}"


def test_no_error_outputs(nb):
    errored = [i for i, c in _code_cells(nb)
               if any(o.get("output_type") == "error" for o in c.get("outputs", []))]
    assert not errored, f"code cells contain error outputs: {errored}"


def test_no_stderr_streams(nb):
    noisy = [i for i, c in _code_cells(nb)
             if any(o.get("output_type") == "stream" and o.get("name") == "stderr"
                    for o in c.get("outputs", []))]
    assert not noisy, f"code cells wrote to stderr (possible failure): {noisy}"


def test_expected_figure_count(nb):
    n_png = sum(
        1
        for _, c in _code_cells(nb)
        for o in c.get("outputs", [])
        if o.get("output_type") in ("execute_result", "display_data")
        and "image/png" in o.get("data", {})
    )
    assert n_png == EXPECTED_FIGURES, f"expected {EXPECTED_FIGURES} figures, found {n_png}"


def test_media_cells_carry_rich_output(nb):
    """Every cell the exporter treats as a media cell must actually hold a rich MIME."""
    media = [i for i, c in _code_cells(nb) if E._cell_has_media(c)]
    assert len(media) == EXPECTED_MEDIA_CELLS, (
        f"expected {EXPECTED_MEDIA_CELLS} media cells, found {len(media)}: {media}")
    for i in media:
        outs = nb.cells[i].get("outputs", [])
        assert any(
            m in o.get("data", {})
            for o in outs
            if o.get("output_type") in ("execute_result", "display_data")
            for m in E.RICH
        ), f"media cell {i} lost its rich output"


def test_skip_list_within_range():
    """The manual table-drop list must reference real media ordinals (1..N)."""
    assert E.SKIP_NEVER_FOCUSED, "skip list unexpectedly empty"
    assert max(E.SKIP_NEVER_FOCUSED) <= EXPECTED_MEDIA_CELLS
    assert min(E.SKIP_NEVER_FOCUSED) >= 1
