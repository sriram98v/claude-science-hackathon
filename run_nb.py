"""In-process notebook executor (sandbox blocks Jupyter kernel sockets).
Runs each code cell through IPython InteractiveShell with the inline backend so
figures, DataFrames, and stdout are captured as real nbformat outputs.
Usage: python run_nb.py <in.ipynb> <out.ipynb> <repo_root> [recompute]
"""
import os, sys, time

def main():
    in_nb, out_nb, repo_root = sys.argv[1], sys.argv[2], sys.argv[3]
    recompute = len(sys.argv) > 4 and sys.argv[4] == "1"
    if recompute:
        os.environ["RECOMPUTE_CAUSAL"] = "1"
        os.environ["RECOMPUTE_CV"] = "1"
    os.chdir(repo_root)
    for m in list(sys.modules):
        if m in ("analysis", "bspline_kan", "causal_helpers",
                 "dag_validation", "second_order_kan"):
            del sys.modules[m]
    import matplotlib
    matplotlib.use("module://matplotlib_inline.backend_inline")
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    import matplotlib_inline
    import nbformat
    from matplotlib.figure import Figure

    shell = InteractiveShell.instance()
    matplotlib_inline.backend_inline.set_matplotlib_formats("png")
    # Register the inline backend's figure-flush hook so bare `ax.plot(...)`
    # cells emit a PNG at cell end (enable_matplotlib is unavailable on the base
    # InteractiveShell in a plain subprocess).
    # post_run_cell passes a result arg in this IPython; the inline backend's
    # flush_figures() takes none -> wrap it in a signature-tolerant shim so the
    # hook doesn't raise TypeError into every cell's captured output.
    def _flush_figures_hook(*args, **kwargs):
        matplotlib_inline.backend_inline.flush_figures()
    shell.events.register("post_run_cell", _flush_figures_hook)

    def exec_cell(src):
        with capture_output() as cap:
            res = shell.run_cell(src, store_history=True)
        outs = []
        if cap.stdout:
            outs.append(nbformat.v4.new_output("stream", name="stdout", text=cap.stdout))
        for o in cap.outputs:
            outs.append(nbformat.v4.new_output("display_data", data=dict(o.data),
                                               metadata=dict(o.metadata or {})))
        if res.result is not None and not isinstance(res.result, Figure):
            data, meta = shell.display_formatter.format(res.result)
            if data:
                outs.append(nbformat.v4.new_output("execute_result", data=data,
                            metadata=meta, execution_count=res.execution_count))
        if cap.stderr:
            outs.append(nbformat.v4.new_output("stream", name="stderr", text=cap.stderr))
        if res.error_in_exec is not None:
            e = res.error_in_exec
            import traceback as tb
            outs.append(nbformat.v4.new_output("error", ename=type(e).__name__,
                        evalue=str(e), traceback=tb.format_exception(type(e), e, e.__traceback__)))
        return outs, res.success

    nb = nbformat.read(in_nb, as_version=4)
    n_code = sum(1 for c in nb.cells if c.cell_type == "code")
    print(f"[runner] {n_code} code cells | recompute={recompute} | repo={repo_root}", flush=True)
    t_start = time.time(); ec = 0; n_err = 0
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        ec += 1; t0 = time.time()
        first = (cell.source.splitlines() or [""])[0][:60]
        print(f"[cell {i:>2} #{ec:>2}] START            | {first}", flush=True)
        _src = "".join(cell.source)
        _src = "\n".join(l for l in _src.splitlines() if not l.strip().startswith("%matplotlib"))
        outs, ok = exec_cell(_src)
        cell.outputs = outs; cell.execution_count = ec
        # A callback exception (e.g. a bad hook) is printed by IPython into the
        # captured stream rather than as an error output; treat any traceback in
        # stream text as a failure too, so it can't pass the tally silently.
        stream_tb = any(o.output_type == "stream"
                        and ("Traceback (most recent call last)" in o.text
                             or "Error in callback" in o.text)
                        for o in outs)
        if stream_tb:
            ok = False
        dt = time.time() - t0; status = "OK " if ok else "ERR"
        if not ok: n_err += 1
        print(f"[cell {i:>2} #{ec:>2}] {status} {dt:6.1f}s | {first}", flush=True)
        if not ok:
            errs = [o for o in outs if o.output_type == "error"]
            if errs:
                print("    " + "\n    ".join(errs[0].evalue.splitlines()[:4]), flush=True)
            elif stream_tb:
                print("    [callback/stream traceback detected in cell output]", flush=True)
        nbformat.write(nb, out_nb)
    print(f"[runner] DONE {ec} cells, {n_err} errors, {time.time()-t_start:.1f}s total", flush=True)
    sys.exit(1 if n_err else 0)

if __name__ == "__main__":
    main()
