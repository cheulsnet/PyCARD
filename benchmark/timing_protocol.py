#!/usr/bin/env python3
"""
PyCARD STAR Protocols timing benchmark.

Changes from v2:
- Granular visualization selection:
    visualize_pie, visualize_prop, visualize_prop_2ct, visualize_cor
- "visualization" remains available as an alias for all four visualization stages.
- CARDfree input handling now exactly mirrors CARDrun.py:
    markerList = [markerList_df[column].tolist() for column in markerList_df.columns]
- CARDfree runs only when explicitly selected with --only cardfree or enabled
  with --include-cardfree during a normal/full run.
- Selective runs prepare only the downstream objects they actually require.
- Live progress and incremental CSV saving are retained.
- PyCARD source code is not modified.

Default benchmark method: 1 warmup + 5 measured runs.
"""

import argparse
import contextlib
import csv
import io
import os
import platform
import statistics
import sys
import time
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import sys
from pathlib import Path

# Allow this script to be executed directly from benchmark/
PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_DIR = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
    
from CARD.dataload import _load_data
from CARD.utilities import create_CARDObject, create_CARDfreeObject
from CARD.CARDprop import CARD_deconvolution
from CARD.CARDrefFree import CARD_refFree
from CARD.CARDrun import _run_CARD_scMapping, _run_CARD_imputation
from CARD.visualization import (
    CARD_visualize_pie,
    CARD_visualize_prop,
    CARD_visualize_prop_2CT,
    CARD_visualize_Cor,
)

ALL_GROUPS = {
    "loading", "object", "deconvolution", "quickstart",
    "scmapping", "imputation", "visualization",
    "visualize_pie", "visualize_prop", "visualize_prop_2ct", "visualize_cor",
    "cardfree"
}

def parse_args():
    p = argparse.ArgumentParser(
        description="PyCARD STAR Protocols timing benchmark"
    )
    p.add_argument("--data-dir", default="./data")
    p.add_argument("--marker-list", default=None)
    p.add_argument("--runs", type=int, default=5)
    p.add_argument("--warmups", type=int, default=1)
    p.add_argument("--output", default=str(BENCHMARK_DIR / "timing_protocol_results.csv"))
    p.add_argument(
        "--only",
        default=None,
        help=("Comma-separated groups: loading,object,deconvolution,quickstart,"
              "scmapping,imputation,visualization,visualize_pie,visualize_prop,"
              "visualize_prop_2ct,visualize_cor,cardfree")
    )
    p.add_argument("--include-cardfree", action="store_true")
    return p.parse_args()

@contextlib.contextmanager
def quiet():
    """Suppress noisy library stdout/stderr while preserving our progress output."""
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield

def selected(group, only):
    return only is None or group in only

def load_inputs(data_dir):
    d = Path(data_dir)
    return (
        _load_data(str(d / "spatial_count.csv"), "csv"),
        _load_data(str(d / "spatial_location.csv"), "csv"),
        _load_data(str(d / "sc_count.csv"), "csv"),
        _load_data(str(d / "sc_meta.csv"), "csv"),
    )

def make_card_object(inputs):
    spatial_count, spatial_location, sc_count, sc_meta = inputs
    return create_CARDObject(
        sc_count=sc_count,
        sc_meta=sc_meta,
        spatial_count=spatial_count,
        spatial_location=spatial_location,
        ct_varname="cellType",
        ct_select=sc_meta["cellType"].unique(),
        sample_varname="sampleInfo",
        minCountGene=100,
        minCountSpot=5,
    )

def append_or_replace_csv(path, row):
    """Persist every completed/failed stage immediately."""
    path = Path(path)
    rows = []
    if path.exists():
        with path.open("r", newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
    rows = [r for r in rows if r.get("stage") != row["stage"]]
    rows.append({k: str(v) for k, v in row.items()})
    fields = [
        "stage", "status", "mean_seconds", "sd_seconds",
        "min_seconds", "max_seconds", "measured_runs",
        "warmup_runs", "error"
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def benchmark(stage, fn, warmups, runs, output, cleanup=None):
    print(f"\n{'='*78}")
    print(f"{stage}")
    print(f"{'='*78}", flush=True)

    try:
        for i in range(warmups):
            print(f"  Warmup {i+1}/{warmups} ... ", end="", flush=True)
            t0 = time.perf_counter()
            with quiet():
                fn()
            if cleanup:
                cleanup()
            dt = time.perf_counter() - t0
            print(f"done ({dt:.3f} s)", flush=True)

        times = []
        for i in range(runs):
            print(f"  Run {i+1}/{runs} ... ", end="", flush=True)
            t0 = time.perf_counter()
            with quiet():
                fn()
            if cleanup:
                cleanup()
            dt = time.perf_counter() - t0
            times.append(dt)
            print(f"done ({dt:.3f} s)", flush=True)

        mean = statistics.mean(times)
        sd = statistics.stdev(times) if len(times) > 1 else 0.0
        row = {
            "stage": stage,
            "status": "OK",
            "mean_seconds": f"{mean:.6f}",
            "sd_seconds": f"{sd:.6f}",
            "min_seconds": f"{min(times):.6f}",
            "max_seconds": f"{max(times):.6f}",
            "measured_runs": runs,
            "warmup_runs": warmups,
            "error": "",
        }
        append_or_replace_csv(output, row)
        print(f"  => Mean {mean:.3f} s | SD {sd:.3f} s | "
              f"Min {min(times):.3f} s | Max {max(times):.3f} s")
        print(f"  Saved -> {output}", flush=True)
        return row

    except Exception as e:
        if cleanup:
            try:
                cleanup()
            except Exception:
                pass
        err = f"{type(e).__name__}: {e}"
        row = {
            "stage": stage,
            "status": "ERROR",
            "mean_seconds": "",
            "sd_seconds": "",
            "min_seconds": "",
            "max_seconds": "",
            "measured_runs": 0,
            "warmup_runs": warmups,
            "error": err,
        }
        append_or_replace_csv(output, row)
        print(f"ERROR", flush=True)
        print(f"  {err}")
        print(f"  Failure recorded; continuing with the next independent stage.")
        print(f"  Saved -> {output}", flush=True)
        return row

def main():
    a = parse_args()
    data_dir = Path(a.data_dir)

    only = None
    if a.only:
        only = {x.strip().lower() for x in a.only.split(",") if x.strip()}
        bad = only - ALL_GROUPS
        if bad:
            raise SystemExit(f"Unknown --only group(s): {', '.join(sorted(bad))}")

    required = [
        data_dir / "spatial_count.csv",
        data_dir / "spatial_location.csv",
        data_dir / "sc_count.csv",
        data_dir / "sc_meta.csv",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required input file(s):\n  " + "\n  ".join(missing))

    print("=" * 78)
    print("PyCARD STAR Protocols Timing Benchmark")
    print("=" * 78)
    print(f"Python:        {platform.python_version()}")
    print(f"Platform:      {platform.platform()}")
    print(f"Data dir:      {data_dir.resolve()}")
    print(f"Warmup runs:   {a.warmups}")
    print(f"Measured runs: {a.runs}")
    print(f"Output CSV:    {Path(a.output).resolve()}")
    if only:
        print(f"Selected:      {', '.join(sorted(only))}")
    print("Results are saved immediately after every stage.", flush=True)

    # Load once for downstream stage setup. This setup load is intentionally
    # outside the timed regions for object/deconvolution/downstream benchmarks.
    print("\nPreparing shared benchmark inputs (not timed for downstream stages) ... ",
          end="", flush=True)
    with quiet():
        inputs = load_inputs(data_dir)
    print("done", flush=True)

    if selected("loading", only):
        benchmark(
            "Input data loading (4 CSV files)",
            lambda: load_inputs(data_dir),
            a.warmups, a.runs, a.output
        )

    if selected("object", only):
        benchmark(
            "Reference-based CARD object creation",
            lambda: make_card_object(inputs),
            a.warmups, a.runs, a.output
        )

    if selected("deconvolution", only):
        # Recreate the object inside each run so each phi condition starts cleanly.
        def deconv_full():
            obj = make_card_object(inputs)
            return CARD_deconvolution(
                CARD_object=obj, phi_grid=None,
                warm_start_phi=False, max_iter=200, epsilon=1e-4
            )
        def deconv_three():
            obj = make_card_object(inputs)
            return CARD_deconvolution(
                CARD_object=obj, phi_grid=[0.1, 0.5, 0.9],
                warm_start_phi=False, max_iter=200, epsilon=1e-4
            )
        def deconv_fixed():
            obj = make_card_object(inputs)
            return CARD_deconvolution(
                CARD_object=obj, phi_grid=[0.9],
                warm_start_phi=False, max_iter=200, epsilon=1e-4
            )

        benchmark("Reference CARD deconvolution - full phi search",
                  deconv_full, a.warmups, a.runs, a.output)
        benchmark("Reference CARD deconvolution - 3-phi grid",
                  deconv_three, a.warmups, a.runs, a.output)
        benchmark("Reference CARD deconvolution - fixed phi=0.9",
                  deconv_fixed, a.warmups, a.runs, a.output)

    if selected("quickstart", only):
        def quickstart():
            inp = load_inputs(data_dir)
            obj = make_card_object(inp)
            return CARD_deconvolution(
                CARD_object=obj, phi_grid=[0.9],
                warm_start_phi=False, max_iter=200, epsilon=1e-4
            )
        benchmark(
            "Standard quickstart total (load + object + fixed phi=0.9)",
            quickstart, a.warmups, a.runs, a.output
        )

    # Prepare one fitted fixed-phi object for downstream workflows.
    visual_groups = {
        "visualization", "visualize_pie", "visualize_prop",
        "visualize_prop_2ct", "visualize_cor"
    }
    need_visual = (
        only is None or
        bool(only & visual_groups)
    )
    need_fitted = (
        selected("scmapping", only) or
        selected("imputation", only) or
        need_visual
    )
    fitted = None
    if need_fitted:
        print("\nPreparing fitted fixed-phi CARD object for downstream stages "
              "(not timed) ... ", end="", flush=True)
        with quiet():
            fitted = CARD_deconvolution(
                CARD_object=make_card_object(inputs),
                phi_grid=[0.9], warm_start_phi=False,
                max_iter=200, epsilon=1e-4
            )
        print("done", flush=True)

    if selected("scmapping", only):
        benchmark(
            "scMapping (numCell=20, ncore=10)",
            lambda: _run_CARD_scMapping(
                fitted, shapeSpot="Square", numCell=20, ncore=10
            ),
            a.warmups, a.runs, a.output
        )

    if selected("imputation", only):
        benchmark(
            "Spatial resolution enhancement (num_grids=2000, in_neighbor=10)",
            lambda: _run_CARD_imputation(
                fitted, num_grids=2000, in_neighbor=10,
                exclude=None, showGrid=False
            ),
            a.warmups, a.runs, a.output
        )

    if need_visual:
        prop = fitted.uns["info_parameters"]["Proportion_CARD"]
        loc = fitted.uns["spatial_location"]
        selected_ct = list(prop.columns[:4])
        selected_ct2 = list(prop.columns[:2])

        # "visualization" means all four. Individual names permit cheap reruns.
        run_all_visual = (only is None or (only is not None and "visualization" in only))

        if run_all_visual or (only is not None and "visualize_pie" in only):
            benchmark(
                "Visualization - CARD_visualize_pie",
                lambda: CARD_visualize_pie(prop, loc),
                a.warmups, a.runs, a.output, cleanup=lambda: plt.close("all")
            )

        if run_all_visual or (only is not None and "visualize_prop" in only):
            benchmark(
                "Visualization - CARD_visualize_prop",
                lambda: CARD_visualize_prop(
                    prop, loc, ct_visualize=selected_ct
                ),
                a.warmups, a.runs, a.output, cleanup=lambda: plt.close("all")
            )

        if run_all_visual or (only is not None and "visualize_prop_2ct" in only):
            benchmark(
                "Visualization - CARD_visualize_prop_2CT",
                lambda: CARD_visualize_prop_2CT(
                    prop, loc, ct2_visualize=selected_ct2
                ),
                a.warmups, a.runs, a.output, cleanup=lambda: plt.close("all")
            )

        if run_all_visual or (only is not None and "visualize_cor" in only):
            benchmark(
                "Visualization - CARD_visualize_Cor",
                lambda: CARD_visualize_Cor(prop),
                a.warmups, a.runs, a.output, cleanup=lambda: plt.close("all")
            )

    run_cardfree = ((only is not None and "cardfree" in only) or
                    (only is None and a.include_cardfree))
    if run_cardfree:
        marker = Path(a.marker_list) if a.marker_list else data_dir / "markerList.csv"
        if not marker.exists():
            row = {
                "stage": "CARDfree setup",
                "status": "ERROR",
                "mean_seconds": "", "sd_seconds": "",
                "min_seconds": "", "max_seconds": "",
                "measured_runs": 0, "warmup_runs": a.warmups,
                "error": f"FileNotFoundError: {marker}",
            }
            append_or_replace_csv(a.output, row)
            print(f"\nCARDfree skipped: marker list not found: {marker}")
        else:
            spatial_count, spatial_location, _, _ = inputs
            marker_df = pd.read_csv(marker)
            # Exactly mirror CARD/CARDrun.py::_run_CARDfree_deconvolution().
            markerList = [
                marker_df[column].tolist()
                for column in marker_df.columns
            ]

            def make_cardfree():
                return create_CARDfreeObject(
                    markerList=markerList,
                    spatial_count=spatial_count,
                    spatial_location=spatial_location,
                    nmfSelect=4,
                    minCountGene=100,
                    minCountSpot=5,
                    marker_ct=None,
                )

            benchmark(
                "CARDfree object creation (nmfSelect=4)",
                make_cardfree, a.warmups, a.runs, a.output
            )

            def cardfree_deconv():
                obj = make_cardfree()
                return CARD_refFree(obj)

            benchmark(
                "CARDfree deconvolution (nmfSelect=4)",
                cardfree_deconv, a.warmups, a.runs, a.output
            )

    print("\n" + "=" * 78)
    print("Timing benchmark finished.")
    print(f"Results: {Path(a.output).resolve()}")
    print("=" * 78)

if __name__ == "__main__":
    main()
