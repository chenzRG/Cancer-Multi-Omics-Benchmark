#!/usr/bin/env python3
"""
run_mlomics.py — Wrapper to run MCluster-VAEs on MLOmics data.

Builds the Hydra config overrides for a given dataset/version and launches main.py.

Usage:
    python run_mlomics.py --dataset ACC --version Top \
        --data_path /path/to/Clustering_datasets --k 2

The script:
1. Resolves filenames using the version suffix (_top / "" / _aligned).
2. Rewrites the survival CSV column names to what MCluster-VAEs expects.
3. Writes a temporary dataset YAML and calls main.py via subprocess.
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil

import pandas as pd

SUFFIX_MAP = {
    "Top":      "_top",
    "Original": "",
    "Aligned":  "_aligned",
}

OMICS = ["mRNA", "miRNA", "CNV", "Methy"]

# MCluster-VAEs CLINICAL_VARS that step_data() tries to use for non-PAN datasets
# We map survival CSV columns to these names so the loader doesn't crash.
SURVIVAL_COL_MAP = {
    "event_observed": "status",
    "survival_times": "days",
}


def build_survival_compat(src_csv: str, dst_csv: str):
    """Write a compatibility survival CSV with MCluster-VAEs column names."""
    df = pd.read_csv(src_csv, index_col=None)
    # Rename columns
    df = df.rename(columns=SURVIVAL_COL_MAP)
    # Set sample as index
    if "sample_name" in df.columns:
        df = df.set_index("sample_name")
    # Add missing CLINICAL_VARS columns as NaN so the loader doesn't fail
    clinical_vars = [
        "status", "days", "gender",
        "age_at_initial_pathologic_diagnosis",
        "pathologic_M", "pathologic_N", "pathologic_T", "pathologic_stage"
    ]
    for col in clinical_vars:
        if col not in df.columns:
            df[col] = float("nan")
    df.to_csv(dst_csv)


def main():
    parser = argparse.ArgumentParser(description="Run MCluster-VAEs on MLOmics data")
    parser.add_argument("--dataset",   required=True,
                        help="Dataset name, e.g. ACC, KIRC, KIRP, LIHC, LUAD, LUSC, PRAD, THCA, THYM")
    parser.add_argument("--version",   default="Top",
                        choices=["Top", "Original", "Aligned"],
                        help="Data version (default: Top)")
    parser.add_argument("--data_path", required=True,
                        help="Root path to Clustering_datasets/ directory")
    parser.add_argument("--k",         type=int, default=2,
                        help="Number of clusters (default: 2)")
    parser.add_argument("--device",    default="cpu",
                        help="PyTorch device (default: cpu)")
    parser.add_argument("--output_dir", default=None,
                        help="Directory to write results (optional; Hydra controls default output)")
    args = parser.parse_args()

    dataset   = args.dataset
    version   = args.version
    data_path = args.data_path
    k         = args.k
    device    = args.device

    suffix = SUFFIX_MAP[version]

    # Resolve data directory
    # Files live in {data_path}/{dataset}/{version}/ e.g. .../ACC/Top/ACC_mRNA_top.csv
    # Survival file lives in {data_path}/{dataset}/{version}/survival_{dataset}.csv
    data_dir = os.path.normpath(os.path.join(data_path, dataset, version))
    if not os.path.isdir(data_dir):
        # Fallback: try {data_path}/{dataset}/ directly (flat layout)
        data_dir_flat = os.path.normpath(os.path.join(data_path, dataset))
        if os.path.isdir(data_dir_flat):
            data_dir = data_dir_flat
        else:
            print(f"ERROR: Data directory not found: {data_dir}", file=sys.stderr)
            sys.exit(1)

    # Check omics files
    csvfn = {}
    for omics in OMICS:
        fname = f"{dataset}_{omics}{suffix}.csv"
        fpath = os.path.join(data_dir, fname)
        if not os.path.isfile(fpath):
            print(f"ERROR: Omics file not found: {fpath}", file=sys.stderr)
            sys.exit(1)
        csvfn[omics] = fname

    # Survival file
    surv_src = os.path.join(data_dir, f"survival_{dataset}.csv")
    if not os.path.isfile(surv_src):
        print(f"WARNING: Survival file not found: {surv_src} — running without clinical data",
              file=sys.stderr)
        surv_compat = None
    else:
        surv_compat = None  # will be written to tmpdir

    # Create temp directory for the compat survival CSV
    tmpdir = tempfile.mkdtemp(prefix="mcluster_vaes_")
    try:
        # Write compatibility survival CSV
        if surv_src and os.path.isfile(surv_src):
            surv_compat_path = os.path.join(tmpdir, f"survival_{dataset}_compat.csv")
            build_survival_compat(surv_src, surv_compat_path)
            surv_compat = surv_compat_path

        # Build Hydra overrides
        # We override individual keys rather than writing a full YAML
        overrides = [
            "dataset=mlomics",
            "train=mlomics",
            "model=middle",
            f"dataset.dat_name={dataset}",
            f"dataset.root={data_dir}",
            f"dataset.version={version}",
            # n_clusters is a dict in the yaml; override just this dataset's entry
            f"dataset.n_clusters.{dataset}={k}",
            f"train.device={device}",
            "dataset.normalize=[mRNA,miRNA,CNV,Methy]",
        ]

        # Override csvfn keys
        for omics, fname in csvfn.items():
            overrides.append(f"dataset.csvfn.{omics}={fname}")

        if args.output_dir:
            os.makedirs(args.output_dir, exist_ok=True)
            overrides.append(f"hydra.run.dir={args.output_dir}/mcluster_vaes")
        else:
            overrides.append("hydra.run.dir=./outputs/mcluster_vaes")

        if surv_compat:
            overrides.append(f"dataset.csvfn.cli={surv_compat}")
        else:
            # No clinical data: set cli to empty string — step_data will skip it
            # (this requires a small patch; alternatively pass a minimal CSV)
            minimal_cli = os.path.join(tmpdir, "minimal_cli.csv")
            # Read sample names from first omics file
            sample_df = pd.read_csv(
                os.path.join(data_dir, next(iter(csvfn.values()))),
                index_col=0, nrows=1
            )
            samples = list(sample_df.columns)
            clinical_vars = [
                "status", "days", "gender",
                "age_at_initial_pathologic_diagnosis",
                "pathologic_M", "pathologic_N", "pathologic_T", "pathologic_stage"
            ]
            cli_df = pd.DataFrame(
                {col: float("nan") for col in clinical_vars},
                index=pd.Index(samples, name="sample_name")
            )
            # Need all samples; read properly
            sample_df_full = pd.read_csv(
                os.path.join(data_dir, next(iter(csvfn.values()))),
                index_col=0
            )
            cli_df = pd.DataFrame(
                {col: float("nan") for col in clinical_vars},
                index=pd.Index(list(sample_df_full.columns), name="sample_name")
            )
            cli_df.to_csv(minimal_cli)
            overrides.append(f"dataset.csvfn.cli={minimal_cli}")

        # Run main.py from the MCluster-VAEs directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cmd = [sys.executable, "main.py"] + overrides
        print("Running: " + " ".join(cmd))
        result = subprocess.run(cmd, cwd=script_dir)
        sys.exit(result.returncode)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
