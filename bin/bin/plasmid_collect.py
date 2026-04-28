#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    header_cols = [
        "#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE",
        "COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
        "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT"
    ]

    dfs = []

    for file_path in args.inputs:
        path = Path(file_path)
        if not path.exists() or path.stat().st_size == 0:
            continue

        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
            if not df.empty:
                dfs.append(df)
        except Exception:
            continue

    if dfs:
        merged = pd.concat(dfs, ignore_index=True)
        merged = merged[~merged.iloc[:, 0].astype(str).eq("#FILE")]
        merged.to_csv(args.output, sep="\t", index=False)
    else:
        pd.DataFrame(columns=header_cols).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()