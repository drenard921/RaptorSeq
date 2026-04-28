# import pandas as pd
# df = pd.read_csv("/Shared/SHL-BUG/workspaces/drenard/nfcore_project/nf-core-raptorseq/results/amrfinder/amrfinder_all.tsv", sep="\t")[["Name","Element symbol","% Identity to reference"]]
# df["% Identity to reference"] = pd.to_numeric(df["% Identity to reference"].astype(str).replace(",", ".", regex=False), errors="coerce")
# mat = df.pivot_table(index="Name", columns="Element symbol", values="% Identity to reference", aggfunc="max", fill_value=0.0)
# mat.sort_index().reindex(sorted(mat.columns), axis=1).to_csv("amrfinder_heatmap_test.tsv", sep="\t")

# # import pandas as pd

# # df = pd.read_csv("/Shared/SHL-BUG/workspaces/drenard/nfcore_project/nf-core-raptorseq/results/amrfinder/amrfinder_all.tsv", sep="\t")[["Name", "Element symbol", "% Identity to reference"]]
# # df["% Identity to reference"] = pd.to_numeric(
# #     df["% Identity to reference"].astype(str).str.replace(",", ".", regex=False),
# #     errors="coerce"
# # )

# # mat = df.pivot_table(
# #     index="Element symbol",
# #     columns="Name",
# #     values="% Identity to reference",
# #     aggfunc="max",
# #     fill_value=0.0
# # )

# # mat = mat.sort_index().reindex(sorted(mat.columns), axis=1)
# # mat.to_csv("amrfinder_heatmap.tsv", sep="\t")

#!/usr/bin/env python3

import argparse
import pandas as pd


def build_amrfinder_heatmap(input_tsv: str, output_tsv: str) -> None:
    df = pd.read_csv(input_tsv, sep="\t")[["Name", "Element symbol", "% Identity to reference"]]

    df["% Identity to reference"] = pd.to_numeric(
        df["% Identity to reference"].astype(str).str.replace(",", ".", regex=False),
        errors="coerce"
    )

    mat = df.pivot_table(
        index="Name",
        columns="Element symbol",
        values="% Identity to reference",
        aggfunc="max",
        fill_value=0.0
    )

    mat = mat.sort_index().reindex(sorted(mat.columns), axis=1)
    mat.to_csv(output_tsv, sep="\t")


def parse_args():
    parser = argparse.ArgumentParser(description="Build AMRFinder heatmap matrix TSV")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to AMRFinder combined TSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output heatmap TSV"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    build_amrfinder_heatmap(args.input, args.output)
    print(f"[INFO] Wrote heatmap TSV to {args.output}")


# Example invocation
# python3 bin/amr_python.py \
#   --input /Shared/SHL-BUG/workspaces/drenard/nfcore_project/nf-core-raptorseq/results/amrfinder/amrfinder_all.tsv \
#   --output amrfinder_heatmap_test.tsv

if __name__ == "__main__":
    main()