"""
1. Expression Data Normalization

Applies a log2(TPM + 1) transformation to the GEUVADIS TPM/RPKM
gene expression matrix.

Input:
    A TSV file in gene x sample format
    (index = gene_id, columns = sample_id)
    Values = TPM or RPKM

Output:
    A log2(TPM + 1)-transformed expression matrix
"""

import argparse
import numpy as np
import pandas as pd


def normalize_expression(input_path: str, output_path: str) -> pd.DataFrame:
    expr = pd.read_csv(input_path, sep="\t", index_col=0)

    # Prevent negative values and handle missing values
    expr = expr.clip(lower=0).fillna(0)

    log2_expr = np.log2(expr + 1)
    log2_expr.to_csv(output_path, sep="\t")
    return log2_expr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TPM/RPKM -> log2(TPM+1) normalization"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input expression matrix (TSV)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the normalized expression matrix"
    )
    args = parser.parse_args()

    result = normalize_expression(args.input, args.output)
    print(
        f"Normalization completed: "
        f"{result.shape[0]} genes x {result.shape[1]} samples "
        f"-> {args.output}"
    )
