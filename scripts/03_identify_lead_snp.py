"""
3. Lead SNP Identification

Identifies the lead SNP showing the strongest association
(minimum P-value) from the cis-eQTL results.

Paper criteria:
P < 1e-50, BETA approximately -29 (rs9270694)
"""

import argparse
import os

import pandas as pd


def identify_lead_snp(
    assoc_path: str,
    p_col: str = "P",
    beta_col: str = "BETA"
) -> pd.Series:
    # Check whether the input file exists
    if not os.path.isfile(assoc_path):
        raise FileNotFoundError(f"Input file not found: {assoc_path}")

    # Load PLINK association results
    df = pd.read_csv(assoc_path, sep=r"\s+")

    # Check required columns
    required_columns = {p_col, beta_col, "SNP"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {', '.join(sorted(missing_columns))}"
        )

    # Remove rows with missing P-values
    df = df.dropna(subset=[p_col])

    if df.empty:
        raise ValueError("No valid SNPs with P-values were found.")

    # Select the SNP with the minimum P-value
    lead_snp = df.loc[df[p_col].idxmin()]

    print(f"Lead SNP: {lead_snp.get('SNP', 'NA')}")
    print(f"  P-value : {lead_snp[p_col]:.3e}")
    print(f"  BETA    : {lead_snp[beta_col]:.3f}")
    print(
        f"  Position: "
        f"chr{lead_snp.get('CHR', 'NA')}:{lead_snp.get('BP', 'NA')}"
    )

    return lead_snp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify the lead SNP from cis-eQTL results"
    )

    parser.add_argument(
        "--assoc",
        required=True,
        help="PLINK .assoc.linear result file"
    )

    args = parser.parse_args()

    identify_lead_snp(args.assoc)

