"""
6. Functional Annotation (RegulomeDB / HaploReg)

Queries RegulomeDB functional rankings and HaploReg
"Motifs changed" annotations for candidate SNPs within
the LD block and merges the results.

Corresponds to Tables 1 and 2.

Notes:
- RegulomeDB provides a public REST API.
- HaploReg does not provide an official API, so bulk annotation
  files are downloaded and queried locally.
- HaploReg bulk data:
  https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php
"""

import argparse
import time

import pandas as pd
import requests


REGULOMEDB_API = "https://www.regulomedb.org/regulome-search/"


def query_regulomedb(rsids: list[str]) -> pd.DataFrame:
    rows = []

    for rsid in rsids:
        try:
            resp = requests.get(
                REGULOMEDB_API,
                params={
                    "regions": rsid,
                    "genome": "GRCh38",
                    "format": "json",
                },
                timeout=15,
            )

            resp.raise_for_status()

            data = resp.json()
            features = data.get("features", [])
            ranking = features[0].get("ranking") if features else None

            rows.append({
                "SNP": rsid,
                "RegulomeDB_rank": ranking,
            })

        except Exception as e:  # noqa: BLE001
            rows.append({
                "SNP": rsid,
                "RegulomeDB_rank": None,
                "error": str(e),
            })

        # Pause between requests to avoid excessive API requests
        time.sleep(0.3)

    return pd.DataFrame(rows)


def annotate_with_haploreg(
    candidate_df: pd.DataFrame,
    haploreg_bulk_path: str
) -> pd.DataFrame:
    """
    Merge Motifs_changed and RFX5-related annotations
    from a HaploReg bulk annotation file.
    """

    haploreg = pd.read_csv(haploreg_bulk_path, sep="\t")

    cols = ["rsID", "Motifs_changed"]
    haploreg = haploreg[
        [c for c in cols if c in haploreg.columns]
    ]

    merged = candidate_df.merge(
        haploreg,
        left_on="SNP",
        right_on="rsID",
        how="left",
    )

    merged["RFX5_related"] = (
        merged["Motifs_changed"]
        .astype(str)
        .str.contains("RFX5", na=False)
    )

    return merged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Functional annotation using RegulomeDB and HaploReg"
    )

    parser.add_argument(
        "--candidates",
        required=True,
        help="LD candidate SNP list (column: SNP)",
    )

    parser.add_argument(
        "--haploreg_bulk",
        required=True,
        help="HaploReg bulk annotation TSV file",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the merged annotation results",
    )

    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates, sep="\t")
    snp_list = candidates["SNP"].dropna().unique().tolist()

    regulome_df = query_regulomedb(snp_list)

    merged = candidates.merge(
        regulome_df,
        on="SNP",
        how="left",
    )

    merged = annotate_with_haploreg(
        merged,
        args.haploreg_bulk,
    )

    merged.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    print(f"Functional annotation completed -> {args.output}")
