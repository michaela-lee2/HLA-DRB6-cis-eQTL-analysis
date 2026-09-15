"""
12. Pathway Enrichment Analysis
   (Enrichr - KEGG_2021_Human, Reactome_Pathways_2024)

Submits the gene list from the PPI network to Enrichr and performs
pathway enrichment analysis.

Statistical significance:
Adjusted P-value (FDR) < 0.05

Enrichr API documentation:
https://maayanlab.cloud/Enrichr/help#api
"""

import argparse
import os
import time

import pandas as pd
import requests


ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_Pathways_2024",
]


def submit_gene_list(
    genes: list[str],
    description: str = "HLA-DRB1_network",
) -> str:
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, description),
    }

    resp = requests.post(
        ENRICHR_ADD_URL,
        files=payload,
        timeout=30,
    )

    resp.raise_for_status()

    return resp.json()["userListId"]


def get_enrichment(
    user_list_id: str,
    library: str,
) -> pd.DataFrame:
    resp = requests.get(
        ENRICHR_ENRICH_URL,
        params={
            "userListId": user_list_id,
            "backgroundType": library,
        },
        timeout=30,
    )

    resp.raise_for_status()

    data = resp.json()[library]

    columns = [
        "rank",
        "term",
        "p_value",
        "z_score",
        "combined_score",
        "overlapping_genes",
        "adjusted_p_value",
        "old_p_value",
        "old_adjusted_p_value",
    ]

    return pd.DataFrame(
        data,
        columns=columns,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform pathway enrichment analysis using Enrichr"
    )

    parser.add_argument(
        "--genes",
        required=True,
        help="Gene list file (one gene per line) or comma-separated gene list",
    )

    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="FDR (adjusted P-value) threshold",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for enrichment results",
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.genes.endswith(".txt") or args.genes.endswith(".tsv"):
        with open(args.genes) as f:
            gene_list = [
                line.strip()
                for line in f
                if line.strip()
            ]
    else:
        gene_list = [
            gene.strip()
            for gene in args.genes.split(",")
            if gene.strip()
        ]

    user_list_id = submit_gene_list(gene_list)

    time.sleep(1)

    for library in LIBRARIES:
        df = get_enrichment(
            user_list_id,
            library,
        )

        sig_df = (
            df[df["adjusted_p_value"] < args.fdr]
            .sort_values("adjusted_p_value")
        )

        out_path = os.path.join(
            args.outdir,
            f"enrichr_{library}.tsv",
        )

        sig_df.to_csv(
            out_path,
            sep="\t",
            index=False,
        )

        print(
            f"[{library}] "
            f"{len(sig_df)} significant pathways "
            f"(FDR < {args.fdr}) -> {out_path}"
        )
