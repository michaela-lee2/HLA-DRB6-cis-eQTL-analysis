```python
"""
11. STRING-Based Protein-Protein Interaction (PPI) Network Construction

Constructs a protein-protein interaction (PPI) network centered on
HLA-DRB1, which is used as a surrogate gene for HLA-DRB6.
Corresponds to Figure 2_A.

STRING API documentation:
https://string-db.org/help/api/
"""

import argparse
from io import StringIO

import pandas as pd
import requests


STRING_API = "https://string-db.org/api"


def get_string_network(
    gene: str,
    species: int = 9606,
    required_score: int = 400,
    limit: int = 50,
) -> pd.DataFrame:
    """
    Retrieve a TSV-formatted protein interaction network
    centered on the specified gene from the STRING API.
    """

    url = f"{STRING_API}/tsv/network"

    params = {
        "identifiers": gene,
        "species": species,  # 9606 = Homo sapiens
        "required_score": required_score,  # STRING confidence score threshold
        "limit": limit,
        "caller_identity": "hla-drb6-eqtl-pipeline",
    }

    resp = requests.get(
        url,
        params=params,
        timeout=30,
    )

    resp.raise_for_status()

    return pd.read_csv(
        StringIO(resp.text),
        sep="\t",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Construct a PPI network using the STRING API"
    )

    parser.add_argument(
        "--gene",
        default="HLA-DRB1",
        help="Central gene (default: HLA-DRB1)",
    )

    parser.add_argument(
        "--score",
        type=int,
        default=400,
        help="STRING confidence score (0-1000)",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the network edge list (TSV)",
    )

    args = parser.parse_args()

    network = get_string_network(
        args.gene,
        required_score=args.score,
    )

    network.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    genes_in_network = sorted(
        set(network["preferredName_A"])
        | set(network["preferredName_B"])
    )

    print(
        f"Number of PPI network nodes: {len(genes_in_network)}"
    )

    print(
        f"Edge list saved: {args.output}"
    )

    print(
        "Network gene list (for Enrichr):",
        ", ".join(genes_in_network),
    )
```
