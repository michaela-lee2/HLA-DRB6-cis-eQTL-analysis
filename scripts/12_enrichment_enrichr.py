"""
12. 기능 농축 분석 (Enrichr - KEGG2021 Human, Reactome Pathways 2024)
PPI 네트워크에 포함된 유전자 목록을 Enrichr에 제출하여 경로 농축을 분석합니다.
통계적 유의성: 조정된 P-value(FDR) < 0.05

Enrichr API 문서: https://maayanlab.cloud/Enrichr/help#api
"""
import argparse
import json
import time
import pandas as pd
import requests

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

LIBRARIES = ["KEGG_2021_Human", "Reactome_Pathways_2024"]


def submit_gene_list(genes: list[str], description: str = "HLA-DRB1_network") -> str:
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, description),
    }
    resp = requests.post(ENRICHR_ADD_URL, files=payload, timeout=30)
    resp.raise_for_status()
    return resp.json()["userListId"]


def get_enrichment(user_list_id: str, library: str) -> pd.DataFrame:
    resp = requests.get(
        ENRICHR_ENRICH_URL,
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()[library]

    columns = [
        "rank", "term", "p_value", "z_score", "combined_score",
        "overlapping_genes", "adjusted_p_value", "old_p_value", "old_adjusted_p_value",
    ]
    return pd.DataFrame(data, columns=columns)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enrichr 경로 농축 분석")
    parser.add_argument("--genes", required=True,
                         help="유전자 목록 파일 (한 줄에 하나) 또는 콤마로 구분된 문자열")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR(adjusted P) 임계값")
    parser.add_argument("--outdir", required=True, help="결과 저장 디렉토리")
    args = parser.parse_args()

    if args.genes.endswith(".txt") or args.genes.endswith(".tsv"):
        with open(args.genes) as f:
            gene_list = [line.strip() for line in f if line.strip()]
    else:
        gene_list = [g.strip() for g in args.genes.split(",")]

    user_list_id = submit_gene_list(gene_list)
    time.sleep(1)

    for library in LIBRARIES:
        df = get_enrichment(user_list_id, library)
        sig_df = df[df["adjusted_p_value"] < args.fdr].sort_values("adjusted_p_value")

        out_path = f"{args.outdir}/enrichr_{library}.tsv"
        sig_df.to_csv(out_path, sep="\t", index=False)
        print(f"[{library}] 유의미한 경로 {len(sig_df)}개 (FDR<{args.fdr}) -> {out_path}")
