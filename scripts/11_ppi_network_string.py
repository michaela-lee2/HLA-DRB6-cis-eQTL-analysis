"""
11. STRING 기반 단백질-단백질 상호작용(PPI) 네트워크 구축
HLA-DRB1을 대리 유전자(surrogate gene)로 사용하여 중심 노드로 하는 PPI 네트워크를 구축합니다.
(그림 2_A에 해당)

STRING API 문서: https://string-db.org/help/api/
"""
import argparse
import pandas as pd
import requests

STRING_API = "https://string-db.org/api"


def get_string_network(gene: str, species: int = 9606, required_score: int = 400,
                        limit: int = 50) -> pd.DataFrame:
    """중심 유전자 주변의 상호작용 네트워크를 TSV 형태로 가져옵니다."""
    url = f"{STRING_API}/tsv/network"
    params = {
        "identifiers": gene,
        "species": species,           # 9606 = Homo sapiens
        "required_score": required_score,  # STRING confidence score threshold
        "limit": limit,
        "caller_identity": "hla-drb6-eqtl-pipeline",
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    from io import StringIO
    return pd.read_csv(StringIO(resp.text), sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STRING API로 PPI 네트워크 구축")
    parser.add_argument("--gene", default="HLA-DRB1", help="중심 유전자 (기본: HLA-DRB1)")
    parser.add_argument("--score", type=int, default=400, help="STRING confidence score (0-1000)")
    parser.add_argument("--output", required=True, help="네트워크 엣지 리스트 출력 경로 (TSV)")
    args = parser.parse_args()

    network = get_string_network(args.gene, required_score=args.score)
    network.to_csv(args.output, sep="\t", index=False)

    genes_in_network = sorted(set(network["preferredName_A"]) | set(network["preferredName_B"]))
    print(f"PPI 네트워크 노드 수: {len(genes_in_network)}")
    print(f"엣지 목록 저장 완료: {args.output}")
    print("네트워크 유전자 목록 (Enrichr 입력용):", ", ".join(genes_in_network))
