"""
6. 기능 주석 (RegulomeDB / HaploReg)
LD 블록 내 후보 SNP들에 대해 RegulomeDB 기능 순위와 HaploReg의
'Motifs changed' 정보를 조회하여 병합합니다. (표1, 표2에 해당)

주의: RegulomeDB는 공개 REST API를 제공합니다.
HaploReg는 공식 API가 없어 다운로드 받은 벌크 주석 파일을 로컬에서 조회하는 방식을 사용합니다.
(HaploReg 벌크 데이터: https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php)
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
                params={"regions": rsid, "genome": "GRCh38", "format": "json"},
                timeout=15,
            )
            resp.raise_for_status()
            data = resp.json()
            features = data.get("features", [])
            ranking = features[0].get("ranking") if features else None
            rows.append({"SNP": rsid, "RegulomeDB_rank": ranking})
        except Exception as e:  # noqa: BLE001
            rows.append({"SNP": rsid, "RegulomeDB_rank": None, "error": str(e)})
        time.sleep(0.3)  # API rate limit 배려
    return pd.DataFrame(rows)


def annotate_with_haploreg(candidate_df: pd.DataFrame, haploreg_bulk_path: str) -> pd.DataFrame:
    """HaploReg 벌크 다운로드 파일에서 Motifs_changed / RFX5 관련 정보를 병합."""
    haploreg = pd.read_csv(haploreg_bulk_path, sep="\t")
    cols = ["rsID", "Motifs_changed"]
    haploreg = haploreg[[c for c in cols if c in haploreg.columns]]
    merged = candidate_df.merge(haploreg, left_on="SNP", right_on="rsID", how="left")
    merged["RFX5_related"] = merged["Motifs_changed"].astype(str).str.contains("RFX5", na=False)
    return merged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RegulomeDB/HaploReg 기능 주석")
    parser.add_argument("--candidates", required=True, help="LD 후보 SNP 목록 (컬럼: SNP)")
    parser.add_argument("--haploreg_bulk", required=True, help="HaploReg 벌크 다운로드 TSV")
    parser.add_argument("--output", required=True, help="병합된 주석 결과 출력 경로")
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates, sep="\t")
    snp_list = candidates["SNP"].dropna().unique().tolist()

    regulome_df = query_regulomedb(snp_list)
    merged = candidates.merge(regulome_df, on="SNP", how="left")
    merged = annotate_with_haploreg(merged, args.haploreg_bulk)

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"기능 주석 완료 -> {args.output}")
