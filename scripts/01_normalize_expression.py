"""
1. 발현 데이터 정규화
GEUVADIS TPM/RPKM 발현 매트릭스에 log2(TPM+1) 변환을 적용합니다.

입력: gene x sample 형태의 TSV (index=gene_id, columns=sample_id), 값=TPM 또는 RPKM
출력: log2(TPM+1) 변환된 매트릭스
"""
import argparse
import numpy as np
import pandas as pd


def normalize_expression(input_path: str, output_path: str) -> pd.DataFrame:
    expr = pd.read_csv(input_path, sep="\t", index_col=0)

    # 음수/결측 방지
    expr = expr.clip(lower=0).fillna(0)

    log2_expr = np.log2(expr + 1)
    log2_expr.to_csv(output_path, sep="\t")
    return log2_expr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TPM/RPKM -> log2(TPM+1) 정규화")
    parser.add_argument("--input", required=True, help="원본 발현 매트릭스 (TSV)")
    parser.add_argument("--output", required=True, help="정규화된 발현 매트릭스 출력 경로")
    args = parser.parse_args()

    result = normalize_expression(args.input, args.output)
    print(f"정규화 완료: {result.shape[0]}개 유전자 x {result.shape[1]}개 샘플 -> {args.output}")
