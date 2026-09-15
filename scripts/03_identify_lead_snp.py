"""
3. 선도 SNP(Lead SNP) 식별
cis-eQTL 결과에서 가장 강력한 연관성(최소 P-value)을 보인 SNP를 선도 SNP로 선택합니다.
논문 기준: P < 1e-50, BETA 약 -29 (rs9270694)
"""
import argparse
import pandas as pd


def identify_lead_snp(assoc_path: str, p_col: str = "P", beta_col: str = "BETA") -> pd.Series:
    df = pd.read_csv(assoc_path, delim_whitespace=True)
    df = df.dropna(subset=[p_col])

    lead_snp = df.loc[df[p_col].idxmin()]

    print(f"선도 SNP: {lead_snp.get('SNP', 'NA')}")
    print(f"  P-value : {lead_snp[p_col]:.3e}")
    print(f"  BETA    : {lead_snp[beta_col]:.3f}")
    print(f"  위치    : chr{lead_snp.get('CHR', 'NA')}:{lead_snp.get('BP', 'NA')}")

    return lead_snp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cis-eQTL 결과에서 선도 SNP 식별")
    parser.add_argument("--assoc", required=True, help="PLINK .assoc.linear 결과 파일")
    args = parser.parse_args()

    identify_lead_snp(args.assoc)
