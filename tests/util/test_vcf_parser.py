from spectre.util.vcf_parser import VCFSNVParser
import pandas as pd


def test_vcf_to_dataframe():
    vcf_parser = VCFSNVParser(as_dev=True)
    vcf_path = "tests/data/vcf_snv/NA12878_GIAB.chr22.vcf.gz"
    df = vcf_parser.vcf_to_dataframe(vcf_path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
