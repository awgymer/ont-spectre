import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, mock_open

from spectre.analysis.analysis import CNVAnalysis
from spectre.util.vcf_parser import VCFSNVParser

@pytest.fixture
def cnv_analysis():
    coverage_file = "dummy_coverage_file"
    coverage_mosdepth_data = MagicMock()
    coverage_mosdepth_data.genome_mean_coverage = 30

    bin_size = 1000
    output_directory = "dummy_output_directory"
    outbed = "dummy_outbed"
    outvcf = "dummy_outvcf"
    genome_info = {
        "chromosomes": ["chr1", "chrX", "chrY"],
        "chr_lengths_by_name": {"chr1": 1000, "chrX": 1000, "chrY": 1000}
    }
    sample_id = "dummy_sample_id"
    snv_file = "dummy_snv_file"
    metadata_ref = {}
    
    return CNVAnalysis(
        coverage_file, coverage_mosdepth_data, bin_size, output_directory, outbed, outvcf, genome_info, 
        sample_id, metadata_ref, snv_file, only_chr_list="chrX,chrY,chr1", as_dev=True
    )

@pytest.mark.parametrize("normalized_cov, expected_ploidy", [
    (1.3, 1),
    (2.0, 2),
    (2.7, 3),
])
def test_detect_chromosome_ploidy(cnv_analysis, normalized_cov, expected_ploidy):
    cnv_analysis.lower_2n_threshold = 1.5
    cnv_analysis.upper_2n_threshold = 2.5

    ploidy = cnv_analysis.detect_chromosome_ploidy(normalized_cov)
    assert ploidy == expected_ploidy

@pytest.mark.parametrize("chrX_copies, chrY_copies", [
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
    (3, 0),
    (3, 1),
    (3, 2),
])
def test_chromosome_ploidy_detection_combinations(monkeypatch, cnv_analysis, chrX_copies, chrY_copies):
    chrX_coverage = 30 * chrX_copies / cnv_analysis.ploidy
    chrY_coverage = 30 * chrY_copies / cnv_analysis.ploidy

    read_data = f'chrX\t0\t1000\t{chrX_coverage}\nchrY\t0\t1000\t{chrY_coverage}\n'
    m = mock_open(read_data=read_data)
    monkeypatch.setattr('builtins.open', m)

    cnv_analysis.coverages_df_diploid = MagicMock()
    cnv_analysis.coverages_df_diploid['coverage_'].median.return_value = 30
    cnv_analysis.coverages_df_diploid['coverage_'].quantile.side_effect = [20, 40]

    cnv_analysis.data_normalization()

    assert cnv_analysis.sex_choromosome_ploidies['chrX'] == chrX_copies
    assert cnv_analysis.sex_choromosome_ploidies['chrY'] == chrY_copies


@pytest.mark.parametrize("snvs, expected_median", [
    # no SNVs: all usable bins with coverage give the diploid coverage
    ([], 31),
    # a heterozygous SNV in a usable bin
    ([("chr1", 500, 0.5)], 30),
    # a heterozygous SNV in a blacklisted bin and one on a chromosome that is not selected
    ([("chr1", 2500, 0.5), ("chr2", 500, 0.5)], 31),
])
def test_diploid_coverage_uses_only_usable_bins(tmp_path, monkeypatch, cnv_analysis, snvs, expected_median):
    coverage_file = tmp_path / "sample.regions.bed"
    coverage_file.write_text(
        "chr1\t0\t1000\t30\n"
        "chr1\t1000\t2000\t32\n"
        "chr1\t2000\t3000\t0.5\n"
        "chr1\t3000\t4000\t0.4\n"
        "chr2\t0\t1000\t0.3\n"
        "chr2\t1000\t2000\t0.2\n"
    )
    cnv_analysis.coverage_full_path = str(coverage_file)
    cnv_analysis.only_chr_list = ["chr1"]
    cnv_analysis.metadata = {"chr1": [("2000", "4000")]}

    vcf_df = pd.DataFrame(snvs, columns=["chrom_", "start_", "af_"]).astype(
        {"chrom_": "string", "start_": "int64", "af_": "float64"})
    monkeypatch.setattr(VCFSNVParser, "vcf_to_dataframe", lambda self, path: vcf_df)

    cnv_analysis.vcf_based_af_bins_annotation()

    assert cnv_analysis.coverages_df["start_"].to_list() == [0, 1000]
    assert cnv_analysis.genome_median() == expected_median
