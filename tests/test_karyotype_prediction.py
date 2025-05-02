import os
import pytest
import urllib.request
import tarfile
import gzip
import shutil

from spectre.main import run_main

# Constants for URLs and file paths
EPI2ME_ARTIFACT_URL = os.environ.get("EPI2ME_ARTIFACT_URL", "default_url")
TEST_DATA_URL = f"{EPI2ME_ARTIFACT_URL}/data/ont-spectre/karyotype_prediction_test_data_v1.tar.gz"
REFERENCE_GENOME_URL = f"{EPI2ME_ARTIFACT_URL}/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bgzf.gz"

def download_file(url, dest):
    """Download a file from a URL to a destination."""
    with urllib.request.urlopen(url) as response, open(dest, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)

def extract_tar_file(tar_path, extract_path):
    """Extract a tar.gz file."""
    with tarfile.open(tar_path, "r:") as tar:
        tar.extractall(path=extract_path)

def decompress_gzip_file(gz_path, out_path):
    """Decompress a .gz file."""
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

@pytest.fixture(scope="session")
def setup_test_environment(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("karyotype_test_data")

    # Download and extract test data into test_dir
    test_data_tar = test_dir / "karyotype_prediction_test_data_v1.tar.gz"
    download_file(TEST_DATA_URL, test_data_tar)
    extract_tar_file(str(test_data_tar), str(test_dir))

    # Download the human hg38 reference genome into test_dir
    reference_genome_gz = test_dir / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bgzf.gz"
    download_file(REFERENCE_GENOME_URL, reference_genome_gz)

    yield test_dir

    shutil.rmtree(test_dir)

@pytest.mark.parametrize("sample_dir", ['GM18501', 'GM18861', 'GM18864', 'NA18310'])
def test_karyotype_prediction(setup_test_environment, sample_dir):
    sample_dir_path = setup_test_environment / "karyotype_prediction_test_data" / sample_dir

    coverage_dir = str(sample_dir_path / "mosdepth")
    snv_path = str(sample_dir_path / "wf_snp.vcf.gz")
    output_dir = str(sample_dir_path / "output_spectre")
    reference_path = str(setup_test_environment / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bgzf.gz")

    args = [
        "CNVCaller",
        "--bin-size", "1000",
        "--coverage", coverage_dir,
        "--snv", snv_path,
        "--sample-id", "sample",
        "--output-dir", output_dir,
        "--reference", reference_path,
        "--metadata", "hg38_metadata",
        "--blacklist", "hg38_blacklist_v1.0"
    ]
    # Call the main function which should use the mocked arguments and exit with 0
    with pytest.raises(SystemExit) as exc_info:
        run_main(args)

    assert exc_info.type == SystemExit
    assert exc_info.value.code == 0

    predicted_file = sample_dir_path / "output_spectre" / "predicted_karyotype.txt"
    with open(predicted_file) as f:
        predicted_karyotype = f.read().strip()

    expected_file = sample_dir_path / "expected_karyotype.txt"
    with open(expected_file) as f:
        expected_karyotype = f.read().strip()

    assert predicted_karyotype == expected_karyotype, f"Karyotype mismatch: {predicted_karyotype} != {expected_karyotype}"
