import numpy as np

from spectre.analysis.call_cnv_coverage import CNVCall
from spectre.analysis.coverage_stats import CoverageData

BIN_SIZE = 1000


def call_candidates(normalized_cov):
    cov_data = CoverageData()
    cov_data.normalized_cov_ploidy = np.array(normalized_cov, dtype=float)
    cov_data.positions = np.arange(len(normalized_cov), dtype=float) * BIN_SIZE
    return CNVCall().cnv_coverage(cov_data, BIN_SIZE, "chr1", "sample", lower_bound=1.5, upper_bound=2.5)


def test_run_followed_by_normal_bin():
    candidates = call_candidates([2.0] * 5 + [3.0] * 20 + [2.0] * 5)

    assert [(c.type, c.start, c.end) for c in candidates] == [("DUP", 5000, 24000)]


def test_run_at_end_of_chromosome():
    candidates = call_candidates([2.0] * 5 + [3.0] * 20)

    assert [(c.type, c.start, c.end) for c in candidates] == [("DUP", 5000, 24000)]


def test_run_followed_only_by_nan_bins():
    candidates = call_candidates([np.nan] * 5 + [2.0] * 5 + [1.0] * 20 + [np.nan] * 5)

    assert [(c.type, c.start, c.end) for c in candidates] == [("DEL", 10000, 29000)]


def test_short_run_at_end_of_chromosome():
    candidates = call_candidates([2.0] * 5 + [3.0] * 9)

    assert candidates == []
