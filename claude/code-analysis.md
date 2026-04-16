# Spectre CNV Caller — Code Analysis

_Generated: 2026-04-16_

## Overview

Spectre is a Python-based Copy Number Variation (CNV) caller for Oxford Nanopore sequencing data. This report documents bugs and performance issues found through a deep static analysis of the source code.

---

## Critical Bugs

### 1. NaN comparison using equality operator
**File:** `spectre/analysis/cnv_candidate.py:92`

```python
# WRONG — np.nan == np.nan is always False per IEEE 754
if median_candidates_coverage == np.inf or median_candidates_coverage == np.nan:
```

Comparing with `np.nan` using `==` always returns `False`. NaN values silently pass through this guard, producing garbage copy-number values from the integer conversion below.

**Fix:** `np.isnan(median_candidates_coverage) or np.isinf(median_candidates_coverage)`

---

### 2. `np.nan` called as a function
**File:** `spectre/analysis/call_cnv_AF.py:69`

```python
return np.nan()  # TypeError: 'float' object is not callable
```

`np.nan` is a constant, not a callable. This raises `TypeError` whenever the AF/CN concordance condition is unmet, crashing SNV-based filtering entirely.

**Fix:** `return np.nan`

---

### 3. Wrong argument passed to `id_generator`
**File:** `spectre/util/cnv_id.py:33`

```python
# n = number of IDs to generate, but it's passed as the 'size' (length) argument
id = cls.id_generator(n, size, chars)
```

The loop counter `n` (total IDs requested) is forwarded as the `size` parameter of `id_generator`. IDs end up far longer than intended and vary in length depending on how many IDs were requested.

**Fix:** `id = cls.id_generator(size=size, chars=chars)`

---

### 4. Bare `except:` silencing all errors
**File:** `spectre/spectreCNVPopulation.py:55,75`

```python
except:
    self.logger.error(f"File does not exist! Provided : {file}")
```

Bare `except:` catches `KeyboardInterrupt`, `SystemExit`, and unexpected exceptions, masking real bugs. A malformed `.spc` file appears to load successfully but with missing data.

**Fix:** `except (FileNotFoundError, json.JSONDecodeError, IOError) as e:`

---

## Major Bugs

### 5. Inverted population-mode flag
**File:** `spectre/util/vcf_parser.py:191`

```python
population_mode = len(df.columns[9:]) < 2  # comment says TRUE when >2 samples — backwards
```

The condition and comment contradict each other. `< 2` means 0–1 samples (single-sample), so `population_mode` is `True` only in single-sample mode — the opposite of what is intended. All multi-sample VCFs are handled as single-sample.

**Fix:** `population_mode = len(df.columns[9:]) >= 2`

---

### 6. Unfinished function that always returns `""`
**File:** `spectre/util/outputWriter.py:194–197`

```python
@staticmethod
def convert_genome_info_to_dictionary(genome_info: dict):
    tmp_genome_inf = dict([(key, []) for key in genome_info.keys()])
    return ""          # always returns empty string; body never executes
    # for info in genome_info;  ← incomplete dead code
```

The function unconditionally returns an empty string. Its parameter is never used and the commented-out loop is syntactically broken. Any caller that depends on the return value silently receives `""`.

---

### 7. Variable shadowing in nested loop
**File:** `spectre/util/outputWriter.py:202–208`

```python
for key, candidates in candidates.items():      # 'candidates' shadows the parameter
    ...
    for key, value in tmp_dict.items():         # 'key' clobbers the chromosome name
```

The outer loop variable `key` (chromosome name) is overwritten by the inner loop. After the first inner iteration, the chromosome key is lost, producing incorrect dictionary writes for every subsequent chromosome.

**Fix:** Rename inner loop variables (`chrom_key` / `field_key`).

---

## Minor Bugs and Code Quality Issues

### 8. Error logged but execution continues with invalid config
**File:** `spectre/main.py:507,512,518`

```python
logger.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
```

The error is logged but the program keeps running with the invalid bin size, producing silent downstream failures.

**Fix:** Follow the log call with `sys.exit(1)`.

---

### 9. Unguarded `IndexError` from `np.where` result
**File:** `spectre/analysis/analysis.py:503`

```python
gap_start = int(np.where(cov_data_chr.positions == gap_in_candidate_start)[0][0])
```

If the position is absent from the coverage array, `np.where` returns an empty array and `[0][0]` raises `IndexError` with no diagnostic information.

**Fix:** Check `len(indices) > 0` before indexing; log a warning and `continue` if empty.

---

### 10. Potential divide-by-zero in allele frequency calculation
**File:** `spectre/analysis/call_cnv_AF.py:98`

```python
af = float(alt_count) / (float(ref_count) + float(alt_count))
```

If both counts are zero (malformed VCF entry), this raises `ZeroDivisionError`.

**Fix:** Guard with `if float(ref_count) + float(alt_count) > 0:` and skip or default otherwise.

---

### 11. Non-Pythonic `__contains__` method calls
**File:** `spectre/spectreCNVPopulation.py:51`

```python
if str(file).__contains__(".spc") or str(file).__contains__(".spc.gz"):
```

Calling `__contains__` directly is non-idiomatic and unnecessary. Minor readability issue; functionally equivalent.

**Fix:** `if file.endswith((".spc", ".spc.gz")):`

---

### 12. Pandas `rolling(step=)` deprecation
**File:** `spectre/analysis/analysis.py:667`

```python
.rolling(self.windows_bins, step=self.windows_bins).agg(...)
```

The `step` keyword was deprecated in newer pandas versions. This will raise warnings or break with pandas >= 2.0.

**Fix:** `rolling(self.windows_bins).agg(...)[::self.windows_bins]`

---

## Performance Issues

### 13. O(n^6) algorithm in population mode overlap detection
**File:** `spectre/spectreCNVPopulation.py:161–189`

The cross-sample candidate comparison uses six nested loops:

```
for sample1 in all_samples:
  for sample2 in all_samples:        # includes self-comparisons
    for chrom1 in sample1.chroms:
      for chrom2 in sample2.chroms:  # cross-chromosome comparisons
        for cnv1 in chrom1.candidates:
          for cnv2 in chrom2.candidates:
```

With even a modest cohort (10 samples × 500 CNVs per chromosome × 24 chromosomes), this is on the order of billions of comparisons. Cross-chromosome pairs can never overlap and are entirely wasted work.

**Fixes (in priority order):**
1. Skip self-comparisons: `if samples_key1 >= samples_key2: continue`
2. Filter to matching chromosomes before the inner loops
3. Replace the linear interval scan with an interval tree (e.g. `ncls` or `intervaltree`) for O(n log n) overlap queries
4. Compute pairwise comparisons only once (symmetric matrix) rather than twice

---

## Summary

| Severity | Count |
|----------|-------|
| Critical | 4 |
| Major | 3 |
| Minor / Code Quality | 5 |
| Performance | 1 |
| **Total** | **13** |

The two highest-priority fixes are the `np.nan()` call (item 2, guaranteed crash) and the inverted `population_mode` flag (item 5, silent wrong-mode execution for all multi-sample VCFs). The O(n^6) loop (item 13) will become a bottleneck the moment population mode is used with a real cohort.

---

## Deep-Dive Pass: Additional Findings

_Extended analysis pass, 2026-04-16_

---

### 14. CNV runs at the end of chromosomes are silently dropped
**File:** `spectre/analysis/call_cnv_coverage.py:42–90`

This is the most impactful newly-found bug. The `cnv_coverage` function emits a candidate **only** when the run is broken by (a) a coverage value returning to the normal range, or (b) a gap/type change. There is no post-loop flush:

```python
for (cov, pos) in zip(...):
    if cov > upper_bound or cov < lower_bound:
        ...  # appends to cnv_pos_cand_list
    else:
        if len(cnv_pos_cand_list) > 1:
            candidate_list.append(cnv_candidates)  # only emit on break
        ...

return candidate_list  # any in-progress run is silently discarded
```

If the final bins of a chromosome are all abnormal coverage (a CNV reaching a telomere), the accumulated `cnv_pos_cand_list` is never pushed and that CNV is lost. This is a systematic bias: subtelomeric CNVs are never reported.

**Fix:** After the loop, flush any remaining run:
```python
if len(cnv_pos_cand_list) >= min_run:
    cnv_candidates = CNVCandidate(sample_origin, self.as_dev)
    cnv_candidates.push_candidates(chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list, current_cnv_type)
    candidate_list.append(cnv_candidates)
return candidate_list
```

---

### 15. First position of every CNV run is silently dropped
**File:** `spectre/analysis/call_cnv_coverage.py:49–51`

When a CNV run begins for the first time on a chromosome (`run_start == 0`), the code records the start position but does **not** append the current bin to the candidate lists:

```python
if run_start == 0:
    current_cnv_type = cnv_type
    run_start = pos        # records start, but does NOT append pos/cov
    # cnv_pos_cand_list stays []
```

The else branch (which does `cnv_pos_cand_list.append(pos)`) is only reached from the second bin onward. Every subsequent run restart (lines 68–72) correctly uses `cnv_pos_cand_list = [pos]`, so only the absolute first run on the chromosome is affected. The result is that the first bin of the first CNV on each chromosome is always missing, biasing start coordinates by one bin width.

**Fix:** Append the position when initialising the run:
```python
if run_start == 0:
    current_cnv_type = cnv_type
    run_start = pos
    cnv_pos_cand_list = [pos]   # add
    cnv_cov_cand_list = [cov]   # add
```

---

### 16. Inverted bounds assignment in coverage plots
**File:** `spectre/plots/plot.py:78`

The call site passes `[lower_bound, upper_bound]` (analysis.py:582), but the unpack reverses the names:

```python
# Caller:  [lower_bound, upper_bound]
[upperb, lowerb] = bounds   # upperb = lower_bound, lowerb = upper_bound
...
self.main_plot.plot(..., np.array([lowerb, lowerb]), ...)  # plots upper value as lower line
self.main_plot.plot(..., np.array([upperb, upperb]), ...)  # plots lower value as upper line
```

Both threshold lines are drawn at the wrong positions. Because both use the same colour, this isn't immediately obvious but makes every generated plot misleading.

**Fix:** `[lowerb, upperb] = bounds`

---

### 17. Incomplete interval overlap check misses containment
**File:** `spectre/spectreCNVPopulation.py:133`

```python
return cnv1.start <= cnv2.start <= cnv1.end or cnv1.start <= cnv2.end <= cnv1.end
```

This correctly detects when an endpoint of `cnv2` falls inside `cnv1`, but silently misses the case where `cnv2` fully contains `cnv1` (e.g. cnv1=[100,200], cnv2=[50,300]). In that scenario both conditions are False and the method returns False, causing population mode to miss valid supporting evidence.

**Fix:**
```python
return not (cnv2.end < cnv1.start or cnv2.start > cnv1.end)
```

---

### 18. `random.sample` raises `ValueError` on small chromosomes
**File:** `spectre/analysis/cnv_metrics.py:138–140`

```python
last_index = len(df_source.index) - 1
cov_window_amount = int(last_index * fraction)   # fraction = 0.1
samples_indices = random.sample(range(0, last_index), cov_window_amount)
```

`random.sample(population, k)` raises `ValueError` when `k > len(population)`. Both conditions can trigger this:
- `last_index < 10` → `cov_window_amount` rounds to 0, then `range(0, last_index)` has fewer than 1 element
- `last_index == 0` (empty coverage file) → `range(0, 0)` is empty, any `k > 0` crashes

**Fix:**
```python
cov_window_amount = max(1, int(last_index * fraction))
cov_window_amount = min(cov_window_amount, last_index)
if last_index == 0:
    return []
```

---

### 19. Unguarded list access after `filter()` in longshot VCF parsing
**File:** `spectre/util/vcf_parser.py:96`

```python
ac = list(filter(lambda x: "AC" in x, info_.split(";")))[0].split("=")[1]
```

If the INFO field contains no `AC=` entry, the filter returns an empty list and `[0]` raises `IndexError`. Any longshot VCF with a non-standard or missing AC field in the INFO column crashes the entire SNV parsing step.

**Fix:**
```python
ac_matches = list(filter(lambda x: x.startswith("AC="), info_.split(";")))
if not ac_matches:
    af = "NA"
else:
    ac = ac_matches[0].split("=")[1]
    ...
```

Note: the original also uses `"AC" in x` which would falsely match tags like `MLEAC=` or `SVAC=`. `x.startswith("AC=")` is safer.

---

### 20. Divide-by-zero in metadata GC% calculation
**File:** `spectre/util/metadata/metadataCollector.py:99`

```python
str((self.__gCnt + self.__cCnt) / self.__nTo * 100)
```

`self.__nTo` is the total base-pair count. If the input FASTA is empty or consists entirely of N regions, `__nTo` remains 0 and this raises `ZeroDivisionError`, crashing the metadata collection step before any analysis begins.

**Fix:** `str(round((self.__gCnt + self.__cCnt) / self.__nTo * 100, 4) if self.__nTo > 0 else 0.0)`

---

### 21. `pop('logger')` without a default raises `KeyError`
**File:** `spectre/util/outputWriter.py:211`

```python
tmp_dict.pop('logger')
```

If a `CNVCandidate` subclass or a manually-constructed candidate lacks a `logger` attribute, `vars()` won't include it and `pop()` raises `KeyError`, crashing intermediate file serialisation.

**Fix:** `tmp_dict.pop('logger', None)`

---

### 22. CNV type `KeyError` when plotting unknown types
**File:** `spectre/plots/plot.py:72`

```python
cnv_color = self.cnv_color[cnv.type]
```

`self.cnv_color` is a dict with only `"DUP"` and `"DEL"` keys. If a candidate has any other type string (e.g. from a future caller change or data corruption), this raises `KeyError` and aborts plotting for the entire chromosome.

**Fix:**
```python
cnv_color = self.cnv_color.get(cnv.type, "#888888")
```

---

## Updated Summary

| Severity | Count |
|----------|-------|
| Critical (guaranteed crash or data loss) | 6 |
| Major (wrong results, silent failures) | 6 |
| Minor / Code Quality | 7 |
| Performance | 1 |
| **Total** | **20** |

**Highest priority new fixes:**
- Item 14 (subtelomeric CNVs silently dropped) — systematic false-negative bias
- Item 15 (first bin of each first run dropped) — systematic coordinate bias
- Item 16 (plot thresholds inverted) — all coverage plots are misleading
- Item 17 (incomplete overlap detection) — population mode misses containment events

---

## Issue #9 Follow-up: `cnv_metrics.py` Cluster

_Reported upstream: nanoporetech/ont-spectre#9_

There are three related bugs in `spectre/analysis/cnv_metrics.py`, all stemming from the same root cause: after `dropna()` creates a subset DataFrame, the code alternates between writing to the subset and reading from the original, so exclusion zones silently have no effect.

---

### 23. `.iloc` with label-based indices → `IndexError` crash
**File:** `spectre/analysis/cnv_metrics.py:229–230`

```python
df_coverage_candidate_no_Nans.iloc[
    excl_zone_blacklist_indices,
    df_coverage_candidate_no_Nans.columns.get_loc("blacklist")] = True
```

`get_exclusion_zones_indices_in_coverage_data` returns **label** indices (via `.index` at line 118). `df_coverage_candidate_no_Nans` was produced by `dropna()`, which preserves original labels but resets positional positions. `.iloc` expects positional integers, so when a label like `2713038` is passed to a DataFrame with only `2703072` rows, Python raises:

```
IndexError: index 2713038 is out of bounds for axis 0 with size 2703072
```

The correct fix was already written and commented out directly above (line 227):
```python
# df_coverage_candidate_no_Nans.loc[excl_zone_blacklist_indices, "blacklist"] = True
```
**Fix:** Uncomment line 227, delete lines 229–230.

---

### 24. Blacklist filter reads from the wrong (unmodified) DataFrame
**File:** `spectre/analysis/cnv_metrics.py:233`

```python
df_coverage_candidate_no_blacklist = df_coverage_candidate_no_Nans.loc[
    df_coverage_candidate.blacklist == False]   # ← original df, not the no-NaN subset
```

`df_coverage_candidate["blacklist"]` was initialised to `False` at line 215 and is **never modified** — all modifications land on `df_coverage_candidate_no_Nans`. The mask `df_coverage_candidate.blacklist == False` is therefore all-`True`, so every row passes the filter and the blacklist has zero effect. All blacklisted regions remain in the random sample pool used for z-score scoring.

**Fix:** `df_coverage_candidate_no_Nans.loc[df_coverage_candidate_no_Nans.blacklist == False]`

---

### 25. Exclusion-zone filter reads from the wrong (unmodified) DataFrame
**File:** `spectre/analysis/cnv_metrics.py:278`

```python
df_coverage_candidate_no_excl_zone = df_coverage_candidate_no_Nans.loc[
    df_coverage_candidate.excl_zone == False]   # ← original df, all False
```

Exactly the same pattern in `__recalculate_exlustion_zone`. Line 275 correctly marks exclusion zones in `df_coverage_candidate_no_Nans`:
```python
df_coverage_candidate_no_Nans.loc[excl_zone_backlist_cnv_indices, "excl_zone"] = True
```
But line 278 filters using `df_coverage_candidate.excl_zone`, which was initialised to `False` at line 261 and never updated. The entire exclusion zone machinery — including previously-called CNV regions — is ignored. The random background sample always includes CNV loci, biasing z-scores toward zero and suppressing significance for real variants.

**Fix:** `df_coverage_candidate_no_Nans.loc[df_coverage_candidate_no_Nans.excl_zone == False]`

---

## Final Summary

| Severity | Count |
|----------|-------|
| Critical (guaranteed crash or data loss) | 7 |
| Major (wrong results, silent failures) | 9 |
| Minor / Code Quality | 7 |
| Performance | 1 |
| **Total** | **24** |

**Highest priority fixes:**
- Items 23–25 (`cnv_metrics.py`) — item 23 is a hard crash; items 24–25 mean blacklist and CNV exclusion zones are silently ignored, corrupting the z-score metric used to score every CNV candidate
- Item 14 (subtelomeric CNVs silently dropped) — systematic false-negative bias
- Item 2 (`np.nan()` crash) — guaranteed `TypeError` on any SNV-filtered run
