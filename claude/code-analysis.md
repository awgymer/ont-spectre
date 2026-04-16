# Spectre CNV Caller — Code Analysis

_Generated: 2026-04-16 | Last verified: 2026-04-16_

## Overview

Spectre is a Python-based Copy Number Variation (CNV) caller for Oxford Nanopore sequencing data. This report documents bugs and performance issues found through deep static analysis of the source code. All items have been verified against the current `master` branch.

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

### 3. Bare `except:` silencing all errors
**File:** `spectre/spectreCNVPopulation.py:55,75`

```python
except:
    self.logger.error(f"File does not exist! Provided : {file}")
```

Bare `except:` catches `KeyboardInterrupt`, `SystemExit`, and unexpected exceptions, masking real bugs. A malformed `.spc` file appears to load successfully but with missing data.

**Fix:** `except (FileNotFoundError, json.JSONDecodeError, IOError) as e:`

---

### 4. Shared `CNVCandidate` object causes all candidates in a chromosome to alias the same instance
**File:** `spectre/spectreCNVPopulation.py:83,91`

```python
def convert_dict_to_candidate_list(self, filename, candidates_dict):
    result = dict([(chrom, []) for chrom in candidates_dict.keys()])
    for chrom, candidates in candidates_dict.items():
        new_candidate = CNVCandidate()          # created ONCE per chromosome
        for candidate in candidates:
            for candidate_key, candidate_value in candidate.items():
                new_candidate.__setattr__(candidate_key, candidate_value)
            new_candidate.reinitialize_candidate_values()
            result[chrom].append(new_candidate) # same object appended every time
    return result
```

`new_candidate` is created once per chromosome but mutated and appended on every inner iteration. All entries in `result[chrom]` are references to the **same object**. The last candidate's attributes overwrite all previous ones, so every element in the list is identical to the last candidate loaded. This corrupts every population-mode run that reads `.spc` files.

**Fix:** Move `new_candidate = CNVCandidate()` inside the inner `for candidate` loop.

---

### 5. VCF INFO dict comprehension crashes on FLAG-type entries
**File:** `spectre/util/vcf_parser.py:185`

```python
info_dict = {i.split('=')[0]: i.split('=')[1] for i in cnv.INFO.split(';')}
```

Standard VCF INFO fields include FLAG-type entries that have no value (e.g. `IMPRECISE;SVTYPE=DEL;END=12345`). `"IMPRECISE".split('=')` returns a one-element list, so `[1]` raises `IndexError`. The comprehension also silently truncates values that themselves contain `=` (e.g. `CIPOS=-10,10` → value becomes `"-10,10"` but `CIEND=` with a negative would drop data the same way only if split produces exactly 2 parts — actually this is fine for standard tags; the FLAG crash is the real issue).

**Fix:**
```python
info_dict = {}
for field in cnv.INFO.split(';'):
    parts = field.split('=', 1)   # maxsplit=1 to preserve values containing '='
    info_dict[parts[0]] = parts[1] if len(parts) == 2 else True
```

---

## Major Bugs

### 6. Inverted population-mode flag
**File:** `spectre/util/vcf_parser.py:191`

```python
population_mode = len(df.columns[9:]) < 2  # TRUE more than 2 samples are present in the VCF
```

The condition and comment contradict each other. `< 2` means 0–1 samples (single-sample), so `population_mode` is `True` only in single-sample mode — the opposite of what is intended. All multi-sample VCFs are handled as single-sample.

**Fix:** `population_mode = len(df.columns[9:]) >= 2`

---

### 7. CNV runs at the end of chromosomes are silently dropped
**File:** `spectre/analysis/call_cnv_coverage.py:42–90`

The `cnv_coverage` loop emits a candidate only when the run is broken by a coverage value returning to the normal range or a gap/type change. There is no post-loop flush:

```python
for (cov, pos) in zip(...):
    if cov > upper_bound or cov < lower_bound:
        ...  # appends to cnv_pos_cand_list
    else:
        if len(cnv_pos_cand_list) > 1:
            candidate_list.append(cnv_candidates)  # only emitted on break
        ...

return candidate_list  # any in-progress run is silently discarded
```

If the final bins of a chromosome are all abnormal coverage, the accumulated `cnv_pos_cand_list` is never pushed. Subtelomeric CNVs are systematically never reported.

**Fix:** After the loop, flush any remaining run:
```python
if len(cnv_pos_cand_list) >= min_run:
    cnv_candidates = CNVCandidate(sample_origin, self.as_dev)
    cnv_candidates.push_candidates(chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list, current_cnv_type)
    candidate_list.append(cnv_candidates)
return candidate_list
```

---

### 8. First position of every chromosome's first CNV run is silently dropped
**File:** `spectre/analysis/call_cnv_coverage.py:49–51`

When a CNV run begins for the first time on a chromosome (`run_start == 0`), the code records the start position but does **not** append the current bin to the candidate lists:

```python
if run_start == 0:
    current_cnv_type = cnv_type
    run_start = pos        # records start, but does NOT append pos/cov
    # cnv_pos_cand_list stays []
```

The else branch (which does `cnv_pos_cand_list.append(pos)`) is only reached from the second bin onward. Every subsequent run restart correctly uses `cnv_pos_cand_list = [pos]`, so only the first run on the chromosome is affected. The result is that the first bin of the first CNV on each chromosome is missing, biasing start coordinates by one bin width.

**Fix:** Append the position when initialising the run:
```python
if run_start == 0:
    current_cnv_type = cnv_type
    run_start = pos
    cnv_pos_cand_list = [pos]
    cnv_cov_cand_list = [cov]
```

---

### 9. Incomplete interval overlap check misses containment
**File:** `spectre/spectreCNVPopulation.py:133`

```python
return cnv1.start <= cnv2.start <= cnv1.end or cnv1.start <= cnv2.end <= cnv1.end
```

This correctly detects when an endpoint of `cnv2` falls inside `cnv1`, but silently misses the case where `cnv2` fully contains `cnv1` (e.g. cnv1=[100,200], cnv2=[50,300]). Population mode misses valid supporting evidence for any CNV that is a strict subset of another.

**Fix:**
```python
return not (cnv2.end < cnv1.start or cnv2.start > cnv1.end)
```

---

### 10. Unfinished function that always returns `""`
**File:** `spectre/util/outputWriter.py:194–197`

```python
@staticmethod
def convert_genome_info_to_dictionary(genome_info: dict):
    tmp_genome_inf = dict([(key, []) for key in genome_info.keys()])
    return ""          # always returns empty string; body never executes
    # for info in genome_info;  ← incomplete dead code
```

The function unconditionally returns an empty string. The commented-out loop is syntactically broken. Any caller depending on the return value silently receives `""`.

---

### 11. Variable shadowing in nested loop
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

### 12. Error logged but execution continues with invalid config
**File:** `spectre/main.py:507,512,518`

```python
logger.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
```

The error is logged but the program keeps running with the invalid bin size, producing silent downstream failures.

**Fix:** Follow the log call with `sys.exit(1)`.

---

### 13. Unguarded `IndexError` from `np.where` result
**File:** `spectre/analysis/analysis.py:503`

```python
gap_start = int(np.where(cov_data_chr.positions == gap_in_candidate_start)[0][0])
```

If the position is absent from the coverage array, `np.where` returns an empty array and `[0][0]` raises `IndexError` with no diagnostic information.

**Fix:** Check `len(indices) > 0` before indexing; log a warning and `continue` if empty.

---

### 14. Potential divide-by-zero in allele frequency calculation
**File:** `spectre/analysis/call_cnv_AF.py:98`

```python
af = float(alt_count) / (float(ref_count) + float(alt_count))
```

If both counts are zero (malformed VCF entry), this raises `ZeroDivisionError`.

**Fix:** Guard with `if float(ref_count) + float(alt_count) > 0:` and skip or default otherwise.

---

### 15. Non-Pythonic `__contains__` method calls
**File:** `spectre/spectreCNVPopulation.py:51`

```python
if str(file).__contains__(".spc") or str(file).__contains__(".spc.gz"):
```

**Fix:** `if file.endswith((".spc", ".spc.gz")):`

---

### 16. Pandas `rolling(step=)` deprecation
**File:** `spectre/analysis/analysis.py:667`

```python
.rolling(self.windows_bins, step=self.windows_bins).agg(...)
```

The `step` keyword was removed in newer pandas versions. Breaks with pandas >= 2.0.

**Fix:** `rolling(self.windows_bins).agg(...)[::self.windows_bins]`

---

### 17. Inverted bounds assignment in coverage plots
**File:** `spectre/plots/plot.py:78`

The call site passes `[lower_bound, upper_bound]` (analysis.py:582), but the unpack reverses the names:

```python
[upperb, lowerb] = bounds   # upperb = lower_bound, lowerb = upper_bound
...
self.main_plot.plot(..., np.array([lowerb, lowerb]), ...)  # plots upper value as lower line
self.main_plot.plot(..., np.array([upperb, upperb]), ...)  # plots lower value as upper line
```

Every generated coverage plot shows the thresholds at the wrong levels.

**Fix:** `[lowerb, upperb] = bounds`

---

### 18. `random.sample` raises `ValueError` on small chromosomes
**File:** `spectre/analysis/cnv_metrics.py:138–140`

```python
last_index = len(df_source.index) - 1
cov_window_amount = int(last_index * fraction)   # fraction = 0.1
samples_indices = random.sample(range(0, last_index), cov_window_amount)
```

`random.sample(population, k)` raises `ValueError` when `k > len(population)`. If `last_index < 10`, `cov_window_amount` rounds to 0 and sampling from an empty or near-empty range crashes.

**Fix:**
```python
if last_index == 0:
    return []
cov_window_amount = max(1, min(int(last_index * fraction), last_index))
```

---

### 19. Unguarded list access after `filter()` in longshot VCF parsing
**File:** `spectre/util/vcf_parser.py:96`

```python
ac = list(filter(lambda x: "AC" in x, info_.split(";")))[0].split("=")[1]
```

If the filter returns an empty list, `[0]` raises `IndexError`. The `"AC" in x` check also falsely matches tags like `MLEAC=` or `SVAC=`.

**Fix:**
```python
ac_matches = list(filter(lambda x: x.startswith("AC="), info_.split(";")))
if not ac_matches:
    af = "NA"
else:
    ac = ac_matches[0].split("=")[1]
```

---

### 20. Divide-by-zero in metadata GC% calculation
**File:** `spectre/util/metadata/metadataCollector.py:99`

```python
str((self.__gCnt + self.__cCnt) / self.__nTo * 100)
```

If the input FASTA is empty or all-N, `__nTo` is 0, raising `ZeroDivisionError`.

**Fix:** `str(round((self.__gCnt + self.__cCnt) / self.__nTo * 100, 4) if self.__nTo > 0 else 0.0)`

---

### 21. `pop('logger')` without a default raises `KeyError`
**File:** `spectre/util/outputWriter.py:211`

```python
tmp_dict.pop('logger')
```

If a candidate object lacks a `logger` attribute, `pop()` raises `KeyError`.

**Fix:** `tmp_dict.pop('logger', None)`

---

### 22. CNV type `KeyError` when plotting unknown types
**File:** `spectre/plots/plot.py:72`

```python
cnv_color = self.cnv_color[cnv.type]
```

`self.cnv_color` only has `"DUP"` and `"DEL"` keys. Any other type string raises `KeyError` and aborts plotting for the entire chromosome.

**Fix:** `cnv_color = self.cnv_color.get(cnv.type, "#888888")`

---

## Performance Issues

### 23. O(n^6) algorithm in population mode overlap detection
**File:** `spectre/spectreCNVPopulation.py:161–189`

The cross-sample candidate comparison uses six nested loops, including self-comparisons and redundant cross-chromosome comparisons:

```
for sample1 in all_samples:
  for sample2 in all_samples:        # includes self-comparisons
    for chrom1 in sample1.chroms:
      for chrom2 in sample2.chroms:  # note: cross-chromosome guarded but still iterated
        for cnv1 in chrom1.candidates:
          for cnv2 in chrom2.candidates:
```

With 10 samples × 500 CNVs × 24 chromosomes, this approaches billions of comparisons.

**Fixes (in priority order):**
1. Skip self-comparisons: `if samples_key1 >= samples_key2: continue`
2. Replace the inner linear scan with an interval tree (`ncls` or `intervaltree`) for O(n log n) overlap queries
3. Compute pairwise comparisons once (symmetric)

---

## Issue #9 Follow-up: `cnv_metrics.py` Cluster

_Reported upstream: nanoporetech/ont-spectre#9_

Three related bugs in `spectre/analysis/cnv_metrics.py`, all from the same root cause: after `dropna()` creates a subset DataFrame, the code writes to the subset but reads from the original, so exclusion zones have no effect.

---

### 24. `.iloc` with label-based indices → `IndexError` crash
**File:** `spectre/analysis/cnv_metrics.py:229–230`

```python
df_coverage_candidate_no_Nans.iloc[
    excl_zone_blacklist_indices,
    df_coverage_candidate_no_Nans.columns.get_loc("blacklist")] = True
```

`get_exclusion_zones_indices_in_coverage_data` returns **label** indices (via `.index`). After `dropna()`, positional and label indices diverge. `.iloc` expects positional integers; passing a large label raises:

```
IndexError: index 2713038 is out of bounds for axis 0 with size 2703072
```

The correct fix was already written and commented out on line 227.

**Fix:** `df_coverage_candidate_no_Nans.loc[excl_zone_blacklist_indices, "blacklist"] = True`

_Note: Fixed on branch `fix/iloc-label-index-cnv-metrics`; not yet merged._

---

### 25. Blacklist filter reads from the wrong (unmodified) DataFrame
**File:** `spectre/analysis/cnv_metrics.py:233`

```python
df_coverage_candidate_no_blacklist = df_coverage_candidate_no_Nans.loc[
    df_coverage_candidate.blacklist == False]   # ← original df, never modified
```

All modifications to the `blacklist` column land on `df_coverage_candidate_no_Nans`. The mask `df_coverage_candidate.blacklist == False` is always all-`True`, so the blacklist never excludes anything from the random sample pool used for z-score scoring.

**Fix:** `df_coverage_candidate_no_Nans.loc[df_coverage_candidate_no_Nans.blacklist == False]`

---

### 26. Exclusion-zone filter reads from the wrong (unmodified) DataFrame
**File:** `spectre/analysis/cnv_metrics.py:278`

```python
df_coverage_candidate_no_excl_zone = df_coverage_candidate_no_Nans.loc[
    df_coverage_candidate.excl_zone == False]   # ← original df, all False
```

Same pattern in `__recalculate_exlustion_zone`. Line 275 correctly marks exclusion zones in `df_coverage_candidate_no_Nans`, but line 278 filters using the original DataFrame where all values are still `False`. CNV and blacklist regions are always included in the background sample, biasing z-scores toward zero and suppressing significance for real variants.

**Fix:** `df_coverage_candidate_no_Nans.loc[df_coverage_candidate_no_Nans.excl_zone == False]`

---

## Retracted Item

### ~~3. Wrong argument passed to `id_generator`~~ — NOT A BUG
**File:** `spectre/util/cnv_id.py:33`

Previously reported: `cls.id_generator(n, size, chars)` was claimed to pass the loop counter `n` as the `size` parameter.

On re-examination: calling an instance method via the class in Python 3 passes `n` as `self`, and `size`/`chars` as the correct second and third parameters. The method body never uses `self`, so the IDs are generated with the correct length. The calling convention is unconventional but functionally correct.

---

## Final Summary

| Severity | Count |
|----------|-------|
| Critical (guaranteed crash or data loss) | 5 |
| Major (wrong results, systematic bias) | 6 |
| Minor / Code Quality | 11 |
| Performance | 1 |
| Retracted | 1 |
| **Total confirmed bugs** | **23** |

**Highest priority fixes:**
1. Item 4 — shared `CNVCandidate` object corrupts every population-mode `.spc` load
2. Item 24 — hard `IndexError` crash during CNV evaluation (issue #9; fix branch exists)
3. Items 25–26 — blacklist and CNV exclusion zones silently ignored, corrupting z-scores
4. Item 2 — guaranteed `TypeError` on any SNV-filtered run
5. Item 7 — subtelomeric CNVs systematically never reported
