# Weighted Vandermonde matrices: code and numerical records

Reproduction materials for **Optimal conditioning and local perturbation criteria
for row-weighted Vandermonde matrices on the unit circle**.

The manuscript uses the complete run `20261028T112716905285Z` (28 October 2026).
The recorded hardware is a Lenovo 12S9CTO1WW with Intel Core i5-13500 and 8 GB
installed RAM. The run used Windows 11 and 64-bit Python 3.12.10.
Dependencies are pinned in [code/requirements.txt](code/requirements.txt).
No external datasets are required.

## Quick start

From the repository root, using Python 3.12:

```sh
cd code
python -m pip install -r requirements.txt
python verify_recorded.py
python reproduce.py --full
```

On Windows, replace `python` with `py -3.12` when selecting that interpreter.
`verify_recorded.py` uses only the standard library and checks the recorded
tables and timing sentences in both manuscripts, raw medians and mathematical
records. A full run also recomputes the mathematical experiments and timings;
run benchmarks on an otherwise idle computer.

To rebuild the three figures from the supplied CSVs without new measurements:

```sh
python reproduce.py --figures-only
```

New outputs are written to `runs/<UTC timestamp>/` at the repository root.
The supplied manuscript, CSVs and figures are preserved. `--output-dir PATH`
must name a new directory outside `code/`. Direct execution of individual
experiment scripts can replace their own outputs; use the top-level runner.
Fresh timings need not equal recorded timings.

## Contents

- [code/](code/): scripts, fixed dependencies, adopted CSVs, environment records
  and the operator-supplied hardware supplement.
- [paper/](paper/): English and Chinese LaTeX/PDF manuscripts and their figures.
  Keep this folder beside `code/`: the verification script reads the manuscripts.
- [reproducibility/](reproducibility/): provenance, file hashes, earlier timing
  records for comparison, and verification results.

See [the experiment README](code/README.md) for the benchmark protocol and
[paper/README.md](paper/README.md) for LaTeX compilation commands.
Table 2 reports preparation times; Figure 3 uses 32 online configurations and
224 raw repetitions. Online timings exclude fixed-operator preparation, except
that Fresh LU includes factorization on every call. BLAS uses one thread.

## Provenance and citation

The adopted data come from the operator's completed full run, not from a new
benchmark during packaging. Only local installation directory prefixes in
public environment records were redacted. Numerical CSVs, manuscript source
and figure PDFs are unchanged. See [PROVENANCE.md](reproducibility/PROVENANCE.md).

When citing these materials, identify the repository URL and the exact release
or commit used, together with the manuscript title. Repository URL, author
metadata and release identifiers have not been invented in this prepared copy.

## 中文说明

本仓库使用2026年10月28日完整运行的数据，中英文论文、表格、图片和CSV相互对应。
请保留 `code/`、`paper/`、`reproducibility/` 的相对目录结构。
在 `code/` 中依次运行安装依赖、`verify_recorded.py` 和 `reproduce.py --full`。
新计算结果写入仓库根目录下独立的 `runs/`，不会覆盖本文记录。
公开环境文件只隐去了本机安装路径前缀，没有改动CPU、版本、线程数或数值数据。
作者信息仍由作者后续填写。
