# Reproduction materials

These files accompany **Optimal conditioning and local perturbation criteria for row-weighted Vandermonde matrices on the unit circle**. No external datasets are used.

## Install and run

Use Python 3.12 (the adopted complete run used 3.12.10) and, from this directory:

```sh
python -m pip install -r requirements.txt
python verify_recorded.py
python reproduce.py --full
```

`verify_recorded.py` uses only the Python standard library. It checks both tables in the manuscript against CSVs, the timing sentence, every median against the seven raw repetitions, solver residuals, and the rational coefficient and effective-block records. The default `reproduce.py` run performs new general-weight, moment-bound, local-slope and quadratic-coefficient calculations, then runs the four-method benchmark and regenerates the three figures. Run the timing stage on an otherwise idle computer. To rerun the original computational accuracy, node-family and input-perturbation experiments as well:

```sh
python reproduce.py --full
```

`reproduce.py` copies the scripts and supplied data into a new `../runs/<UTC timestamp>/repro` directory and runs there. It preserves the supplied CSVs, manuscript figures and text, and records source hashes and completion in `run_manifest.json`. Choose another new output location with `--output-dir PATH`; existing output directories are rejected. Direct execution of an individual experiment script still replaces its own local outputs, so use the top-level runner when preserving the paper's records matters.

The supplied CSVs are the recorded runs used for the manuscript. Timing values vary with hardware, libraries and load. The original experiments use seed 20261023; structural checks and revised timings use 20261024. Directions and second-order checks are deterministic.

To regenerate all three manuscript figures from the supplied CSVs, including Figure 3, without collecting new timings:

```sh
python reproduce.py --figures-only --output-dir ../recorded_figure_rebuild
```

The three rebuilt PDFs appear in that output directory's `figures/` folder. Figure 3 reads `numerical_experiments/data/batch_summary.csv`; the same records are checked against `batch_raw.csv` by `verify_recorded.py`. The runner also reconstructs Table 2 rows from `preparation.csv`.

To rerun all mathematical checks and rebuild figures while retaining the recorded timings:

```sh
python reproduce.py --theory-only --output-dir ../theory_verification
```

## Files and manuscript locations

| Directory / script | Content | Manuscript |
|---|---|---|
| `theory_extensions/verify_extensions.py` | General weights, equality configurations, Kadec-range observations, signed-square reconstruction and higher powers | Sections 3–6 |
| `moment_sensitivity/moment_sensitivity.py` | 48 nonzero moment-bound cases, 6 exact-grid cases and precision rechecks | Figure 2 |
| `moment_sensitivity/directional_checks.py` | 45 configurations, 360 signed perturbation matrices, 80/120-digit checks | Theorem 3.6, Corollary 3.7, Figure 1 |
| `numerical_experiments/verify_general_second_order.py` | 21 rational coupled-bound checks; 22 nonconstant flat directions, 88 full matrices, 8 constant directions and two 120-digit replays | Proposition 3.4, Theorem 3.6 |
| `numerical_experiments/verify_alternating_example.py` | Exact two-by-two blocks for W=diag(1,2,1,2) and 12 full-matrix checks at 80/120 digits | General-weight example after Theorem 3.6 |
| `numerical_experiments/verify_second_order.py` | Exact rational effective blocks for degrees 2–64; 42 independent 80/120-digit matrix calculations | Proposition 3.8 |
| `numerical_experiments/benchmark.py` | FFT, cached dense transform, cached LU and fresh LU, including separate preparation measurements | Table 2, Figure 3 |
| `numerical_experiments/plot_from_data.py` | Regenerates Figures 1 and 2 directly from CSVs with distinguishable grayscale styles | Figures 1–2 |
| `mathematics/experiments.py` | Original accuracy, node-family and input-error calculations | Table 1 and Section 6.3 |

Every numerical check is finite; the manuscript supplies the proofs. Endpoint Kadec cases are observations only and do not establish any dimension-uniform endpoint theorem.

## Adopted manuscript records

All four main experiment modules use the data and environment records from the
operator's completed run `20261028T112716905285Z` (28 October 2026).
Table 2, Figure 3 and the manuscript's inline timing values use this same run.
The 224 raw online repetitions and all seven repetitions of each of the 16
preparation stages are preserved. Table 1 and all non-timing numerical results
are unchanged because they agree exactly with the earlier records.
All three manuscript figure PDFs are copied from this completed run.

The accompanying `../article.tex` snapshot and this code directory form the supplementary materials. The complete primary source and figure files are supplied separately. All data in the four main modules use the adopted run; old timing runs are not included here.

## Benchmark protocol

All four methods receive identical Fortran-contiguous complex128 right-hand sides, with `d = 8, 32, 128, 512` and batch sizes 1 and 128. NumPy and SciPy BLAS thread pools are both limited to one thread. All methods are warmed up. Seven repetitions are interleaved in seeded random order; short calls use calibrated inner loops. The reported online time is the median per batch call, excluding input generation and fixed-operator preparation.

`Fresh LU` factors the fixed rounded matrix and solves on every call. `Cached LU` uses the same SciPy factorization/solve routines and options, but reuses the factors. Neither overwrites the input. Finite-input scans are disabled consistently for both LU variants. The FFT includes coefficient normalization and final scaling. The dense transform is constructed directly from the Fourier formula.

`preparation.csv` measures weight construction, matrix construction, dense-operator construction and LU factorization separately. A cached-LU setup includes matrix construction plus factorization; these stages should be added when an end-to-end preparation cost is wanted. Online loops, each raw repetition, solver agreement/residual checks and the complete runtime/thread-pool versions are retained in `batch_raw.csv`, `solver_checks.csv` and `environment.json`.

The dense transform and FFT implement the ideal-grid formula, whereas LU solves the independently rounded explicit matrix. The saved agreement and dense-residual checks quantify that distinction.


## Environment and provenance

The adopted run used a Lenovo 12S9CTO1WW computer, Intel Core i5-13500
(reported base frequency 2.50 GHz), 8 GB installed RAM, Windows 11 build 26200,
64-bit Python 3.12.10, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2,
mpmath 1.4.1 and threadpoolctl 3.7.0. The saved NumPy and SciPy BLAS records
identify OpenBLAS 0.3.34.106.0 and 0.3.31.dev, respectively, each with one thread.
The OS-reported memory is 7.66 GiB (rounded TotalPhysicalMemory / 2**30),
not free memory measured during the benchmark. `hardware_supplement.json`
records the CPU/model/RAM information supplied by the operator after the run,
and the provenance of the threadpoolctl version. The public environment.json files redact only personal installation directory
prefixes. The numerical CSVs and software-version fields are unchanged. OpenBLAS's Haswell kernel target is not a CPU model.

All non-timing CSVs match the earlier run. Only `numerical_experiments/data/batch_summary.csv`,
`numerical_experiments/data/batch_raw.csv`, `numerical_experiments/data/preparation.csv` and
`mathematics/data/batch_timings.csv` differ. Figure 3 uses 32 online algorithm
configurations; Table 2 uses 16 preparation-stage records. The timing changes
do not change the stated distinction between small and large matrix sizes.

On Windows, use `py -3.12 -m pip install -r requirements.txt`,
`py -3.12 verify_recorded.py` and `py -3.12 reproduce.py --full`.
The equivalent `pip install -r requirements.txt`, `python verify_recorded.py`
and `python reproduce.py --full` commands work when pip and python refer to
that environment. No fresh benchmark was collected while adopting this run.
