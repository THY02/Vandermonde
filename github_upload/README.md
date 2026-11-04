# Code and numerical results

Python programs and numerical records accompanying *Optimal conditioning and local perturbation criteria for row-weighted Vandermonde matrices on the unit circle*.

Use Python 3.12 (recorded run: 3.12.10). From this directory:

```sh
cd code
python -m pip install -r requirements.txt
python verify_recorded.py
python reproduce.py --full
```

On Windows, `py -3.12` can replace `python`. New outputs go to `runs/<UTC timestamp>/` beside `code/`; supplied files are preserved. `--output-dir PATH` must name a new directory outside `code/`. Run benchmarks on an otherwise idle computer. To rebuild all three figures without new timings, use `python reproduce.py --figures-only`.

All main numerical records come from the complete run `20261028T112716905285Z`. Raw timing repetitions and preparation measurements are included. Environment files retain the software and thread-pool information; only personal installation path prefixes are redacted. Hardware details supplied after the run are recorded in `code/hardware_supplement.json`.

`article.tex` is a reference manuscript snapshot used to check recorded tables and inline values. The final manuscript is supplied separately and may differ in editorial text. The included article.tex is a manuscript snapshot used by the verification script; it is not a standalone LaTeX compilation package. No external datasets are required.

The four methods receive identical complex128 inputs, use one BLAS thread, are warmed up and are timed in seven repetitions. Online timings are medians and exclude fixed-operator preparation; Fresh LU includes factorization at each call. Preparation costs are recorded separately. Small matrix sizes can favor cached dense multiplication over FFT.

## Additional structural checks

The structural checks comprise 130 general-weight configurations (`L=4,5,8,9,16,17`,
minimum weight in every row, complex phase and equal-modulus weights) and 240
triangular-weight node perturbations, of which 200 satisfy the strict range of the Kadec bound in the manuscript and 40 are endpoint observations at `alpha=1/2`; all pass the
applicable checks. The real-square construction is verified on 64 coefficient-convolution
reconstructions with maximum normalized error about `6.10e-14`. The remaining
node-family, input-perturbation and higher-power checks are kept in the reproduction
materials.

The moment_sensitivity figure PDFs are the original recorded exports. The separately supplied manuscript uses later figure exports; the underlying recorded data are retained here.
