# ELA revision checks and Figure 3

Use the parent `reproduce.py` runner to keep the supplied manuscript data unchanged:

```sh
cd ..
python verify_recorded.py
python reproduce.py --figures-only --output-dir ../recorded_figure_rebuild
python reproduce.py --theory-only --output-dir ../theory_verification
python reproduce.py --full --output-dir ../full_verification
```

Both output locations must be new. Install the pinned versions in `../requirements.txt` first. No external datasets are needed.

| Script / file | Purpose |
|---|---|
| `verify_general_second_order.py` | General effective coefficients and positive lower bounds; 88 full-matrix checks, two 120-digit replays, constant directions and exact coupled-bound checks |
| `verify_second_order.py` | Construct the effective matrices with rational arithmetic for degrees 2–64; independently compute 42 full-matrix eigenvalue cases at 80/120 digits |
| `benchmark.py` | Four-method timing, preparation costs and residual checks |
| `benchmark.py --plot-only` | Rebuild Figure 3 and Table 2 rows from existing CSVs without timing |
| `plot_from_data.py` | Rebuild Figures 1–2 from the moment and directional CSVs |
| `data/batch_raw.csv` | All 224 individual online repetitions |
| `data/batch_summary.csv` | The 32 online summaries used by Figure 3 |
| `data/preparation.csv` | All seven repetitions for each of 16 preparation stages |
| `data/solver_checks.csv` | Agreement with FFT and explicit-matrix residuals |
| `data/environment.json` | Runtime, library, hardware and thread-pool record for the published timing data |

Directly running scripts in this directory replaces their local outputs. The parent runner copies this directory before execution. Fresh elapsed times are expected to differ from the manuscript. For the complete protocol and interpretation, see `../README.md`.
