# Run and publication provenance

The adopted run is `20261028T112716905285Z`, executed with
`py -3.12 reproduce.py --full`. Its original run manifest reports
`mode: full`, `completed: true`, and `source_preserved: true`.
The current paper's Table 2, Figure 3 and inline timings use this run.
The original machine-local run archive is retained by the operator.

This public package omits old manuscript drafts, internal revision reports,
compilation logs, exploratory programs and duplicated run directories.
`previous_records/` contains the four earlier timing CSVs, their benchmark
environment records and earlier Figure 3 solely for comparison.
All other numerical CSVs of the main experiments agree with the earlier run.
`csv_comparison.json` is the earlier comparison report; its exploratory entries
are historical and do not imply those background experiments were rerun here.

Public copies of environment.json redact only the installation directory
prefix in BLAS `filepath` fields. The CPU identifier, BLAS version, kernel
target, thread count and library basename remain intact. Historical hashes
and current public-copy hashes are distinguished in publication_manifest.json.
The per-module sha256_manifest.json files were refreshed to describe this
public package; original module-manifest hashes are retained in publication_manifest.json.

CPU product name, computer model and memory were supplied by the operator
after the run, as documented in code/hardware_supplement.json. Device IDs,
product IDs and machine names are not included. No performance samples were
collected during preparation of this package.

The manuscript PDFs and LaTeX sources are unchanged from the adopted release.
The full text is included because verify_recorded.py checks manuscript values.
No repository URL, DOI, author identity, open-source license or published
release has been assigned by the packaging process. Existing third-party
template notices remain intact.
