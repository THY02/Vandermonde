"""Reproducible high-precision sweep for the weighted Fourier-moment bound.

Run: python moment_sensitivity.py
The node perturbation is deterministic; no timing measurements are made.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys

import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
EPSILONS = ["1e-12", "1e-10", "1e-8", "1e-6", "1e-4", "1e-3", "1e-2", "1e-1"]
SIZES = [5, 9, 17]
FAMILIES = ["interior_minimum", "triangular"]
DIGITS = 80


def evaluate(L, family, epsilon, digits=DIGITS):
    with mp.workdps(digits):
        j = L // 2 if family == "interior_minimum" else 0
        if family == "interior_minimum":
            weights = [mp.mpf(2) + mp.mpf(n) / L for n in range(L)]
            weights[j] = mp.mpf(1)
        else:
            weights = [mp.mpf(min(n + 1, L - n)) for n in range(L)]
        eps = mp.mpf(epsilon)
        direction = [mp.mpf(L * (-1)**k - 1) / (L + 1) for k in range(L)]
        assert abs(mp.fsum(direction)) < mp.mpf(10)**(-digits + 5)
        assert max(map(abs, direction)) == 1
        angles = [2 * mp.pi * (k + eps * direction[k]) / L for k in range(L)]
        matrix = mp.matrix([[weights[n] * mp.exp(mp.j * n * t)
                             for t in angles] for n in range(L)])
        gram = matrix * matrix.H / L
        eig = mp.eighe(gram, eigvals_only=True)
        lo, hi = eig[0], eig[L - 1]
        a, b = min(weights), max(weights)
        eta = hi / lo / (b / a)**2 - 1
        lhs = mp.fsum(weights[n]**2 * abs(mp.fsum(
            mp.exp(mp.j * (n - j) * t) for t in angles) / L)**2
            for n in range(L) if n != j)
        tolerance = mp.mpf(10)**(-digits + 15)
        if epsilon == "0":
            assert abs(eta) < tolerance
            assert lhs < tolerance
            eta = mp.mpf(0)
        else:
            assert eta > 0
        spectral_bound = (a*a - lo) * (hi - a*a) / (a*a)
        loose_rhs = eta / (1 + eta) * (b*b * (1 + eta) - a*a)
        threshold = 1 - a*a/(b*b)
        if eta <= threshold:
            rhs = (b*b-a*a)*eta/(1+eta)
            branch = "endpoint"
        else:
            rhs = (b*b*(1+eta)-a*a)**2/(4*b*b*(1+eta))
            branch = "vertex"
        assert lhs <= spectral_bound + tolerance
        assert spectral_bound <= rhs + tolerance
        assert rhs <= loose_rhs + tolerance
        result = dict(L=L, family=family, minimum_row=j, epsilon=epsilon,
                      eta=eta, moment_lhs=lhs, spectral_variance_bound=spectral_bound,
                      moment_rhs=rhs, lhs_over_rhs=lhs/rhs if rhs > 0 else mp.mpf(0),
                      a=a,b=b,branch=branch,previous_bound=loose_rhs)
        return {key: mp.nstr(value, 55) if isinstance(value, mp.mpf) else value
                for key, value in result.items()}


def scientific_tex(value):
    value = float(value)
    coefficient, exponent = f"{value:.3e}".split("e")
    return rf"${coefficient}\times10^{{{int(exponent)}}}$"


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    rows = [evaluate(L, family, epsilon)
            for family in FAMILIES for L in SIZES for epsilon in EPSILONS]
    zeros = [evaluate(L, family, "0") for family in FAMILIES for L in SIZES]
    precision_checks = []
    with mp.workdps(120):
        for family in FAMILIES:
            original = next(r for r in rows if r["family"] == family and
                            r["L"] == 17 and r["epsilon"] == "1e-12")
            repeated = evaluate(17, family, "1e-12", 120)
            errors = {key: abs(mp.mpf(original[key]) / mp.mpf(repeated[key]) - 1)
                      for key in ["eta", "moment_lhs", "moment_rhs"]}
            assert max(errors.values()) < mp.mpf("1e-45")
            precision_checks.append(dict(family=family, L=17, epsilon="1e-12",
                maximum_relative_change=mp.nstr(max(errors.values()), 8)))
    with (DATA / "moment_sweep.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (DATA / "zero_perturbation_checks.json").write_text(
        json.dumps(zeros, indent=2), encoding="utf-8")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), layout="constrained")
    colors = {5: "#2465a8", 9: "#d26923", 17: "#228266"}
    for ax, family, title in zip(axes, FAMILIES,
            ["(a) Unique interior minimum", "(b) Triangular weights"]):
        for L in SIZES:
            group = [r for r in rows if r["family"] == family and r["L"] == L]
            ax.loglog([float(r["eta"]) for r in group],
                      [float(r["lhs_over_rhs"]) for r in group],
                      "o-", ms=4.5, lw=1.3, color=colors[L], label=f"L = {L}")
        ax.axhline(1, color="#555555", ls="--", lw=1, label="Bound (ratio = 1)")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"Condition-number excess $\eta$")
        ax.set_ylabel("Moment sum / upper bound")
        ax.grid(True, which="major", alpha=.2)
        ax.legend(fontsize=8.5, frameon=False, loc="best")
    fig.savefig(FIGURES / "moment_bound_ratio.pdf")
    fig.savefig(FIGURES / "moment_bound_ratio.png", dpi=180)
    plt.close(fig)
    table_rows = []
    for family in FAMILIES:
        label = r"内部极小行" if family == "interior_minimum" else r"三角权重"
        for L in [9, 17]:
            row = next(r for r in rows if r["family"] == family and r["L"] == L
                       and r["epsilon"] == "1e-6")
            table_rows.append(" & ".join([label, str(L), scientific_tex(row["eta"]),
                scientific_tex(row["moment_lhs"]), scientific_tex(row["moment_rhs"]),
                scientific_tex(row["lhs_over_rhs"])]) + r" \\")
    (DATA / "moment_rows.tex").write_text("\n".join(table_rows) + "\n", encoding="utf-8")
    ranges = {family: dict(min_ratio=min(float(r["lhs_over_rhs"]) for r in rows if r["family"] == family),
                          max_ratio=max(float(r["lhs_over_rhs"]) for r in rows if r["family"] == family))
              for family in FAMILIES}
    report = dict(passed=True, precision_digits=DIGITS, positive_perturbation_cases=len(rows),
                  zero_perturbation_cases=len(zeros), precision_rechecks=precision_checks,
                  ratio_ranges=ranges, sizes=SIZES, epsilon_values=EPSILONS,
                  scope="A deterministic finite sweep; no claim of sharpness or universal asymptotic rate.")
    (DATA / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (DATA / "environment.json").write_text(json.dumps(dict(
        python=sys.version, mpmath=mp.__version__, matplotlib=matplotlib.__version__,
        platform=platform.platform()), indent=2), encoding="utf-8")
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(ROOT.rglob("*")) if p.is_file()
              and p.name != "sha256_manifest.json" and "__pycache__" not in p.parts}
    (DATA / "sha256_manifest.json").write_text(json.dumps(hashes, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
