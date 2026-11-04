"""High-precision checks of the local Fourier criterion and one-sided slopes."""
from pathlib import Path
import csv
import json
import sys
import platform
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
DIGITS = 80


def direction(L, kind):
    if kind == "rotation":
        return [mp.mpf(1)] * L
    if kind == "projected":
        raw = [mp.mpf((-1)**k) for k in range(L)]
        coefficient = mp.fsum(raw[k]*mp.exp(2j*mp.pi*k/L) for k in range(L))/L
        return [raw[k]-2*mp.re(coefficient*mp.exp(-2j*mp.pi*k/L)) for k in range(L)]
    mode = int(kind[3:])
    return [mp.cos(2*mp.pi*mode*k/L) for k in range(L)]


def weight_vector(L, family):
    if family == "triangular":
        magnitude = [mp.mpf(min(n+1, L-n)) for n in range(L)]
    elif family == "unique_extrema":
        magnitude = [2+mp.mpf(n)/L for n in range(L)]
        magnitude[L//2] = mp.mpf(1)
    elif family == "repeated_extrema":
        magnitude = [mp.mpf(2)]*L
        for n in [1, 3]: magnitude[n] = mp.mpf(1)
        for n in [0, L-1]: magnitude[n] = mp.mpf(3)
    else:
        magnitude = [mp.mpf(1)]*L
    # Complex phases test the claimed invariance as well.
    return [magnitude[n]*mp.exp(mp.j*n/7) for n in range(L)], magnitude


def slopes(L, magnitude, h):
    a, b = min(magnitude), max(magnitude)
    indices = [[n for n in range(L) if magnitude[n] == extreme] for extreme in [a, b]]
    blocks = []
    for group in indices:
        block = mp.matrix([[mp.j*(n-m)*mp.fsum(h[k]*mp.exp(2j*mp.pi*(n-m)*k/L)
            for k in range(L))/L if n != m else mp.mpf(0)
            for m in group] for n in group])
        blocks.append(mp.eighe(block, eigvals_only=True))
    lower, upper = blocks
    plus = (upper[len(upper)-1] - lower[0])/2
    minus = (lower[len(lower)-1] - upper[0])/2
    return plus, minus


def excess(L, weights, magnitude, h, t):
    angles = [2*mp.pi*k/L+t*h[k] for k in range(L)]
    matrix = mp.matrix([[weights[n]*mp.exp(mp.j*n*angle) for angle in angles] for n in range(L)])
    eig = mp.eighe(matrix*matrix.H/L, eigvals_only=True)
    return mp.sqrt(eig[L-1]/eig[0])/(max(magnitude)/min(magnitude))-1


def main():
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    rows, cases, summary_rows = [], [], []
    with mp.workdps(DIGITS):
        configurations = [(L,"triangular") for L in [3,5,9,17]]
        configurations += [(L,"unique_extrema") for L in [5,8,9]]
        configurations += [(L,"repeated_extrema") for L in [8,9]]
        configurations += [(5,"equal_modulus")]
        for L,family in configurations:
            kinds = ["rotation", "cos1", "projected"]
            if L >= 5: kinds.append("cos2")
            if L >= 8: kinds.append("cos3")
            for kind in kinds:
                h = direction(L,kind)
                weights, magnitude = weight_vector(L,family)
                plus, minus = slopes(L,magnitude,h)
                flat = max(abs(plus),abs(minus)) < mp.mpf("1e-65")
                rotation = max(h)-min(h) < mp.mpf("1e-65")
                if family == "triangular":
                    hat1 = mp.fsum(h[k]*mp.exp(2j*mp.pi*k/L) for k in range(L))/L
                    expected = (L-1)/2*abs(hat1)
                    assert abs(plus-expected) < mp.mpf("1e-65")
                    assert abs(minus-expected) < mp.mpf("1e-65")
                all_scaled = []
                for sign in [-1,1]:
                    predicted = plus if sign == 1 else minus
                    for t_text in ["1e-2","1e-4","1e-6","1e-8"]:
                        t_abs = mp.mpf(t_text)
                        actual = excess(L,weights,magnitude,h,sign*t_abs)
                        assert actual > -mp.mpf("1e-65")
                        actual = max(mp.mpf(0),actual)
                        slope = actual/t_abs
                        if rotation:
                            assert actual < mp.mpf("1e-65")
                        if t_text == "1e-8":
                            assert abs(slope-predicted) < mp.mpf("2e-5")*max(1,abs(predicted))
                        scaled = actual/(t_abs*t_abs) if flat else actual/t_abs
                        all_scaled.append((sign,t_text,scaled))
                        rows.append(dict(L=L,family=family,direction=kind,sign=sign,t=t_text,
                            predicted_slope=mp.nstr(predicted,45),condition_excess=mp.nstr(actual,45),
                            measured_slope=mp.nstr(slope,45),quadratic_scale=mp.nstr(actual/t_abs**2,45),
                            flat=flat,rotation=rotation))
                if not rotation:
                    for sign in [-1,1]:
                        v6 = next(v for s,t,v in all_scaled if s==sign and t=="1e-6")
                        v8 = next(v for s,t,v in all_scaled if s==sign and t=="1e-8")
                        assert abs(v6/v8-1) < mp.mpf("0.002"), (L,family,kind,sign)
                cases.append(dict(L=L,family=family,direction=kind,flat=flat,rotation=rotation))
        # Independent higher-precision replay of one linear and one quadratic case.
        replay = []
        for kind in ["cos1","cos2"]:
            old = next(r for r in rows if r["L"]==17 and r["family"]=="triangular"
                       and r["direction"]==kind and r["sign"]==1 and r["t"]=="1e-8")
            with mp.workdps(120):
                w,mag = weight_vector(17,"triangular")
                new = excess(17,w,mag,direction(17,kind),mp.mpf("1e-8"))
                change = abs(mp.mpf(old["condition_excess"])/new-1)
                assert change < mp.mpf("1e-40")
                replay.append(dict(L=17,direction=kind,relative_change=mp.nstr(change,8)))
    with (DATA/"directional_sweep.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,ax=plt.subplots(figsize=(7.2,3.8),layout="constrained")
    for kind,label,color in [("cos1","First harmonic: linear","#2465a8"),
                             ("cos2","Second harmonic: quadratic","#d26923"),
                             ("projected","Projected direction: quadratic","#228266")]:
        selected=[r for r in rows if r["L"]==17 and r["family"]=="triangular"
                  and r["direction"]==kind and r["sign"]==1]
        ax.loglog([float(r["t"]) for r in selected],[float(r["condition_excess"]) for r in selected],
                  "o-",color=color,label=label,ms=5)
    times=[1e-8,1e-2]
    ax.loglog(times,[4*t for t in times],"--",color="#2465a8",alpha=.6,label=r"Predicted $4|t|$ ($d=8$)")
    ax.set_xlabel(r"Perturbation parameter $|t|$")
    ax.set_ylabel(r"$\kappa_2(DV(t))/(d+1)-1$")
    ax.legend(frameon=False,fontsize=9)
    ax.grid(True,which="major",alpha=.2)
    fig.savefig(FIGURES/"directional_orders.pdf")
    fig.savefig(FIGURES/"directional_orders.png",dpi=180)
    plt.close(fig)
    report=dict(passed=True,precision_digits=DIGITS,direction_configurations=len(cases),
        evaluated_matrices=len(rows),flat_directions=sum(c["flat"] for c in cases),
        rotation_directions=sum(c["rotation"] for c in cases),precision_rechecks=replay,
        tested_families=sorted(set(c["family"] for c in cases)),
        scope="Finite checks of one-sided slopes and convergence orders; proof is in the manuscript.")
    (DATA/"directional_validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
