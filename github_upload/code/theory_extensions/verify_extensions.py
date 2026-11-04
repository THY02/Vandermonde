"""Numerical falsification checks for manuscript_2, independent of old outputs.

Run: python verify_extensions.py
Dependencies: numpy, mpmath, threadpoolctl. All data are generated locally.
These finite checks complement the proofs; they do not prove uniqueness.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import numpy as np
import mpmath as mp
from threadpoolctl import threadpool_limits, threadpool_info

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
DATA.mkdir(exist_ok=True)
SEED = 20261024
rng = np.random.default_rng(SEED)
TOL = 2e-10


def save_csv(name, rows):
    with (DATA / name).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check_nodes(w, theta, j):
    L = len(w)
    V = np.exp(1j * np.arange(L)[:, None] * theta[None, :])
    M = w[:, None] * V
    singular = np.linalg.svd(M, compute_uv=False)
    cond = float(singular[0] / singular[-1])
    a, b = np.min(abs(w)), np.max(abs(w))
    rho = b / a
    raw_eta = (cond / rho)**2 - 1
    assert raw_eta >= -TOL
    eta = max(0., raw_eta)
    moments = np.exp(1j * (np.arange(L) - j)[:, None] * theta).mean(axis=1)
    terms = abs(w)**2 * abs(moments)**2
    terms[j] = 0
    lhs = float(terms.sum())
    rhs = float(eta / (1 + eta) * (b*b*(1+eta) - a*a))
    assert lhs <= rhs + TOL * b*b, (L, j, lhs, rhs)
    return cond, rho, eta, lhs, rhs, singular


def weighted_checks():
    rows = []
    for L in [4, 5, 8, 9, 16, 17]:
        for j in range(L):
            mag = rng.uniform(1.2, 3., L)
            mag[j] = 1.
            w = mag * np.exp(1j*rng.uniform(-np.pi, np.pi, L))
            base = .37 + 2*np.pi*np.arange(L)/L
            for kind, theta in [('roots', base), ('jitter', base + rng.uniform(-.4, .4, L)*np.pi/L)]:
                cond, rho, eta, lhs, rhs, ss = check_nodes(w, theta, j)
                if kind == 'roots':
                    assert np.allclose(ss, np.sort(abs(w))[::-1]*np.sqrt(L), atol=TOL, rtol=TOL)
                    assert abs(cond/rho-1) < TOL
                rows.append(dict(L=L, minimum_row=j, kind=kind, condition=cond,
                                 rho=rho, eta=eta, moment_lhs=lhs, moment_rhs=rhs))
        w = np.exp(1j*rng.uniform(-np.pi, np.pi, L))
        for kind, theta in [('equal_modulus_roots', base), ('equal_modulus_jitter', base+rng.uniform(-.4,.4,L)*np.pi/L)]:
            j = int(np.argmin(abs(w)))
            cond, rho, eta, lhs, rhs, _ = check_nodes(w, theta, j)
            rows.append(dict(L=L, minimum_row=j, kind=kind, condition=cond,
                             rho=rho, eta=eta, moment_lhs=lhs, moment_rhs=rhs))
    save_csv('general_weights.csv', rows)
    return rows


def perturbation_checks():
    rows = []
    for d in [2, 8, 32, 64]:
        L = 2*d+1
        w = np.minimum(np.arange(L)+1, L-np.arange(L)).astype(float)
        base = 2*np.pi*np.arange(L)/L
        for alpha in [0., .10, .25, .40, .49, .50]:
            phi = 1 - np.cos(np.pi*alpha/2) + np.sin(np.pi*alpha/2)
            bound = (1+phi)/(1-phi) if alpha < .5 else None
            for trial in range(10):
                delta = rng.uniform(-alpha, alpha, L)*np.pi/L
                cond, rho, eta, lhs, rhs, _ = check_nodes(w, base+delta, 0)
                if bound is not None:
                    assert cond/rho <= bound + TOL
                rows.append(dict(d=d,L=L,alpha=alpha,trial=trial,actual_alpha=float(np.max(abs(delta))*L/np.pi),
                                 condition=cond,ratio=cond/rho,kadec_ratio_bound=bound,
                                 eta=eta,moment_lhs=lhs,moment_rhs=rhs))
    save_csv('node_perturbations.csv', rows)
    lines = []
    for alpha in [0., .10, .25, .40, .49, .50]:
        group = [r for r in rows if r['alpha']==alpha]
        bound = group[0]['kadec_ratio_bound']
        value = f'{bound:.4f}' if bound is not None else r'\text{不适用}'
        lines.append(f"{alpha:.2f} & {max(r['ratio'] for r in group):.4f} & {value}" + r' \\')
    (DATA/'perturbation_rows.tex').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return rows


def real_checks():
    rows=[]
    for d in [1,2,8,32]:
        L=2*d+1
        w=np.minimum(np.arange(L)+1,L-np.arange(L)).astype(float)
        inputs=[('zero',np.zeros(L)),('negative_root_only',-w.copy()),('positive_root_only',w.copy())]
        inputs += [(f'sparse_{j}',np.eye(L)[j]) for j in [0,d,L-1]]
        inputs += [(f'random_{j}',rng.normal(size=L)) for j in range(10)]
        for kind,A in inputs:
            half=np.fft.rfft(A/w)/L
            full=np.concatenate([half,np.conj(half[1:][::-1])])
            assert np.allclose(full,np.fft.fft(A/w)/L,atol=TOL,rtol=TOL)
            q0=np.ones(d+1)
            out=half[0].real*np.convolve(q0,q0)
            for k in range(1,d+1):
                coeff=np.sqrt(complex(half[k]))*np.exp(2j*np.pi*k*np.arange(d+1)/L)
                U,T=coeff.real,coeff.imag
                out += 2*(np.convolve(U,U)-np.convolve(T,T))
            error=float(np.linalg.norm(out-A)/max(1.,np.linalg.norm(A)))
            assert error < TOL, (d,kind,error)
            rows.append(dict(d=d,kind=kind,scaled_coefficient_error=error))
    save_csv('real_signed_squares.csv',rows)
    return rows


def high_precision_checks():
    """Resolve the moment inequality near equality without double subtraction."""
    rows=[]
    with mp.workdps(80):
        for L in [4,5,8,9]:
            j=L//2
            for eps_string in ['0','1e-12','1e-6','0.1']:
                eps=mp.mpf(eps_string)
                w=[mp.mpf(2)+mp.mpf(n)/L for n in range(L)]
                w[j]=mp.mpf(1)
                theta=[2*mp.pi*k/L+eps*(-1)**k/L for k in range(L)]
                M=mp.matrix([[w[n]*mp.exp(mp.j*n*t) for t in theta] for n in range(L)])
                eig=mp.eighe(M*M.H,eigvals_only=True)
                rho=max(w)/min(w)
                raw_eta=eig[L-1]/eig[0]/rho**2-1
                assert raw_eta > -mp.mpf('1e-65')
                eta=max(mp.mpf(0),raw_eta)
                lhs=mp.fsum(w[n]**2*abs(mp.fsum(mp.exp(mp.j*(n-j)*t) for t in theta)/L)**2 for n in range(L) if n!=j)
                rhs=eta/(1+eta)*(max(w)**2*(1+eta)-1)
                assert lhs <= rhs + mp.mpf('1e-65')
                rows.append(dict(L=L,minimum_row=j,epsilon=eps_string,eta=mp.nstr(eta,30),
                                 moment_lhs=mp.nstr(lhs,30),moment_rhs=mp.nstr(rhs,30)))
    save_csv('high_precision_moments.csv',rows)
    return rows


def exact_power_checks():
    from fractions import Fraction
    def mul(x,y):
        z=[0]*(len(x)+len(y)-1)
        for i,a in enumerate(x):
            for j,b in enumerate(y): z[i+j]+=a*b
        return z
    count=0
    for d in range(1,5):
        for r in range(1,5):
            weights=[1]
            for _ in range(r): weights=mul(weights,[1]*(d+1))
            for beta in [Fraction(-2),Fraction(0),Fraction(1,3),Fraction(2)]:
                coeff=[1]
                for _ in range(r): coeff=mul(coeff,[beta**j for j in range(d+1)])
                assert coeff==[c*beta**n for n,c in enumerate(weights)]
                assert min(weights)>0
                count+=1
    # Zero coefficient invalidates the full-space basis even for distinct roots.
    V=np.exp(2j*np.pi*np.arange(5)[:,None]*np.arange(5)[None,:]/5)
    assert np.linalg.matrix_rank(np.array([1,2,0,2,1])[:,None]*V)==4
    return count


def main():
    with threadpool_limits(limits=1):
        general=weighted_checks()
        perturb=perturbation_checks()
        real=real_checks()
        precise=high_precision_checks()
        exact=exact_power_checks()
        env=dict(python=sys.version,numpy=np.__version__,mpmath=mp.__version__,platform=platform.platform(),
                 seed=SEED,threadpools=threadpool_info())
    report=dict(passed=True,general_weight_cases=len(general),node_perturbation_cases=len(perturb),
                strict_kadec_cases=sum(r['alpha']<.5 for r in perturb),endpoint_observations=sum(r['alpha']==.5 for r in perturb),
                real_reconstruction_cases=len(real),high_precision_moment_cases=len(precise),
                exact_power_expansions=exact,max_real_scaled_error=max(r['scaled_coefficient_error'] for r in real),
                double_check_tolerance=TOL,high_precision_digits=80,
                scope='Finite falsification checks complement proofs; no novelty certification.')
    (DATA/'environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    (DATA/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.iterdir()) if p.name!='sha256_manifest.json'}
    hashes['verify_extensions.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (DATA/'sha256_manifest.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
