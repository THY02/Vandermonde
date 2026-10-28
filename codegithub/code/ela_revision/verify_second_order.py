"""Exact rational block checks plus independent 80/120-digit eigenvalue checks."""
from pathlib import Path
from fractions import Fraction as Q
import csv, json, sys
import mpmath as mp

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'

def explicit(d):
    return Q(21,8) if d==2 else Q((d-1)**2,4*d)+Q(9,16)+Q((2*d-1)**2,6)

def exact_blocks(d):
    L=2*d+1; w=[min(n+1,L-n) for n in range(L)]
    def hat(r): return Q(1,2) if r%L in [2,L-2] else Q(0)
    def hat2(r): return Q(1,2) if r%L==0 else (Q(1,4) if r%L in [4,L-4] else Q(0))
    def ei(n,m): return (n-m)*w[n]*w[m]*hat(n-m)  # E = i * ei
    def f(n,m): return -Q((n-m)**2*w[n]*w[m],2)*hat2(n-m)
    def block(J,mu):
        return [[f(n,m)-sum((ei(n,q)*ei(q,m)/Q(mu-w[q]**2)
                            for q in range(L) if q not in J),Q(0))
                 for m in J] for n in J]
    top=block([d],(d+1)**2)[0][0]/(d+1)**2
    bottom=block([0,L-1],1)
    assert bottom[0][0]==bottom[1][1] and bottom[0][1]==bottom[1][0]
    coeff=(top-bottom[0][0]+abs(bottom[0][1]))/2
    assert top==Q((d-1)**2,2*d)
    assert bottom[0][0]==-Q(9,8)-Q((2*d-1)**2,3)
    assert bottom[0][1]==(-Q(7,8) if d==2 else Q(0))
    assert coeff==explicit(d)
    return dict(d=d,coefficient=str(coeff),top_normalized=str(top),
                bottom_diagonal=str(bottom[0][0]),bottom_offdiagonal=str(bottom[0][1]))

def measured(d,t,digits):
    with mp.workdps(digits):
        L=2*d+1; w=[min(n+1,L-n) for n in range(L)]
        h=[mp.cos(4*mp.pi*k/L) for k in range(L)]
        # Independent full matrix; no effective-block formula used here.
        V=mp.matrix([[w[n]*mp.exp(mp.j*n*(2*mp.pi*k/L+t*h[k]))
                      for k in range(L)] for n in range(L)])
        eig=mp.eighe(V*V.H/L,eigvals_only=True)
        observed=(mp.sqrt(eig[L-1]/eig[0])/(d+1)-1)/t**2
        C=explicit(d); expected=mp.mpf(C.numerator)/C.denominator
        error=abs(observed-expected)
        assert error < 1000*max(1,expected)*abs(t), (d,t,error)
        return dict(d=d,L=L,t=mp.nstr(t,12),digits=digits,expected=mp.nstr(expected,45),
                    observed=mp.nstr(observed,45),absolute_error=mp.nstr(error,20))

def main():
    DATA.mkdir(exist_ok=True)
    exact=[exact_blocks(d) for d in range(2,65)]
    rows=[]
    for d in [2,3,4,5,8,12]:
        for sign in [-1,1]:
            for ts in ['1e-4','1e-6','1e-8']:
                with mp.workdps(80): rows.append(measured(d,sign*mp.mpf(ts),80))
        print(f'Checked second-order coefficient d={d}',flush=True)
    for d in [2,3,8]:
        for sign in [-1,1]:
            with mp.workdps(120): rows.append(measured(d,sign*mp.mpf('1e-10'),120))
    for name,data in [('second_order_exact.csv',exact),('second_order_high_precision.csv',rows)]:
        with (DATA/name).open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    summary=dict(exact_degrees=len(exact),high_precision_matrices=len(rows),python=sys.version,
                 mpmath=mp.__version__,assertions_passed=True,
                 max_error_at_1e8=max(float(r['absolute_error']) for r in rows if abs(float(r['t']))==1e-8))
    (DATA/'second_order_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
