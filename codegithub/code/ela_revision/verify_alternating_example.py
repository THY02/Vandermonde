"""Exact effective blocks and full-matrix checks for W=diag(1,2,1,2)."""
from pathlib import Path
from fractions import Fraction as Q
import csv,json,sys,platform
import mpmath as mp

ROOT=Path(__file__).resolve().parent
def exact_block(J,mu):
    w=[1,2,1,2]
    # h=(1,0,-1,0): hat(h)_r=1/2 for odd r;
    # hat(h^2)_r=1/2 for even r. E is i times the rational matrix below.
    e=lambda n,m: Q((n-m)*w[n]*w[m],2) if (n-m)%2 else Q(0)
    f=lambda n,m: -Q((n-m)**2*w[n]*w[m],4) if (n-m)%2==0 else Q(0)
    return [[f(n,m)-sum((e(n,q)*e(q,m)/Q(mu-w[q]**2)
                        for q in range(4) if q not in J),Q(0)) for m in J] for n in J]

def main():
    bottom=exact_block([0,2],1); top=exact_block([1,3],4)
    assert bottom==[[Q(-10,3),Q(-5,3)],[Q(-5,3),Q(-2,3)]]
    assert top==[[Q(2,3),Q(-10,3)],[Q(-10,3),Q(10,3)]]
    assert sum(bottom[i][i] for i in range(2))==-4
    assert bottom[0][0]*bottom[1][1]-bottom[0][1]**2==Q(-5,9)
    assert sum(top[i][i] for i in range(2))==4
    assert top[0][0]*top[1][1]-top[0][1]**2==Q(-80,9)
    rows=[]
    for digits,times in [(80,['1e-2','1e-4','1e-6','1e-8']),(120,['1e-8','1e-10'])]:
        with mp.workdps(digits):
            c=mp.mpf(5)/4+mp.sqrt(29)/12+mp.sqrt(41)/6
            for sign in [-1,1]:
                for ts in times:
                    t=sign*mp.mpf(ts); h=[1,0,-1,0]; w=[1,2,1,2]
                    A=mp.matrix([[w[n]*mp.exp(mp.j*n*(2*mp.pi*k/4+t*h[k]))
                                  for k in range(4)] for n in range(4)])
                    ev=mp.eighe(A*A.H/4,eigvals_only=True)
                    observed=(mp.sqrt(ev[3]/ev[0])/2-1)/(t*t)
                    err=abs(observed-c)
                    assert err<100*t*t
                    assert observed>mp.mpf(5)/6
                    rows.append(dict(digits=digits,t=mp.nstr(t,15),coefficient=mp.nstr(c,60),
                                     observed=mp.nstr(observed,60),absolute_error=mp.nstr(err,25)))
    out=ROOT/'data'; out.mkdir(exist_ok=True)
    with (out/'alternating_example.csv').open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    report=dict(passed=True,L=4,weights=[1,2,1,2],direction=[1,0,-1,0],
        bottom_block=[[str(x) for x in r] for r in bottom],top_block=[[str(x) for x in r] for r in top],
        coefficient='5/4+sqrt(29)/12+sqrt(41)/6',coefficient_decimal=rows[0]['coefficient'],
        diagonal_lower_bound='25/12',moment_lower_bound='5/3',finite_t_lower_bound='5/6',
        full_matrix_checks=len(rows),precisions=[80,120],python=sys.version,mpmath=mp.__version__,
        platform=platform.platform())
    (out/'alternating_example_summary.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
