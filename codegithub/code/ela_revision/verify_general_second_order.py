"""Check the coupled moment bound and general flat-direction coefficient.

No timings. Uses exact rational optimization checks and independent 80-digit
full Gram eigenvalues, with two selected 120-digit replays. Requires the supplied
directional_sweep.csv, or run directional_checks.py first. Results are finite
verification evidence, not substitutes for the proofs.
"""
from pathlib import Path
from fractions import Fraction as Q
import csv, json, platform, sys
import mpmath as mp

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'moment_sensitivity'))
from directional_checks import direction, weight_vector, excess

def coupled(A,B,eta):
    if eta<=1-A/B: return (B-A)*eta/(1+eta)
    return (B*(1+eta)-A)**2/(4*B*(1+eta))

def rational_checks():
    checked=0
    for A,B in [(Q(1),Q(1)),(Q(1),Q(4)),(Q(4),Q(9))]:
        for eta in [Q(0),Q(1,100),Q(1,4),Q(3,4),Q(1),Q(2),Q(10)]:
            lo=1/(1+eta); vertex=(B*(1+eta)+A)/(2*B*(1+eta))
            optimum=min(Q(1),max(lo,vertex))
            f=lambda x:(1-x)*(B*(1+eta)*x-A)
            value=coupled(A,B,eta)
            assert f(optimum)==value
            assert value<=eta*(B*(1+eta)-A)/(1+eta)
            for i in range(101): assert f(lo+(1-lo)*Q(i,100))<=value
            if A==B: assert value==A*eta*eta/(4*(1+eta))
            checked+=1
    # A real two-node Vandermonde example catches the incorrect plus sign in
    # the proposed transition: beta=(-19+10*i*sqrt(62))/81, weights=(1,2).
    assert Q(19**2+100*62,81**2)==1
    assert Q(5,9)*Q(40,9)==4-Q(124,81)
    assert Q(40,9)/Q(5,9)/4-1==1
    assert Q(3,2)<Q(124,81)<=coupled(Q(1),Q(4),Q(1))==Q(49,32)
    return checked

def coefficient(L,mag,h):
    a,b=min(mag),max(mag)
    assert b>a
    lower=[n for n in range(L) if mag[n]==a]
    upper=[n for n in range(L) if mag[n]==b]
    hats={r:mp.fsum(h[k]*mp.exp(2*mp.pi*mp.j*r*k/L) for k in range(L))/L
          for r in range(1-L,L)}
    squares={r:mp.fsum(h[k]**2*mp.exp(2*mp.pi*mp.j*r*k/L) for k in range(L))/L
             for r in range(1-L,L)}
    E=mp.matrix([[mp.j*(n-m)*mag[n]*mag[m]*hats[n-m] for m in range(L)] for n in range(L)])
    F=mp.matrix([[-mp.mpf('.5')*(n-m)**2*mag[n]*mag[m]*squares[n-m]
                  for m in range(L)] for n in range(L)])
    def block(J,mu):
        assert max(abs(E[n,m]) for n in J for m in J)<mp.mpf('1e-65')
        K=mp.matrix([[F[n,m]+mp.fsum(E[n,q]*E[q,m]/(mu-mag[q]**2)
                     for q in range(L) if q not in J) for m in J] for n in J])
        assert mp.norm(K-K.H)<mp.mpf('1e-65')
        return K
    low,high=block(lower,a*a),block(upper,b*b)
    c2=(mp.eighe(high,eigvals_only=True)[len(upper)-1]/(b*b)
        -mp.eighe(low,eigvals_only=True)[0]/(a*a))/2
    diagonal=(max(mp.re(high[i,i])/(b*b) for i in range(len(upper)))
              -min(mp.re(low[j,j])/(a*a) for j in range(len(lower))))/2
    Cmax=max(mp.fsum(mag[n]**2*(n-j)**2*abs(hats[n-j])**2
                     for n in range(L) if n!=j) for j in lower)
    moment=Cmax/(2*(b*b-a*a))
    tol=mp.mpf('1e-65')
    assert c2+tol>=diagonal and diagonal+tol>=moment
    return c2,diagonal,moment

def main():
    exact=rational_checks()
    with (ROOT.parent/'moment_sensitivity/data/directional_sweep.csv').open(encoding='utf8') as f:
        supplied=list(csv.DictReader(f))
    cases=sorted({(int(r['L']),r['family'],r['direction']) for r in supplied
                  if r['flat']=='True' and r['rotation']=='False'})
    assert len(cases)==22
    results=[]; rotations=0; replays=[]
    with mp.workdps(80):
        for L,family,kind in cases:
            w,mag=weight_vector(L,family); h=direction(L,kind)
            c2,diagonal,moment=coefficient(L,mag,h)
            assert moment>0
            for sign in [-1,1]:
                for t_text in ['1e-5','1e-7']:
                    t=sign*mp.mpf(t_text)
                    observed=excess(L,w,mag,h,t)/(t*t)
                    error=abs(observed-c2)/max(1,abs(c2))
                    assert error<mp.mpf('0.005')
                    if t_text=='1e-7': assert error<mp.mpf('5e-5')
                    assert observed>=moment/2
                    results.append(dict(L=L,family=family,direction=kind,sign=sign,t=t_text,
                        c2=mp.nstr(c2,55),diagonal_lower=mp.nstr(diagonal,55),
                        moment_lower=mp.nstr(moment,55),observed=mp.nstr(observed,55),
                        scaled_error=mp.nstr(error,20)))
        for L,family in sorted({(L,family) for L,family,_ in cases}):
            w,mag=weight_vector(L,family)
            c2,diagonal,moment=coefficient(L,mag,[mp.mpf(1)]*L)
            assert max(abs(c2),abs(diagonal),abs(moment))<mp.mpf('1e-65')
            rotations+=1
        for target in [(17,'triangular','cos2'),(9,'unique_extrema','cos1')]:
            old=next(r for r in results if (r['L'],r['family'],r['direction'])==target
                     and r['sign']==1 and r['t']=='1e-7')
            with mp.workdps(120):
                L,family,kind=target; w,mag=weight_vector(L,family); h=direction(L,kind)
                c2,_,_=coefficient(L,mag,h)
                measured=excess(L,w,mag,h,mp.mpf('1e-7'))/mp.mpf('1e-14')
                change=max(abs(mp.mpf(old['c2'])-c2),abs(mp.mpf(old['observed'])-measured))
                assert change<mp.mpf('1e-48')
                replays.append(dict(L=L,family=family,direction=kind,absolute_change=mp.nstr(change,15)))
    DATA=ROOT/'data'; DATA.mkdir(exist_ok=True)
    with (DATA/'general_second_order.csv').open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0])); writer.writeheader(); writer.writerows(results)
    report=dict(passed=True,precision_digits=80,coupled_rational_cases=exact,
                transition='1-a^2/b^2',incorrect_transition_counterexample_checked=True,
                nonconstant_flat_directions=len(cases),full_matrix_checks=len(results),
                constant_direction_checks=rotations,precision_rechecks=replays,
                maximum_scaled_error_at_1e7=max(float(r['scaled_error']) for r in results if r['t']=='1e-7'),
                scope='Fixed-dimension finite checks; no uniform remainder assertion.',
                environment=dict(python=sys.version,mpmath=mp.__version__,platform=platform.platform()))
    (DATA/'general_second_order_summary.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
