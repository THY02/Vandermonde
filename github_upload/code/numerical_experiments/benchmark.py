"""Fair fixed-matrix batch benchmark. Run: python benchmark.py.

All four methods use identical complex128 inputs. Preparation is measured separately.
Every warm online repetition is retained in CSV; no recorded time is invented.
"""
from pathlib import Path
import argparse, csv, json, platform, sys, time, statistics
import numpy as np
import scipy
from scipy.linalg import lu_factor, lu_solve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits, threadpool_info

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
FIG = ROOT / 'figures'
SEED = 20261024
REPEATS = 7

def write_csv(name, rows):
    with (DATA/name).open('w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def measure(fn):
    fn()
    start = time.perf_counter(); fn(); pilot = time.perf_counter()-start
    loops = max(1, min(2000, int(.02/max(pilot, 1e-9))))
    values = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        for _ in range(loops): fn()
        values.append((time.perf_counter()-start)/loops)
    return values, loops

def main():
    DATA.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    summary, raw, prep, checks = [], [], [], []
    with threadpool_limits(limits=1):
        for d in [8,32,128,512]:
            L = 2*d+1
            def make_weights():
                n=np.arange(L); return np.minimum(n+1,L-n).astype(float)
            N=make_weights()
            phase=2*np.pi*np.arange(L)/L
            def make_matrix():
                beta=np.exp(1j*phase)
                return np.asfortranarray(N[:,None]*beta[None,:]**np.arange(L)[:,None])
            def make_operator():
                return np.exp(-1j*np.outer(np.arange(L),phase))/L/N[None,:]
            M=make_matrix(); R=make_operator()
            factors=lu_factor(M,check_finite=False,overwrite_a=False)
            setup={}
            for name,fn in [('weights',make_weights),('matrix',make_matrix),
                            ('dense_operator',make_operator),
                            ('lu_factor',lambda:lu_factor(M,check_finite=False,overwrite_a=False))]:
                values,loops=measure(fn); setup[name]=statistics.median(values)
                prep.append(dict(d=d,L=L,stage=name,median_seconds=setup[name],loops=loops,
                                 **{f'repeat_{i+1}':v for i,v in enumerate(values)}))
            for B in [1,128]:
                A=np.asfortranarray(rng.normal(size=(L,B))+1j*rng.normal(size=(L,B)))
                methods={
                    'FFT':lambda:np.fft.fft(A/N[:,None],axis=0)/L,
                    'Cached dense':lambda:R@A,
                    'Cached LU':lambda:lu_solve(factors,A,check_finite=False,overwrite_b=False),
                    'Fresh LU':lambda:lu_solve(lu_factor(M,check_finite=False,overwrite_a=False),
                                              A,check_finite=False,overwrite_b=False)}
                reference=methods['FFT']()
                for name,fn in methods.items():
                    X=fn()
                    relative=float(np.linalg.norm(X-reference)/np.linalg.norm(reference))
                    residual=float(np.linalg.norm(M@X-A)/np.linalg.norm(A))
                    assert relative<1e-9 and residual<1e-8, (d,B,name,relative,residual)
                    checks.append(dict(d=d,L=L,batch=B,method=name,
                                       relative_to_fft=relative,dense_residual=residual))
                # Warm all methods, then interleave their repetitions in seeded random order.
                loops={}
                for name,fn in methods.items():
                    fn(); start=time.perf_counter(); fn(); elapsed=time.perf_counter()-start
                    loops[name]=max(1,min(2000,int(.02/max(elapsed,1e-9))))
                observed={name:[] for name in methods}
                for rep in range(REPEATS):
                    for name in rng.permutation(list(methods)):
                        start=time.perf_counter()
                        for _ in range(loops[name]): methods[name]()
                        elapsed=(time.perf_counter()-start)/loops[name]
                        observed[name].append(elapsed)
                        raw.append(dict(d=d,L=L,batch=B,method=name,repeat=rep+1,
                                        loops=loops[name],seconds_per_batch=elapsed))
                for name,values in observed.items():
                    summary.append(dict(d=d,L=L,batch=B,method=name,
                                        median_seconds=statistics.median(values),
                                        min_seconds=min(values),max_seconds=max(values)))
            print(f'Completed d={d}',flush=True)
        env=dict(python=sys.version,platform=platform.platform(),processor=platform.processor(),
                 numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
                 seed=SEED,repetitions=REPEATS,threadpools=threadpool_info(),
                 clock='time.perf_counter',dtype='complex128',rhs_layout='Fortran contiguous',
                 preparation='weights; matrix; dense operator; LU factorization measured separately',
                 online='input generation and fixed-matrix setup excluded; Fresh LU includes factorization',
                 finite_checks=False,overwrite_inputs=False)
    for name,rows in [('batch_summary.csv',summary),('batch_raw.csv',raw),
                      ('preparation.csv',prep),('solver_checks.csv',checks)]: write_csv(name,rows)
    (DATA/'environment.json').write_text(json.dumps(env,indent=2),encoding='utf8')
    render_recorded_outputs(summary, prep)
    print(json.dumps({r['method']:r['median_seconds']*1000 for r in summary if r['d']==512 and r['batch']==128},indent=2))


def render_recorded_outputs(summary, prep):
    """Render Figure 3 and Table 2 rows without making fresh measurements."""
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(9,3.3))
    styles=['o-','s--','^-.','D:']
    for ax,B in zip(axes,[1,128]):
        for (name,style) in zip(['FFT','Cached dense','Cached LU','Fresh LU'],styles):
            rows=[r for r in summary if r['batch']==B and r['method']==name]
            ax.loglog([r['L'] for r in rows],[r['median_seconds']*1000 for r in rows],
                      style,label=name,linewidth=1.5,markersize=4)
        ax.set(xlabel='Number of nodes L',ylabel='Milliseconds per batch',title=f'Batch B = {B}')
        ax.grid(True,alpha=.2); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG/'batch_timings.pdf',bbox_inches='tight')
    fig.savefig(FIG/'batch_timings.png',dpi=180,bbox_inches='tight'); plt.close(fig)
    lines=[]
    for d in [8,32,128,512]:
        stage={r['stage']:r['median_seconds']*1000 for r in prep if r['d']==d}
        lines.append(f"{d} & "+' & '.join(f'{stage[k]:.3g}' for k in ['weights','matrix','dense_operator','lu_factor'])+r' \\')
    (DATA/'preparation_rows.tex').write_text('\n'.join(lines),encoding='utf8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot-only',action='store_true',
                        help='Rebuild Figure 3 and Table 2 rows from existing CSVs; do not time.')
    args=parser.parse_args()
    if args.plot_only:
        def read_numeric(name):
            with (DATA/name).open(encoding='utf8',newline='') as f:
                return [{k:(v if k in {'method','stage'} else float(v)) for k,v in row.items()}
                        for row in csv.DictReader(f)]
        render_recorded_outputs(read_numeric('batch_summary.csv'),read_numeric('preparation.csv'))
        print('Rebuilt Figure 3 and Table 2 rows from recorded CSVs; no timings measured.')
    else:
        main()
