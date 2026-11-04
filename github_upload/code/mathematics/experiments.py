"""Reproduce all experiments, figures and TeX tables. Run from this directory.

Requires numpy, matplotlib, mpmath, threadpoolctl. No downloaded datasets.
All timings are measured locally; reruns will have different timing values.
"""
from pathlib import Path
import os, sys, csv, json, platform, time, statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
from threadpoolctl import threadpool_limits, threadpool_info

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; FIG=ROOT/'figures'
DATA.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
SEED=20261023
rng=np.random.default_rng(SEED)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                     'axes.grid':True,'grid.alpha':.22,'savefig.bbox':'tight'})

def weights(d):
    n=np.arange(2*d+1)
    return np.minimum(n+1,2*d-n+1).astype(float)

def matrix(d,beta):
    return weights(d)[:,None]*beta[None,:]**np.arange(2*d+1)[:,None]

def decompose(A,d):
    """Standard exact roots-of-unity convention; axis zero holds coefficients."""
    N=weights(d)
    return np.fft.fft(A/N.reshape((-1,)+(1,)*(A.ndim-1)),axis=0)/(2*d+1)

def save_csv(name,rows):
    with (DATA/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def savefig(name):
    plt.savefig(FIG/f'{name}.pdf');plt.savefig(FIG/f'{name}.png',dpi=190);plt.close()

def plot_nodes(rows):
    limit=1/np.finfo(float).eps
    names={'roots':'Roots of unity','rotated':'Rotated roots','jittered':'Jittered roots',
           'random_circle':'Random circle','quarter_arc':'Quarter arc','real_grid':'Real grid'}
    for fam,label in names.items():
        ds=sorted(set(r['d'] for r in rows))
        med=np.array([statistics.median(r['condition_estimate'] for r in rows if r['d']==d and r['family']==fam) for d in ds])
        line,=plt.loglog([2*d+1 for d in ds],np.minimum(med,limit),'o-',label=label)
        bad=med>=limit
        if np.any(bad):plt.scatter(np.array([2*d+1 for d in ds])[bad],np.full(sum(bad),limit),marker='^',s=70,color=line.get_color(),zorder=4)
    plt.xlabel('Number of nodes L = 2d + 1');plt.ylabel('Estimated spectral condition number')
    plt.axhline(limit,color='black',ls=':',label='Resolution limit (censored above)')
    plt.legend(fontsize=8,ncol=2);savefig('node_conditioning')

def condition_table(rows):
    lines=[]
    for d in [2,4,8,16,32,64]:
        vals=[statistics.median(r['condition_estimate'] for r in rows if r['d']==d and r['family']==fam) for fam in ['roots','jittered','random_circle','quarter_arc','real_grid']]
        def fmt(v):return r'\text{失辨}' if v>=1/np.finfo(float).eps else f'{v:.2e}'
        lines.append(f"{d} & "+' & '.join(fmt(v) for v in vals)+' \\\\')
    (DATA/'condition_rows.tex').write_text('\n'.join(lines),encoding='utf-8')

def plot_noise(rows):
    for d in [2,8,32,128]:
        ss=[r for r in rows if r['d']==d and r['kind']=='sharp']
        x=[r['epsilon'] for r in ss]; y=[r['relative_weight_change'] for r in ss]
        line,=plt.loglog(x,y,'o',label=f'd={d}, measured')
        plt.loglog(x,[(d+1)*v for v in x],'--',color=line.get_color(),label=f'd={d}, bound')
    plt.xlabel('Relative input perturbation');plt.ylabel('Relative change in decomposition weights')
    plt.legend(fontsize=8,ncol=2);savefig('sharp_perturbation')

def normmp(x): return mp.sqrt(mp.fsum(abs(z)**2 for z in x))
def tomp(z): return mp.mpc(float(np.real(z)),float(np.imag(z)))

def mp_ideal_reference(A,d,dps=80):
    """Independent O(L^2) DFT at ideal roots, high precision, no FFT roundtrip."""
    with mp.workdps(dps):
        L=2*d+1; av=[tomp(x)/int(n) for x,n in zip(A,weights(d))]
        out=[]
        for k in range(L):
            q=mp.exp(-2*mp.pi*mp.j*k/L); p=mp.mpc(1);s=mp.mpc(0)
            for x in av: s+=x*p;p*=q
            out.append(s/L)
        return out

def mp_error(x,ref):
    with mp.workdps(80): return float(normmp([tomp(z)-r for z,r in zip(x,ref)])/normmp(ref))

def node_experiment():
    rows=[]; families=['roots','rotated','jittered','random_circle','quarter_arc','real_grid']
    for d in [2,4,8,16,32,64]:
        L=2*d+1;base=2*np.pi*np.arange(L)/L
        for trial in range(10):
            nodes={
                'roots':np.exp(1j*base),
                'rotated':np.exp(1j*(base+.37)),
                'jittered':np.exp(1j*(base+rng.uniform(-.25,.25,L)*2*np.pi/L)),
                'random_circle':np.exp(1j*rng.uniform(0,2*np.pi,L)),
                'quarter_arc':np.exp(1j*base/4),
                'real_grid':np.linspace(-1,1,L).astype(complex)}
            for fam in families:
                M=matrix(d,nodes[fam]);s=np.linalg.svd(M,compute_uv=False)
                cond=float(s[0]/s[-1]) if s[-1]>0 else float('inf')
                # Label values near/beyond 1/u as unresolved numerical estimates.
                rows.append(dict(d=d,L=L,trial=trial,family=fam,condition_estimate=cond,
                                 beyond_double_resolution=int(cond*np.finfo(float).eps>=1)))
    save_csv('node_conditions.csv',rows)
    # Reference solutions for rounded M and A; same A across node families in each trial.
    reference=[]
    for d in [4,8,16]:
        L=2*d+1;base=2*np.pi*np.arange(L)/L
        for trial in range(3):
            A=rng.normal(size=L)+1j*rng.normal(size=L)
            for fam,beta in [('roots',np.exp(1j*base)),
                             ('jittered',np.exp(1j*(base+rng.uniform(-.25,.25,L)*2*np.pi/L))),
                             ('real_grid',np.linspace(-1,1,L).astype(complex))]:
                M=matrix(d,beta); x=np.linalg.solve(M,A)
                with mp.workdps(100):
                    ref=mp.lu_solve(mp.matrix([[tomp(z) for z in row] for row in M]),mp.matrix([tomp(z) for z in A]))
                    ferr=float(normmp([tomp(z)-r for z,r in zip(x,ref)])/normmp(ref))
                res=float(np.linalg.norm(M@x-A)/np.linalg.norm(A))
                reference.append(dict(d=d,trial=trial,family=fam,forward_error=ferr,coefficient_residual=res))
    save_csv('rounded_system_reference.csv',reference)
    plot_nodes(rows)
    return rows,reference

def accuracy_experiment():
    rows=[]
    for d in [2,8,32,128]:
        L=2*d+1;N=weights(d)
        for fam in ['random_complex','sparse','wide_scale']:
            for trial in range(3):
                if fam=='random_complex':A=rng.normal(size=L)+1j*rng.normal(size=L)
                elif fam=='sparse':
                    A=np.zeros(L,dtype=complex);idx=rng.choice(L,min(4,L),replace=False)
                    A[idx]=rng.normal(size=len(idx))+1j*rng.normal(size=len(idx))
                else:A=10.**rng.uniform(-8,8,L)*np.exp(1j*rng.uniform(0,2*np.pi,L))
                x=decompose(A,d);ref=mp_ideal_reference(A,d)
                forward=mp_error(x,ref)
                beta=np.exp(2j*np.pi*np.arange(L)/L)
                dense=float(np.linalg.norm(matrix(d,beta)@x-A)/np.linalg.norm(A))
                roundtrip=float(np.linalg.norm(N*(L*np.fft.ifft(x))-A)/np.linalg.norm(A))
                rows.append(dict(d=d,L=L,family=fam,trial=trial,weight_forward_error=forward,
                                 dense_coefficient_residual=dense,fft_roundtrip_residual=roundtrip))
    save_csv('fft_accuracy.csv',rows)
    for fam in ['random_complex','sparse','wide_scale']:
        ds=[2,8,32,128]
        yy=[max(r['weight_forward_error'] for r in rows if r['d']==d and r['family']==fam) for d in ds]
        plt.loglog([2*d+1 for d in ds],yy,'o-',label=fam)
    plt.xlabel('L = 2d + 1');plt.ylabel('Maximum relative weight error (80-digit reference)')
    plt.legend();savefig('fft_accuracy')
    return rows

def perturbation_experiment():
    rows=[]
    for d in [2,8,32,128]:
        L=2*d+1
        for eps in [1e-12,1e-10,1e-8,1e-6,1e-4]:
            A=np.zeros(L);A[d]=1;delta=np.zeros(L);delta[0]=eps
            x=decompose(A,d);xp=decompose(A+delta,d)
            rel=float(np.linalg.norm(xp-x)/np.linalg.norm(x));amp=rel/eps
            rows.append(dict(d=d,epsilon=eps,kind='sharp',trial=0,relative_weight_change=rel,
                             amplification=amp,theoretical_upper=d+1))
            for trial in range(20):
                A=rng.normal(size=L)+1j*rng.normal(size=L)
                delta=rng.normal(size=L)+1j*rng.normal(size=L)
                delta*=eps*np.linalg.norm(A)/np.linalg.norm(delta)
                x=decompose(A,d);xp=decompose(A+delta,d)
                rel=float(np.linalg.norm(xp-x)/np.linalg.norm(x))
                rows.append(dict(d=d,epsilon=eps,kind='random',trial=trial,relative_weight_change=rel,
                                 amplification=rel/eps,theoretical_upper=d+1))
    save_csv('input_perturbations.csv',rows)
    plot_noise(rows)
    return rows

def measured_time(fn):
    fn();t=time.perf_counter();fn();single=time.perf_counter()-t
    loops=max(1,min(1000,int(.025/max(single,1e-9))))
    values=[]
    for _ in range(7):
        t=time.perf_counter()
        for j in range(loops): ans=fn()
        values.append((time.perf_counter()-t)/loops)
    return statistics.median(values),min(values),max(values),loops

def batch_experiment():
    rows=[]
    for d in [8,32,128,512]:
        L=2*d+1;N=weights(d);phase=2*np.pi*np.arange(L)/L
        t=time.perf_counter();beta=np.exp(1j*phase);M=matrix(d,beta);matrix_setup=time.perf_counter()-t
        t=time.perf_counter();R=np.exp(-1j*np.outer(np.arange(L),phase))/L/N[None,:];inverse_setup=time.perf_counter()-t
        for B in [1,128]:
            A=rng.normal(size=(L,B))+1j*rng.normal(size=(L,B))
            functions={'fft':lambda:decompose(A,d),'cached_dense':lambda:R@A,'direct_solve':lambda:np.linalg.solve(M,A)}
            for name,fn in functions.items():
                med,mn,mx,loops=measured_time(fn)
                rows.append(dict(d=d,L=L,batch=B,method=name,median_seconds=med,min_seconds=mn,
                                 max_seconds=mx,loops=loops,matrix_setup_seconds=matrix_setup,
                                 cached_operator_setup_seconds=inverse_setup,cached_operator_bytes=R.nbytes))
    save_csv('batch_timings.csv',rows)
    fig,axes=plt.subplots(1,2,figsize=(10,3.5))
    for ax,B in zip(axes,[1,128]):
        for name in ['fft','cached_dense','direct_solve']:
            ss=[r for r in rows if r['batch']==B and r['method']==name]
            ax.loglog([r['L'] for r in ss],[r['median_seconds'] for r in ss],'o-',label=name)
        ax.set_title(f'Batch B = {B}');ax.set_xlabel('L = 2d + 1');ax.set_ylabel('Seconds per batch');ax.legend(fontsize=8)
    fig.tight_layout();savefig('batch_timings')
    return rows

def tables_and_summary(node,ref,acc,noise,timing):
    def sci(x):return f'{x:.2e}'
    # TeX fragments deliberately contain data only, no hand-entered numerical claims.
    lines=[]
    for d in [2,8,32,128]:
        ss=[r for r in acc if r['d']==d]
        lines.append(f"{d} & {2*d+1} & {sci(max(r['weight_forward_error'] for r in ss))} & {sci(max(r['dense_coefficient_residual'] for r in ss))} & {sci(max(r['fft_roundtrip_residual'] for r in ss))} \\\\")
    (DATA/'accuracy_rows.tex').write_text('\n'.join(lines),encoding='utf-8')
    lines=[]
    for d in [8,32,128,512]:
        s={r['method']:r['median_seconds'] for r in timing if r['d']==d and r['batch']==128}
        lines.append(f"{d} & {2*d+1} & {sci(s['fft']*1e3)} & {sci(s['cached_dense']*1e3)} & {sci(s['direct_solve']*1e3)} \\\\")
    (DATA/'timing_rows.tex').write_text('\n'.join(lines),encoding='utf-8')
    condition_table(node)
    summary={
      'seed':SEED,'node_trials':len(node),'rounded_system_references':len(ref),'accuracy_trials':len(acc),
      'perturbation_trials':len(noise),'timing_configurations':len(timing),
      'max_weight_forward_error':max(r['weight_forward_error'] for r in acc),
      'max_dense_residual':max(r['dense_coefficient_residual'] for r in acc),
      'max_fft_roundtrip_residual':max(r['fft_roundtrip_residual'] for r in acc),
      'max_root_condition_relative_discrepancy':max(abs(r['condition_estimate']/(r['d']+1)-1) for r in node if r['family'] in ['roots','rotated']),
      'sharp_amplification_max_relative_discrepancy':max(abs(r['amplification']/(r['d']+1)-1) for r in noise if r['kind']=='sharp'),
      'random_amplification_max_fraction_of_bound':max(r['amplification']/(r['d']+1) for r in noise if r['kind']=='random'),
      'timing_at_d512_B128':{r['method']:r['median_seconds'] for r in timing if r['d']==512 and r['batch']==128}}
    (DATA/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2),flush=True)
    assert summary['max_weight_forward_error']<1e-12
    assert summary['max_root_condition_relative_discrepancy']<1e-10
    assert summary['random_amplification_max_fraction_of_bound']<1.001

def main():
    with threadpool_limits(limits=1):
        meta={'seed':SEED,'python':sys.version,'platform':platform.platform(),'processor':platform.processor(),
              'numpy':np.__version__,'matplotlib':matplotlib.__version__,'mpmath':mp.__version__,
              'threadpools':threadpool_info(),'timing_clock':'time.perf_counter','repetitions':7,
              'reference_digits':{'ideal_roots':80,'rounded_linear_system':100}}
        (DATA/'environment.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
        print('1/4 node selection and independent linear-system references',flush=True);node,ref=node_experiment()
        print('2/4 FFT accuracy against ideal-root high-precision DFT',flush=True);acc=accuracy_experiment()
        print('3/4 input perturbations',flush=True);noise=perturbation_experiment()
        print('4/4 measured single-thread batch timings',flush=True);timing=batch_experiment()
        tables_and_summary(node,ref,acc,noise,timing)

if __name__=='__main__':
    if '--plots-only' in sys.argv:
        def load_csv(name,text_fields):
            with (DATA/name).open(encoding='utf-8') as f:
                return [{k:(v if k in text_fields else float(v)) for k,v in row.items()} for row in csv.DictReader(f)]
        nodes=load_csv('node_conditions.csv',{'family'})
        plot_nodes(nodes);condition_table(nodes)
        plot_noise(load_csv('input_perturbations.csv',{'kind'}))
    else:main()
