"""Check recorded tables, raw timing medians, solver checks and exact coefficients.

This reads the supplied evidence without replacing it. It uses only the standard
library. Fresh mathematical computations and benchmark reruns use reproduce.py.
"""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext
import argparse, csv, json, math, statistics

ROOT=Path(__file__).resolve().parent


def rows(relative):
    with (ROOT/relative).open(encoding='utf8',newline='') as f:
        return list(csv.DictReader(f))


def check():
    source=ROOT.parent/'paper/article.tex'
    if not source.exists(): source=ROOT.parent/'article.tex'
    manuscript=source.read_text(encoding='utf8')
    accuracy=rows('mathematics/data/fft_accuracy.csv')
    table1=[]
    for d in [2,8,32,128]:
        group=[r for r in accuracy if int(r['d'])==d]
        assert len(group)==9
        values=[max(float(r[k]) for r in group) for k in
                ['weight_forward_error','dense_coefficient_residual','fft_roundtrip_residual']]
        table1.append(f'{d} & {2*d+1} & '+' & '.join(f'{v:.2e}' for v in values)+r' \\')
    assert all(line in manuscript for line in table1)

    summary=rows('numerical_experiments/data/batch_summary.csv')
    raw=rows('numerical_experiments/data/batch_raw.csv')
    prep=rows('numerical_experiments/data/preparation.csv')
    solves=rows('numerical_experiments/data/solver_checks.csv')
    assert len(summary)==32 and len(raw)==224 and len(prep)==16 and len(solves)==32
    methods=['FFT','Cached dense','Cached LU','Fresh LU']
    assert {r['method'] for r in summary}==set(methods)
    for r in summary:
        group=[x for x in raw if all(x[k]==r[k] for k in ['d','batch','method'])]
        assert {int(x['repeat']) for x in group}==set(range(1,8))
        values=[float(x['seconds_per_batch']) for x in group]
        assert float(r['median_seconds'])==statistics.median(values)
        assert float(r['min_seconds'])==min(values) and float(r['max_seconds'])==max(values)
    for r in prep:
        assert float(r['median_seconds'])==statistics.median(float(r[f'repeat_{i}']) for i in range(1,8))
    table2=[]
    for d in [8,32,128,512]:
        stage={r['stage']:float(r['median_seconds'])*1000 for r in prep if int(r['d'])==d}
        table2.append(f'{d} & '+' & '.join(f'{stage[k]:.3g}' for k in
                      ['weights','matrix','dense_operator','lu_factor'])+r' \\')
    assert all(line in manuscript for line in table2)
    assert '\n'.join(table2)==(ROOT/'numerical_experiments/data/preparation_rows.tex').read_text().strip()
    times={r['method']:float(r['median_seconds'])*1000 for r in summary
           if r['d']=='512' and r['batch']=='128'}
    time_text=', '.join(f'{times[m]:.3g}' for m in methods[:3])+f" and {times[methods[3]]:.3g} ms"
    assert time_text in manuscript
    chinese_source=source.parent/'chinese/article.tex'
    if chinese_source.exists():
        chinese=chinese_source.read_text(encoding='utf8')
        assert all(line in chinese for line in table1+table2)
        chinese_times='、'.join(f'${times[m]:.3g}$' for m in methods[:3])+f' 和 ${times[methods[3]]:.3g}$ ms'
        assert chinese_times in chinese
    assert max(float(r['relative_to_fft']) for r in solves)<1e-9
    assert max(float(r['dense_residual']) for r in solves)<1e-8

    exact=rows('numerical_experiments/data/second_order_exact.csv')
    assert [int(r['d']) for r in exact]==list(range(2,65))
    for r in exact:
        d=int(r['d'])
        expected=Q(21,8) if d==2 else Q(2*d*d,3)-Q(5*d,12)+Q(11,48)+Q(1,4*d)
        assert Q(r['coefficient'])==expected
        assert Q(r['bottom_diagonal'])==-Q(9,8)-Q((2*d-1)**2,3)
        assert Q(r['bottom_offdiagonal'])==(-Q(7,8) if d==2 else Q(0))
    high_precision=rows('numerical_experiments/data/second_order_high_precision.csv')
    assert len(high_precision)==42
    assert all(math.isfinite(float(r['observed'])) for r in high_precision)
    moments=rows('moment_sensitivity/data/moment_sweep.csv')
    assert len(moments)==48
    with localcontext() as ctx:
        ctx.prec=80
        for r in moments:
            a,b,eta=[D(r[k]) for k in ['a','b','eta']]
            threshold=1-a*a/(b*b)
            if eta<=threshold:
                bound=(b*b-a*a)*eta/(1+eta); branch='endpoint'
            else:
                bound=(b*b*(1+eta)-a*a)**2/(4*b*b*(1+eta)); branch='vertex'
            assert r['branch']==branch
            assert abs(bound-D(r['moment_rhs']))<D('1e-45')*max(1,bound)
            assert D(r['moment_lhs'])<=bound+D('1e-45')
            assert bound<=D(r['previous_bound'])+D('1e-45')
            assert abs(D(r['lhs_over_rhs'])-D(r['moment_lhs'])/bound)<D('1e-45')
    coefficients=rows('numerical_experiments/data/general_second_order.csv')
    assert len(coefficients)==88
    for r in coefficients:
        c,diag,moment=[float(r[k]) for k in ['c2','diagonal_lower','moment_lower']]
        assert c+1e-12*max(1,c)>=diag>=moment-1e-12*max(1,moment)>0
        assert float(r['observed'])>=moment/2
    general=json.loads((ROOT/'numerical_experiments/data/general_second_order_summary.json').read_text())
    assert general['passed'] and general['nonconstant_flat_directions']==22
    assert '$0.334$' in manuscript and r'\label{eq:general-second}' in manuscript
    example=json.loads((ROOT/'numerical_experiments/data/alternating_example_summary.json').read_text())
    assert example['passed'] and example['full_matrix_checks']==12
    assert example['bottom_block']==[['-10/3','-5/3'],['-5/3','-2/3']]
    assert example['top_block']==[['2/3','-10/3'],['-10/3','10/3']]
    return dict(passed=True,accuracy_cases=len(accuracy),table1_matches=True,table2_matches=True,
                timing_sentence_matches=True,benchmark_configurations=len(summary),
                raw_repetitions=len(raw),raw_medians_match=True,exact_degrees=len(exact),
                expanded_coefficient_identity=True,high_precision_records=len(high_precision),
                coupled_moment_records=len(moments),general_second_order_records=len(coefficients),
                general_coefficient_lower_bounds_checked=True,
                alternating_example_checks=example['full_matrix_checks'],
                figure3_data='numerical_experiments/data/batch_summary.csv',
                figure_rebuild_command='python reproduce.py --figures-only',
                max_solver_difference_from_fft=max(float(r['relative_to_fft']) for r in solves),
                max_dense_solver_residual=max(float(r['dense_residual']) for r in solves))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    result=check()
    text=json.dumps(result,indent=2)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(text,encoding='utf8')
    print(text)
