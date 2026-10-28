"""Reproduce in a new output directory, preserving the supplied manuscript records."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parent

def run(root, relative, *arguments):
    print('\nRunning '+relative,flush=True)
    script=root/relative
    subprocess.run([sys.executable,str(script),*arguments],cwd=script.parent,check=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true')
    parser.add_argument('--theory-only',action='store_true',
                        help='Rerun mathematical checks and rebuild figures, retaining recorded timings.')
    parser.add_argument('--figures-only',action='store_true',
                        help='Rebuild all three figures from recorded CSVs without remeasurement.')
    parser.add_argument('--output-dir',type=Path,
                        help='New directory for the copied scripts, outputs and run manifest.')
    args=parser.parse_args()
    if sum([args.full,args.figures_only,args.theory_only])>1:
        parser.error('--full, --figures-only and --theory-only are mutually exclusive')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=(args.output_dir or ROOT.parent/'runs'/stamp).resolve()
    if output==ROOT or ROOT in output.parents:
        parser.error('--output-dir must be outside the source repro directory')
    if output.exists():
        parser.error('--output-dir must be new; existing outputs are never replaced')
    def snapshot():
        return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    before=snapshot()
    output.mkdir(parents=True)
    copied=output/'repro'
    shutil.copytree(ROOT,copied,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    report={'mode':'figures-only' if args.figures_only else ('theory-only' if args.theory_only else ('full' if args.full else 'default')),
            'source':str(ROOT),'output':str(output),'source_sha256':before,'completed':False}
    report_path=output/'run_manifest.json'
    report_path.write_text(json.dumps(report,indent=2),encoding='utf8')
    if not args.figures_only:
        if args.full: run(copied,'mathematics/experiments.py')
        for script in ['theory_extensions/verify_extensions.py',
                       'moment_sensitivity/moment_sensitivity.py',
                       'moment_sensitivity/directional_checks.py',
                       'ela_revision/verify_second_order.py',
                       'ela_revision/verify_general_second_order.py',
                       'ela_revision/verify_alternating_example.py']:
            run(copied,script)
        # Timing is sequential, after all mathematical computations have completed.
        run(copied,'ela_revision/benchmark.py',*(['--plot-only'] if args.theory_only else []))
    else:
        run(copied,'ela_revision/benchmark.py','--plot-only')
    run(copied,'ela_revision/plot_from_data.py')
    shutil.copy2(copied/'ela_revision/figures/batch_timings.pdf',output/'figures/batch_timings.pdf')
    report['source_preserved']=before==snapshot()
    if not report['source_preserved']:
        raise RuntimeError('Source reproduction files changed during the run')
    report['completed']=True
    report_path.write_text(json.dumps(report,indent=2),encoding='utf8')
    print('\nFinished. Outputs: '+str(output))
    print('Supplied data, manuscript figures and inline timing values were preserved.')
    if not args.figures_only and not args.theory_only:
        print('Fresh timings describe this rerun only. See README.md before updating the paper.')

if __name__=='__main__': main()
