"""Rebuild Figures 1 and 2 from recorded CSVs, with distinct grayscale styles."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT.parent/'figures'

def read(name):
    with (ROOT/'moment_sensitivity/data'/name).open(encoding='utf8',newline='') as f:
        return list(csv.DictReader(f))

def main():
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    rows=read('directional_sweep.csv')
    fig,ax=plt.subplots(figsize=(7.2,3.8),layout='constrained')
    for kind,label,color,style in [
        ('cos1','First harmonic: linear','#2465a8','o-'),
        ('cos2','Second harmonic: quadratic','#d26923','s--'),
        ('projected','Projected direction: quadratic','#228266','^-.')]:
        group=[r for r in rows if r['L']=='17' and r['family']=='triangular'
               and r['direction']==kind and r['sign']=='1']
        ax.loglog([float(r['t']) for r in group],
                  [float(r['condition_excess']) for r in group],style,
                  color=color,label=label,ms=5,lw=1.4)
    ts=[1e-8,1e-2]
    ax.loglog(ts,[4*t for t in ts],':',color='black',alpha=.7,
              label=r'Predicted $4|t|$ ($d=8$)')
    ax.set(xlabel=r'Perturbation parameter $|t|$',ylabel=r'$\kappa_2(DV(t))/(d+1)-1$')
    ax.legend(frameon=False,fontsize=9); ax.grid(True,which='major',alpha=.2)
    fig.savefig(FIG/'directional_orders.pdf'); plt.close(fig)

    rows=read('moment_sweep.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
    for ax,family,title in zip(axes,['interior_minimum','triangular'],
            ['(a) Unique interior minimum','(b) Triangular weights']):
        for L,color,style in [(5,'#2465a8','o-'),(9,'#d26923','s--'),(17,'#228266','^-.')]:
            group=[r for r in rows if r['family']==family and r['L']==str(L)]
            ax.loglog([float(r['eta']) for r in group],[float(r['lhs_over_rhs']) for r in group],
                      style,ms=4.5,lw=1.3,color=color,label=f'L = {L}')
        ax.axhline(1,color='#555555',ls=':',lw=1.2,label='Bound (ratio = 1)')
        ax.set(title=title,xlabel=r'Condition-number excess $\eta$',ylabel='Moment sum / upper bound')
        ax.grid(True,which='major',alpha=.2); ax.legend(fontsize=8.5,frameon=False,loc='best')
    fig.savefig(FIG/'moment_bound_ratio.pdf'); plt.close(fig)
    print('Rebuilt Figures 1 and 2 from raw data; line and marker styles distinguish series in grayscale.')

if __name__=='__main__': main()
