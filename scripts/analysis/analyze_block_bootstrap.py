#!/usr/bin/env python3
"""Autocorrelation and block-bootstrap feasibility diagnostics.

Uses the project's original C++ AutoCorr implementation without modifying it.
Does not run bootstrap or multihistogram. ASCII source only.
"""
import argparse
import csv
import math
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

OBS = ('energy', 'abs_m', 'abs_ec', 'abs_d', 'm_squared', 'ec_squared', 'd_squared')
DRIVER = r'''
#include <vector>
using namespace std;
#include "Autocorr.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 5) return 2;
    ifstream in(argv[1]);
    if (!in) { cerr << "Cannot open input\n"; return 2; }
    const double vol = pow(stod(argv[2]), stoi(argv[3]));
    const bool square = string(argv[4]) == "square";
    vector<vector<double>> series(7);
    string line;
    int lineno=0;
    while (getline(in,line)) {
        ++lineno;
        size_t p=line.find_first_not_of(" \t\r");
        if (p==string::npos || line[p]=='#') continue;
        istringstream row(line);
        double x[10];
        for (int k=0;k<10;k++) {
            if (!(row >> x[k]) || !isfinite(x[k])) {
                cerr << "Invalid data at line " << lineno << "\n";
                return 3;
            }
        }
        const double ec=(x[2]-x[3])/2.0;
        const double d=x[0]+(square ? 2.0 : 1.0)*ec;
        const double values[7]={x[1]/vol, fabs(x[0])/vol, fabs(ec)/vol,
             fabs(d)/vol, pow(x[0]/vol,2), pow(ec/vol,2), pow(d/vol,2)};
        for (int k=0;k<7;k++) series[k].push_back(values[k]);
    }
    const char* names[7]={"energy","abs_m","abs_ec","abs_d",
                           "m_squared","ec_squared","d_squared"};
    cout << setprecision(17);
    for (int k=0;k<7;k++) {
        const auto &v=series[k];
        double tau=0,err=0;
        string status="ok";
        if (v.size()<4) status="too_short";
        else {
            bool constant=true;
            for (double a:v) if(a!=v[0]) { constant=false; break; }
            if (constant) status="constant";
            else {
                AutoCorr(v,tau,err);
                if (!isfinite(tau)||!isfinite(err)||tau<=0) status="failed";
            }
        }
        cout << names[k] << "," << v.size() << ",";
        if(status=="ok") cout << tau << "," << err;
        else cout << "nan,nan";
        cout << "," << status << "\n";
    }
}
'''

def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def percentile(values,p):
    return float(np.quantile(values,p)) if values else math.nan

def block_diagnostics(tau, N, factor, min_blocks):
    if not math.isfinite(tau) or N<=0:
        return math.nan, math.nan, math.nan, 'invalid_tau'
    # tau >= 0.5 is a conservative reference floor for positive correlations.
    # It is used for block sizing only; raw estimates remain in the output.
    length=max(1,math.ceil(factor*max(0.5,tau)))
    blocks=N//length
    unused=N%length
    return length, blocks, unused, ('feasible' if blocks>=min_blocks else 'too_few_blocks')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--autocorr-source', type=Path, required=True)
    p.add_argument('--autocorr-header-dir', type=Path)
    p.add_argument('--lattice-type', choices=('square','triangular'), required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--sizes', nargs='+', type=int)
    p.add_argument('--dimension', type=int, default=2)
    p.add_argument('--compiler', default='g++')
    p.add_argument('--block-factor', type=float, default=10.0,
                   help='Proposed nonoverlapping block length = ceil(factor * max(tau, 0.5)).')
    p.add_argument('--min-blocks', type=int, default=20)
    p.add_argument('--no-plots', action='store_true')
    a=p.parse_args()
    if a.dimension!=2: p.error('EC/D definitions implemented only for 2D.')
    if a.block_factor<=0 or not math.isfinite(a.block_factor): p.error('block-factor must be positive.')
    if a.min_blocks<2: p.error('min-blocks must be >= 2.')
    if not a.input_dir.is_dir(): p.error('Input directory does not exist.')
    if not a.autocorr_source.is_file(): p.error('Autocorr source does not exist.')
    files=[]
    for path in sorted(a.input_dir.rglob('conf-*.dat')):
        m=re.fullmatch(r'conf-([0-9.eE+\-]+)-(\d+)\.dat',path.name)
        if not m: continue
        L=None
        for parent in path.parents:
            sm=re.fullmatch(r'(?:size_)?(\d+)',parent.name)
            if sm: L=int(sm[1]); break
        if L is None:
            raise ValueError(f'Cannot infer L from {path}; expected size_L or L directory.')
        if a.sizes and L not in a.sizes: continue
        files.append((L,float(m[1]),int(m[2]),path.resolve()))
    if not files: p.error('No matching conf files found.')
    out=a.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    with tempfile.TemporaryDirectory(prefix='autocorr_') as tmp:
        driver=Path(tmp)/'driver.cpp'
        driver.write_text(DRIVER,encoding='ascii')
        exe=Path(tmp)/'autocorr_driver'
        header=(a.autocorr_header_dir or a.autocorr_source.parent).resolve()
        subprocess.run([a.compiler,'-O2','-std=c++11','-I',str(header),str(driver),
                        str(a.autocorr_source.resolve()),'-o',str(exe)],check=True)
        for j,(L,T,run,path) in enumerate(files,1):
            proc=subprocess.run([str(exe),str(path),str(L),str(a.dimension),a.lattice_type],
                                capture_output=True,text=True)
            if proc.returncode:
                raise RuntimeError(f'{path}: {proc.stderr}')
            for rec in csv.reader(proc.stdout.splitlines()):
                if len(rec)!=5: raise RuntimeError(f'Unexpected driver output for {path}: {rec}')
                obs,n,tau,err,status=rec
                n=int(n);tau=float(tau);err=float(err)
                rows.append(dict(lattice_type=a.lattice_type,L=L,T=T,run=run,observable=obs,
                                 N=n,tau_int=tau,tau_err=err,
                                 variance_inflation=2*tau if status=='ok' else math.nan,
                                 N_eff_mean=n/(2*tau) if status=='ok' else math.nan,
                                 status=status,file=str(path)))
            if j%50==0 or j==len(files): print(f'Processed {j}/{len(files)} files',flush=True)
    allfields=list(rows[0])
    write_csv(out/'autocorr_all.csv',rows,allfields)
    write_csv(out/'autocorr_invalid.csv',[r for r in rows if r['status']!='ok'],allfields)

    # Exactly one common block length per time series, selected from all seven observables.
    grouped=defaultdict(list)
    for r in rows: grouped[(r['L'],r['T'],r['run'],r['file'])].append(r)
    blockrows=[]
    for (L,T,run,path),group in sorted(grouped.items()):
        valid=[r for r in group if r['status']=='ok']
        invalid=[r['observable'] for r in group if r['status']!='ok']
        Ns={r['N'] for r in group}
        if len(Ns)!=1: raise RuntimeError(f'Inconsistent N for {path}')
        N=Ns.pop()
        worst=max(valid,key=lambda r:r['tau_int']) if valid else None
        tau=worst['tau_int'] if worst else math.nan
        ell,nb,unused,feas=block_diagnostics(tau,N,a.block_factor,a.min_blocks)
        # Incomplete coverage means the proposed block is only provisional.
        if invalid: feas='incomplete_autocorr' if valid else 'no_valid_autocorr'
        tau_limit=N/(a.block_factor*a.min_blocks)
        blockrows.append(dict(lattice_type=a.lattice_type,L=L,T=T,run=run,N=N,
                              valid_observables=len(valid),invalid_observables=';'.join(invalid),
                              tau_max=tau,tau_max_observable=worst['observable'] if worst else '',
                              tau_err_at_max=worst['tau_err'] if worst else math.nan,
                              block_factor=a.block_factor,block_length=ell,n_blocks=nb,
                              unused_measurements=unused,min_blocks=a.min_blocks,
                              max_tau_for_min_blocks=tau_limit,status=feas,file=path))
    blockfields=list(blockrows[0])
    write_csv(out/'block_sizes_by_series.csv',blockrows,blockfields)
    write_csv(out/'block_problem_cases.csv',
              [r for r in blockrows if r['status']!='feasible'],blockfields)

    summary=[]
    byL=defaultdict(list)
    for r in blockrows: byL[r['L']].append(r)
    for L,g in sorted(byL.items()):
        valid=[r for r in g if r['valid_observables']==len(OBS)]
        t=[r['tau_max'] for r in valid]
        lengths=[r['block_length'] for r in valid]
        blocks=[r['n_blocks'] for r in valid]
        statuscounts={name:sum(r['status']==name for r in g)
                      for name in ('feasible','too_few_blocks','incomplete_autocorr','no_valid_autocorr')}
        summary.append(dict(L=L,n_series=len(g),n_complete=len(valid),**statuscounts,
                            tau_max_median=percentile(t,.5),tau_max_q95=percentile(t,.95),
                            tau_max_overall=max(t) if t else math.nan,
                            block_length_median=percentile(lengths,.5),
                            block_length_q95=percentile(lengths,.95),
                            block_length_max=max(lengths) if lengths else math.nan,
                            n_blocks_min=min(blocks) if blocks else math.nan,
                            n_blocks_median=percentile(blocks,.5)))
    write_csv(out/'block_summary_by_L.csv',summary,list(summary[0]))
    autocorr_summary=[]
    byObs=defaultdict(list)
    for r in rows: byObs[(r['L'],r['observable'])].append(r)
    for (L,obs),g in sorted(byObs.items()):
        good=[r for r in g if r['status']=='ok']
        values=[r['tau_int'] for r in good]
        worst=max(good,key=lambda r:r['tau_int']) if good else None
        autocorr_summary.append(dict(L=L,observable=obs,files=len(g),valid=len(good),
            invalid=len(g)-len(good),T_min=min(r['T'] for r in g),
            T_max=max(r['T'] for r in g),tau_min=min(values) if values else math.nan,
            tau_median=percentile(values,.5),tau_q95=percentile(values,.95),
            tau_max=max(values) if values else math.nan,
            T_at_tau_max=worst['T'] if worst else math.nan,
            run_at_tau_max=worst['run'] if worst else '',
            tau_err_at_max=worst['tau_err'] if worst else math.nan,
            N_eff_min=min(r['N_eff_mean'] for r in good) if good else math.nan,
            count_tau_above_1=sum(v>1 for v in values)))
    write_csv(out/'autocorr_summary_by_L.csv',autocorr_summary,list(autocorr_summary[0]))

    total=len(blockrows)
    notes=[
        'Autocorrelation and block-bootstrap feasibility diagnostics',
        f'Lattice: {a.lattice_type}', f'Input files: {total}',
        f'Block rule: length = ceil({a.block_factor:g} * max(0.5, tau_max))',
        f'Minimum nonoverlapping blocks: {a.min_blocks}',
        'The maximum autocorrelation is taken over energy, |M|, |EC_sym|, |D| and three squared observables.',
        'Tau units: intervals between saved measurements (ndecorr already applied by simulation).',
        'AutoCorr convention: tau_int = 0.5 + sum rho(t).',
        'Factor 2 belongs in N_eff = N/(2*tau_int), NOT an additional factor in the 10*tau block rule.',
        'Proposed lengths assume nonoverlapping blocks; N % length observations would be unused.',
        'A status of incomplete_autocorr marks any failed observable; do not trust its proposed block length.',
        'The tau>=0.5 floor is only for block sizing; actual tau estimates are not changed.',
        'These are provisional block sizes; test bootstrap stability against block length.',
        'No bootstrap or multihistogram calculation was performed.',
        '',
        f'Feasible: {sum(r["status"]=="feasible" for r in blockrows)}/{total}',
        f'Too few blocks: {sum(r["status"]=="too_few_blocks" for r in blockrows)}/{total}',
        f'Incomplete autocorr: {sum(r["status"]=="incomplete_autocorr" for r in blockrows)}/{total}',
        f'No valid autocorr: {sum(r["status"]=="no_valid_autocorr" for r in blockrows)}/{total}',
        '', 'L series complete feasible too_few incomplete no_valid tau_med tau_q95 ell_med blocks_min'
    ]
    for r in summary:
        notes.append(f'{r["L"]} {r["n_series"]} {r["n_complete"]} {r["feasible"]} '
          f'{r["too_few_blocks"]} {r["incomplete_autocorr"]} {r["no_valid_autocorr"]} '
          f'{r["tau_max_median"]:.4g} {r["tau_max_q95"]:.4g} '
          f'{r["block_length_median"]:.4g} {r["n_blocks_min"]}')
    (out/'block_summary.txt').write_text('\n'.join(notes)+'\n',encoding='ascii')

    if not a.no_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print('Matplotlib unavailable; saved CSV and text only.')
        else:
            for obs in OBS:
                fig,ax=plt.subplots(figsize=(9,5))
                for L in sorted({r['L'] for r in rows}):
                    g=sorted((r for r in rows if r['L']==L and r['observable']==obs and r['status']=='ok'),
                             key=lambda r:(r['T'],r['run']))
                    if g:
                        ax.errorbar([r['T'] for r in g],[r['tau_int'] for r in g],
                                    yerr=[r['tau_err'] for r in g],fmt='o',ms=3,
                                    alpha=.6,capsize=2,label=f'L={L}')
                ax.axhline(.5,color='black',ls='--',lw=1,label='Independent reference')
                ax.set(xlabel='Temperature T',ylabel='Integrated autocorrelation time',title=obs)
                ax.legend(fontsize=8,ncol=2)
                ax.grid(alpha=.2);fig.tight_layout()
                fig.savefig(out/f'autocorr_vs_T_{obs}.png',dpi=180);plt.close(fig)
            for key,ylabel,title,filename in (
                ('tau_max','Maximum estimated tau_int','Maximum autocorrelation vs L','tau_max_vs_L.png'),
                ('block_length','Proposed block length','Block length vs L','block_length_vs_L.png'),
                ('n_blocks','Number of nonoverlapping blocks','Available blocks vs L','n_blocks_vs_L.png')):
                fig,ax=plt.subplots(figsize=(9,5))
                data=[r for r in blockrows if r['valid_observables']==len(OBS)]
                if data:
                    sizes=sorted({r['L'] for r in data})
                    ax.boxplot([[r[key] for r in data if r['L']==L] for L in sizes],
                               positions=np.arange(len(sizes)),showfliers=True)
                    ax.set_xticks(np.arange(len(sizes)),[str(L) for L in sizes])
                    if key=='n_blocks': ax.axhline(a.min_blocks,color='red',ls='--',label='Minimum blocks');ax.legend()
                ax.set(xlabel='L',ylabel=ylabel,title=title)
                ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/filename,dpi=180);plt.close(fig)
            fig,ax=plt.subplots(figsize=(8,5))
            vals=[int(r['n_blocks']) for r in blockrows if r['valid_observables']==len(OBS)]
            if vals:
                ax.hist(vals,bins=min(40,max(5,len(set(vals)))))
                ax.axvline(a.min_blocks,color='red',ls='--',label=f'Minimum {a.min_blocks}')
                ax.legend()
            ax.set(xlabel='Number of nonoverlapping blocks',ylabel='Series count',title='Block availability')
            ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'n_blocks_histogram.png',dpi=180);plt.close(fig)
    print(f'Results saved to {out}')

if __name__=='__main__':
    main()
