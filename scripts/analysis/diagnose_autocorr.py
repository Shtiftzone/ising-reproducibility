#!/usr/bin/env python3
"""Read existing conf-T-run.dat files; reuse the supplied C++ AutoCorr function.
No bootstrap or multihistogram reweighting is performed. ASCII source only.
"""
import argparse
import csv
import math
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import numpy as np

DRIVER = r'''
#include <vector>
using namespace std;
#include "Autocorr.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
int main(int argc, char** argv) {
    if (argc != 5) return 2;
    ifstream in(argv[1]);
    if (!in) { cerr << "Cannot open input\n"; return 2; }
    double vol = pow(stod(argv[2]), stoi(argv[3]));
    bool square = string(argv[4]) == "square";
    vector<vector<double>> series(7);
    string line; int lineno=0;
    while (getline(in,line)) {
        ++lineno;
        if (line.find_first_not_of(" \t\r") == string::npos) continue;
        if (line[line.find_first_not_of(" \t\r")] == '#') continue;
        istringstream row(line); double x[10];
        for (int k=0;k<10;k++) {
            if (!(row >> x[k]) || !isfinite(x[k])) {
                cerr << "Invalid data at line " << lineno << "\n"; return 3;
            }
        }
        double q=(x[2]-x[3])/2.0;
        double d=x[0]+(square ? 2.0 : 1.0)*q;
        double values[7]={x[1]/vol,abs(x[0])/vol,abs(q)/vol,abs(d)/vol,
                          pow(x[0]/vol,2),pow(q/vol,2),pow(d/vol,2)};
        for (int k=0;k<7;k++) series[k].push_back(values[k]);
    }
    const char* names[7]={"energy","abs_m","abs_ec","abs_d","m_squared","ec_squared","d_squared"};
    cout << setprecision(17);
    for(int k=0;k<7;k++) {
        const auto &v=series[k]; double tau=0,err=0;
        string status="ok";
        if(v.size()<4) status="too_short";
        else {
            bool constant=true;
            for(double a:v) if(a!=v[0]) { constant=false; break; }
            if(constant) status="constant";
            else { AutoCorr(v,tau,err);
                if(!isfinite(tau)||!isfinite(err)||tau<=0) status="failed";
            }
        }
        cout << names[k] << "," << v.size() << ",";
        if(status=="ok") cout << tau << "," << err;
        else cout << "nan,nan";
        cout << "," << status << "\n";
    }
}
'''

def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True,help='Recursively scan conf-T-run.dat files below this directory.')
    p.add_argument('--autocorr-source',type=Path,required=True,help='Existing Autocorr.cpp (or .cc/.c) defining AutoCorr.')
    p.add_argument('--autocorr-header-dir',type=Path,help='Directory containing Autocorr.h; default: source directory.')
    p.add_argument('--lattice-type',choices=['square','triangular'],required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--sizes',nargs='+',type=int,help='Optional size filter.')
    p.add_argument('--dimension',type=int,default=2)
    p.add_argument('--compiler',default='g++')
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args()
    if a.dimension!=2: p.error('These EC/D definitions are for two-dimensional lattices.')
    if not a.input_dir.is_dir(): p.error('Input directory does not exist.')
    if not a.autocorr_source.is_file(): p.error('AutoCorr source does not exist.')
    files=[]
    for path in sorted(a.input_dir.rglob('conf-*.dat')):
        m=re.fullmatch(r'conf-([0-9.eE+\-]+)-(\d+)\.dat',path.name)
        if not m: continue
        L=None
        for parent in path.parents:
            sm=re.fullmatch(r'(?:size_)?(\d+)',parent.name)
            if sm: L=int(sm[1]); break
        if L is None: raise ValueError(f'Cannot infer L from {path}; expected size_L or L directory.')
        if a.sizes and L not in a.sizes: continue
        files.append((L,float(m[1]),int(m[2]),path.resolve()))
    if not files: p.error('No matching files found. Expected size_L/conf-T-run.dat or L/conf-T-run.dat.')
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    rows=[]
    with tempfile.TemporaryDirectory(prefix='autocorr_') as tmp:
        driver=Path(tmp)/'driver.cpp'; driver.write_text(DRIVER)
        exe=Path(tmp)/'autocorr_driver'
        header=(a.autocorr_header_dir or a.autocorr_source.parent).resolve()
        subprocess.run([a.compiler,'-O2','-std=c++11','-I',str(header),str(driver),str(a.autocorr_source.resolve()),'-o',str(exe)],check=True)
        for j,(L,T,run,path) in enumerate(files,1):
            result=subprocess.run([str(exe),str(path),str(L),str(a.dimension),a.lattice_type],capture_output=True,text=True)
            if result.returncode:
                raise RuntimeError(f'{path}: {result.stderr}')
            for rec in csv.reader(result.stdout.splitlines()):
                obs,n,tau,err,status=rec
                n=int(n); tau=float(tau); err=float(err)
                rows.append(dict(lattice_type=a.lattice_type,L=L,T=T,run=run,observable=obs,N=n,tau_int=tau,tau_err=err,
                    variance_inflation=2*tau if status=='ok' else math.nan,
                    N_eff_mean=n/(2*tau) if status=='ok' else math.nan,status=status,file=str(path)))
            if j%50==0 or j==len(files): print(f'Processed {j}/{len(files)} files',flush=True)
    write_csv(out/'autocorr_all.csv',rows)
    groups=defaultdict(list)
    for r in rows: groups[(r['L'],r['observable'])].append(r)
    summary=[]
    for (L,obs),g in sorted(groups.items()):
        good=[r for r in g if r['status']=='ok']; vals=sorted(r['tau_int'] for r in good)
        med=(vals[len(vals)//2] if len(vals)%2 else (vals[len(vals)//2-1]+vals[len(vals)//2])/2) if vals else math.nan
        worst=max(good,key=lambda r:r['tau_int']) if good else None
        summary.append(dict(L=L,observable=obs,files=len(g),valid=len(good),invalid=len(g)-len(good),
            T_min=min(r['T'] for r in g),T_max=max(r['T'] for r in g),tau_min=min(vals) if vals else math.nan,
            tau_median=med, tau_q95=float(np.quantile(vals, 0.95)) if vals else math.nan, tau_max=max(vals) if vals else math.nan,T_at_tau_max=worst['T'] if worst else math.nan,
            run_at_tau_max=worst['run'] if worst else '',tau_err_at_max=worst['tau_err'] if worst else math.nan,
            N_eff_min=min(r['N_eff_mean'] for r in good) if good else math.nan,
            count_tau_above_1=sum(v>1 for v in vals)))
    write_csv(out/'autocorr_summary_by_L.csv',summary)
    write_csv(out/'autocorr_invalid.csv',[r for r in rows if r['status']!='ok'])
    lines=['AutoCorr diagnostics',f'Lattice: {a.lattice_type}',f'Input files: {len(files)}',
        'Tau units: intervals between stored measurements.',
        'Independent samples: tau_int = 0.5 in the supplied estimator convention.',
        'N_eff_mean = N/(2*tau_int) is an approximate diagnostic for each series mean.',
        'Values above N are retained when estimated correlations are negative.',
        'No independence decision is made from a fixed tau threshold.',
        'Short series and the supplied window estimator can miss slow correlations.',
        'These diagnostics do not by themselves validate errors of reweighted peaks.',
        'Runs and temperatures are analyzed separately, never concatenated.',
        'Signed D: square M+2*EC_sym; triangular M+EC_sym.',
        'All series are normalized per spin; this does not change autocorrelation.',
        'Quadratic series are included because susceptibility uses second moments.',
        '', 'L observable valid/files tau_median tau_max T_at_max']
    for r in summary:
        lines.append(f"{r['L']} {r['observable']} {r['valid']}/{r['files']} {r['tau_median']:.4g} {r['tau_max']:.4g} {r['T_at_tau_max']:.8g}")
    (out/'autocorr_summary.txt').write_text('\n'.join(lines)+'\n')
    if not a.no_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print('Matplotlib unavailable; CSV and text summaries were saved.')
        else:
            observables=list(dict.fromkeys(r['observable'] for r in rows))
            for obs in observables:
                fig,ax=plt.subplots(figsize=(9,5))
                for L in sorted({r['L'] for r in rows}):
                    g=sorted([r for r in rows if r['L']==L and r['observable']==obs and r['status']=='ok'],key=lambda r:(r['T'],r['run']))
                    if g: ax.errorbar([r['T'] for r in g],[r['tau_int'] for r in g],yerr=[r['tau_err'] for r in g],fmt='o',ms=3,alpha=.65,capsize=2,label=f'L={L}')
                ax.axhline(.5,color='black',ls='--',lw=1,label='Independent-sample reference')
                ax.set(xlabel='Temperature T',ylabel='Integrated autocorrelation time',title=obs)
                ax.legend(fontsize=8,ncol=2); ax.grid(alpha=.2); fig.tight_layout()
                fig.savefig(out/f'autocorr_vs_T_{obs}.png',dpi=180); plt.close(fig)
    print(f'Results saved to {out}')

if __name__=='__main__': main()
