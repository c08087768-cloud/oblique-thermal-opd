"""Build active manuscript figures from the approved optical-parameter run.

Usage: --run-directory outputs/optical_parameter_update/<run> --part main/all.
All input arrays stay in the isolated run; this edits presentation assets only.
Existing plotting helpers are reused. Figure5 compares a thermal approximation,
whereas Figure2 compares independent numerical methods and has two columns.
Numerical gates are in the run summaries, not inferred from visual similarity.
Added 2026-11-03. Rebuild only after all corresponding source checks pass.
"""
from pathlib import Path
import argparse
import csv
import json
import hashlib
from dataclasses import asdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from thermal_opd.plotting import save_temperature_opd_figure
from thermal_opd.parameters import PhysicalParameters

ROOT=Path(__file__).resolve().parents[2]
PAPER=ROOT/'manuscript'
COLORS={0:'#33658a',30:'#c06b25',60:'#39846a'}


def load(run,angle):
    folder=run/f'sa_theta{angle}'
    summary=json.loads((folder/'summary.json').read_text())
    if summary['physical_parameters_SI'] != asdict(PhysicalParameters()):
        raise ValueError('Saved optical fields do not match the active material inputs')
    return summary,np.load(folder/'fields.npz',allow_pickle=False)


def build_main(run):
    summary,fields=load(run,60)
    save_temperature_opd_figure(PAPER/'figures/representative_steady_60deg.png',
        fields['temperature_r_m'],fields['temperature_phi_rad'],fields['steady_temperature_K'],
        fields['opd_x_m'],fields['opd_y_m'],fields['mask'],fields['steady_opd_m'],
        title='Steady response at 60 degrees',pv_m=summary['steady_sampled_pv_nm']*1e-9,
        continuous_pv_m=summary['continuous_extrema']['steady']['pv'],figsize=(7.4,6.7))
    fields.close()
    fig,axes=plt.subplots(1,2,figsize=(9.,3.1),layout='constrained')
    for angle in (0,30,60):
        s,_f=load(run,angle);_f.close()
        with (PAPER/f'data/theta{angle}_history.csv').open(encoding='utf-8') as stream:
            thermal=[{k:float(v)for k,v in row.items()}for row in csv.DictReader(stream)]
        thermal=[row for row in thermal if 0<row['time_s']<=75.000001]
        axes[0].semilogx([v['time_s']for v in thermal],[v['center_front_K']for v in thermal],color=COLORS[angle],label=f'{angle}°')
        optical=[v for v in s['history']if v['time_s']>0]
        axes[1].semilogx([v['time_s']for v in optical],[v['sampled_pv_nm']for v in optical],color=COLORS[angle])
        obs=[v for v in optical if v['time_s']in (.001,.0025,.01,.05,.1,.625,2.5,20.,50.,75.)]
        axes[1].scatter([v['time_s']for v in obs],[v['sampled_pv_nm']for v in obs],s=16,facecolor='white',edgecolor=COLORS[angle],zorder=3)
        limit=json.loads((run/f'sa_limits_theta{angle}.json').read_text())['steady_sampled_pv_nm']
        axes[1].axhline(limit,color=COLORS[angle],ls='--',lw=.85,alpha=.75)
    for ax in axes:
        ax.set_xlabel('Time (s)');ax.set_xlim(.0008,95);ax.grid(alpha=.2)
    axes[0].set_ylabel('Center-front temperature rise (K)');axes[1].set_ylabel('OPD-PV (nm)')
    axes[0].set_title('(a) Surface temperature',loc='left');axes[1].set_title('(b) Optical-path variation',loc='left')
    axes[0].legend(frameon=False,title='External angle')
    fig.savefig(PAPER/'figures/transient_evolution.pdf');plt.close(fig)
    # The same pupil/rays are retained when m>0 temperature terms are omitted.
    s,f=load(run,60);j=int(np.flatnonzero(f['times_s']==75.)[0])
    full,zero=f['opd_m'][j]*1e9,f['opd_m0_m'][j]*1e9
    mask=f['mask'];x=f['opd_x_m']*1e3;y=f['opd_y_m']*1e3
    delta=full-zero;delta[mask]-=delta[mask].mean()
    fig,axes=plt.subplots(2,2,figsize=(9.2,6.8),layout='constrained')
    lo=min(np.nanmin(full),np.nanmin(zero));hi=max(np.nanmax(full),np.nanmax(zero))
    for ax,values,label in [(axes[0,0],full,'(a) Full OPD'),(axes[0,1],zero,'(b) OPD from $m=0$ temperature')]:
        im=ax.pcolormesh(x,y,np.ma.masked_invalid(values),shading='auto',cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
        ax.set_title(label,loc='left',fontsize=10);fig.colorbar(im,ax=ax,label='OPD (nm)',shrink=.8)
        ax.plot(0,0,'k+',ms=5,mew=.9)
    amp=float(np.nanmax(abs(delta)))
    im=axes[1,0].pcolormesh(x,y,np.ma.masked_invalid(delta),shading='auto',cmap='coolwarm',vmin=-amp,vmax=amp,rasterized=True)
    axes[1,0].set_title('(c) Effect of omitted thermal components',loc='left',fontsize=10)
    fig.colorbar(im,ax=axes[1,0],label='Piston-removed difference (nm)',shrink=.8)
    for ax in (axes[0,0],axes[0,1],axes[1,0]):
        ax.set_aspect('equal');ax.set_xlabel('$x_0$ (mm)');ax.set_ylabel('$y_0$ (mm)')
        ax.set(xlim=(-11.8,11.8),ylim=(-11.8,11.8))
    ax=axes[1,1]
    for angle in (0,30,60):
        data,_f=load(run,angle);_f.close();rows=[v for v in data['history']if v['time_s']>0]
        ax.semilogx([v['time_s']for v in rows],[100*v['rms_over_full_pv']for v in rows],color=COLORS[angle],label=f'{angle}°')
    ax.set_title('(d) Normalized RMS difference',loc='left',fontsize=10)
    ax.set(xlabel='Time (s)',ylabel='$E_{\\mathrm{nonax}}$ (%)',xlim=(.0008,95),ylim=(-.2,9))
    ax.grid(alpha=.2);ax.legend(frameon=False,title='External angle',loc='upper right')
    for ax in fig.axes:
        ax.tick_params(labelsize=12)
        ax.xaxis.label.set_size(12);ax.yaxis.label.set_size(12)
        ax.set_title(ax.get_title(loc='left'),loc='left',fontsize=12)
    axes[1,1].get_legend().set_title('External angle',prop={'size':11})
    fig.savefig(PAPER/'figures/azimuthal_omission.pdf');plt.close(fig);f.close()


def build_steady_comparison(run):
    with np.load(run/'steady_fem_theta60.npz',allow_pickle=False) as f:
        fig,axes=plt.subplots(2,3,figsize=(10,6.4),layout='constrained')
        for row,kind in enumerate(('temperature','opd')):
            mask=f['temperature_mask'] if row==0 else f['mask']
            x=f['temperature_x_m']if row==0 else f['opd_x_m']
            y=f['temperature_y_m']if row==0 else f['opd_y_m']
            raw=[f['temperature_sa_K'][0],f['temperature_fem_K'][0]] if row==0 else [f['opd_sa_m']*1e9,f['opd_fem_m']*1e9]
            sa=np.full(mask.shape,np.nan);fe=sa.copy();sa[mask],fe[mask]=raw
            diff=fe-sa
            if row==1:diff[mask]-=diff[mask].mean()
            amp=np.nanmax(abs(diff));lo=min(np.nanmin(sa),np.nanmin(fe));hi=max(np.nanmax(sa),np.nanmax(fe))
            for col,(values,title)in enumerate(zip((sa,fe,diff),('Semi-analytical','FEM','FEM - SA'))):
                ax=axes[row,col];im=ax.pcolormesh(x*1e3,y*1e3,np.ma.masked_invalid(values),shading='auto',rasterized=True,
                    cmap='coolwarm'if col==2 else ('inferno'if row==0 else 'viridis'),vmin=-amp if col==2 else lo,vmax=amp if col==2 else hi)
                ax.set_aspect('equal');ax.set_title(title,fontsize=10);ax.set_xlabel('$x_0$ (mm)'if row else '$x$ (mm)')
                if col==0:ax.set_ylabel('$y_0$ (mm)'if row else '$y$ (mm)')
                fig.colorbar(im,ax=ax,shrink=.78,label='nm'if row else 'K')
        fig.savefig(PAPER/'figures/steady_fem_comparison.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-directory',type=Path)
    parser.add_argument('--part',choices=('main','all'),default='all')
    args=parser.parse_args()
    run=args.run_directory
    if run is None:
        manifest=json.loads((PAPER/'optical_revision_manifest.json').read_text())
        run=ROOT/manifest['run_directory']
        for name,expected in manifest.get('new_data_sha256',{}).items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
                raise ValueError('Recorded optical source file changed: '+name)
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False}):
        build_main(run)
        if args.part=='all':
            from build_transient_comparison import make_transient_comparison
            from build_supplement_figures import plot_steady_sweeps,plot_fem_refinement
            make_transient_comparison()
            plot_steady_sweeps(json.loads((run/'supplement_sweeps.json').read_text()))
            plot_fem_refinement()
            build_steady_comparison(run)
    print('FIGURES_COMPLETE',args.part,flush=True)


if __name__=='__main__':
    main()
