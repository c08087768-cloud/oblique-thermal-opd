"""Replot saved diagnostics and copy saved maps for Supplement 1; no solver runs.

Existing build_artifacts.py and verification plotters were inspected first.
Their entry points also write main-text assets or create verification runs.
This standalone presentation wrapper reads only the hashed saved inventory.
The max-over-metrics reduction matches the existing refinement plotter.
"""
from pathlib import Path
import hashlib
import json
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / 'outputs/experiment_audit/20261017T121521368371Z_f107298c'
OUT = ROOT / 'manuscript/figures'
WORK = ROOT / 'tmp/supplement_figures_20261027'
INDEX = json.loads((AUDIT / 'summary.json').read_text(encoding='utf-8'))
STAGES = {s['name']: s for s in INDEX['stages'] if s.get('run_directory')}
INPUTS, OUTPUTS, CURVES = {}, {}, {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved(stage, relative):
    name = STAGES[stage]['run_directory'] + '/' + relative
    path = AUDIT / 'data/snapshot' / name
    assert sha(path) == INDEX['files'][name], name
    INPUTS[name] = sha(path)
    return path


def read(stage, relative='summary.json'):
    return json.loads(saved(stage, relative).read_text(encoding='utf-8'))


def max_metric_percent(errors):
    return 100 * max(float(x) for values in errors.values()
                     for x in np.atleast_1d(values))


def plot_steady_sweeps(data=None):
    data = read('sweeps') if data is None else data
    assert data.get('incidence_angle_deg', 60) == 60
    fig, axes = plt.subplots(3, 2, figsize=(7, 6.6), layout='constrained')
    labels = [
        ('azimuthal_order', '(a) Azimuthal truncation', 'Maximum angular order, $M$'),
        ('radial_modes', '(b) Radial truncation', 'Radial modes per order, $N$'),
        ('projection_radial', '(c) Radial source projection', 'Radial quadrature points'),
        ('projection_angular', '(d) Angular source projection', 'Azimuthal quadrature points'),
        ('opd_sampling', '(e) Entrance-pupil sampling', 'Cartesian points per direction'),
        ('ray_quadrature', '(f) Ray integration', 'Quadrature points per ray'),
    ]
    for ax, (key, title, xlabel) in zip(axes.flat, labels):
        group = data['sweeps'][key]
        rows = group['cases']
        x = np.array([r['varied_value'] for r in rows])
        pv = np.array([r['opd_pv_m'] for r in rows]) * 1e9
        if key == 'ray_quadrature':
            y = np.array([r['relative_pv_change_vs_reference'] for r in rows])
            ax.set_ylabel('Relative PV change')
        else:
            y = pv
            ax.set_ylabel('Sampled OPD-PV (nm)')
        CURVES[key] = {'x': x.tolist(), 'y': y.tolist(),
                       'fixed_settings': group['reference_settings'],
                       'varied_setting': group['varied_setting']}
        ax.plot(x, y, 'o-', color='#246387', ms=4, lw=1.3)
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(key == 'ray_quadrature')
        formatter.set_powerlimits((-3, 4))
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        if key in ('azimuthal_order', 'radial_modes', 'ray_quadrature'):
            ax.set_xticks(x)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.set_title(title, loc='left', fontsize=10)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=.22)
    path = OUT / 'supp_steady_resolution.pdf'
    fig.savefig(path, metadata={'Title':'Saved steady resolution sweeps at 60 degrees'})
    plt.close(fig)
    OUTPUTS[path.name] = sha(path)


def plot_fem_refinement():
    verification = read('transient_fem_convergence', 'verification.json')
    assert verification['status'] == 'passed'
    checks = {c['name']: c for c in verification['checks']}
    fig, axes = plt.subplots(3, 1, figsize=(7, 6.2), sharex=True, sharey=True,
                             layout='constrained')
    styles = [('xy','Transverse mesh','#236a92','o'),
              ('z','Axial mesh','#db8230','s'),
              ('time','Time step','#2a8561','^'),
              ('joint','Joint refinement','#92519b','D')]
    for ax, angle, panel in zip(axes, (0,30,60), ('a','b','c')):
        for axis, label, color, marker in styles:
            check = checks[f'theta{angle}_{axis}_convergence']
            m = check['metrics']
            times = np.asarray(m['times_s'])
            used = np.asarray(m['used_for_gate_by_time'], dtype=bool)
            # Thermal physics is unchanged by the optical-input update. Only
            # its temperature refinement evidence is reused for current plots.
            values = np.array([100*max(float(x) for key,values in e.items()
                if key.startswith('temperature_') for x in np.atleast_1d(values))
                for e in m['errors_by_time']])
            assert np.all(np.isfinite(values)) and np.all(values > 0)
            assert np.array_equal(used, times >= .0025)
            target = float(check['threshold']) * 100
            assert np.all(values[used] <= target)
            CURVES[f'fem_{angle}_{axis}'] = {'times_s':times.tolist(),
                'max_metric_percent':values.tolist(), 'used_for_refinement':used.tolist(),
                'target_percent':target}
            ax.loglog(times[used], values[used], '-', color=color, lw=1.25,
                      marker=marker, ms=3.7, label=label)
            ax.scatter(times[~used], values[~used], marker=marker, s=27,
                       facecolor='none', edgecolor=color, zorder=4)
        ax.axhline(.5, color='#555555', lw=.9, ls='--', label='Spatial / joint target: 0.5%')
        ax.axhline(.1, color='#555555', lw=.9, ls=':', label='Time-step target: 0.1%')
        ax.set_title(f'({panel}) External angle {angle}'+r'$^\circ$', loc='left', fontsize=10)
        ax.set_ylabel('')
        ax.set_ylim(7e-6, .9)
        ax.grid(which='major', alpha=.24)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.supylabel('Largest temperature change (%)', fontsize=10)
    fig.legend(handles, labels, loc='outside upper center', ncol=2,
               frameon=False, fontsize=8)
    axes[-1].set_xlabel('Time after step turn-on (s)')
    path = OUT / 'supp_transient_fem_refinement.pdf'
    fig.savefig(path, metadata={'Title':'Saved transient FEM refinement comparisons'})
    plt.close(fig)
    OUTPUTS[path.name] = sha(path)


def copy_saved_maps():
    summary = read('transient_fields_theta60')
    for time, name in ((.625, 'supp_transient_60deg_0p625s.png'),
                       (75., 'supp_transient_60deg_75s.png')):
        frame = next(s for s in summary['snapshots'] if s['time_s'] == time)
        path = saved('transient_fields_theta60', 'figures/' + frame['figure'])
        target = OUT / name
        shutil.copyfile(path, target)
        assert sha(path) == sha(target)
        OUTPUTS[target.name] = sha(target)
        CURVES[f'snapshot_{time:g}s'] = {'time_s':time,
            'sampled_pv_nm':1e9*frame['sampled_pv_m'],
            'continuous_pv_nm':1e9*frame['continuous_pv_m']}


def main():
    # Current optical inputs: never silently rebuild the old-wavelength assets.
    from build_optical_revision import build_steady_comparison
    manifest=json.loads((ROOT/'manuscript/optical_revision_manifest.json').read_text())
    run=ROOT/manifest['run_directory']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    plot_steady_sweeps(json.loads((run/'supplement_sweeps.json').read_text()))
    plot_fem_refinement()
    build_steady_comparison(run)
    print('Updated three supplementary figures from current optical/thermal evidence.')


if __name__ == '__main__':
    main()
