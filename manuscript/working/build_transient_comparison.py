"""Replot the saved 60-degree OPD frames as a two-method, two-time comparison.

Run this file to rebuild only figures/transient_fem_comparison.pdf. The saved
export must match its manifest hash. No interpolation or solver is used; the
original common-pupil mask, raw OPD values and per-row shared scales are kept.
The return value preserves the parent assembler's source/diagnostic records.
Check the saved-data hash, four panel arrays and row scales before delivery.
Updated 2026-11-03: remove only the displayed residual column at author request.
"""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def make_transient_comparison(out=None, source_directory=None):
    out = Path(out) if out is not None else Path(__file__).resolve().parents[1]
    data_dir = Path(source_directory) if source_directory is not None else out / 'data/optical_update_20261103'
    path = data_dir / 'transient_60deg_comparison_frames.npz'
    manifest = json.loads((data_dir / 'transient_60deg_export_manifest.json').read_text())
    source_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert source_hash == manifest['export_sha256'] and manifest['source_files_unchanged']
    numbers = {}
    with np.load(path, allow_pickle=False) as data, plt.rc_context({
        'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
        'axes.spines.top': False, 'axes.spines.right': False,
    }):
        x, y, pm, cm = [data[k] for k in (
            'opd_x_m', 'opd_y_m', 'opd_physical_mask', 'opd_fem_common_mask')]
        assert np.array_equal(data['times_s'], [0.625, 75.0])
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), constrained_layout=True)
        for row, time in enumerate(data['times_s']):
            sa = np.full(pm.shape, np.nan)
            fe = sa.copy()
            sa[pm] = data['opd_sa_m'][row] * 1e9
            fe[cm] = data['opd_fem_m'][row] * 1e9
            sa[~cm] = np.nan
            lo = min(np.nanmin(sa), np.nanmin(fe))
            hi = max(np.nanmax(sa), np.nanmax(fe))
            for col, (field, label) in enumerate(zip((sa, fe), ('Semianalytical', 'FEM'))):
                ax = axes[row, col]
                mesh = ax.pcolormesh(x*1e3, y*1e3, np.ma.masked_invalid(field),
                                     cmap='viridis', vmin=lo, vmax=hi,
                                     shading='auto', rasterized=True)
                ax.set_aspect('equal')
                ax.set(xlim=(-11.8, 11.8), ylim=(-11.8, 11.8), xlabel='$x_0$ (mm)')
                if col == 0:
                    ax.set_ylabel('$y_0$ (mm)')
                ax.set_title(f'{label}\n$t={time:g}$ s', fontsize=10)
                fig.colorbar(mesh, ax=ax, shrink=.78, pad=.025, label='nm')
            # Retain the existing read-only diagnostic; it is no longer a panel.
            diff = fe - sa
            diff[cm] -= diff[cm].mean()
            denom = np.ptp(sa[cm])
            assert denom > 0
            numbers[f'export_60deg_{time:g}s'] = {
                'opd_pv_sa_nm': denom, 'opd_pv_fem_nm': np.ptp(fe[cm]),
                'piston_removed_linf_over_sa_pv': np.nanmax(abs(diff)) / denom,
            }
        fig.savefig(out / 'figures/transient_fem_comparison.pdf')
        plt.close(fig)
    return {path.relative_to(out.parent).as_posix(): source_hash}, numbers


if __name__ == '__main__':
    make_transient_comparison()
