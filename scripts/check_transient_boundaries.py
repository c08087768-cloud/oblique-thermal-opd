"""Add/recompute Robin checks for a saved transient run, then validate the run.

Reuses saved coefficients and roots. Does not rerun projection, OPD or FEM,
and does not replace the original fields/summary/provenance.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from thermal_opd.run_io import provenance
from thermal_opd.run_validation import validate_saved_run
from thermal_opd.transient_validation import model_from_saved_transient,write_transient_boundary_check

side_z_points = 65
relative_tolerance,absolute_tolerance_W_per_m2 = 1e-10,1e-9
repository_root = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',required=True,type=Path)
    args = parser.parse_args()
    directory = args.results.resolve()
    summary = json.loads((directory/'summary.json').read_text(encoding='utf-8'))
    if summary['status']!='COMPLETE':
        raise ValueError('Need a completed transient run')
    print(f'Checking saved transient boundaries: {directory}',flush=True)
    with np.load(directory/'fields.npz',allow_pickle=False) as saved:
        model = model_from_saved_transient(summary,saved)
        keys = ('temperature_r_m','temperature_phi_rad') if summary['kind']=='transient_fine_fields' else ('r_m','phi_rad')
        write_transient_boundary_check(directory,model,saved['times_s'],saved[keys[0]],saved[keys[1]],
            np.linspace(0,model.parameters.thickness,side_z_points),angle_deg=summary['angle_deg'],
            relative_tolerance=relative_tolerance,absolute_tolerance=absolute_tolerance_W_per_m2,
            provenance=provenance(repository_root,__file__))
    try:
        report = validate_saved_run(directory)
    except (ValueError,OSError,KeyError,IndexError,TypeError) as error:
        report = dict(status='FAIL',error=str(error),run_directory=str(directory))
        (directory/'check_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        raise SystemExit(f'FAIL: {error}') from error
    (directory/'check_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS: {directory}\nBoundary table: {directory/"boundary_residuals.csv"}',flush=True)


if __name__=='__main__':
    main()
