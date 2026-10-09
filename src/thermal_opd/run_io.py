"""Unique run directories and reproducibility metadata for command-line studies."""

from datetime import datetime, timezone
from dataclasses import asdict,is_dataclass
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import re
from uuid import uuid4


def create_run(repository_root, prefix, output_dir=None):
    root = Path(repository_root)
    if not re.fullmatch(r'[A-Za-z0-9_-][A-Za-z0-9_.-]*',prefix):
        raise ValueError('Run label must be a simple name; decimal angle labels are allowed')
    experiment_type=prefix
    for start,category in (('steady_refracted','steady_case'),('transient_fields','transient_case'),
                           ('transient_smoke','fem_comparison'),('transient_full','fem_comparison'),
                           ('fem_smoke','fem_comparison'),('fem_full','fem_comparison'),
                           ('transient_truncation','convergence'),('truncation','convergence')):
        if prefix.startswith(start):
            experiment_type=category
            break
    now=datetime.now(timezone.utc)
    run_id=f'{now:%Y%m%dT%H%M%S%fZ}_{uuid4().hex[:8]}'
    path = Path(output_dir).resolve() if output_dir else root/'outputs'/experiment_type/run_id
    path.mkdir(parents=True, exist_ok=False)
    (path/'data').mkdir();(path/'figures').mkdir()
    (path/'run.json').write_text(json.dumps(dict(schema_version=2,experiment_type=experiment_type,
        run_label=prefix,created_at_utc=now.isoformat(),run_id=run_id),indent=2)+'\n',encoding='utf-8')
    return path


def data_directory(directory):
    """Numeric artifacts for new runs; preserve reading legacy flat runs."""
    path=Path(directory)
    return path/'data' if (path/'run.json').is_file() else path


def figure_directory(directory):
    path=Path(directory)
    return path/'figures' if (path/'run.json').is_file() else path


def parse_output_directory(description=None):
    """Small shared option for inspection scripts; no second execution path."""
    import argparse
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('--output-dir',type=Path)
    return parser.parse_args().output_dir


def script_defaults(namespace):
    """Snapshot public scalar/dataclass settings; no configuration framework."""
    def plain(value):
        if is_dataclass(value) and not isinstance(value,type):
            return asdict(value)
        if isinstance(value,(str,int,float,bool,type(None))):
            return value
        if isinstance(value,(tuple,list)):
            return [plain(v) for v in value]
        if isinstance(value,dict):
            return {str(k):plain(v) for k,v in value.items()}
        raise TypeError
    settings={}
    for key,value in namespace.items():
        if key.startswith('_') or key in ('parameters','repository_root'):
            continue
        try:
            settings[key]=plain(value)
        except TypeError:
            continue
    return settings


def provenance(repository_root, script):
    root = Path(repository_root)
    files = [Path(script), *sorted((root/'src/thermal_opd').glob('*.py'))]
    versions = {}
    for name in ('numpy', 'scipy', 'scikit-fem', 'matplotlib'):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    revision = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=False).stdout.strip()
    changed = subprocess.run(['git', '-C', str(root), 'status', '--porcelain'],
                             capture_output=True, text=True, check=False).stdout.strip()
    return {'python': platform.python_version(), 'platform': platform.platform(),
            'packages': versions, 'git_revision': revision or None,
            'working_tree_dirty': bool(changed),
            'source_sha256': {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}


def write_summary(directory, summary, *, defaults=None):
    # allow_nan=False prevents apparently complete JSON containing NaN/Inf.
    if defaults is not None:
        summary={**summary,'script_defaults':defaults}
    path = Path(directory)/'summary.json'
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not (Path(directory)/'run.json').is_file():
        return
    configuration_keys={'physical_parameters_SI','numerical_settings','modal_settings','reference_settings','reference',
        'enlarged_reference','history','angle_deg','angles_deg','incidence_angle_deg','incidence_angles_deg','axial_counts',
        'selected_axial_modes','times_s','precision_scope','aperture_radius_ratio','tolerance','reference_stability_tolerance',
        'provenance','config','kind','mode','script_defaults','selection_directory','verification_references'}
    config={key:value for key,value in summary.items() if key in configuration_keys}
    (Path(directory)/'config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if summary.get('status')=='COMPLETE':
        metrics={key:value for key,value in summary.items() if key not in configuration_keys and key not in ('source_sha256','artifacts','limits')}
        (Path(directory)/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        verification_path=Path(directory)/'verification.json'
        if not verification_path.exists():
            from thermal_opd.verification import verification_report,verification_check
            report=verification_report([verification_check('scientific_verification','not_run',
                scope='Run completed; execute or attach scoped independent scientific checks separately')],config=config,run_directory=directory)
            verification_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
