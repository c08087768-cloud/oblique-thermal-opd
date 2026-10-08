"""Unique run directories and reproducibility metadata for command-line studies."""

from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
import hashlib
import json
import platform
import subprocess


def create_run(repository_root, prefix, output_dir=None):
    root = Path(repository_root)
    run_id = f'{prefix}_{datetime.now(timezone.utc):%Y%m%dT%H%M%S%fZ}'
    path = Path(output_dir).resolve() if output_dir else root/'results'/run_id
    path.mkdir(parents=True, exist_ok=False)
    return path


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


def write_summary(directory, summary):
    # allow_nan=False prevents apparently complete JSON containing NaN/Inf.
    path = Path(directory)/'summary.json'
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
