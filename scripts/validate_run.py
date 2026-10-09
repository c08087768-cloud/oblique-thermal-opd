"""Validate a saved truncation/FEM run; writes PASS only after reloading fields."""

import argparse
import json
from pathlib import Path
from thermal_opd.run_validation import validate_saved_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True, type=Path)
    parser.add_argument('--require-scientific-pass',action='store_true',help='Exit nonzero if any scientific check failed or was not run')
    args = parser.parse_args()
    try:
        report = validate_saved_run(args.results)
    except (ValueError, OSError, KeyError, StopIteration, IndexError, TypeError) as error:
        report = {'status': 'FAIL', 'error': str(error), 'run_directory': str(args.results.resolve())}
        if args.results.is_dir():
            (args.results/'check_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        raise SystemExit(f'FAIL: {error}') from error
    (args.results/'check_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    (args.results/'verification.json').write_text(json.dumps(report['verification'],indent=2)+'\n',encoding='utf-8')
    science=report['verification']['status']
    print(f"PASS (saved-file integrity): {args.results}; scientific status={science}")
    if args.require_scientific_pass and science!='passed':
        raise SystemExit('Scientific verification is failed or incomplete; inspect verification.json')


if __name__ == '__main__':
    main()
