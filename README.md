# Thermo-optic response of cylindrical optics under oblique Gaussian heating

[简体中文](README.zh-CN.md) · [Archived Chinese guide](README.legacy.zh-CN.md)

This repository computes the steady and transient temperature rise in a finite cylinder subject to oblique Gaussian surface heating and Robin heat exchange. The semi-analytical solver expands the thermal field in Fourier–Bessel modes. A separately assembled finite-element model (FEM) on a prismatic mesh integrates the Gaussian load directly for comparison. The optical path difference (OPD) is the first-order refractive-index contribution integrated along fixed, straight refracted rays that reach the rear face.

## Requirements and installation

Use Python 3.10 or later and Git. The core solver uses NumPy and SciPy. The installation below also provides Matplotlib for plots, scikit-fem for finite-element runs, and pytest for software tests. Run the commands from the repository root. [Dependency specification](pyproject.toml)

**Linux / WSL (Bash)**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,fem]"
python -m pip check
```

**Windows (PowerShell)**

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev,fem]"
.\.venv\Scripts\python.exe -m pip check
```

On Windows, use `.\.venv\Scripts\python.exe` in place of `python` in the commands below.

## Repository structure

| Path | Contents |
| --- | --- |
| [`src/thermal_opd/`](src/thermal_opd/) | Heat-source projection, modal and FEM solvers, OPD integration, and validation |
| [`scripts/`](scripts/) | Single cases, convergence studies, independent benchmarks, diagnostics, and replotting |
| [`run.py`](run.py), [`run_all.py`](run_all.py) | Interactive experiment launcher and batch workflow |
| [`tests/`](tests/) | Software regression tests |
| [`manuscript/`](manuscript/) | LaTeX sources and selected derived data used in the manuscript |

## Reproducing the numerical experiments

These commands start new runs. Replaying the exact archived manuscript comparisons requires the run directories listed in the [optical revision manifest](manuscript/optical_revision_manifest.json).

### Steady example at 60°

Choose an unused output directory name if rerunning the example.

```bash
python scripts/cases/run_steady_case.py --angle 60 --output-dir outputs/steady_case/example-60deg
python scripts/verification/validate_run.py --results outputs/steady_case/example-60deg
```

The script saves the temperature and OPD arrays in `data/fields.npz` and a plot in `figures/steady_case.png` under the run directory. The sampled peak-to-valley OPD is stored in metres as `opd_pv_m` in `summary.json`. The validator reports saved-file integrity and scientific status separately. [Saved-run validator](scripts/verification/validate_run.py)

The nominal cylinder has a radius of 12.7 mm and a thickness of 5 mm. The Gaussian beam's 1/e² intensity radius normal to propagation is 0.30 times the cylinder radius. The assigned absorbed power is 110 mW on the full plane. The optical constants use a reference wavelength of 1064 nm and a temperature of 300 K. The example retains azimuthal orders 0–6 and eight radial modes per order. It uses 400 radial by 32 angular projection samples, 16 ray quadrature points, and a 201 × 201 entrance-plane grid for OPD with an aperture radius of 0.9 times the cylinder radius. [Physical inputs](src/thermal_opd/parameters.py) · [Steady case settings](scripts/cases/run_steady_case.py)

### Batch workflow

```bash
python run_all.py --mode smoke --angles 60
python run_all.py --mode full --angles 0 30 60 --plan
python run_all.py --mode full --angles 0 30 60
```

The smoke batch runs the software tests, the independent uniform-flux FEM benchmark, and a small transient comparison with boundary checks. The full batch adds truncation selection, detailed fields, and FEM refinement studies. Use `--plan` to inspect the stages without starting calculations. The full run can be computationally demanding and may stop if a prerequisite fails. Its `summary.json` reports the aggregate scientific status. For a full batch, exit code `2` means execution finished with incomplete or failed scientific evidence. [Batch implementation](run_all.py)

The experiment drivers generate the heat loads, meshes, and sampling grids from their inputs. The workflows above require no random seed or external input dataset. Individual experiments record settings and provenance in `config.json` and scientific checks in `verification.json`. [Run metadata](src/thermal_opd/run_io.py)

## Recorded results and data

The 60° semi-analytical case gives a sampled peak-to-valley OPD of **135.07 nm** over entrance-plane points whose rays reach the rear face. The steady example reproduced the tabulated value with a saved-file integrity status of `PASS` and an aggregate scientific status of `not_run`. Separate finite-reference checks report differences for 0°, 30°, and 60°. [Steady OPD table](manuscript/data/supplement_details/pupil_pv.csv) · [Semi-analytical checks](manuscript/data/supplement_details/sa_checks.csv)

The repository includes selected manuscript tables and derived optical data. Generated run directories under `outputs/` are excluded from Git. For the updated optical model, the manifest records semi-analytical reference checks and cross-method comparisons with FEM. A new FEM refinement check for the updated optical quantities falls outside that evidence. [Derived optical data](manuscript/data/optical_update_20261103/) · [Optical revision manifest](manuscript/optical_revision_manifest.json)
