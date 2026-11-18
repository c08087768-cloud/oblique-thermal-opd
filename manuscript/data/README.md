# Manuscript data inputs

`theta0_history.csv`, `theta30_history.csv`, and `theta60_history.csv` are existing historical solver exports, now versioned without changing their bytes. Each contains 131 time rows. The optical rebuild reads their positive times through 75 s; the figure builder reads the time and front-center temperature columns. Their historical `sampled_PV_nm` column predates the matched optical parameters and must not be used as the current paper OPD result. Current values are in `supplement_details/pupil_pv.csv` and the optical revision manifest.

The optical rebuild is `scripts/cases/recompute_paper_optics.py --run-directory outputs/optical_parameter_update/<new-run> --angle 60 --stage sa`; create the new parent run directory first. Select 0 or 30 for the other cases. The SA stage requires these versioned CSV inputs and installed project dependencies; it does not require a previous outputs archive. Other stages that replay archived FEM solutions and build the complete paper still require the archives named in `manuscript/optical_revision_manifest.json`.

The CSV time grids were checked against the completed optical run `20261102T211531Z`. This is a provenance and input-availability repair, not a new numerical validation of the historical columns.

| Input | SHA-256 |
| --- | --- |
| theta0_history.csv | 27820d4f84da40abc6a841c2754ed552e9974fbe92e4ed035b9f1e564ff52530 |
| theta30_history.csv | 8be0f746da813c29ee8b8b000f787dbf71ce2ef9b38d8089f4fc46144f13b170 |
| theta60_history.csv | 2694b40fb828dd92744e6a3b05cd2e580268aec90fc9ea3444fc005b88ebeb02 |
