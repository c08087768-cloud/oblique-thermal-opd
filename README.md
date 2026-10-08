# Oblique thermal OPD

有限圆柱基底上的斜入射 Gaussian 表面加热：Fourier–Bessel 半解析稳态温度、沿基底内折射光线的热折射 OPD，以及独立三维 FEM 对照。内部统一使用 SI 单位，前表面 `z=0`，后表面 `z=L`，入射平面为 y-z。

本仓库服务于同一篇论文及其补充材料。`manuscript/`、`Notes/` 和生成结果已被 Git 忽略。已实现稳态、零初始温升的瞬态、方位角贡献分析和独立 FEM；阶跃瞬态截断已在指定角度/时刻按0.1%目标选择，FEM最终时空收敛仍需加密验证。

## 安装

在仓库目录并激活自己的 Python 环境后运行：

```bash
python -m pip install -e ".[dev,fem]"
python -m pytest -q -p no:cacheprovider
```

基础模型需要 NumPy/SciPy；FEM 使用 scikit-fem 和 Matplotlib，不需要 Gmsh、FEniCS 或外部 FEM 文件。Windows/ Linux 使用同一套 Python 入口。

## 参数位置

- `src/thermal_opd/parameters.py`：几何、材料、热损失、功率、光学参数。
- 运行脚本顶部：M、N、投影网格、温度/OPD 采样、积分点、FEM 网格和误差目标。没有另设 YAML/TOML 算例配置。
- `dn_dT` 保留论文原始表值；背景折射率是标注来源的名义值。改变材料、光斑或角度范围后须重新选择数值参数。
- 瞬态新增密度 2203 kg/m³、比热 772 J/(kg K)，来自 [Heraeus 石英玻璃板资料](https://www.heraeus-conamic.com/media/Media/Documents/Products_and_Solutions/SEMI/EN/Quartz_Glass_Plates_EN.pdf)。比热为 0–100°C 平均值，采用常数近似；这两项并非原论文给出的参数。

## 当前计算流程

```bash
# 1. 0/30/60 度的参考稳定性、截断组合与 OPD 采样选择
python scripts/select_truncation.py --output-dir results/truncation_my_run
python scripts/validate_run.py --results results/truncation_my_run

# 2. 本地小批量：0/60 度，三档小 FEM 网格
python scripts/compare_fem.py --mode smoke --selection results/truncation_my_run --output-dir results/fem_smoke_my_run
python scripts/validate_run.py --results results/fem_smoke_my_run
```

`--output-dir` 指定的目录必须尚未存在；不指定则自动创建带 UTC 时间的唯一结果目录。每次运行使用新的名称。截断脚本完成后在 `summary.json` 查看 `selected_settings` 和 `selected_figure_xy`；它检查有限参考的数值变化，目标 0.1%，不等同于物理精度或 FEM 误差。

FEM 独立积分原始 Gaussian 热源，前/后/侧面使用相同 Robin 散热。网格是内接多边形圆盘的三角单元与厚度线单元张量积，形成六节点三维棱柱。光线在穿过三角边和厚度层时分段，每段使用二点 Gauss 积分；不在基底外外推。所有网格和半解析解使用记录下来的共同采样域，PV 也在同一共同孔径比较。几何偏差、因多边形而排除的光线数和每级网格变化均写入结果。

`check_report.json` 的 `PASS` 表示文件、指标、能量和运行检查通过。**smoke 的 PASS 不代表最终 FEM 网格已经收敛。** `summary.json` 的 `mesh_convergence_established` 单独报告 full 网格变化是否通过目标；这一判断仍基于有限网格。

## 之后在服务器统一重跑

进入服务器仓库并激活环境后，运行下面的 Bash 块。它只执行当前环境中的程序；不自动 SSH、部署、提交或推送。截断和 FEM 都会重新计算，日志位于结果目录外。

```bash
set -euo pipefail
python -m pip install -e ".[dev,fem]"
python -m pytest -q -p no:cacheprovider
mkdir -p results
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
TRUNCATION="results/truncation_${RUN_ID}"
FEM_RUN="results/fem_full_${RUN_ID}"
python -u scripts/select_truncation.py --output-dir "$TRUNCATION" 2>&1 | tee "${TRUNCATION}.run.log"
python scripts/validate_run.py --results "$TRUNCATION"
python -u scripts/compare_fem.py --mode full --selection "$TRUNCATION" --output-dir "$FEM_RUN" 2>&1 | tee "${FEM_RUN}.run.log"
python scripts/validate_run.py --results "$FEM_RUN"
```

full 使用 0/30/60 度和三档更细网格，计算量显著高于本地 smoke。未实测 full，不给运行时间承诺。完整结果需检查 `comparison.csv`、`summary.json`、`check_report.json` 和各角度误差图。

## 其他入口

- `run_steady_case.py`：单个完整稳态温度/折射 OPD 算例、边界残差和能量平衡。
- `check_convergence.py`：60 度六种设置的逐项诊断，原脚本不做联合参数选择。
- `inspect_heat_source.py`、`inspect_steady_temperature.py`、`inspect_opd.py`：单环节可视化与检查。
- `inspect_azimuthal_modes.py`：固定半径的热源角向分解。
- `run_transient_case.py`：瞬态温度/折射 OPD、完整径向投影后的各方位角阶次贡献、轴向/横向/光线诊断和独立瞬态 FEM。

## 瞬态与方位角贡献

```bash
python scripts/run_transient_case.py --mode smoke --angle 60 --history step --output-dir results/transient_my_run
python scripts/validate_run.py --results results/transient_my_run
# 相同模型可切换补充材料 S4 的功率爬升过程
python scripts/run_transient_case.py --mode smoke --history linear --ramp-time 5
python scripts/run_transient_case.py --mode smoke --history exponential --ramp-time 5
```

数值设置仍在脚本顶部。观察时间为 0、0.625、2.5、20、50、75 s；半解析另按最慢衰减率检查长期稳态极限。`full` 使用脚本中的较细网格与时间步，按角度分别运行；full 名称本身不表示已经收敛。

`transient.py` 复用稳态热流投影、径向根、稳态轴向响应、温度空间求和和折射光线积分。新增轴向齐次 Robin 扩散模态，每个三维模式的衰减率为 `alpha*(lambda**2+beta**2)`。通过“稳态加载项减去各扩散模态的时间修正”加速谱求和。有限 J 对阶跃的任意接近零的正时间不保证精度，需要在最早观察时刻检查；代码在 t=0 返回严格零初值。

瞬态 FEM 与稳态共用空间矩阵和直接 Gaussian 表面加载；新增 `rho*cp` 一致质量矩阵及后向 Euler。网格与时间步分别改变，逐步记录输入、储能变化和散热。早期后表面温升可能极小，结果同时保留绝对温差和相对误差。

查看 `modal_contributions.csv` 和 `transient_modes.png` 可比较 m=0/2/4/6 的温度和 OPD，CSV 也包含稳态行。完整源、各阶场、时间、径向根及投影系数保存在 `fields.npz`；FEM 节点和时间步记录在 `fem_run*.npz`。温度贡献份额指完整热圆盘上各阶的 L2 范数平方份额，不是热量份额。OPD 是每个温度阶次沿折射光线的积分贡献，不能直接称为入口平面上的同阶 Fourier 谐波；各阶 PV 不可相加。

上述smoke/full FEM诊断输出 `transient_truncation_established=false` 与 `fem_convergence_established=false`：它提供有限参考和加密诊断，不自行执行联合截断选择，也未建立最终 FEM 时空收敛。精细fields模式读取下述完成的联合选择，单独记录截断验证范围。

精细瞬态另有与稳态相同的0.1%联合截断流程：

```bash
python scripts/select_transient_truncation.py --output-dir results/transient_selection_my_run
python scripts/validate_run.py --results results/transient_selection_my_run
python scripts/run_transient_case.py --mode fields --angle 60 --selection results/transient_selection_my_run --output-dir results/transient_fields_my_run
python scripts/validate_run.py --results results/transient_fields_my_run
```

选择覆盖0/30/60度、阶跃加载和0.625/2.5/20/50/75s，参考加大变化目标0.02%，候选变化目标0.1%，并加密共同孔径和对照连续数值PV。`fields`模式只有在选择完成、物理参数/角度/加载匹配后出图；温度采样161×180、OPD采样201×201，复用稳态四联图函数。结果在`results/<run_id>`，每时刻四联图、演化图、早期厚度图、历史曲线及模态贡献图在`outputs/figures/<run_id>`。密集历史中的中间时间不逐点声明独立参考认证。

FEM smoke/full也可传入`--selection`使用已选横向参数。精细半解析与FEM时空加密分别验证，`fields`不会自动运行大FEM批量。

瞬态 `fields/smoke/full` 还逐快照检查前面 `hT-kT_z-g(t)q_modal`、后面 `hT+kT_z` 和侧面 `hT+kT_r`，使用解析导数。`boundary_residuals.csv` 保存各时刻三面的最大绝对热流残差（W/m²）和除以满功率 Gaussian 峰值后的比例；比例值不是百分数。`boundary_checks.npz` 保存温度、导数、加载和有符号残差，`boundary_summary.json` 记录门槛与检查范围。门槛是 `1e-9 W/m² + 1e-10 × 满功率峰值热流`，与截断收敛的 0.1% 分开。

`validate_run.py` 从保存的温度/导数/加载重新计算残差，并从原始模态系数重建模型核对边界温度和实际快照，检查 CSV/JSON 一致性。新运行必须包含三面检查；没有这些文件的既有运行明确标为 `transient_boundary_check.status=NOT_RUN`。已有瞬态场结果可单独补检查：

```bash
python scripts/check_transient_boundaries.py --results results/<已有瞬态场结果目录>
python scripts/validate_run.py --results results/<已有瞬态场结果目录>
```

补检查复用保存的系数、根、物理参数、J 和功率历史，只新增边界检查文件并更新 `check_report.json`，保留原始场和计算记录。检查针对截断热源 `q_modal`；原始 Gaussian 的投影误差另算。边界成立不代表瞬态 PDE/截断或 FEM 已收敛，FEM 的整体能量检查也不等同于逐点边界通量检查。此处覆盖保存的场快照，密集历史曲线的中间采样时刻未另做三面扫描。

三维 FEM 的节点解和共同采样场保存为 NPZ；单位在键名中标注。绘图和数据均直接由运行程序生成，不手工修改数值或把不同网格的结果改名成另一组参数。
