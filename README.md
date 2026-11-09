# Oblique thermal OPD

研究有限圆柱基底在斜入射 Gaussian 表面加热下的稳态、瞬态温升和热折射 OPD。半解析方法采用 Fourier–Bessel 空间展开，FEM 独立积分原始热源。内部使用 SI 单位，前表面 z=0、后表面 z=L，入射平面为 y-z。温度接口返回相对环境的**温升**。

第一次使用，从 [配置环境](#1-配置-python-环境) 和 [运行稳态算例](#2-运行第一个稳态算例) 开始。已有结果可以直接[复核或重绘](#4-查看复核和重绘结果)。

## 1. 配置 Python 环境

从仓库根目录运行命令，也就是能看到 `pyproject.toml` 和 `run.py` 的目录。项目要求 **Python 3.10 或更新版本**。若还没装 Python，可从 [Python 官网](https://www.python.org/downloads/) 安装。

### Windows（PowerShell）

如果在 VS Code 里操作，先用「文件 → 打开文件夹」打开仓库，再打开「终端 → 新建终端」。已有仓库可以直接进入其所在目录。下面以 `E:\oblique-thermal-opd` 为例。

```powershell
Set-Location 'E:\oblique-thermal-opd'
py --version
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev,fem]"
.\.venv\Scripts\python.exe -m pip check
```

如果 `py` 找不到、但 `python --version` 正常，把前两处 `py` 换成 `python`。上面的命令直接调用虚拟环境里的 Python，适用于 PowerShell 限制运行激活脚本的情况。也可以在每次打开新终端后执行 `.\.venv\Scripts\Activate.ps1`，激活后就可以用短命令 `python`。 [Python 虚拟环境文档](https://docs.python.org/3/library/venv.html#how-venvs-work)

在 VS Code 里可用 `Ctrl+Shift+P` 打开命令面板，运行 **Python: Select Interpreter**，选项目里的 `.venv\Scripts\python.exe`。

### Linux（Bash，包括 WSL）

```bash
python3 --version
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,fem]"
python -m pip check
```

每次开新终端后重新激活 `.venv`。Linux 和 Windows 各自创建自己的虚拟环境。若 Linux 提示缺少 `venv` 或 `ensurepip`，使用发行版的包管理器安装与 Python 版本对应的 venv 组件。

### 安装了哪些依赖

| 用途 | 安装命令 |
|---|---|
| 跑算例、画图、运行 FEM 和软件测试 | `python -m pip install -e ".[dev,fem]"` |
| 跑半解析算例、画图或运行无需 FEM 的测试 | `python -m pip install -e ".[dev]"` |
| 只从 Python 调用半解析模型 | `python -m pip install -e .` |

基础依赖是 NumPy 和 SciPy。`dev` 增加 pytest 和 Matplotlib，`fem` 增加 scikit-fem 和 Matplotlib。FEM 使用代码生成的网格，无需另装 Gmsh 或 FEniCS。依赖范围见 [pyproject.toml](pyproject.toml)。`-e` 让 Python 直接使用当前仓库的源码。 [pip 可编辑安装说明](https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs)

若遇到 `No module named 'thermal_opd'`，先用 `python -c "import sys; print(sys.executable)"` 检查当前解释器，再用同一个解释器从仓库根目录执行安装。若 Windows 控制台输出中文时报编码错误，可以用 `python -X utf8 run.py --list`。

如果 pip 配置的镜像找不到构建依赖，可以在安装命令里加上 `--index-url https://pypi.org/simple`，为这次安装指定官方 PyPI。

## 2. 运行第一个稳态算例

从仓库根目录运行下面一条命令。Windows 即使没有激活虚拟环境，也可以复制执行。

```powershell
.\.venv\Scripts\python.exe scripts/cases/run_steady_case.py --angle 60
```

Linux 或已激活虚拟环境的终端使用以下命令。后面的示例也用 `python`，Windows 未激活时把它换成 `.\.venv\Scripts\python.exe`。

```bash
python scripts/cases/run_steady_case.py --angle 60
```

算例计算温升和折射光线的 OPD。运行结束后，终端会打印 `Results:` 和 `Figure:`，对应本次运行目录和结果图。打开 `outputs/steady_case/<本次运行ID>/figures/steady_case.png` 看图，回到本次运行目录打开 `summary.json` 看摘要。第一次可以只跑这一条命令，它使用脚本自带的物理和数值设置。

需要换角度时，把 `--angle 60` 改成 `--angle 30` 等，单位是度。用脚本的 `--help` 查看可以从命令行传入的参数。几何、材料、热源等名义参数在 [parameters.py](src/thermal_opd/parameters.py)。稳态的模态截断、投影网格和光线积分点数在 [run_steady_case.py](scripts/cases/run_steady_case.py) 顶部。物理参数改变后需要重新检查数值收敛。

### 用菜单运行同一个算例

```bash
python -X utf8 run.py --list
python -X utf8 run.py
```

输入 `steady`，看到配置预览后输入 `y` 开始，直接回车取消。也可以用以下命令预览设置，或跳过菜单直接运行。

```bash
python -X utf8 run.py --experiment steady --show-config --angle 60
python -X utf8 run.py --experiment steady --yes --angle 60
```

预览分开显示脚本默认值和传入的命令行参数。运行后的 `config.json` 记录实际设置。

## 3. 选择其他实验

脚本按用途分在 `scripts/cases/`、`scripts/convergence/`、`scripts/verification/`、`scripts/diagnostics/` 和 `scripts/tools/`。用 `run.py --list` 可以看到当前可用的实验名。

`run.py` 每次启动一个实验。`run_all.py` 按依赖顺序运行一批实验并复核保存数据，分别记录执行状态和科学验证状态。

以下均可通过 `run.py` 的对应名称选择，也可独立运行。默认自动创建唯一运行目录；指定 `--output-dir` 时必须采用本仓库 `outputs/<实验类型>/<自定义唯一ID>` 结构，且目录尚未存在。例如稳态使用 `--output-dir outputs/steady_case/trial-001`，瞬态 smoke 使用 `outputs/fem_comparison/trial-001`。不接受旧 `results/`、共享 `outputs/figures/` 或归档目录作为新输出位置。

| 交互名称 | scripts/ 中的入口 | 用途 |
|---|---|---|
| steady | cases/run_steady_case.py | 单个稳态温度/OPD、三面边界及能量 |
| transient | cases/run_transient_case.py | fields 精细场，或 smoke/full 瞬态 FEM 诊断 |
| convergence | convergence/select_steady_truncation.py | 稳态联合截断选择 |
| transient-convergence | convergence/select_transient_truncation.py | 瞬态联合截断选择 |
| semi-validation | verification/verify_semi_analytical.py | 独立有限差分 PDE、能量、初值附近、投影及已知极限 |
| fem-validation | verification/verify_fem.py | 均匀热流/绝热侧面的一维独立基准 |
| fem-comparison | verification/compare_steady_fem.py | 稳态半解析与 FEM 同条件对照 |
| transient-fem-convergence | convergence/check_transient_fem_convergence.py | 运行瞬态 P2 FEM 候选与加密参考 |
| transient-fem-validation | verification/verify_transient_fem_convergence.py | 按清单验收已保存的瞬态 FEM 节点数据 |
| validate | verification/validate_run.py | 重读保存的数据与报告 |
| boundaries | verification/check_transient_boundaries.py | 给既有瞬态场补三面边界证据 |
| source / azimuthal | diagnostics/inspect_heat_source.py / diagnostics/inspect_azimuthal_modes.py | 热源投影、单圆周方位角分解 |
| temperature / opd | diagnostics/inspect_steady_temperature.py / diagnostics/inspect_opd.py | 单环节检查与图像 |
| sweeps | convergence/check_steady_convergence.py | 六种数值设置逐项诊断，不进行联合选择 |
| replot / regression | tools/replot_case.py / tools/compare_refactor_baseline.py | 已有场重绘、重构基线比较 |

### 瞬态先跑小算例

```bash
python scripts/cases/run_transient_case.py --mode smoke --angle 60 --history step
```

需要安装 FEM 依赖。输出在 `outputs/fem_comparison/<运行ID>/`。可以用 `--history linear --ramp-time 5` 或 `--history exponential --ramp-time 5` 改为爬升加载。

### 先选参，再生成精细结果或做对照

先运行相应的选参脚本，从终端输出或 `outputs/convergence/` 找到这次选参目录。把下面引号中的占位目录换成实际路径。

```bash
python scripts/convergence/select_steady_truncation.py
python scripts/verification/validate_run.py --results "outputs/convergence/<稳态选择ID>"
python scripts/verification/compare_steady_fem.py --mode smoke --selection "outputs/convergence/<稳态选择ID>"

python scripts/convergence/select_transient_truncation.py
python scripts/cases/run_transient_case.py --mode fields --angle 60 --selection "outputs/convergence/<瞬态选择ID>"
```

`fields` 读取已完成且参数匹配的瞬态选择。默认温度网格为 161×180，OPD 网格为 201×201，会生成逐时刻温度和 OPD 图、稳态极限、演化、早期厚度、历史和模态贡献图，输出在 `outputs/transient_case/<运行ID>/`。正式瞬态 FEM `--mode full` 还需要通过独立复核的 `--fem-convergence` 证据目录，见后面的 P2 正式流程。

## 4. 查看、复核和重绘结果

```text
outputs/<experiment_type>/<UTC时间戳_唯一后缀>/
├── run.json             # 目录模式、运行标识与创建时间
├── config.json          # 物理/数值设置、脚本默认值及版本记录
├── summary.json         # 运行状态、关键结果与关联记录
├── metrics.json         # 当前运行状态及已产生的指标
├── verification.json    # 科学检查项目、阈值、范围与状态
├── data/                # NPZ、CSV 等数值证据
└── figures/             # 图像
```

自动目录采用 UTC 时间和唯一后缀，同一秒多次运行也不覆盖。创建运行时即建立以上文件和两个目录，初始运行状态为 `RUNNING`、科学状态为 `not_run`；失败或中止仍保留这些位置，指标只反映实际已保存的内容。稳态场归入 `steady_case`，精细瞬态场归入 `transient_case`，稳态/瞬态截断归入 `convergence`，FEM 对照归入 `fem_comparison`；其他检查按各自实验名称分类。

`config.json` 包含可获得的 Git revision、工作区是否有修改、源文件哈希和依赖版本；未提交的修改不能仅靠 commit hash 复现，应一并保留代码及哈希。运行完成本身不表示科学验证通过。新摘要中的文件和关联运行路径相对本次运行目录保存，跨运行引用如 `../../convergence/<ID>`；读取时以摘要所在运行目录解析，不依赖当前终端目录。不同 Windows 盘符之间无法表示相对路径时保留绝对路径。

旧结果按原仓库相对路径归档：`results/` → `outputs/legacy/results/`，旧共享图像 `outputs/figures/` → `outputs/legacy/outputs/figures/`。原文件内容、早期归档、失败及中止结果完整保留；归档数据维持原结构，不强行拆成新布局。实际移动状态见审计目录的 `data/migration_manifest.csv`；被其他程序占用而尚未移动的文件保留在原处，单独列于 `data/unmigrated_files.json`。

按本次整理时的最新约定，不提供旧路径自动跳转，也不修复归档元数据中的历史跨目录引用。旧结果可按迁移清单定位，但不保证仍能完整校验；需要新结果时修改代码后另行运行。原有平铺数组读取分支保留，`validate_run.py`、`replot_case.py` 可以处理其已支持且依赖完整的格式；使用时显式传入实际目录。旧 normal-path OPD、未完成结果等不恢复旧算法或伪装为已通过。审计清单、逐项哈希、引用例外及前后目录树位于 `outputs/migration_audit/<审计ID>/data/`。

```bash
python scripts/verification/validate_run.py --results "outputs/steady_case/<运行ID>"
python scripts/tools/replot_case.py --results "outputs/steady_case/<运行ID>"
python scripts/verification/check_transient_boundaries.py --results "outputs/transient_case/<运行ID>"
```

重绘只读取保存的 fields.npz、history.npz 和摘要，图像进入新的 replot 运行；不重算模态或修改原始场。边界补检查使用保存的根、系数、J 和功率历史；新增边界证据及检查报告。

### 软件测试和科学验证

安装 `.[dev,fem]` 后，可以运行全套软件测试。科学验证按物理条件、采样范围和门槛逐项记录 **passed / failed / not_run**。任一检查失败时，汇总为 `failed`。其余情况下，有未执行项则汇总为 `not_run`。

```bash
python -m pytest -q -p no:cacheprovider
python scripts/verification/verify_semi_analytical.py --angle 60
python scripts/verification/verify_fem.py
python scripts/verification/validate_run.py --results "outputs/steady_case/<运行ID>" --require-scientific-pass
```

`validate_run.py` 的 `PASS (saved-file integrity)` 表示已重新读取并核对保存的数组、指标和适用的运行检查；同一行另打印科学汇总状态。科学检查存在失败时，证据文件仍可以完整且一致。严格标志会在科学状态为 failed/not_run 时返回非零退出码。

统一入口按保存的 `kind` 调用稳态或瞬态验证模块。`scientific_verification` 是可包含两类物理检查的通用证据格式，其标量数组、阈值和状态由 `validation_common.py` 直接复核；旧报告不含模型子类型也可读取。新代码统一采用 validation 命名；`verification.json`、`verification_fields.npz` 和既有 JSON/NPZ 字段继续沿用原数据格式。状态记录与汇总只整理已执行检查，不由文件 I/O 判断物理结果。稳态/瞬态共享空间差分及网格几何检查，FEM 求解与载荷积分仍保持独立。

半解析独立实验包括内部差分 PDE、三面 Robin、独立体积储能/表面散热积分、零热源、初始和长期极限、固定 M/N 的投影积分加密、早期 J 加密以及短时深厚基底解析极限。热边界检查针对截断热源；它不能代替原始 Gaussian 投影误差和截断收敛。

FEM 独立实验采用均匀前面热流、绝热侧面：稳态对照精确线性轴向解，瞬态对照独立一维有限体积且精确时间演化的参考；参考也单独加密。网格和时间步分别改变，检查全部零初值节点、离散方程及逐步能量。该局部基准使用 2% 门槛，一般 Gaussian 问题的最终网格/时间和逐点边界通量仍标为 not_run。方法间接近不能代替两种方法各自验证。

## 5. Python 接口与源码导览

```python
import numpy as np
from thermal_opd import ThermalModel, ModalResolution, PhysicalParameters, PowerHistory

p = PhysicalParameters()
resolution = ModalResolution(M=6, N=8, projection_r=400,
                             projection_phi=32, ray_points=32)
model = ThermalModel(p, angle_deg=60, resolution=resolution)  # 仅描述问题
prepared = model.prepare()                                 # 径向根与热源投影
steady = prepared.solve_steady()
transient = prepared.solve_transient(num_axial_modes=48,
                                     history=PowerHistory("step"))

T_steady = steady.temperature(r=0., phi=0., z=0.)
T_transient = transient.temperature(0., 0., 0., time=2.5)
dT_dr, dT_dz = transient.temperature_gradient(0., 0., 0., time=2.5)
opd = steady.opd(np.array([0., .001]), np.array([0., -.001]))
opd_history = transient.opd_series([0.], [0.], times=[0., .625, 2.5])
```

这里的数值设置用于演示，精度须按具体角度、时间和参数另行检查。`temperature_gradient` 返回径向与轴向导数，不包含角向分量。

`ThermalModel` 固定物理参数、角度和数值设置，构造时不计算模态。`PreparedModel` 以只读映射及数组保存共同数据，两个解复用同一份准备结果；瞬态专用轴向模态只在创建瞬态解时计算。显式再次调用 `prepare()` 会重新准备，没有隐藏的全局缓存。解对象只求值，普通求解不会自动绘图、写文件或执行昂贵验证。

`TransientSolution` 直接定义在 `transient.py`，保存功率历史、轴向模态、衰减率和幅值，并执行瞬态响应与物理量求值。共同数据统一从 `solution.prepared` 读取：`parameters`、`coefficients` 和 `roots`。`SteadySolution` 位于 `solutions.py`，调用独立稳态数值函数；两个解都提供 `temperature`、`temperature_gradient`、`opd`、`prepare_opd`，瞬态另有 `opd_series`。

```text
ThermalModel.prepare()
  └─ heat_source + eigenmodes + projection → PreparedModel
       ├─ solve_steady()    → SteadySolution    → steady 轴向响应
       └─ solve_transient() → TransientSolution → transient 轴向/时间响应
            两者共享 modal_evaluation 空间求和、梯度及 opd 光线积分

实验脚本 → 模型求值 → metrics / plotting → run_io
独立验证脚本 → verification + validation + convergence → 验证证据
FEM 实验 → fem 独立组装/求解 → 独立参考或方法间对照
```

| 模块 | 职责 |
|---|---|
| parameters.py | 几何、材料、环境及光学物理参数 |
| model.py、solutions.py | 问题/共同准备、稳态求值接口；瞬态解直接位于 transient.py |
| heat_source.py、azimuthal.py、eigenmodes.py、projection.py | 热源、圆周分解、径向根与投影 |
| steady.py、transient.py、modal_evaluation.py | 各自轴向算法、功率历史、共同空间模态求和 |
| opd.py | Snell 折射、透射掩膜、分块光线积分与预备回调 |
| fem.py | 棱柱 FEM 的共享空间组装、稳态/瞬态及分段光线积分 |
| metrics.py、plotting.py | PV、去 piston RMS、绘图和读取已有场重绘 |
| steady_validation.py | 稳态边界、能量、独立 PDE、截断/FEM 和保存结果复核 |
| transient_validation.py | 瞬态边界、初值/时间/长期极限、独立参考、截断/FEM 和保存结果复核 |
| validation_common.py、convergence.py | 共用差分、网格几何、误差与保存证据核对；场误差和数值极值 |
| run_io.py、experiments.py | 文件读写、路径、报告状态元数据、实验入口目录 |

`scripts/` 下是可直接运行的实验和诊断，`src/thermal_opd/` 下是共享模型与求解器，`tests/` 下是软件测试，`manuscript/` 下是论文与补充材料。

旧瞬态求解器类及其转发包装已删除，统一使用 `prepared.solve_transient(...)`。底层 `solve_modal_source` 保留实际源准备职责；共享空间函数从 `modal_evaluation.py` 导入。没有通用 Experiment 继承体系，也没有额外 YAML/TOML 算例配置。

### 物理及数值约定

- 温度场由截断系数直接求值；热源 q 用于投影、加载与验证，不作为温度的插值来源。
- 瞬态保留稳态加载项与轴向扩散模态修正，衰减率为 alpha*(lambda²+beta²)；支持阶跃、线性爬升和指数爬升。
- OPD 始终沿基底内折射光线积分：x=x0、y=y0+z*tan(theta_t)、ds=dz/cos(theta_t)。只保留后表面透射光线，原始 OPD 保存 piston。采用未加热基底中的直线参考光线，未包含热致弯曲。
- PV 为共同有效孔径内的极差；连续数值极值与网格采样 PV 分别记录。RMS 为等间距 Cartesian 有效采样点去 piston 后的 RMS，不用于不等权积分网格。
- 方位角温度份额是各阶 L2 范数平方份额，非热量份额；各温度阶次沿折射光线的 OPD 贡献也不等于入口平面同阶 Fourier 谐波，各阶 PV 不可相加。
- parameters.py 保留 dn_dT=0.96e-5/K 的原文表值，其标注波长为 643.8 nm，而名义光学波长为 1064 nm。密度 2203 kg/m³、比热 772 J/(kg K) 的来源在源码中，比热是 0–100°C 平均值，采用常数近似。改变材料、波长、几何或光斑后需核对物理参数并重新检查数值收敛。

## 6. 批量实验与恢复

批量入口为 `run_all.py`。先用 `--plan` 查看命令和依赖顺序，这一步不启动计算，也不创建运行目录。实际运行批量流程需要安装 `.[dev,fem]`。

```bash
python run_all.py --mode full --plan
python run_all.py --mode smoke --angles 60
python run_all.py --mode full
```

默认使用 smoke。它运行软件测试、FEM 独立基准、一个小规模瞬态 FEM 算例及其边界复核。默认角度集合为 0°、30°、60°，smoke 的瞬态算例取集合中的最大角度。full 串行运行截断选择、独立验证、稳态与精细瞬态场、正式 FEM 对照、诊断和重绘，包含计算量较大的收敛实验。

中断后可用相同代码、模式和角度恢复。将引号内的 `<批次ID>` 替换成已有目录名。

```bash
python run_all.py --mode full --resume "outputs/experiment_suite/<批次ID>"
```

批次摘要和阶段清单位于 `outputs/experiment_suite/<ID>/summary.json` 和 `data/stages.csv`，日志位于 `data/logs/`。单项实验结果仍写到对应的输出分类。稳态脚本的固定设置如果与新选择不同，批次会停下来并提示更新设置。

### 代码版本接续与证据

后续步骤引用本批次生成或经显式版本接续核验继承的选择目录。

同一残差公式只提高累加算术精度、原容差及保存温度 float64 不变时，可显式创建新版本批次，保留旧批次和已经完成的科学数据：

```bash
python run_all.py --mode full --angles 0 30 60 \
  --continue-from outputs/experiment_suite/<旧批次ID> \
  --continuation-evidence <版本接续证据.json> \
  --fem-reuse-index <已完成FEM节点索引.json>
```

证据 JSON 必须包含 `kind="experiment_suite_code_continuation"`、`version=1`、`purpose="same_formula_residual_accumulation_precision"`、旧 `summary.json` 的 `parent_suite_sha256`、旧 `code_sha256` 字典的规范 JSON SHA-256（`sort_keys=True,separators=(',',':')`）作为 `parent_code_manifest_sha256`，以及 `parent_source_snapshot` 目录。程序核对快照内 `run_all.py` 的出生哈希，执行其计划构建函数（不运行主入口），要求计划、模式、角度与当前完全相同。

`changes` 数组逐项列出所有变更代码的 `{path,before_sha256,after_sha256,reason}`，新增或删除的一端用 `null`；不能漏列或额外忽略哈希。`numerical_policy` 必须为 `{"residual_formula_unchanged":true,"residual_tolerance_unchanged":true,"saved_temperature_dtype":"float64"}`。`independent_test_reports` 数组列出 `{data_file,sha256,status:"passed"}`，引用实际 JUnit XML 回归报告；程序核对文件哈希并要求有执行通过的测试，且失败及错误数均为零。这些是软件回归证据，不替代科学独立验证。证据内相对路径以证据 JSON 所在目录为基准。

只有原状态为 `COMPLETE` 且全部已记录稳定产物、阶段日志核验成功的阶段会继承；其原科学 `failed/not_run`、产物哈希和出生代码清单保持可追溯。失败阶段从新批次重新执行，FEM 完整节点只通过显式复用索引采用。新批次保存当前真实代码清单及父批次/证据版本链；之后 `--resume` 仍要求当前代码、已完成产物、继承日志和接续证据哈希一致。

计算或保存文件复核失败会停止依赖流程。恢复时会核对已完成产物的哈希，并为失败尝试创建新目录。独立科学检查失败但文件完整时，批次仍可继续收集其他实验。退出码 `0` 表示批次成功结束，`1` 表示执行或数据复核失败，`2` 表示 full 执行完成但科学证据失败或不完整，`130` 表示中断。smoke 退出码为 `0` 时也应查看科学状态。

## 进阶参考与历史记录

### FEM 接口与收敛证据

FEM 的网格配置名称不代表收敛结论，原稿二次四面体网格的 fine 结论也不能套用于当前线性棱柱网格。一般 Gaussian 问题需分别检查横向、厚度方向和联合网格加密；瞬态还需独立检查时间步，并覆盖实际使用的早期时刻。温度场和 OPD 场均按既有误差定义及阈值判定。

`assemble_thermal_fem`、`solve_steady_fem`、`solve_transient_fem` 可通过 `axial_nodes_m` 指定恰好 `axial_layers+1` 个严格递增节点，端点必须为 `0` 和实际厚度；省略时仍采用原均匀网格。瞬态还可传入 `time_step_by_interval_s`，在观测时间包含 `t=0` 时，为相邻观测区间分别指定不超过 `time_step` 的正步长上限。每个区间仍等分，默认使用原隐式欧拉和一致容量矩阵；早期与后期可使用不同步长，但都必须独立加密验证。默认调用的矩阵、解和逐步诊断已与修改前代码做完全一致的回归检查。

`solve_transient_fem(..., time_integrator='bdf2')` 是可选的二阶时间积分试验路径，默认仍为 `'backward_euler'`。它复用同一独立FEM装配，以隐式欧拉启动；按实际相邻步长计算BDF2系数，步长增长比超过1.5时重新启动。此路径尚不代表正式时空收敛已通过，也没有接入全量主流程。已知时间响应、独立一维有限体积参考和故意损坏的保存证据用于检验新方法。

BDF2保存的 `storage_rate_W`、`heat_loss_W`、`balance_absolute_error_W` 保留原有能量增量率、端点热损失及两者的差额含义；新增 `bdf_storage_rate_W`、`bdf_balance_absolute_error_W`、三个 `bdf_alpha*` 系数和 `time_order` 记录该方法的离散导数与方程平衡。保存结果复核按 `time_integrator` 选择对应时间导数，保持原1e-7能量/弱残差及0.1%时间收敛门槛，独立复算节点能量、Robin积分及系数的多项式精确性。候选与加密参考必须使用相同积分方法，不能将不同方法的比较标作单独时间步加密。

瞬态 CG 可选择 `linear_preconditioner='axial_lines'`，逐条反演同一横向节点上的轴向块，改善前表面加密网格的求解速度；默认仍为 `'jacobi'`。它不改变组装矩阵、载荷、隐式欧拉方法或 CG 容差，也不使用半解析系数。独立矩阵分块对照和直接求解器回归检查用于验证该可选预条件器。

瞬态 FEM 保存数据验收入口：

```bash
python scripts/verification/verify_transient_fem_convergence.py --manifest <验证清单.json>
python scripts/verification/validate_run.py --results <生成的验收目录> --require-scientific-pass
```

清单只指定已有节点结果及其配对，不重新运行热求解。以下展示单个角度的配对；正式清单应包含需要验收的全部角度。`source_run` 相对清单文件所在目录解析，`case_index` 对应源目录的 `data/fem_run<序号>.npz`；历史平铺目录和“终值数组+初值数组”格式也支持读取。

```json
{
  "space_tolerance": 0.005,
  "joint_tolerance": 0.005,
  "time_tolerance": 0.001,
  "diagnostic_only_times_s": [0.001],
  "cases": [{
    "angle_deg": 0,
    "candidate": {"source_run": "candidate-run", "case_index": 0},
    "xy_reference": {"source_run": "xy-refined-run", "case_index": 0},
    "z_reference": {"source_run": "z-refined-run", "case_index": 0},
    "time_reference": {"source_run": "time-refined-run", "case_index": 0},
    "joint_reference": {"source_run": "joint-refined-run", "case_index": 0}
  }]
}
```

默认验收角度为0/30/60°，时刻为0、.001、.0025、.01、.05、.1、.625、2.5、20、50、75 s；温度使用201×201完整公共热圆盘、前/中/后截面及随时间变化的扩散深度点，OPD使用独立101×101公共折射光瞳。角度、时刻和取样点数可以在清单中显式指定，此时结论仅限该范围。按2026-10-13用户新决定，瞬态FEM的横向/轴向网格与联合加密温度、OPD误差门槛为0.5%，时间步门槛仍为0.1%；1e-7方程/能量门槛、严格零初值和边界标准保持原值。t<.625 s保留原global温度范数，之后保留原逐截面范数。旧报告和显式使用旧单一 `tolerance=0.001` 的清单仍按原0.1%标准复核；新清单分别记录 `space_tolerance`、`joint_tolerance`、`time_tolerance`，不自动改写旧判定。稳态与半解析标准不变。

按2026-10-16用户决定，**0.001 s仍逐次记录温度/OPD及各方向加密差异，作为诊断时刻，不参与瞬态FEM空间、时间和联合收敛的通过判定**。新验收CLI默认写入 `diagnostic_only_times_s: [0.001]`；其他时刻及所有时间步上的初值、能量、Robin和弱残差检查保持。仅有这一诊断帧不能证明收敛，记为 `not_run`。历史保存报告缺少此字段时，旧结果读取入口仍按其原规则重放；对历史源数据执行一次新CLI复核则明确记录新政策，并保留原报告。不能借此排除其他时刻。

验证会检查物理参数/功率历史一致、实际网格和步长是否构成所声明的单独或联合加密，并重新计算温度/OPD误差、Gaussian载荷、节点储能与三面Robin散热积分。原散热诊断一致性门槛 `rtol=1e-9, atol=1e-12 W` 继续保留；方程残差取自逐步归档数据，范围明确标为FEM弱形式，不能代替逐点法向通量检查。缺少参考、角度或时间覆盖记为 `not_run`；失败的源实验不会被标为最终通过。只有全部已声明验收门槛通过，才设置 `fem_convergence_established=true`。固定载荷积分阶次的收敛、半解析截断和跨方法一致性仍需各自证据。

新报告保存到 `outputs/fem_convergence_validation/<唯一ID>/`，包含 `verification.json`、`check_report.json`、`data/comparison_fields.npz`、`data/convergence_by_time.csv` 和 `figures/convergence_by_time.png`。CSV中误差为比例（0.001代表0.1%、0.005代表0.5%），`threshold_fraction` 保存每项使用的门槛，图中以百分比显示相应目标线。复核会再次读取源节点文件、验证SHA-256并重算指标；需要保留被引用的源实验目录。

CSV新增 `used_for_gate` 与 `time_role` 标记，1ms行仍保留原数值但标为 `False/diagnostic_only`；`data/diagnostic_0p001s.json` 单独收集该帧各角度/各方向对比，图中用空心点显示。`field_values_within_target` 只表示数值是否落在对照限值内，不代表该诊断行构成验收门槛。若输入缺少该帧，诊断证据完整性记为 `not_run` 并指出来源，不能用空文件冒充已记录；这与该帧数值超限是两种情况。

2026-10-13 的当前 P1 棱柱稳态加密实验已对 0/30/60°分别检查横向、轴向和联合加密，并从保存节点做 201×201 取样复核。候选 `(disk_refinements=7, axial_layers=48)` 相对联合参考 `(8,96)` 在该稠密复核中的最大指标分别为 0.0317%、0.0239%、0.0156%，均小于 0.1%。因此稳态 FEM `full` 的三档配置更新为 `(6,48),(7,48),(8,96)`；最后一档为加密参考。此结论限于记录的名义参数、角度、温度截面和公共折射光瞳，不构成瞬态空间/时间收敛结论。瞬态仍严格保留原来的温度归一化、初值、边界、能量和 OPD 判据；不能用整体温度尺度替代原逐截面门槛来宣布通过。

### P2与镜像提速后的正式流程（2026-10-17）

生产FEM已增加可选P2场与精确双镜像线性求解，几何仍是仿射多边形棱柱。直接调用求解器的默认P1/BE保留；正式实验脚本显式选择P2，瞬态使用BDF2。镜像求解重构完整温度节点，能量、Robin和弱残差仍在原空间检查；P2折射OPD在每段使用三点Gauss积分。未知量减少不等于预先证明固定加速倍率。

```bash
python run_all.py --mode full --plan
python run_all.py --mode full
```

正式流程新增 `scripts/convergence/check_transient_fem_convergence.py`：独立运行候选、XY、Z、时间、联合五种配置，保存每例节点和逐步证据，再调用现有保存结果验收入口。当前待验证候选是P2 XY4/双面加密Z分辨率64，参考分别为XY5、Z128和时间步减半；配置与2176步区间计划在脚本明确，不把它们提前称为合格参数。每例先实测短样本并受预算控制，超预算保留结果并退出，不能盲目延长。

已有完整结果可通过 `--fem-reuse-index <JSON>` 显式交给 `run_all.py`；程序核对配置、源码范围与数据SHA-256。已审查的历史P2原型只能走显式来源桥接，保留出生源码与原摘要，不复制大节点、不伪称当前源码重算。普通同源码复用规则不放宽。

`run_transient_case.py --mode full` 现在必须提供 `--fem-convergence <验收目录>`。只有独立重读passed且 `fem_convergence_established=true` 才能采用候选节点；CLI核对物性、角度、功率历史、时间覆盖和适用源码。正式算例温度采用201×201极坐标采样、OPD采用101×101笛卡尔采样，证书另有201×201笛卡尔公共热圆盘；两种采样含义分别记录。smoke和旧保存结果读取保留。

稳态正式FEM改用P2候选 `(4,32)`，XY参考 `(5,32)`、Z参考 `(4,64)`、联合参考 `(5,64)`，温度201/OPD101及载荷10→14阶检查。是否收敛仍由服务器的实际结果与独立保存复核决定，稳态0.1%标准不变。P2保存场积分使用解析形函数权重，不能按P1顶点平均/梯形积分读取。

本次接入和软件测试不等于三角度实验已完成。1ms诊断政策、瞬态FEM网格/联合0.5%、时间0.1%、所有物理门槛均保持；半解析选参没有新增1ms豁免。COMPLETE、文件完整性PASS与科学passed继续分开记录。下列日期较早的P1配置和进度保留为历史证据。

### 当时的证据与待解决问题（2026-10-13，历史记录）

- 最近一次结构重构验收：218项完整快速软件测试通过，38项轻量入口检查通过；34个固定基线数组完全一致。下列118/141项及97数组比较是较早阶段的历史记录。
- 2026-10-13增加批次入口后，227项快速测试通过；真实smoke的测试、FEM独立基准、瞬态小规模FEM和边界复核完成，断点恢复核对哈希后跳过全部完成步骤。全量34步尚未实际执行；其科学验收限制继续保留。
- 本轮已完成瞬态解包装清理、稳态/瞬态科学验证模块整理、统一输出路径和17个实验脚本的分类迁移。所有新实验输出默认位于 outputs/<类型>/<唯一ID>/；论文、Notes 和生成结果不随 Git 提交。
- 全量批次前仍需补齐 FEM 实验验收：瞬态 full 模式目前记录网格/时间步变化，但没有最终收敛接受判定。半解析独立验证现已通过 `--selection` 绑定实际选出的 M/N/J；不带该选项仍是固定参数的诊断算例。不能把 COMPLETE 或保存文件完整性 PASS 当成全部科学检查通过。
- 2026-10-09进一步移除瞬态包装：本次修改前后174个数组完全一致（0/60/-30°、三种历史、J12/J48、0至2500 s的温度/梯度及三种OPD求值路径），最大差0；此前34数组基线也保持差异0。141项软件测试仍全部通过。
- 重构前 118 项测试，重构后 141 项测试通过。固定回归基线的 34 个温度/梯度/系数/OPD 数组完全一致；实际旧/新稳态及瞬态 CLI 的 97 个数组中 86 个完全一致，其余 OPD 最大绝对差 2.12e-22 m。这证明所测计算的重构一致性，不是物理精度认证。
- 既有稳态选择为 M6/N8、投影400×32、光线16点；既有瞬态选择为 M10/N20、投影1200×128、光线32点、J12。瞬态选择仅覆盖名义参数、0/30/60°、阶跃及 0.625/2.5/20/50/75 s，0.1% 为相对有限参考的数值变化目标。
- **额外早期检查发现限制：**60°、M6/N8、J12 在 t=0.001 s 的前表面温度相对 J512 的峰值归一化最大差约 **1.849（184.9%）**，J256/J512 稳定性差约 2.25e-12。报告保留 failed。精确 t=0 零初值和很小的边界残差均不能证明任意接近零的正时刻精度；这项检查不属于此前 t>=0.625 s 的选择范围。本次不修改瞬态公式来消除该失败。
- 本地小 FEM 与独立一维基准已运行；一般 Gaussian 算例最终 FEM 时空收敛、逐点 FEM 边界通量及最终方法间精度接受仍未建立。正式联合选择和全量 FEM 本次没有重新大批量执行。

本地重构基线在 `outputs/refactor_baseline/20261008T134605332619Z`，包含配置、源文件审查、原始数组和 CLI 比较记录，可运行：

```bash
python scripts/tools/compare_refactor_baseline.py --baseline outputs/refactor_baseline/20261008T134605332619Z
```

完整审查、分阶段记录与本地证据目录见 `Notes/仓库重构审查与实施计划.md`；Notes 和生成结果不随 Git 提交。

### 瞬态选参的科学验收（2026-10-13）

`select_transient_truncation.py` 现在同时搜索预设的横向/积分配置与 J=12/48/128/256，正时间检查包含 0.001/0.0025/0.01/0.05/0.1/0.625/2.5/20/50/75 s。两组参考的 J 提高为512/768，参考稳定性仍要求0.02%，候选温度、梯度及OPD误差仍要求0.1%。参考不稳定或没有候选满足门槛时，明确失败，不自动放宽阈值。

新范围内 t<0.625 s 的温度误差使用参考场最大截面范数，避免未加热深度的零分母；旧时刻保留逐截面的误差定义。新增梯度误差为采样处 `(dT/dr,dT/dz)` 差向量的最大范数除以参考梯度向量最大范数；两分量原始数据分别保留。温度和梯度还采样随时间移动的前表面热扩散层。OPD去活塞、参考PV归一化及连续数值PV规则保持不变。

候选还必须通过：t=0零温升/零OPD、全部保留基函数的径向和轴向Robin关系、三类组装边界残差、独立差分PDE（包括早期热扩散层）、体积储能/表面散热能量平衡、固定模态的投影积分加密及长期稳态极限。边界基函数关系不依赖时间，数值残差覆盖零时刻、检查时刻和到长期极限的对数采样；它针对截断热流，不替代原始Gaussian模态截断误差。固定有限J仍不证明任意接近0+的精度。

新选参结果额外保存 `data/acceptance.npz`，`validate_run.py` 重建所选模型、复算物理门槛和证据。旧选择数据仍按原有字段与误差规则读取，不补造新验收证据。精细瞬态算例对全部快照和全部输出历史时刻保存边界检查；历史格式兼容保留。

全量批次顺序改为先选参，再运行带 `--selection` 的半解析独立验证、精细场和FEM对照；FEM侧的半解析J与观察时刻也取自选择结果。修改代码后必须新建批次，旧代码批次不能直接resume；仅可按上文明确限定的版本接续规则继承核验通过的完成阶段。正式扩展范围的联合选择和最终FEM收敛尚需服务器实际运行。
