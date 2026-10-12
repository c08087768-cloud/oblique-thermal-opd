# Oblique thermal OPD

研究有限圆柱基底在斜入射 Gaussian 表面加热下的稳态、瞬态温升和热折射 OPD。半解析方法采用 Fourier–Bessel 空间展开，FEM 独立积分原始热源。内部使用 SI 单位，前表面 z=0、后表面 z=L，入射平面为 y-z。温度接口返回相对环境的**温升**。

本仓库服务于同一篇论文及其补充材料。`manuscript/`、`Notes/`、`results/`、`outputs/` 已被 Git 忽略。重构保留现有物理方程、字典模态数据和瞬态轴向算法；一般 Gaussian 算例的最终 FEM 时空收敛仍待完成。

## 安装与开始

在仓库目录、自己的 Python 环境中运行：

```bash
python -m pip install -e ".[dev,fem]"
python -m pytest -q -p no:cacheprovider
python run.py
```

服务器批次入口为 `run_all.py`，默认 smoke，不会默认启动全量：

```bash
python run_all.py --mode full --plan
python run_all.py --mode smoke --angles 60
python run_all.py --mode full
# 中断后用同一代码、同一 mode/angles 恢复：
python run_all.py --mode full --resume outputs/experiment_suite/<批次ID>
```

默认全量计划包含34步：快速测试、FEM独立基准、三个角度的半解析验证、两类截断选择、稳态/精细瞬态场、正式稳态/瞬态FEM、五项诊断、瞬态边界复核和保存结果重绘。串行运行；不需要旧 results、legacy 或开发 baseline。smoke 仅运行测试、FEM独立基准、60°瞬态小规模FEM及其边界复核，不执行正式截断、精细场或大网格。

批次摘要和日志清单位于 `outputs/experiment_suite/<ID>/summary.json`、`data/stages.csv`、`data/logs/`；各实验仍写入自己的 outputs 分类目录。后续步骤只引用本批次生成并复核的选择目录。稳态算例的固定默认设置如果与新选择不同，批次会停止，要求先更新算例设置，不会默默使用不一致的配置。

计算或保存文件复核失败会停止依赖流程；恢复时已完成的数值文件、图像及配置会核对哈希，失败尝试另建目录，保留原文件。科学报告中的 failed/not_run 不改写；独立科学检查失败但文件完整时仍继续收集其它实验。退出码0表示 smoke控制路径完成，1表示执行/数据复核失败，2表示全量执行完成但科学证据失败或不完整，130表示中断。现有瞬态full的收敛验收和早期失败仍可能使全量返回2，不能直接把跑完当成论文验收。

基础模型需要 NumPy/SciPy；绘图需要 Matplotlib，FEM 另需 scikit-fem。无需 Gmsh、FEniCS 或外部网格文件。只使用模型接口可以安装 `python -m pip install -e .`。Windows 与 Linux 使用同一套入口。

`python run.py` 提供实验选择、脚本默认设置与命令行参数预览、运行确认。取消不会启动计算。直接脚本和交互入口调用同一份实验逻辑。

```bash
python run.py --list
python run.py --experiment steady --show-config
python run.py --experiment steady --yes --angle 60
# 等价的直接入口
python scripts/cases/run_steady_case.py --angle 60
```

预览显示脚本默认值及传入的命令行参数；运行后 `config.json` 保存实际设置。M、N、投影网格、采样点、观察时刻和 FEM 网格等仍在对应 Python 脚本顶部定义，支持的 CLI 参数用 `--help` 查看。

## 模型接口与架构

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

旧瞬态求解器类及其转发包装已删除，统一使用 `prepared.solve_transient(...)`。底层 `solve_modal_source` 保留实际源准备职责；共享空间函数从 `modal_evaluation.py` 导入。没有通用 Experiment 继承体系，也没有额外 YAML/TOML 算例配置。

## 物理及数值约定

- 温度场由截断系数直接求值；热源 q 用于投影、加载与验证，不作为温度的插值来源。
- 瞬态保留稳态加载项与轴向扩散模态修正，衰减率为 alpha*(lambda²+beta²)；支持阶跃、线性爬升和指数爬升。
- OPD 始终沿基底内折射光线积分：x=x0、y=y0+z*tan(theta_t)、ds=dz/cos(theta_t)。只保留后表面透射光线，原始 OPD 保存 piston。采用未加热基底中的直线参考光线，未包含热致弯曲。
- PV 为共同有效孔径内的极差；连续数值极值与网格采样 PV 分别记录。RMS 为等间距 Cartesian 有效采样点去 piston 后的 RMS，不用于不等权积分网格。
- 方位角温度份额是各阶 L2 范数平方份额，非热量份额；各温度阶次沿折射光线的 OPD 贡献也不等于入口平面同阶 Fourier 谐波，各阶 PV 不可相加。
- parameters.py 保留 dn_dT=0.96e-5/K 的原文表值，其标注波长为 643.8 nm，而名义光学波长为 1064 nm。密度 2203 kg/m³、比热 772 J/(kg K) 的来源在源码中，比热是 0–100°C 平均值，采用常数近似。改变材料、波长、几何或光斑后需核对物理参数并重新检查数值收敛。

## 实验入口

脚本按 cases、convergence、verification、diagnostics、tools 分类；统一入口的实验名称保持不变。verify_* 是本次指定的脚本名称，核心科学验证模块仍使用 validation 命名。

`run.py` 每次启动一个实验；`run_all.py` 按依赖顺序运行现有全量实验并复核保存数据。批次执行完成与科学验证通过分别记录；它不会补造尚未实现的收敛证据，也不会消除原检查中的失败。

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
| validate | verification/validate_run.py | 重读保存的数据与报告 |
| boundaries | verification/check_transient_boundaries.py | 给既有瞬态场补三面边界证据 |
| source / azimuthal | diagnostics/inspect_heat_source.py / diagnostics/inspect_azimuthal_modes.py | 热源投影、单圆周方位角分解 |
| temperature / opd | diagnostics/inspect_steady_temperature.py / diagnostics/inspect_opd.py | 单环节检查与图像 |
| sweeps | convergence/check_steady_convergence.py | 六种数值设置逐项诊断，不进行联合选择 |
| replot / regression | tools/replot_case.py / tools/compare_refactor_baseline.py | 已有场重绘、重构基线比较 |

典型流程（`<...>` 换成实际目录）：

```bash
python scripts/convergence/select_steady_truncation.py
python scripts/verification/validate_run.py --results <稳态选择目录>
python scripts/verification/compare_steady_fem.py --mode smoke --selection <稳态选择目录>

python scripts/convergence/select_transient_truncation.py
python scripts/cases/run_transient_case.py --mode fields --angle 60 --selection <瞬态选择目录>
python scripts/cases/run_transient_case.py --mode smoke --angle 60 --history step
python scripts/cases/run_transient_case.py --mode smoke --history linear --ramp-time 5
python scripts/cases/run_transient_case.py --mode smoke --history exponential --ramp-time 5
```

`fields` 读取已完成且参数匹配的选择；默认温度网格 161×180、OPD 网格 201×201，包含逐时刻四联图、稳态极限、演化、早期厚度、历史和模态贡献图。`smoke/full` 的命名不构成精度证明；完整 FEM 批量须另行执行并评估收敛。

## 输出与结果复用

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
python scripts/verification/validate_run.py --results <已有结果目录>
python scripts/verification/check_transient_boundaries.py --results <已有瞬态结果目录>
python scripts/tools/replot_case.py --results <稳态或精细瞬态场目录>
```

重绘只读取保存的 fields.npz、history.npz 和摘要，图像进入新的 replot 运行；不重算模态或修改原始场。边界补检查使用保存的根、系数、J 和功率历史；新增边界证据及检查报告。

## 科学验证与软件测试

`tests/` 是快速软件测试。科学验证按明确的物理条件、采样范围和门槛逐项记录 **passed / failed / not_run**。任何失败使汇总状态为 failed；没有失败但仍有未执行项目时汇总为 not_run。不把未执行项目默认通过。

```bash
python scripts/verification/verify_semi_analytical.py --angle 60
python scripts/verification/verify_fem.py
python scripts/verification/validate_run.py --results <验证实验目录>
# 自动化需要严格接受科学状态时使用：
python scripts/verification/validate_run.py --results <结果目录> --require-scientific-pass
```

`validate_run.py` 的 `PASS (saved-file integrity)` 表示已重新读取并核对保存的数组、指标和适用的运行检查；同一行另打印科学汇总状态。科学检查存在失败时，证据文件仍可以完整且一致。严格标志会在科学状态为 failed/not_run 时返回非零退出码。

统一入口按保存的 `kind` 调用稳态或瞬态验证模块。`scientific_verification` 是可包含两类物理检查的通用证据格式，其标量数组、阈值和状态由 `validation_common.py` 直接复核；旧报告不含模型子类型也可读取。新代码统一采用 validation 命名；`verification.json`、`verification_fields.npz` 和既有 JSON/NPZ 字段继续沿用原数据格式。状态记录与汇总只整理已执行检查，不由文件 I/O 判断物理结果。稳态/瞬态共享空间差分及网格几何检查，FEM 求解与载荷积分仍保持独立。

半解析独立实验包括内部差分 PDE、三面 Robin、独立体积储能/表面散热积分、零热源、初始和长期极限、固定 M/N 的投影积分加密、早期 J 加密以及短时深厚基底解析极限。热边界检查针对截断热源；它不能代替原始 Gaussian 投影误差和截断收敛。

FEM 独立实验采用均匀前面热流、绝热侧面：稳态对照精确线性轴向解，瞬态对照独立一维有限体积且精确时间演化的参考；参考也单独加密。网格和时间步分别改变，检查全部零初值节点、离散方程及逐步能量。该局部基准使用 2% 门槛，一般 Gaussian 问题的最终网格/时间和逐点边界通量仍标为 not_run。方法间接近不能代替两种方法各自验证。

## 当前证据与待解决问题（2026-10-13）

- 最近一次结构重构验收：218项完整快速软件测试通过，38项轻量入口检查通过；34个固定基线数组完全一致。下列118/141项及97数组比较是较早阶段的历史记录。
- 2026-10-13增加批次入口后，227项快速测试通过；真实smoke的测试、FEM独立基准、瞬态小规模FEM和边界复核完成，断点恢复核对哈希后跳过全部完成步骤。全量34步尚未实际执行；其科学验收限制继续保留。
- 本轮已完成瞬态解包装清理、稳态/瞬态科学验证模块整理、统一输出路径和17个实验脚本的分类迁移。所有新实验输出默认位于 outputs/<类型>/<唯一ID>/；论文、Notes 和生成结果不随 Git 提交。
- 全量批次前仍需补齐实验验收：瞬态 full 模式目前记录网格/时间步变化，但没有最终收敛接受判定；半解析独立验证固定使用 M6/N8、J12，尚未绑定精细瞬态选出的 M10/N20 配置。不能把 COMPLETE 或保存文件完整性 PASS 当成全部科学检查通过。
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
