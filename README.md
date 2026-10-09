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

基础模型需要 NumPy/SciPy；绘图需要 Matplotlib，FEM 另需 scikit-fem。无需 Gmsh、FEniCS 或外部网格文件。只使用模型接口可以安装 `python -m pip install -e .`。Windows 与 Linux 使用同一套入口。

`python run.py` 提供实验选择、脚本默认设置与命令行参数预览、运行确认。取消不会启动计算。直接脚本和交互入口调用同一份实验逻辑。

```bash
python run.py --list
python run.py --experiment steady --show-config
python run.py --experiment steady --yes --angle 60
# 等价的直接入口
python scripts/run_steady_case.py --angle 60
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
| model.py、solutions.py | 四个公共对象、准备流程和统一求值接口 |
| heat_source.py、azimuthal.py、eigenmodes.py、projection.py | 热源、圆周分解、径向根与投影 |
| steady.py、transient.py、modal_evaluation.py | 各自轴向算法、功率历史、共同空间模态求和 |
| opd.py | Snell 折射、透射掩膜、分块光线积分与预备回调 |
| fem.py | 棱柱 FEM 的共享空间组装、稳态/瞬态及分段光线积分 |
| metrics.py、plotting.py | PV、去 piston RMS、绘图和读取已有场重绘 |
| validation.py、transient_validation.py | 三面边界、能量、边界证据保存与模型恢复 |
| convergence.py、verification.py | 场误差、数值极值、独立参考/PDE检查与科学状态 |
| run_validation.py、run_io.py、experiments.py | 重读结果检查、输出元数据、实验入口目录 |

底层 `solve_modal_source`、`TransientModel` 等原有接口仍可用。没有通用 Experiment 继承体系，也没有额外 YAML/TOML 算例配置。

## 物理及数值约定

- 温度场由截断系数直接求值；热源 q 用于投影、加载与验证，不作为温度的插值来源。
- 瞬态保留稳态加载项与轴向扩散模态修正，衰减率为 alpha*(lambda²+beta²)；支持阶跃、线性爬升和指数爬升。
- OPD 始终沿基底内折射光线积分：x=x0、y=y0+z*tan(theta_t)、ds=dz/cos(theta_t)。只保留后表面透射光线，原始 OPD 保存 piston。采用未加热基底中的直线参考光线，未包含热致弯曲。
- PV 为共同有效孔径内的极差；连续数值极值与网格采样 PV 分别记录。RMS 为等间距 Cartesian 有效采样点去 piston 后的 RMS，不用于不等权积分网格。
- 方位角温度份额是各阶 L2 范数平方份额，非热量份额；各温度阶次沿折射光线的 OPD 贡献也不等于入口平面同阶 Fourier 谐波，各阶 PV 不可相加。
- parameters.py 保留 dn_dT=0.96e-5/K 的原文表值，其标注波长为 643.8 nm，而名义光学波长为 1064 nm。密度 2203 kg/m³、比热 772 J/(kg K) 的来源在源码中，比热是 0–100°C 平均值，采用常数近似。改变材料、波长、几何或光斑后需核对物理参数并重新检查数值收敛。

## 实验入口

以下均可通过 `run.py` 的对应名称选择，也可独立运行。指定 `--output-dir` 时目录必须尚未存在。

| 交互名称 | scripts/ 中的入口 | 用途 |
|---|---|---|
| steady | run_steady_case.py | 单个稳态温度/OPD、三面边界及能量 |
| transient | run_transient_case.py | fields 精细场，或 smoke/full 瞬态 FEM 诊断 |
| convergence | select_truncation.py | 稳态联合截断选择 |
| transient-convergence | select_transient_truncation.py | 瞬态联合截断选择 |
| semi-verification | verify_semi_analytical.py | 独立有限差分 PDE、能量、初值附近、投影及已知极限 |
| fem-verification | verify_fem.py | 均匀热流/绝热侧面的一维独立基准 |
| fem-comparison | compare_fem.py | 稳态半解析与 FEM 同条件对照 |
| validate | validate_run.py | 重读保存的数据与报告 |
| boundaries | check_transient_boundaries.py | 给既有瞬态场补三面边界证据 |
| source / azimuthal | inspect_heat_source.py / inspect_azimuthal_modes.py | 热源投影、单圆周方位角分解 |
| temperature / opd | inspect_steady_temperature.py / inspect_opd.py | 单环节检查与图像 |
| sweeps | check_convergence.py | 六种数值设置逐项诊断，不进行联合选择 |
| replot / regression | replot_case.py / compare_refactor_baseline.py | 已有场重绘、重构基线比较 |

典型流程（`<...>` 换成实际目录）：

```bash
python scripts/select_truncation.py
python scripts/validate_run.py --results <稳态选择目录>
python scripts/compare_fem.py --mode smoke --selection <稳态选择目录>

python scripts/select_transient_truncation.py
python scripts/run_transient_case.py --mode fields --angle 60 --selection <瞬态选择目录>
python scripts/run_transient_case.py --mode smoke --angle 60 --history step
python scripts/run_transient_case.py --mode smoke --history linear --ramp-time 5
python scripts/run_transient_case.py --mode smoke --history exponential --ramp-time 5
```

`fields` 读取已完成且参数匹配的选择；默认温度网格 161×180、OPD 网格 201×201，包含逐时刻四联图、稳态极限、演化、早期厚度、历史和模态贡献图。`smoke/full` 的命名不构成精度证明；完整 FEM 批量须另行执行并评估收敛。

## 输出与结果复用

```text
outputs/<experiment_type>/<UTC时间戳_唯一后缀>/
├── run.json             # 目录模式、运行标识与创建时间
├── config.json          # 物理/数值设置、脚本默认值及版本记录
├── summary.json         # 运行状态、关键结果与关联记录
├── metrics.json         # 完成运行的指标
├── verification.json    # 科学检查项目、阈值、范围与状态
├── data/                # NPZ、CSV 等数值证据
└── figures/             # 图像
```

自动目录采用 UTC 时间和唯一后缀，同一秒多次运行也不覆盖。`config.json` 包含可获得的 Git revision、工作区是否有修改、源文件哈希和依赖版本；未提交的修改不能仅靠 commit hash 复现，应一并保留代码及哈希。普通完成运行的科学状态初始为 `not_run`，运行完成本身不表示验证通过。

旧 `results/<run>/` 平铺数据和旧图像位置继续可读，不自动搬动、删除或覆盖已有场。读取旧选择可以生成采用新布局的新实验。

```bash
python scripts/validate_run.py --results <已有结果目录>
python scripts/check_transient_boundaries.py --results <已有瞬态结果目录>
python scripts/replot_case.py --results <稳态或精细瞬态场目录>
```

重绘只读取保存的 fields.npz、history.npz 和摘要，图像进入新的 replot 运行；不重算模态或修改原始场。边界补检查使用保存的根、系数、J 和功率历史；新增边界证据及检查报告。

## 科学验证与软件测试

`tests/` 是快速软件测试。科学验证按明确的物理条件、采样范围和门槛逐项记录 **passed / failed / not_run**。任何失败使汇总状态为 failed；没有失败但仍有未执行项目时汇总为 not_run。不把未执行项目默认通过。

```bash
python scripts/verify_semi_analytical.py --angle 60
python scripts/verify_fem.py
python scripts/validate_run.py --results <验证实验目录>
# 自动化需要严格接受科学状态时使用：
python scripts/validate_run.py --results <结果目录> --require-scientific-pass
```

`validate_run.py` 的 `PASS (saved-file integrity)` 表示已重新读取并核对保存的数组、指标和适用的运行检查；同一行另打印科学汇总状态。科学检查存在失败时，证据文件仍可以完整且一致。严格标志会在科学状态为 failed/not_run 时返回非零退出码。

半解析独立实验包括内部差分 PDE、三面 Robin、独立体积储能/表面散热积分、零热源、初始和长期极限、固定 M/N 的投影积分加密、早期 J 加密以及短时深厚基底解析极限。热边界检查针对截断热源；它不能代替原始 Gaussian 投影误差和截断收敛。

FEM 独立实验采用均匀前面热流、绝热侧面：稳态对照精确线性轴向解，瞬态对照独立一维有限体积且精确时间演化的参考；参考也单独加密。网格和时间步分别改变，检查全部零初值节点、离散方程及逐步能量。该局部基准使用 2% 门槛，一般 Gaussian 问题的最终网格/时间和逐点边界通量仍标为 not_run。方法间接近不能代替两种方法各自验证。

## 当前证据与待解决问题（2026-10-08）

- 重构前 118 项测试，重构后 141 项测试通过。固定回归基线的 34 个温度/梯度/系数/OPD 数组完全一致；实际旧/新稳态及瞬态 CLI 的 97 个数组中 86 个完全一致，其余 OPD 最大绝对差 2.12e-22 m。这证明所测计算的重构一致性，不是物理精度认证。
- 既有稳态选择为 M6/N8、投影400×32、光线16点；既有瞬态选择为 M10/N20、投影1200×128、光线32点、J12。瞬态选择仅覆盖名义参数、0/30/60°、阶跃及 0.625/2.5/20/50/75 s，0.1% 为相对有限参考的数值变化目标。
- **额外早期检查发现限制：**60°、M6/N8、J12 在 t=0.001 s 的前表面温度相对 J512 的峰值归一化最大差约 **1.849（184.9%）**，J256/J512 稳定性差约 2.25e-12。报告保留 failed。精确 t=0 零初值和很小的边界残差均不能证明任意接近零的正时刻精度；这项检查不属于此前 t>=0.625 s 的选择范围。本次不修改瞬态公式来消除该失败。
- 本地小 FEM 与独立一维基准已运行；一般 Gaussian 算例最终 FEM 时空收敛、逐点 FEM 边界通量及最终方法间精度接受仍未建立。正式联合选择和全量 FEM 本次没有重新大批量执行。

本地重构基线在 `outputs/refactor_baseline/20261008T134605332619Z`，包含配置、源文件审查、原始数组和 CLI 比较记录，可运行：

```bash
python scripts/compare_refactor_baseline.py --baseline outputs/refactor_baseline/20261008T134605332619Z
```

完整审查、分阶段记录与本地证据目录见 `Notes/仓库重构审查与实施计划.md`；Notes 和生成结果不随 Git 提交。
