# 斜入射 Gaussian 加热下圆柱光学元件的热光响应

[English](README.md)

本仓库计算有限圆柱在斜入射 Gaussian 表面加热和 Robin 换热条件下的稳态与瞬态温升。半解析求解器用 Fourier–Bessel 模态展开温度场。有限元模型（FEM）在棱柱网格上单独组装，直接积分 Gaussian 载荷，用于方法间对照。这里的光程差（OPD）是一阶折射率变化的贡献，沿能够到达后表面的固定直线折射光线积分。

## 环境要求与安装

需要 Python 3.10 或更高版本及 Git。核心求解器使用 NumPy 和 SciPy。下列安装命令还会安装用于绘图的 Matplotlib、用于有限元计算的 scikit-fem 和用于软件测试的 pytest。请在仓库根目录执行。 [依赖配置](pyproject.toml)

**Linux / WSL（Bash）**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,fem]"
python -m pip check
```

**Windows（PowerShell）**

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev,fem]"
.\.venv\Scripts\python.exe -m pip check
```

在 Windows 上，执行下文命令时将 `python` 替换为 `.\.venv\Scripts\python.exe`。

## 仓库结构

| 路径 | 内容 |
| --- | --- |
| [`src/thermal_opd/`](src/thermal_opd/) | 热源投影、模态与有限元求解、OPD 积分和验证 |
| [`scripts/`](scripts/) | 单个算例、收敛研究、独立基准、诊断和重绘 |
| [`run.py`](run.py)、[`run_all.py`](run_all.py) | 交互式实验入口和批量流程 |
| [`tests/`](tests/) | 软件回归测试 |
| [`data/reference/`](data/reference/) | 历史采样输入、部分数值表格、派生光学数组和来源记录 |

## 复现数值实验

下列命令根据物理参数和 Gaussian 热源模型生成新的数值结果，无需准备外部实验数据集。运行结果保存在 `outputs/` 下，该目录不纳入 Git。重放部分历史光学对照属于另一套流程，还需要相应的归档运行文件，详见[参考数据说明](data/reference/README.md)。

### 60° 稳态算例

重复运行时，请选一个尚未使用的输出目录名。

```bash
python scripts/cases/run_steady_case.py --angle 60 --output-dir outputs/steady_case/example-60deg
python scripts/verification/validate_run.py --results outputs/steady_case/example-60deg
```

脚本在运行目录下将温升和 OPD 数组保存到 `data/fields.npz`，将图保存到 `figures/steady_case.png`。`summary.json` 的 `opd_pv_m` 字段记录采样峰谷值，单位为米。复核脚本分别报告保存文件的完整性和科学状态。 [保存结果复核](scripts/verification/validate_run.py)

名义圆柱半径为 12.7 mm，厚度为 5 mm。Gaussian 光束在垂直于传播方向的平面上定义的 1/e² 强度半径为圆柱半径的 0.30 倍。整个平面上的设定吸收功率为 110 mW。折射率和折射率温度导数均由 Corning 7980 的同一 Sellmeier 拟合在 1064 nm、300 K 下求得。算例保留 0、2、4、6 阶方位模态，每阶保留八个径向模态。投影使用 400 个径向采样点和 32 个角向采样点，每条光线使用 16 个积分点。OPD 在入口平面的 201 × 201 网格上求值，孔径半径为圆柱半径的 0.9 倍。 [物理参数](src/thermal_opd/parameters.py) · [稳态算例设置](scripts/cases/run_steady_case.py) · [Sellmeier 系数，表 3](https://ntrs.nasa.gov/api/citations/20070018851/downloads/20070018851.pdf#page=6)

### 稳态数值设置

可用 `python run.py --experiment numerical-settings --yes --angles 0 30 60` 单独复核稳态后处理设置。
程序读取生产算例的默认值，分别改变极值搜索种子、OPD 分批大小、侧壁采样和能量积分点数。
变化量与耗时保存在 `outputs/numerical_settings_audit/<运行ID>/data/comparison.csv`，
验收结果保存在 `verification.json`。
最密种子网格用于数值对照，分批耗时只作诊断。
每个测试的分批大小都须通过逐点 OPD 一致性检查。

### 批量流程

```bash
python run_all.py --mode smoke --angles 60
python run_all.py --mode full --angles 0 30 60 --plan
python run_all.py --mode full --angles 0 30 60
```

smoke 批次运行软件测试、均匀热流下的独立有限元基准，以及带边界检查的小规模瞬态对照。full 批次还包括截断选参、精细场和有限元加密研究。用 `--plan` 可以在启动计算前查看各阶段。全量运行的计算开销较大，前置条件不满足时可能中止。批次目录的 `summary.json` 汇总科学状态。全量批次的退出码 `2` 表示执行结束，但科学证据不完整或未通过。 [批量流程实现](run_all.py)

实验脚本根据输入生成热载荷、网格和采样点。上述流程无需设置随机种子或准备外部输入数据集。单项实验在 `config.json` 中记录设置和来源信息，在 `verification.json` 中记录科学检查。 [运行元数据](src/thermal_opd/run_io.py)

模态热源准备要求有效表面换热系数为正。全绝热情况需要尚未实现的轴对称径向零模态，程序会明确拒绝这类输入。

## 数值参考结果与数据

60° 半解析算例在入口平面上对能到达后表面的光线取样，得到 OPD 峰谷值 **135.07 nm**。稳态算例复现了表中数值，保存文件完整性为 `PASS`，科学汇总状态为 `not_run`。另有有限参考检查记录了 0°、30° 和 60° 下的差异。 [稳态 OPD 表](data/reference/supplement_details/pupil_pv.csv) · [半解析检查表](data/reference/supplement_details/sa_checks.csv)

仓库在 [`data/reference/`](data/reference/) 中保留了部分数值参考表、派生光学数组和来源记录。已记录的对照包括半解析有限参考检查、与有限元的方法间比较、折射光路下的有限元光学加密检查，以及稳态前表面热载荷积分检查。其中加密和积分检查复用了归档温度解，量化的是数值敏感性，不代表实验验证或标准加载器支持的完整科学验收。[有限元光学加密证据](data/reference/supplement_details/fem_optical_refinement_provenance.json) · [派生光学数据](data/reference/optical_update_20261103/)

三个 `data/reference/theta*_history.csv` 文件保留了可选光学重放脚本所用的历史时间网格和温度绘图值，其中旧光学列不代表当前统一材料参数下的结果。详见[输入范围和重放要求](data/reference/README.md)。
