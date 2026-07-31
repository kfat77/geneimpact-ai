# GeneImpact AI

状态：测试通过 | 许可证：MIT | 版本：1.1.0 | Python：>= 3.11

GeneImpact AI 是一套面向动物基因组编辑研究的证据感知型预测与评估工具，覆盖从 sgRNA 设计、编辑效率预测、脱靶分析到可审计报告生成的完整流程。

## 1. 项目简介

### 1.1 背景

CRISPR 基因编辑的预测工具常将模型输出呈现为不透明的"通过/失败"结论，忽视结果的适用范围、不确定性来源与数据血缘。在科研场景中，这种呈现方式难以支撑可复现的研究决策，也不利于伦理与动物福利审查。

### 1.2 解决的问题

- 预测结果缺乏适用性边界与不确定性量化，无法直接支撑研究决策。
- 脱靶检测在大规模参考序列上计算开销高。
- 参考基因组序列依赖研究者手动下载与整理。
- 缺乏将预测输出、生物证据与数据来源绑定为可审计记录的能力。

### 1.3 核心功能

- sgRNA 设计：支持 SpCas9（NGG）、SaCas9（NNGRRT）、Cas12a（TTTV）的 PAM 扫描与候选设计。
- 编辑效率预测：38 特征 RuleSet2-Enhanced 模型（位置权重矩阵 + 热力学稳定性 + 二核苷酸组成 + 物种特异性逻辑校准），置信度 0.55–0.85。
- 脱靶检测：基于 k-mer 种子延伸与哈希索引的加速算法，支持位置加权错配评分与特异性指标。
- 基因组序列获取：通过 Ensembl REST API（NCBI Datasets 兜底）自动下载，本地 SHA-256 缓存。
- 端到端流水线：设计、预测、脱靶检测与报告生成一体化。
- 交互式 Web 应用：基于 Flask 的 REST API 与单页前端。
- 证据评分框架：将预测输出映射至多维有界证据评分。

### 1.4 近期更新（v1.1.0）

| 方向 | v1.0.0 | v1.1.0 | 改进 |
|------|--------|--------|------|
| 效率模型 | 启发式（置信度约 0.35） | 38 特征 RuleSet2-Enhanced，物种校准 | 置信度 0.55–0.85 |
| 脱靶搜索 | 纯 Python 字符串匹配 | k-mer 种子延伸 + 哈希索引 | 100 kb 以上参考序列约 10 倍加速 |
| 基因组下载 | 手动下载 FASTA | Ensembl/NCBI 自动下载 + SHA-256 缓存 | 无需手动步骤 |
| Web 应用 | CLI + 静态 HTML | Flask REST API + 交互式单页应用 | 完整浏览器交互体验 |

在线演示：https://kfat77.github.io/geneimpact-ai/

## 2. 技术栈

- 编程语言：Python 3.11 及以上。
- 构建与打包：hatchling。
- 核心依赖：openpyxl、xlrd（Excel 数据读写）、flask>=3.0（Web 应用 REST API 与单页前端）。
- 标准库模块：urllib（HTTP 下载）、gzip、hashlib、json、dataclasses、pathlib。
- 测试框架：pytest。
- 部署：GitHub Pages 与 Vercel 托管静态演示页面（配置见 vercel.json）。
- 可视化：自包含 HTML 与 SVG 图表，无前端构建步骤。

## 3. 安装与使用说明

### 3.1 环境依赖

- Python 3.11 或更高版本。
- pip 包管理器。

### 3.2 安装

```bash
git clone https://github.com/kfat77/geneimpact-ai.git
cd geneimpact-ai

python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install --upgrade pip
pip install -e .                  # 安装含 Flask 在内的全部运行时依赖
pip install -e ".[dev]"           # 如需运行测试，追加 dev 可选依赖
```

### 3.3 验证安装

```bash
pytest
```

### 3.4 基本用法（命令行）

```bash
# 设计 sgRNA
geneimpact design-sgrna --input target.fasta --species mouse

# 预测编辑效率
geneimpact predict --guide GAGTCTGCTGACAGAGCTC --species mouse

# 检测脱靶位点
geneimpact offtarget --guide GAGTCTGCTGACAGAGCTC --reference ref.fasta

# 运行端到端流水线
geneimpact pipeline --input target.fasta --species mouse --output report.json

# 自动下载参考序列
geneimpact download-genome --species mouse --chrom 1 --start 100000 --end 200000

# 启动 Web 应用
geneimpact webapp --host 0.0.0.0 --port 5000
```

### 3.5 Python API 用法

```python
from geneimpact.advanced_models import score_ruleset2
from geneimpact.genome_downloader import download_sequence
from geneimpact.fast_offtarget import build_seed_index, fast_find_offtargets

# 效率预测
result = score_ruleset2("GAGTCTGCTGACAGAGCTC", species="mouse")
print(result.calibrated_score, result.confidence)

# 下载序列并检索脱靶
seq = download_sequence("mouse", "1", 100000, 200000)
index = build_seed_index({"chr1": seq}, k=12)
offtargets = fast_find_offtargets("GAGTCTGCTGACAGAGCTCG", index, max_mismatches=4)
```

## 4. 项目结构概览

```
geneimpact-ai/
├── src/geneimpact/          # 核心 Python 包
│   ├── cli.py               # 命令行入口与子命令定义
│   ├── pipeline.py          # 端到端分析流水线
│   ├── sgrna_design.py      # sgRNA 设计与 PAM 扫描
│   ├── offtarget.py         # 脱靶检测（暴力比对）
│   ├── fast_offtarget.py    # 脱靶检测（k-mer 加速）
│   ├── efficiency.py        # 效率预测接口
│   ├── advanced_models.py   # 38 特征 RuleSet2-Enhanced 模型
│   ├── genome_downloader.py # Ensembl/NCBI 序列自动下载
│   ├── webapp.py            # Flask Web 应用
│   ├── webapp_static/       # Web 应用前端资源
│   ├── evidence.py          # 多维证据评分框架
│   ├── visualization.py     # HTML/SVG 报告生成
│   ├── crisprscan.py, crispritz.py, behive.py ...  # 预测器与验证适配器
│   └── species.py, calibration.py, readiness.py ... # 物种注册与校准
├── tests/                   # pytest 测试套件
├── demo/                    # 静态演示页面
├── docs/                    # 方法学、模型卡、治理策略文档
├── examples/                # 示例请求与数据
├── pyproject.toml           # 项目元数据与依赖
└── vercel.json              # 静态部署配置
```

## 5. 配置说明

本项目不读取环境变量，也不依赖独立的配置文件。运行参数通过命令行参数传入。

### 5.1 download-genome 命令

| 参数 | 说明 | 默认值 |
|------|------|--------|
| --species | 物种键（mouse、rat、zebrafish、fruit_fly、macaque 等） | 必填 |
| --chrom | 染色体标识 | 必填 |
| --start / --end | 区域起止坐标（1-based） | 必填 |
| --source | 数据源，ensembl 或 ncbi | ensembl |
| --cache-dir | 本地缓存目录 | ./genome_cache |
| --force-refresh | 忽略缓存，强制重新下载 | 关闭 |

### 5.2 webapp 命令

| 参数 | 说明 | 默认值 |
|------|------|--------|
| --host | 监听地址 | 127.0.0.1 |
| --port | 监听端口 | 5000 |
| --debug | 启用调试模式 | 关闭 |

### 5.3 缓存与部署

- 下载的序列按 SHA-256 校验存入 --cache-dir（默认 ./genome_cache），命中缓存时跳过网络请求。
- vercel.json 定义静态部署：构建命令为空，输出目录为 demo，适用于纯静态演示托管。

## 6. 贡献指南

### 6.1 代码规范

- 遵循 PEP 8 风格，使用 4 空格缩进。
- 公共接口提供类型注解。
- 数据载体优先使用 dataclass。
- 新增或修改行为必须附带 pytest 测试。
- 新增预测器或适配器须记录生物域、证据来源、假设与局限。

### 6.2 Pull Request 流程

1. Fork 仓库并创建主题分支。
2. 在分支上实现变更并补充测试。
3. 本地运行 pytest 确保全部通过。
4. 更新相关文档（方法学、模型卡或治理策略）。
5. 提交 Pull Request，附简明的技术与科学依据说明。

### 6.3 提交信息格式

采用约定式提交（Conventional Commits）：

```
<type>: <subject>
```

- type：feat（新功能）、fix（缺陷修复）、docs（文档）、test（测试）、ci（持续集成）、refactor（重构）、chore（杂项）。
- subject：祈使语气，简洁描述，不超过 72 字符，不以句号结尾。
- 示例：feat: add k-mer seed-and-extend off-target search

### 6.4 数据合规

禁止提交原始测序文件、受限研究数据、设施记录、动物个体标识或未经授权可再分发的第三方模型资产。

## 7. 许可证

本项目以 MIT 许可证发布，详见 LICENSE 文件。第三方数据集、出版物、服务与模型集成可能适用各自条款，再分发或投入运行前请查阅 docs/third-party-notices.md。

---

免责声明：GeneImpact AI 为科研决策支持软件，不授权动物基因组编辑、不构成安全性认定、不替代实验验证，亦不替代伦理、生物安全、兽医、动物福利或监管审查。
