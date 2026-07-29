# 猕猴基因编辑预测：一手证据、可复现接口与能力边界

> 核验日期：2026-07-29
>
> 物种：恒河猴 *Macaca mulatta*；食蟹猴 *Macaca fascicularis*
>
> 编辑系统：CRISPR/Cas9 核酸酶、碱基编辑、prime editing
>
> 目标：判断哪些公开的一手数据足以实现可复现的胚胎/体内编辑结果、修复谱、脱靶或表型风险外部验证接口，而不是把小样本动物实验包装成通用预测器。

## 结论先行

1. **目前不能建立或宣称“适用于猴的通用基因编辑效果/安全概率预测器”。** 两个物种均缺少同一实验协议下、覆盖数百至数千条 guide、多个基因和足够动物/胚胎的直接标签。现有研究的有效样本单位常只有 1–11 个靶点、2–5 只动物或少量胚胎，且递送、胚胎阶段、测序和群体来源高度异质。
2. **食蟹猴最强可实现切口是胚胎碱基编辑的“外部迁移基准”，不是物种内训练集。** Zhang 等 2020 提供 11 个 on-target loci、30 个目标碱基/上下文记录、逐胚胎 clone counts，涵盖 BE3、BE4-Gam、ABE7.10、SaKKH-BE3 和 multiplex；文章、补充与原始深测序/WGS 均公开。它足以检验一个在人/鼠等数据上训练的模型，能否在锁定协议和 legacy assembly 上保持排序、编辑窗与旁观者/indel 方向的一致性，但 11 loci 远不足以拟合或校准食蟹猴概率。[论文与补充](https://pmc.ncbi.nlm.nih.gov/articles/PMC7214463/)、[SRA 深测序 PRJNA561611](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA561611)、[WGS PRJNA505503](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA505503)。
3. **恒河猴最强碱基编辑风险证据不是跨靶点活性数据，而是 ABEmax 的匹配 SCNT 胚胎 DNA/RNA 安全实验。** 11 个 ABEmax SCNT 样本、11 个未编辑 SCNT 对照和 4 个供体样本做了 >60× WGS，另有 5+5 囊胚 RNA-seq；研究未发现明显 DNA SNV 增量，却发现大规模 A-to-G RNA 脱靶。它能作为“编辑器本体的转录组危害”外部基准，但靶点是外源 GFP，不能训练 endogenous guide activity。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9307242/)、[GSA CRA004092](https://ngdc.cncb.ac.cn/gsa/browse/CRA004092)、[分析代码 Zenodo 6605030](https://doi.org/10.5281/zenodo.6605030)。
4. **体内/长期效果只能做情境化验证。** 食蟹猴 PCSK9 肝脏 LNP-ABE 有两项独立研究，能检验剂量—肝编辑—蛋白—LDL 的链路；恒河猴 HBB ABE8e 编辑 HSC 自体移植能检验 200 天的编辑持久性、血系选择和候选脱靶。它们都只有一个生物靶点，不得转化为跨基因效应概率。
5. **Cas9 大缺失与嵌合风险必须作为“案例型 hazard benchmark”单独输出。** Mauritian 食蟹猴 CCR5 blastomere WGS、恒河猴 MCPH1 trio-WGS、PINK1 大缺失和 MYO7A 胚胎/活产均显示，仅用短 PCR 或少量候选位点会漏掉重要事件；但阳性案例太少，不能校准绝对发生率。
6. **截至核验日，没有找到同行评议且公开原始标签、可用于恒河猴或食蟹猴胚胎/体内 prime-editing 训练或校准的数据集。** 2023 年 Prime Medicine 会议稿报告健康食蟹猴肝脏 SLC37A4 surrogate 编辑，但未披露可复现实验级原始数据、样本数、assembly 或开放许可；它只能登记为 `corporate_nonreproducible_evidence`，不能进入模型或 benchmark。[原始会议稿](https://primemedicine.com/wp-content/uploads/2023/10/ESGCT-2023-Prime-editing-precisely-corrects-prevalent-pathogenic-mutations-causing-Glycogen-Storage-Disease-Type-1b-GSD1b.pdf)。

建议的平台能力状态：

| 物种 | 通用预测 | 当前可交付的最强状态 |
|---|---|---|
| *M. fascicularis* | `insufficient_public_data` | `bounded_external_benchmark`：胚胎碱基编辑迁移；PCSK9 肝脏 PK/PD；CCR5 大缺失/脱靶案例 |
| *M. mulatta* | `insufficient_public_data` | `bounded_external_benchmark`：ABE RNA hazard；Cas9 胚胎嵌合/大缺失；HSC 长期编辑选择 |
| 两种猕猴 prime editing | `no_reproducible_public_labels_found` | 只允许证据登记，不生成效能或安全预测 |

## 物种、群体与参考基因组不可合并

### 物种与群体

- *M. mulatta* 与 *M. fascicularis* 必须是两个顶层物种键，不能用 `macaque` 直接汇总。
- Mauritian cynomolgus macaque（MCM）是明确的地理/遗传群体；CCR5 blastomere 数据只能标为 `Macaca_fascicularis / Mauritian`。
- 大量中国或北美灵长类中心研究没有披露 Indian/Chinese/Vietnamese 等祖源。未披露时必须保存 `population="not_reported"`，不得根据机构反推。
- 多篇实验用了同一动物、同一胚胎批次或后续表型研究中的同一 founder。划分训练/验证集时必须按 `animal_id/family_id/embryo_batch/publication_lineage` 分组，不能按测序文件随机拆分。

### Assembly 规则

历史编辑研究主要使用：

- 食蟹猴：`Macaca_fascicularis_5.0 / GCF_000364345.1`（常写作 `macFas5` 或 `M_fascicularis_5.0`）。
- 恒河猴：`Mmul_8.0.1 / rheMac8`，较新研究也使用 `Mmul_10 / rheMac10`。

截至核验日，NCBI 当前参考分别为：

- 恒河猴 `T2T-MMU8v2.0 / GCF_049350105.2`。[NCBI 官方注释报告](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Macaca_mulatta/GCF_049350105.2-RS_2025_08/)
- 食蟹猴 `T2T-MFA8v1.1 / GCF_037993035.2`。[NCBI 官方注释报告](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Macaca_fascicularis/GCF_037993035.2-RS_2025_03/)

因此接口必须同时保存：

```text
species
population
source_assembly_accession
source_contig
source_coordinate
source_ref_alt
target_assembly_accession
liftover_status
target_sequence_verified
```

不得静默 liftOver。只有目标序列、PAM、编辑窗碱基和等位基因在新 assembly 上逐项复核后，才能给出新坐标；否则输出 `unresolved_legacy_coordinate`。

## 食蟹猴一手证据

### A. 最强：胚胎碱基编辑序列级迁移基准

**Zhang et al., Nature Communications 2020, “Multiplex precise base editing in cynomolgus monkeys”**

- 物种/群体：*M. fascicularis*；1 只 11 岁雄性供精，4 只 5–9 岁雌性供卵；祖源未报告。
- 阶段与递送：ICSI 后 10–12 h 的合子胞质注射；base-editor mRNA 100 ng/µL、每条 sgRNA 50 ng/µL；培养 ≥3 天或移植 8-cell 至 blastocyst。
- 编辑器/标签：
  - BE3，FAH exon 4：16 胚胎，11 编辑，10 个出现期望 C>T，3 个有 indel。
  - ABE7.10，APP：9 胚胎，6 编辑，均为 A>G，未检出 indel。
  - multiplex BE3/ABE/SaKKH-BE3：FAH 多 exon、HBB、TP53、EMX1、FANCF、BRCA1；逐胚胎 clone counts、目标碱基比例、非目标转换和 indel。
  - 合计可整理为 11 loci、30 个目标碱基/sequence-context records；不能把同一 sgRNA 编辑窗内多个碱基当作独立 guide。
- 体内/安全：67 个带原核胚胎注射 BE3+FAH sgRNA，56 个发育并移植至 5 个代孕，3 个妊娠/胎儿；一个 FAH W78Y 胎儿与双亲做 30–40× WGS，得到 66 个 de novo SNV、5 个 indel，无法区分自然 de novo 与 BE3 非向导依赖事件。
- Assembly：`GCF_000364345.1 / Macaca_fascicularis_5.0`。
- 文件：补充表含 guide/候选脱靶/逐胚胎 clone 结果；SRA `PRJNA561611` 为 deep sequencing，`PRJNA505503` 为 WGS；正文与 Source Data 可重建图表。
- 许可：文章和补充为 CC BY 4.0；SRA 为公开归档，但仓库可访问不等同于额外的数据库版权许可，发布衍生数据时仍应保留 provenance。
- 可校准性：**不可做食蟹猴概率校准**。可做冻结的外部迁移测试：target-base rank、editing-window profile、bystander composition、indel presence 和 embryo-level uncertainty。
- 训练重叠：同一胚胎内 clone、同一 guide 的多个碱基、同一 multiplex 注射批次都必须同组；FAH fetus WGS 与胚胎 FAH-E4 不是独立靶点。

来源：[论文/补充/Associated Data](https://pmc.ncbi.nlm.nih.gov/articles/PMC7214463/)、[PRJNA561611](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA561611)、[PRJNA505503](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA505503)。

### B. 最强体内碱基编辑情境验证：PCSK9 肝脏

**Rothgangl et al., Nature Biotechnology 2021**

- 物种/群体：*M. fascicularis*，约 2 岁雄性；祖源未报告。
- 系统：LNP 共递送 1-methoxyuridine-modified ABEmax mRNA 与化学修饰 sgRNA_hP01，靶向 PCSK9 intron 1 splice donor。
- 设计：4 组、每组 3 只；0.75 或 1.5 mg/kg，总 RNA 单次或间隔两周重复静脉输注；Day 29 终点。
- 标签：每只 6 个肝活检；高剂量单次约 27.6% A>G、重复约 24.1%，indel 最高 0.27%；同时有血清 PCSK9、LDL、ALT/AST、细胞因子、组织分布和抗 Cas9/TadA 抗体。
- 脱靶：CIRCLE-seq、CHANGE-seq、iGUIDE，另对 macaque genomic DNA 做 CHANGE-seq；深测 top candidates。结论只适用于该 guide、候选发现流程与检测限。
- 文件：GEO `GSE168365` 含 targeted amplicon/RNA-seq，ENA/SRA `PRJEB41832` 含 WGS；文章补充含 guide、候选位点和组织结果。
- 许可：文章 CC BY 4.0。
- 边界：可用于 `delivery × dose × organ × time` 外部验证；只有一个 locus，不能检验跨序列活性。

来源：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8352781/)、[GSE168365](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168365)、[PRJEB41832](https://www.ebi.ac.uk/ena/browser/view/PRJEB41832)。

**Musunuru et al., Nature 2021**

- 物种：*M. fascicularis*；群体未报告。
- 系统：LNP 递送 ABE8.8 mRNA 与 PCSK9-1 gRNA，靶向同一类 exon 1 splice-donor 机制。
- 设计：短期有 1.0 mg/kg 的 3 只处理动物；长期 3.0 mg/kg 为 4 只处理、2 只 PBS，对 PCSK9/LDL 随访 238 天；另有 non-targeting gRNA 和 PBS 组。
- 标签：肝编辑、血清 PCSK9/LDL、肝功能、药代、组织分布、Digenome-seq/候选脱靶和源数据 XLSX；报告约 90% PCSK9 和约 60% LDL 降低并至少维持 8 个月。
- 原始数据：SRA `PRJNA716270`；论文明确使用 `macFas5 / Macaca_fascicularis_5.0`。
- 许可：论文是订阅内容，Source Data/补充和 SRA 可公开访问，但未授予 CC 文章再利用许可；不要将可下载等同于可任意再发布。
- 训练重叠：与 Rothgangl 研究是不同 editor/guide/LNP 组合，但都是 PCSK9 splice disruption，不能计作两个独立 sequence loci。它们适合做“同一生物靶点、不同实现”的机制一致性检查。

来源：[论文与 Data Availability/Source Data](https://www.nature.com/articles/s41586-021-03534-y)、[PRJNA716270](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA716270)。

### C. Cas9 胚胎大缺失/脱靶案例：Mauritian CCR5

**Schmidt et al., Frontiers in Genome Editing 2023**

- 群体：Mauritian cynomolgus macaque；双亲为 12 岁雌性和 6 岁雄性。
- 系统：one-cell IVF embryo 注射 Cas9 RNP + 2 条 CCR5 exon 2 sgRNA。
- 样本：2 个停滞胚胎（6-cell 与 9-cell）；WGS 覆盖 embryo 4 的 6 个 blastomere、embryo 5 的 8/9 个 blastomere，加双亲约 30× WGS。单细胞 WGA 后覆盖高度不均（0.81–77.77×）。
- Assembly：`M_fascicularis_5.0`。
- 标签：短/长 PCR、on-target indel/大缺失、93 个 2–3 mismatch 候选位点、de novo SNV/SV；SFMBT2 和 LIPC 的候选脱靶由 Sanger 进一步检查。
- 关键限制：没有未操作 blastomere 对照；WGA bias、allele dropout、短读长 SV calling 和不完整旧 assembly 都会制造假阳性/假阴性。论文自己也不能确认所有 SV 均由 Cas9 引起。
- 文件/许可：WGS `PRJNA880597`；补充表；CC BY 4.0。
- 可实现性：做 `hazard_case` 和 pipeline regression fixture；**不能**用 2 个胚胎估计大缺失或脱靶概率。

来源：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9877282/)、[PRJNA880597](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA880597)。

### D. 表型风险只能作为单基因 evidence card

**STXBP1 R292H base-edited monkeys**

- A3A (hA3A-BE3-Y130F) mRNA 100 ng/µL + sgRNA 50 ng/µL 注射 one-cell embryos；13 个 blastocyst 中 12 个有期望单碱基编辑。
- 有 2 个活产编辑猴、3 个 WT 对照用于 on-target/indel；主要表型分析最终是极小动物数，并有多组织 mosaicism、EEG、行为和单细胞 RNA-seq。
- WGS 使用 `Macaca_fascicularis_5.0`，Cas-OFFinder/Cas-OT 为候选位点流程；文章没有给出可核验的公共 SRA/GEO accession。
- 许可为 CC BY-NC-ND 4.0。它只能证明特定 STXBP1 R292H/L291L 编辑组合与严重神经表型之间的实验关联，不能变成通用“表型安全分数”。

来源：[论文与补充](https://pmc.ncbi.nlm.nih.gov/articles/PMC9171284/)。

## 恒河猴一手证据

### A. 最强编辑器本体风险基准：ABEmax SCNT 胚胎

**Kang et al., Science Advances 2022**

- 物种/来源：*M. mulatta*；12 岁 GFP transgenic monkey 的成纤维细胞作 SCNT 核供体。所有实验胚胎共享该供体遗传背景，这是优点，也是禁止随机拆分的原因。
- 系统：SCNT activation 后 6–8 h 注射 ABEmax mRNA 100 ng/µL + 1 条 GFP sgRNA 50 ng/µL。
- 规模：87 个 SCNT、207 个 SCNT-ABE；68 个对照胚胎移植至 18 个 surrogate，154 个编辑胚胎移植至 31 个 surrogate；一个编辑活产，12 h 后死亡，作者将其归因为 SCNT 发育/胎盘问题而非 ABE，但该因果不能由单例证实。
- 直接标签：40 个 SCNT-ABE 胚胎中的 on-target/旁观者 clone 结果；活产 5 个组织均为目标 A4/A9 100%。
- WGS：11 SCNT、11 SCNT-ABE、4 donor-cell samples，平均 >60×；9+8 blastocysts、一个 aborted control 的 2 个组织和 edited newborn 的 3 个组织；`Mmul_8.0.1`。
- RNA-seq：5 SCNT + 5 SCNT-ABE blastocysts。ABE 组 A>G RNA edits 平均约 7,888，对照约 2,297；>75% 的 exonic edits 为 nonsynonymous。DNA SNV/indel 总数无显著增加，但 large insertion 信号与检测/统计能力仍有限。
- 文件：GSA `CRA004092`；Zenodo `6605030` 与 GitHub `daishaoxing/OA-SCNT` 提供代码。
- 许可：文章为 CC BY-NC；Zenodo/GitHub 代码必须读取各自 license 文件后再决定再分发，不能由论文许可自动推断。
- 可实现性：`editor_global_rna_hazard` 外部基准；不能训练 endogenous target efficiency，也不能把“未见明显 DNA 增量”写成“无 DNA 脱靶”。

来源：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9307242/)、[CRA004092](https://ngdc.cncb.ac.cn/gsa/browse/CRA004092)、[代码快照](https://doi.org/10.5281/zenodo.6605030)。

### B. 最强 Cas9 trio-WGS 阴性基准：MCPH1

**Luo et al., Nature Communications 2019**

- 系统：Cas9 mRNA 20 ng/µL + 两条 sgRNA 各 10 ng/µL，分别靶向 MCPH1 exon 2 与 exon 4，注射 IVF zygote。
- 胚胎/动物：15 个测试胚胎中 13 个 KO-positive；另注射 30 个 zygote，24 个发育并移植至 6 个 surrogate，产生 5 个处理后代，其中 4 个检出编辑。
- WGS：5 个处理后代 + 3 个 WT parents，约 46×；死亡个体有 brain/liver/muscle 多组织；靶区另有 PacBio。SpeedSeq 的 4,807 个候选和 Cas-OFFinder ≤7 mismatch 候选未检出与 guide 相关突变。
- Assembly：`Mmul_8.0.1 / rheMac8`；文中对一个既有食蟹猴 OCT4 trio 进行了再分析，后者不是新的独立动物。
- 文件/许可：`PRJNA588331`，Source Data 与 Supplementary Data；CC BY 4.0。
- 可实现性：高深度 family-aware WGS regression benchmark。几乎全是候选脱靶阴性，且只有 2 guides，不能拟合 specificity probability；阴性结论受短读长、mosaic VAF 和候选定义限制。

来源：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892871/)、[PRJNA588331](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA588331)。

### C. 胚胎递送、嵌合与 TE biopsy 不可靠性：MYO7A

**Ryu et al., Scientific Reports 2022**

- 3 条 MYO7A exon 3 guides 混合，比较 Cas9 mRNA + hybrid gRNA 与 Cas9 protein + synthetic sgRNA。
- mRNA/hyb-gRNA：125 oocytes、53 fertilized、36 analyzable、14 edited（38.9%）。
- RNP：409 oocytes、234 fertilized、211 analyzable、161 edited（76.3%），但 embryo arrest 更高。
- 6 个 blastocysts、5 次 transfer 产生 1 个活产；TE biopsy 预测约 92.1% mutant，而活产多个组织只有约 40–50%；单细胞 13/22 为 homozygous insertion、9/22 WT。
- 只检测 9 个候选脱靶；没有 SRA/GEO，参考 assembly 未明确报告。文章/补充为 CC BY 4.0。
- 可实现性：可测试“递送形式—编辑率—发育停滞”的方向和 `TE biopsy != whole embryo` 警告；三条 guide 作为混合物递送，不能拆出单 guide efficiency。

来源：[论文与补充](https://pmc.ncbi.nlm.nih.gov/articles/PMC9203743/)。

### D. 多 guide 大缺失/表型 hazard：PINK1

**Yang et al., Cell Research 2019**

- one-cell rhesus embryos，Cas9 + PINK1 exon 2/exon 4 两条 guides；胚胎编辑率约 61.5%。
- 87 个 embryo 移植至 28 个 surrogate，11 个 live births，8 个 mutants。
- 两靶点之间出现约 7,237 bp 大缺失，并观察严重神经退行性表型；3 个 mutant 做 WGS、分析 2,189 个 ≤5 mismatch candidate sites。
- 文章与补充为 CC BY 4.0；未给可核验的公共原始 WGS accession，群体与 assembly 披露不足。
- 可实现性：大缺失和复杂表型的 `hazard_case`；不能估计发生率或将表型归因泛化到其他基因。

来源：[论文与补充](https://pmc.ncbi.nlm.nih.gov/articles/PMC6461954/)。

### E. 最强长期碱基编辑选择/候选脱靶基准：HBB HSC 自体移植

**Radtke et al., Science Translational Medicine 2025**

- 3 只健康 juvenile rhesus macaques；群体/祖源未报告。
- CD34+CD90+ HSC 离体电转 ABE8e-NRCH mRNA + 为恒河猴修订的 HBB guide，TBI 后自体回输；追踪至约 200 天。
- 由于猴不携带 sickle allele，以 A9 synonymous bystander 作为 surrogate：infusion product >60%，长期外周血约 20–30%，各血系和骨髓约 25–35%。
- 罕见标签：附近 C substitution <0.1%，单个 indel <0.15%、总 indel <1%。
- 脱靶：从旧 `rheMac3` candidates liftOver 至 `Mmul_10/rheMac10`；4,324 个 ≤4 mismatch candidates 中选择 1,000 个设计 rhAMPseq，实际覆盖 951，发现 8 个显著 A>G sites，并随时间持续。
- 文件：SRA `PRJNA1036686`；代码 `KiemLab-RIS/Makassar_project`，Zenodo `15190299`。
- 许可：文章 CC BY 4.0；代码按仓库/Zenodo 独立许可处理。
- 可实现性：`longitudinal_edit_persistence`、`lineage_selection` 和候选脱靶外部验证。951 个位点是相似性预选，不是全基因组无偏 negatives；仅 1 guide、3 animals，不得做跨 guide 校准。

来源：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12490786/)、[PRJNA1036686](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1036686)、[代码](https://github.com/KiemLab-RIS/Makassar_project)、[Zenodo](https://doi.org/10.5281/zenodo.15190299)。

### F. 其他可登记但不能升级能力的研究

| 研究 | 可用信息 | 为什么不是训练集 |
|---|---|---|
| HTT CAG expansion HDR embryos (2024) | 438 oocytes→274 zygotes→105 analyzable embryos；16 HDR、23 NHEJ；2 guides + 76-CAG ssDNA；文章/补充 CC BY | 单一 locus，原始 reads 未归档，guide 2 还有候选脱靶阳性；只适合作为 HDR/NHEJ case benchmark。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11119628/) |
| DMD rhesus WGS (2018) | 2 个 4 岁 edited monkeys + 1 WT sibling，约 60× WGS，17 个候选/847 homologous sites | 无公共 raw WGS；非 trio；样本极小；文章非 CC。只登记阴性安全方法证据。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6066302/) |
| Rhesus airway ABE RNP (2023) | ABE8e-Cas9 RNP + CCR5 guide 气管内递送；少量 juvenile monkeys，气道上皮最高约 5.3%；`PRJNA1043615` | 1 guide、每条件约 1–2 animals，只能验证组织递送。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10698009/)、[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1043615) |
| Fetal PCSK9 HMEJ/PET (2026) | 4 个 rhesus fetuses，第二孕期肝内双 AAV9，SaCas9 + PCSK9 guide/HMEJ donor；有 indel、正确与异常 junction、PET 时序 | 1 guide、无公共原始测序、非开放文章许可；只能做胎儿递送和 junction hazard case。[论文](https://journals.sagepub.com/doi/10.1177/10430342261438941) |
| Multiplex CD33+HBG ABE HSC transplant (2025) | 2 只雄性 rhesus，2 guides，随访约 500–600 天；SRA `PRJNA1222463` 含大量 longitudinal samples | animal n=2；部分 off-target 分析来自 human CD34，不能自动归入猴标签；文章 CC BY-NC-ND。[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1222463)、[论文 DOI](https://doi.org/10.1038/s41467-025-59713-2) |

## Prime editing：明确的缺口

Prime Medicine 2023 ESGCT 幻灯片称：

- 健康成年食蟹猴；
- universal liver LNP，单次静脉输注；
- prime-editor mRNA + surrogate pegRNA + nicking gRNA；
- 在 SLC37A4 p.L348 codon 安装单碱基替换；
- whole-liver precise editing 最高约 50%，按“肝细胞占肝脏 60%”推算为最高约 83% hepatocytes；
- “unintended edits”只定义为目标上下游 300 bp 内 SNV/indel；
- 有 naive 与此前接受同一 LNP、另一未知靶点的 non-naive 动物。

但材料未给出：

- 各组精确 animal n、性别、年龄、群体；
- pegRNA/ngRNA 序列、剂量、assembly；
- per-animal measurements、FASTQ/BAM、BioProject/GEO/GSA；
- genome-wide discovery 的可执行代码和阈值；
- 可再利用许可或同行评议方法全文。

因此最高状态是：

```yaml
technology: prime_editing
species: Macaca_fascicularis
status: corporate_nonreproducible_evidence
training_eligible: false
external_benchmark_eligible: false
claim_allowed: "public presentation reports feasibility"
claim_forbidden: "validated efficacy, calibrated safety, or species predictor"
```

来源：[Prime Medicine ESGCT 2023 原始幻灯片](https://primemedicine.com/wp-content/uploads/2023/10/ESGCT-2023-Prime-editing-precisely-corrects-prevalent-pathogenic-mutations-causing-Glycogen-Storage-Disease-Type-1b-GSD1b.pdf)。

## 推荐实现的 5 个受限接口

### 1. `cyno_embryo_base_edit_transfer_v1`

用途：冻结 Zhang 2020 标签，评估非猕猴训练模型在食蟹猴胚胎上的外部迁移。

最小输入：

```text
assembly_accession, target_sequence, pam, editor, guide_sequence,
target_base_position, embryo_batch, injection_protocol
```

最小输出：

```text
observed_target_fraction_by_embryo
observed_bystander_spectrum_by_embryo
observed_indel_presence
model_rank_metrics
calibration_status = "not_calibratable_11_loci"
```

拆分：按 guide/gene/injection batch 分组；clone 不能跨折，编辑窗内多个碱基不能跨折。

### 2. `macaque_editor_global_hazard_v1`

用途：以恒河猴 SCNT ABEmax 5+5 RNA-seq 和 11+11 WGS 检查 RNA A>G burden、motif/genic distribution 与 DNA event burden。

必须同时输出：

```text
editor = ABEmax
target = exogenous_GFP
developmental_context = SCNT_blastocyst
matched_donor_background = true
generalization = "editor/context only; not guide sequence"
```

### 3. `primate_nuclease_structural_hazard_cases_v1`

用途：用 CCR5 MCM blastomeres、MCPH1 trio、PINK1 和 MYO7A 作为大缺失、嵌合、TE-biopsy mismatch 与 candidate-off-target pipeline 的回归案例。

输出不是概率，而是：

```text
hazard_detected
detection_assay
assay_blind_spots
evidence_strength
causality_status
```

### 4. `primate_in_vivo_context_validation_v1`

用途：PCSK9 cynomolgus liver、rhesus airway CCR5 和 fetal PCSK9 仅按原研究的 route/dose/organ/time 验证。

禁止 sequence-level pooling。不同递送平台、组织和发育阶段必须各自成为 context key。

### 5. `rhesus_longitudinal_selection_v1`

用途：HBB ABE8e-NRCH HSC 自体移植，评估输入产品到外周血/骨髓/血系的编辑持久性、旁观者和候选脱靶选择。

可报告纵向轨迹和检测限；不能从 3 只动物推断其他 guide 的克隆优势或肿瘤风险。

## 统一标签与不确定性规范

每条记录至少包含：

```text
species
population
assembly_accession
gene
locus_id
guide_id
editor_or_nuclease
delivery
developmental_stage
animal_id
family_id
embryo_id
embryo_batch
tissue
timepoint
assay
coverage_or_read_count
detection_limit
raw_accession
article_license
data_license
code_license
outcome_type
observed_value
causality_status
training_eligibility
benchmark_eligibility
```

强制不确定性：

- `not_detected` 不能改写为 `absent`。
- 候选位点阴性不能改写为 genome-wide negative。
- WGS 阴性必须绑定 depth、mosaic VAF 检出能力、SV caller 和 assembly。
- embryo clone counts 用 beta-binomial/Dirichlet-multinomial 等层级观测模型时，胚胎才是生物重复，clone/read 不是。
- 只有 1 locus 的 dose-response 可拟合该实验曲线，但不得生成跨 locus calibration curve。
- 表型必须标记 founder mosaicism、旁观者、large deletion、背景和是否有同窝/家系对照。

## 许可与再利用规则

1. 文章、补充、原始测序和代码是四种不同资产，许可必须分别记录。
2. CC BY 文章可在署名条件下重用文章/所含补充；CC BY-NC、CC BY-NC-ND 不适合默认商业再分发或改作训练包。
3. NCBI SRA/ENA/GSA “公开可下载”不自动等同于 CC0；保留 accession、submitter、论文和使用条款。
4. GitHub 仓库没有明确 LICENSE 时，默认只能读和复现实验，不能假定可复制进产品。
5. Nature/SAGE 等订阅文章的公开 Source Data 或作者稿，不把整篇文章转成开放许可。

## 最终能力声明

当前可以诚实交付：

> 系统提供针对恒河猴与食蟹猴的 assembly-aware guide 设计、既有跨物种模型的受限外部迁移评估、公开一手实验的证据卡，以及针对特定编辑器/递送/组织/发育阶段的 hazard 与 longitudinal benchmark。结果显示证据覆盖、检测限和不可外推项，不替代动物实验与独立测序验证。

当前不能声称：

- “预测任意猴基因编辑成功率”；
- “预测多基因协作后的猴整体表型”；
- “证明某条 guide 无脱靶或安全”；
- “预测胚胎发育、活产或疾病风险的绝对概率”；
- “已支持猴 prime editing 效果预测”；
- “把食蟹猴结果直接迁移为恒河猴结果”。

要把状态升级为可校准 predictor，最低需要一个前瞻性公共项目：每个物种至少数百条 guide、多基因、固定 editor/delivery/stage，逐胚胎或逐动物的 amplicon repair spectrum、long-read/SV、unbiased off-target、RNA off-target、发育/组织/表型标签，公开原始 reads、样本表、代码、明确许可，并预注册按 gene/animal/family/batch 隔离的外部测试集。
