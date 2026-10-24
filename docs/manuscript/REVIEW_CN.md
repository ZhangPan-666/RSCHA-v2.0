# 稿件全文核查记录

修订：`manuscript-20261024-r2`，核查日期：2026-10-24。

当前稿件为 **Stable Computation of Revised Spherical Cap Harmonics**，38 页、5 张图、12 张表、50 条参考文献。此次核查涵盖正文、公式的内部一致性、图表及标签、主要统计与现有记录的对应、参考文献元数据、复现材料和 PDF 版面。没有重新执行完整 C++ 数值验收、吞吐测试或区域反演，因此本记录不是新的数值精度认证。

## 已完成的修改

| 问题 | 本次处理 |
|---|---|
| 第 2 节缺少几何示意，参数与坐标轴难以对应 | 插入作者确认的四子图定稿，作为图 1；补充 O、z_g、局部球冠轴、P、r1/r2、S1/S2、theta/theta0/phi 的说明。图 1 在 PDF 第 3 页，处于第 2 节内。 |
| 算法图与正文衔接不足 | 第 2 节补充算法概览，第 4 节及其根、范数、残差压缩段落引用图 1 对应子图。图注明确三类函数的输出尺度及 Mehler 运行时不进行 Riccati ODE 积分。 |
| 插图后旧图号及说明可能失配 | LaTeX 标签自动顺延原四图为图 2–5；更新根目录及稿件、复现说明。原图文件名保留，以维持历史数据来源和哈希记录。 |
| 同一记号用于不同物理/数值量 | 反演阻尼改用 lambda_d，BiCGSTAB 停止容差改用 tau_stop；公式与说明同步，避免同前文谱参数 lambda、tau 混淆。 |
| 总强度异常算子未明确说明其线性化 | 补明其为异常磁场在预设背景场方向上的线性投影，与当前代码和原有审计说明一致。 |
| “resolves shorter-wavelength regional structure” 容易被理解为已经证明独立空间分辨率 | 改为“can represent shorter-wavelength regional structure”，明确独立观测仍是确认新增波长是否可分辨的必要证据。 |
| 一条文献保留在线预发表页码 | 按期刊正式目录，将 Zhang et al. (2025), DOI 10.26464/epp2025057 的页码由 1–11 修正为 686–696。 |
| 数据代码声明未直接指向当前稿件与绘图代码 | 补入稿件修订标识及 GitHub docs/manuscript 入口；软件仍为 v2.0.1。 |
| 绘图代码和图件分散在对话工作区 | 工程根目录新增 Geometry_Algorithm_Figure，保存可运行代码、依赖、图注、PDF/PNG/SVG 和清单；正文图件为定稿 PDF 的逐字节副本。 |

论文题目、作者次序、已确认基金、第一作者单位 1/2/3、既有四幅数据图、C++ 数值源码与历史验收数据保持原样。替换前稿件位于 Archive/Before_GeometryIntegration_20261024_101912。

## 核查与验证结果

- 108 个标签、81 次交叉引用：没有重复标签或缺失目标；算法引用改为标签引用。
- 50 条参考文献全部被引用：没有缺失引用键、重复文献键或未引用条目。BibTeX 无警告。
- 五张图、十二张表及其正文引用、图注/表注已检查。全部图件存在，图 1 PDF 的字体嵌入，图中无栅格图层；绘图代码重绘 PNG 与作者确认的定稿逐像素一致。
- pdflatex + bibtex + pdflatex + pdflatex 编译成功。最终日志没有未定义引用、重复标签、Overfull/Underfull 版面警告。38 页均已渲染并逐页检查版面；图 1 页另作放大检查，没有发现明显裁切、重叠、空白页或错位。
- 主要数值表中的精度、分位数、吞吐、灵敏度、根与范数数据，与历史 acceptance.json 在印刷舍入精度内一致；共 18 项统计/图件一致性检查通过。此次没有把历史记录重新绑定为当前源码验收。
- 新运行六个已归档高精度参考请求的 Python 再生，并运行参考程序原有 24 项控制，全部通过。仅为参考脚本可运行性检查，不等于重新完成 40,995 请求的 C++ 精度比较。
- 五条 Highlights 为 60–71 字符，均不超过 85 字符。
- 全部 39 个 DOI 条目进行注册元数据查询：32 项成功，7 项接口查询未成功，保留原始查询记录；“查询失败”不等于“文献不存在”。对未成功的关键算法文献另查作者/机构或期刊来源。Blakely et al. (2005) 的注册记录只给起始页 445，不据此错误地删除原文页码 445–448。
- MHGA 的 clip 写法与原论文的饱和函数等价，pseudo-vertical derivative 的说明与原文一致；当前代码中的 Hz/den−pi/3 与正文表达对应。应用于 Bz 且没有 RTP/MNLM 的处理差异已在正文说明。

## 仍须明确保留的限制与待办

1. **AI 使用声明继续暂缓。** 作者在本次核查中明确选择暂缓，故没有写入正文和图注。submission/AI_Declaration_Draft.txt 仅是草稿。Elsevier 当前政策要求披露此类稿件和图示辅助；投稿前应由作者完成最终说明及责任确认，本报告不把它标记为已完成。
2. 两个区域应用缺少实际运行日志和完整执行清单，无法从当前保存文件确认实际迭代次数、停止原因或证明收敛。正文仅将 100 写为上限，明确残差不能替代收敛证据；完整重复区域反演仍需要处理后的输入及处理历史。
3. 新增高次数谱功率与地质对应不能替代留出观测验证，也不证明相对于全球参考模型的独立精度优势。现稿明确区分数值基函数误差、代数自洽、拟合残差与空间分辨率。
4. 本轮对关键恒等式、符号、尺度和实现描述作内部核查，并与现有记录对照；没有提供全参数域的数学误差证明。对整个连续参数域或全部 50 篇文献内容的独立验证不在此次检查覆盖之内。
5. JCP 投稿界面的最新必填项、全体作者同意、资助及单位最终表单仍需投稿时核对；本次没有登录或提交期刊系统。

## 主要外部核对来源

- 正式页码：Earth and Planetary Physics, 2025, 9(3)，[期刊目录](https://www.eppcgs.org/article/2025/3)，Zhang et al. 文章为 686–696。
- MHGA 原文：[期刊文章页](https://journal.geo.sav.sk/cgg/article/view/514)，公式 (10)–(11) 及 pseudo-vertical derivative 说明。
- Gil et al. (2012)：[CWI 机构记录](https://ir.cwi.nl/pub/18881)，[作者出版列表](https://personales.unican.es/segurajj/publirk.html)。
- AI 披露：[Elsevier 官方政策](https://www.elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals)。按作者指示列为待办。

## 同目录的证据文件

- structure_audit.json：标签、引用、图件和 Highlights 检查。
- statistics_audit.json：18 项统计/图件一致性检查。
- bibliography_registry_audit.json：39 个 DOI 的注册元数据查询结果。
- reference_regeneration_summary.json、reference_self_check.json：本次六请求及 24 控制结果。
- build_qa.json：页数、图表数、参考文献数及编译日志检查。
- MANUSCRIPT_SHA256.json 位于源码包根目录；逐文件 SHA-256 记录源码包内容。历史 release.json 与原数值验收记录不作覆盖。
