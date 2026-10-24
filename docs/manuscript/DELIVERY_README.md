# R-SCHA 论文交付说明

更新时间：2026-10-24

## 当前交付版本

- 稿件修订标识：`manuscript-20261024-r2`（软件代码目录仍为 v2.0；当前发布编号为 v1.0）
- 论文：`Zhang et al., 2026, Manuscript_RSCHAFunctions`
- 当前 PDF：38 页、5 张图、50 条参考文献
- 正文源文件：`Zhang et al., 2026, Manuscript_RSCHAFunctions.tex`
- 参考文献：`references.bib`
- 编译参考文献：同名 `.bbl`
- 正式 PDF：同名 `.pdf`
- 本次核查定稿 PDF：`Zhang et al., 2026, Manuscript_RSCHAFunctions - manuscript-20261024-r2.pdf`；旧同名 PDF 当前被桌面预览锁定时，以该修订文件和源码包为准。
- 投稿亮点：`JCP_Highlights.txt`

## 目录内容

- `figures/`：正文图件和高清图件预览。
- `reproducibility/`：数值验证请求、参考生成脚本、谱线数据、反演参数和应用残差审计。
- `submission/`：JCP 投稿信草稿、声明文本和投稿准备说明。
- `Archive/`：正式修改前的备份。

## 2026-10-24 第二次修订：几何图与全文核查

- 定稿几何与算法图插入第 2 节，当前为图 1；原四张图的显示编号顺延，保留原始文件名。
- 工程根目录新增 `Geometry_Algorithm_Figure/`，包含可运行绘图代码、依赖、PDF/PNG/SVG、图注与校验清单；代码已重绘验证，PNG 与作者确认的图件逐像素一致。
- 补写局部坐标、球锥壳边界和算法概览；相关算法章节引用图 1 的对应子图。
- 阻尼记号采用 `lambda_d`，停止容差采用 `tau_stop`，避免与谱参数重名；补明总强度异常的线性投影算子。
- 将“解析出更短波长结构”改为“可表示更短波长结构”，保留独立分辨率验证的限制；修正 Zhang et al. (2025) 的正式页码为 686–696。
- 检查标签、交叉引用、50 条参考文献、全部图表与主要统计；六点 Python 参考再生和原有 24 项控制检查通过。没有重新执行完整 C++ 精度实验、性能测试或区域反演。
- 核查报告位于 `reproducibility/editorial_audit_20261024/`。AI 使用声明按作者于本次核查中的明确要求继续暂缓，仅保留草稿和待办。
- 当前 PDF 和独立源码包同步到 GitHub `docs/manuscript/`；当前软件 release 为 v1.0。

## 2026-10-24 第一次修订（历史稿，37 页）

- 论文题目保持 **Stable Computation of Revised Spherical Cap Harmonics**。
- 第 5 节及其参考误差评价小节的 “error measures” 改为 “error evaluation”。
- 基金资助更新为国家自然科学基金（42574077、42250103、42674085）及湖北省自然科学基金（2024AFE003）；尚未获批的实验室开放基金不列入。
- 第一作者 Pan Zhang 关联单位 1、2、3，保留现单位并补列地大学院及湖北省重点实验室。
- 投稿声明与说明同步更新；图件、数值源码和历史验收结果未改动。
- 已按 `pdflatex + bibtex + pdflatex + pdflatex` 重新编译，仍为 37 页、50 条参考文献；已检查首页、两处标题、基金及相邻页，最终编译日志无未定义引用或版面溢出警告。
- 当前稿件及可独立编译的源文件压缩包同步至 GitHub 的 `docs/manuscript/`，修订标识为 `manuscript-20261024`；当前软件 release 为 v1.0。

## 代码与下载

R-SCHA v1.0 的 OpenMP 和 MPI 源代码、查表文件、构建脚本、验证汇总和数学基础手册位于：

- GitHub 项目：https://github.com/ZhangPan-666/RSCHA-v2.0
- v1.0 发布版：https://github.com/ZhangPan-666/RSCHA-v2.0/releases/tag/v1.0

## 编译

在本目录依次运行：

```powershell
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
bibtex "Zhang et al., 2026, Manuscript_RSCHAFunctions"
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
```

## 后续修改与 GitHub 同步要求

**后续任何修改都必须同步更新 GitHub。** 适用范围包括正文、参考文献、C++ 代码、OpenMP/MPI 版本、图件、复现脚本、验证数据、投稿附件和 README。

每次修改后按以下顺序执行：

1. 在本地完成修改并重新编译、检查 PDF 和必要的数值结果。
2. 更新版本号、变更说明和源文件压缩包。
3. 将对应文件提交并推送到 GitHub；代码或发布包有实质变化时创建新的 release/tag。
4. 同步修改论文中的 **Data and code availability**，使下载地址和版本号与 GitHub 当前内容一致。
5. 在本目录保留最终 PDF、TeX、BibTeX、图件和复现材料，并在 `Archive/` 保存替换前版本。

不能只更新本地文件而不更新 GitHub；也不能让论文写的版本号、下载地址和 GitHub 实际发布内容不一致。

## 投稿前人工核对

上传前请核对作者姓名、单位、邮箱、ORCID、通讯作者、作者贡献、经费声明、利益冲突声明、全体作者同意投稿以及期刊系统的文件类型要求。

