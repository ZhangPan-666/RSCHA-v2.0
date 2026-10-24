# 数学手册全面复核报告（2026-10-24）

## 结果

本次对 `Mathematical Foundations of Revised Spherical Cap Harmonic Analysis` 的 24 个章节、图注、交叉引用、参考文献和最终 PDF 进行了重新核查。新版图件已全部改用 `figures/vector/*.pdf`，正文已重新编译为 167 页。

- LaTeX 编译成功；无未定义引用、无缺失文献键。
- 结构检查：1,222 个标签，重复标签 0，缺失引用 0，缺失图件 0，环境不平衡 0。
- 编译日志没有 overfull/underfull box，也没有 citation/reference warning。
- 逐页渲染 167 页并检查缩略图；几何图、索引图、三维基函数图、导数图、坐标旋转图和参考文献页均正常。

## 主要数学修正

1. **Mehler 展开统一为矩形索引。** 正文默认采用独立的 `1≤p≤P_max` 与 `0≤m≤M_Meh`。系数数目改为 `P_max(2M_Meh+1)`；`P_max=M_Meh=4` 时为 36 个实模式。旧的 `m≤p` 只保留为 legacy restricted model，并明确给出计数。
2. **三角截断完备性表述修正。** 在第 13 节加入 `Ψ_* = Ψ^c_{1,2}` 的明确遗漏模式。它与三角族及常数族在能量空间正交，因此三角族到它的误差有严格正下界；原先的“是否完备”开放表述改为“全空间不稠密，研究其闭包内的逼近率”。
3. **Legendre 递推式修正。** 第 9 节第一条递推改为 degree recurrence；第二条的分母改为 `sqrt(1-x^2)`。Schmidt 半归一化版本同步修正，并限制 order-lowering 公式的适用范围为 `m≥1`。
4. **Wigner-D 约定统一。** 按正文给出的 defining identity 重新核对系数转置和 Euler 首尾相位，并在映射章节明确 active/passive 矩阵约定。
5. **边界与径向细节修正。** 重复根 `μ=1/4` 明确只在退化的一表面构造中出现，不是非退化球壳的径向本征值；非整数外部支在原点的行为改成有限球壳内的陈述；极轴处说明了球坐标/NED 分量依赖经线，需先变换到固定 Cartesian frame 再取唯一极限。
6. **替代边界条件的陈述修正。** 侧面 Neumann、径向 Dirichlet 的互补分解补上 Green 恒等式证明，说明其也具有严格梯度正交性，但谱、权重和边缘条件需要独立推导。
7. **正文图件更新。** 所有 17 张正文图改用最新矢量 PDF；高密度 Legendre/Mehler 图册改为横向整页，常数族三张图整合在同一页，便于印刷阅读。

## 独立数值核验

`reproducibility/verify_handbook.py` 使用 45 位工作精度，独立检查 Legendre 递推与导数、根处范数、球谐投影边界公式、Mehler 双线性/参数范数、常数族 Lerch 范数、径向 Neumann 本征函数与闭式重叠、两边界系数解、笛卡尔场/Hessian、NED 六分量旋转、Wigner-D 旋转、体能量和固定半径表面 Gram 块。共 450 组比较通过；最大缩放误差约 `1.52×10^-13`，来自数值积分的体能量块。结果保存在 `reproducibility/handbook_math_validation.json`。

Legendre 递推和特殊函数的分支/递推约定参照 NIST DLMF 第 14.3 和 14.10 节：[DLMF 14.3](https://dlmf.nist.gov/14.3)、[DLMF 14.10](https://dlmf.nist.gov/14.10)。

## 仍保留的事项

- AI 使用声明按作者此前指示继续暂缓，未加入正文；它是投稿前的编辑待办。
- 第 23 节仍保留 R-SCHA2D 组合延拓空间、边缘奇异性、固定表面谱对角化和非球面边界等研究问题。
- 本次数值核验覆盖代表性参数和有限阶公式，不等同于对任意高阶、任意 cap aperture 或所有 C++ 实现路径的形式证明。
