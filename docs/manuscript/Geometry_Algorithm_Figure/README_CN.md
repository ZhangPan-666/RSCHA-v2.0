# R-SCHA 几何与算法图

本目录保存作者确认的四子图定稿、绘图代码及图注。几何图相对初始方位绕地理轴顺时针旋转 12 度；外球面采用稀疏白色虚网格，线宽 1.2 pt、不透明度 0.80。

## 文件

- `draw_geometry_algorithm.py`：可独立运行的 Python 绘图代码，不读取外部观测数据。
- `images/RSCHA_Geometry_Algorithm.pdf`：论文使用的矢量图，字体已嵌入。
- `images/RSCHA_Geometry_Algorithm.svg`：可编辑矢量版。
- `images/RSCHA_Geometry_Algorithm.png`：完整预览；同目录的 `_geometry.png` 为几何子图预览。
- `caption.tex`：与论文对应的 LaTeX 图注。
- `requirements.txt`、`manifest.json`：已验证的依赖版本和文件校验记录。

## 运行

使用 Python 3.12 或更高版本，在本目录执行：

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python draw_geometry_algorithm.py
```

默认输出到本目录的 `images/`；也可将输出目录作为第一个参数传入。PDF 元数据中的生成时间会变化，重绘后的文件哈希可能不同；几何及版式参数由代码固定。

## 正文对应

该图插入论文第 2 节，标签为 `fig:geometry_algorithm`，为当前图 1。正文使用的副本是 `../论文初稿/figures/Fig1_Geometry_Algorithm.pdf`。原有四张图的显示编号顺延，保留原始文件名以维持既有数值结果和图件来源记录。

以上路径描述研究工程的目录布局；在独立稿件源码压缩包中，正文位于解压根目录，本目录作为 `Geometry_Algorithm_Figure/` 一并提供，图件代码仍可独立运行。

图 (b) 概括三类函数的稳定计算；图 (c) 区分离线残差压缩与在线恢复；图 (d) 展示 Legendre 次数灵敏度、边界根与帽域归一化的联系。图中 Mehler 的 Riccati 标签表示归一化恒等式，运行时不进行 Riccati ODE 积分。示意图不按比例绘制。
