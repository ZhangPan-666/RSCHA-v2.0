# Work3：SH 到 R-SCHA 的旋转与解析投影

按[仓库说明](../../README.md)准备查找表并编译。程序可以从源代码目录运行，也可将 `RSCHA_CR_TABLE_DIR` 设为包含 `cr.bin` 的绝对路径。

## Work3：双源 SH 正投影

```text
Work3_LocalizeSHModel.exe [input_directory] [output_directory]
```

可省略两个目录参数，输入默认当前目录、输出默认输入目录。输入保留 Parameters_Modelling.dat、Parameters_SHModel.dat，新增可选 Parameters_SHModel_External.dat；缺省外部源为零。程序直接写输出，建议每次指定不同输出目录。

SH 输入为 n=0..N 完整三角表，按 n/m 排序，四列 n,m,g,h。采用实 Schmidt 半归一化、无 Condon–Shortley 相位，R_E=6371200 米。磁场入口 n_min=1，n=0 为占位行；m=0 正弦项不参与实势函数。两份源可以有不同最大度数，补零到共同 N。

建模参数九行顺序：帽心经度、纬度、半角（度）、内部 Legendre 最大 k、外部 Legendre 最大 k、Mehler 径向上限 p、Mehler 与 Constant 共用的周向上限 m、上球面半径、下球面半径（米）。

可选 `Parameters_SideProjection.dat` 含两行非负整数，依次覆盖建模参数的第 6、7 行。Mehler 的每个余弦或正弦系数块有 `P_Max*(M_Max+1)` 项，按径向阶优先排列。Constant 的每块有 `M_Max+1` 项。实际使用的截断写入 `Projection_Metadata.txt`。

Results_RSCHAModel.dat 保留原九行头部和四族 G/H 排列，可由当前 Work2 直接读取；另输出 Projection_Metadata.txt 和 Results_RSCHASpectrum.dat。该谱仅包含两族 Legendre。

新实现用逐度 Wigner 被动旋转，包含 ii/ie/ei/ee 四个源交叉块和两族侧边界。Mehler、常数族现在包含 R_E/a，因此相对旧版这两族按此比例改变。近重合度数用帽区积分求极限。旧 LocalizeSH2RSCHA.h 内部源 API 已委托新模块，要求完整三角表、n_Start=1，传入根缓存应对应相同半角和截断。
