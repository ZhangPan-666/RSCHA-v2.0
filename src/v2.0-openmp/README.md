# v2.0：直接sbatch提交

本目录是独立部署目录，包含全部项目头文件、SHRSCHA子目录、cr.bin和其校验文件。整目录内容上传至Linux共享文件夹RSCHAC++后可独立使用。修改run.pbs开头的用户操作区；脚本会自动加载桌面样例中的Intel编译器环境：

```bash
WORK=1                         # 1=反演；2=VBT正演；3=SH转R-SCHA
DATA_DIR="data/Work1"           # 测区输入/输出目录
```

只执行选中的一个Work。平时只改WORK这一项；资料目录不变时DATA_DIR无需修改。测区二改为data/Work2，所有dat直接放在代码目录时改为"."。DATA_DIR相对提交目录解释，也可以写绝对Linux路径，路径有空格会正确处理。

从放有run.pbs的目录直接提交：

```bash
sbatch run.pbs
```

队列cpu，作业时限无限（Slurm --time=0，受分区及账号限制）。默认1节点、1进程、45个OpenMP线程。脚本会自动核对cr.bin、检查必要输入、只编译选中的程序、运行并将输出留在DATA_DIR。文件名run.pbs是历史名称，脚本由Slurm sbatch提交。

Work1 defaults: IncludeTMI=true, IncludeVector=false, IncludeDiffTMI=true, IncludeDiffVector=true, matching the current near-surface TMI and satellite endpoint-difference data. Enable IncludeVector only when all three single-point vector input files are supplied. The submission script checks every required file and then runs a single-process --validate-only pass before starting the full inversion. This pass reads and validates inputs without constructing a design matrix or calculating basis functions.

Work2读取Results_RSCHAModel.dat和Parameters_Coordinates.dat，输出Results_ForwardResults.dat；Work3读取Parameters_Modelling.dat和Parameters_SHModel.dat，外部SH及侧边配置文件可选，输出Results_RSCHAModel.dat及谱和元数据。Work1初值、权重和先验文件必须准备。两测区请分别用data/Work1及data/Work2存放资料和结果。

### Work2 OpenMP 正演进度条（2026-10-16）

OpenMP 版 Work2 原先只有阶段日志，现在在 V/B/T kernel 正演阶段增加坐标完成进度。每完成一个坐标，由原子计数器累加；约每阶段25次进入命名的 OpenMP critical 区域输出一行，避免线程输出交错及进度倒退。串行回退分支也会逐点计数。

```text
[Forward progress] [##############--------------------------] 166074/474496 (35.00%); 120.0 points/s; ETA 42.8 min
```

显示内容包括完成数、总数、百分比、平均点速率和估计剩余时间；每次输出均换行并立即刷新，适合查看 Slurm 输出文件。进度只覆盖 kernel 正演，后续坐标旋转和写文件仍使用原有阶段日志。若长时间没有坐标完成，也不会产生新的进度行；这不是单点计算内部的定时心跳。

Work2 保留动态调度、数值计算及13列输出格式。坐标筛选改为一次遍历原地压缩，保持原顺序和筛选阈值；合并重复的基函数初始化，移除多余临时变量及未使用的头文件，结果写入使用缓冲换行。OpenMP 进度回调合并为一个，使用单调时钟计时。

已用本机 Intel oneAPI C++14/OpenMP 编译，并用小型合成数据完成修改前后回归：OpenMP 使用 2、4、8 线程，MPI 使用 1、3、6 进程，覆盖 3、41、96 个有效坐标、内外阶数大小的三种关系、无 Mehler 项、无效坐标筛选及空任务进程。共 7 组对比的 13 列结果文件逐字节一致，全部进度完成到 100%。集群 Intel 2019 和完整 1800 进程作业仍需实际提交验证。

把更新后的 Work2_RSCHAForward.cpp 上传至原源码目录，保持此文件名，设置 WORK=2 后重新提交 run.pbs；make 会根据源码修改重新编译。已经启动的旧二进制不会自动获得此功能。

Intel或GCC编译器通过Makefile自动检测，MPI需加载匹配的MPI模块；模块名称按集群提供的环境。常规MPI启动使用mpirun，集群要求srun时可将脚本下方MPI_LAUNCHER默认值改为srun并配置集群PMI/PMIx。cr.bin在本目录，不依赖另一版本或根目录CR-Atlas/CRTab。

既有实验验收与 release 属于旧源码快照，保持原记录。当前 Work2 的本机合成回归见进度条说明；上传源码不包含 Windows 可执行文件。

## 根据桌面集群样例修正编译器初始化

先加载/opt/intel/bin/compilervars.sh intel64，再检测C++编译器并传给make。MPI版也会在存在时加载桌面Makefile指定的Intel MPI 2019.7.217的intel64/bin/mpivars.sh。桌面Makefile原来是Fortran配置，mpiifort不能用来编译本项目C++；对应MPI C++使用mpiicpc（Classic）或mpiicpx（oneAPI），OpenMP版使用icpc/icpx，另支持已有GNU环境。

无需在操作区添加加载命令，仍然选择WORK=1/2/3和DATA_DIR后sbatch run.pbs。新增日志会打印加载路径和实际C++编译器版本。make clean不再因为未加载编译器而在解析Makefile时退出。保留set -eo pipefail、本地查表、cpu无限时限和原并行资源。

桌面参考文件原样保留；本机没有编译或计算，集群实际编译仍须实际运行验证。

## 旧集群标准库兼容

构建采用C++14，不再使用std::clamp、std::filesystem和std::make_unique。新RSCHACompatibility.h提供等同比较逻辑的clamp，以及本项目所需的路径、目录和文件状态接口；Windows使用宽字符系统接口，Linux使用getcwd/stat/mkdir。数值容差、根求解策略、查表数据与计算公式未改。旧验收仅绑定旧源哈希，不在本次重绑定。集群实际编译需重新提交验证。
