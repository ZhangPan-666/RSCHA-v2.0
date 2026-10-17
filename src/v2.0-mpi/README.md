# v2.0mpi：直接sbatch提交

本目录是独立部署目录，包含全部项目头文件、SHRSCHA子目录、cr.bin和其校验文件。整目录内容上传至Linux共享文件夹RSCHAC++后可独立使用。修改run.pbs开头的用户操作区；脚本会自动加载桌面样例中的Intel编译器及MPI环境：

```bash
WORK=1                         # 1=反演；2=VBT正演；3=SH转R-SCHA
DATA_DIR="data/Work1"           # 测区输入/输出目录
```

只执行选中的一个Work。平时只改WORK这一项；资料目录不变时DATA_DIR无需修改。测区二改为data/Work2，所有dat直接放在代码目录时改为"."。DATA_DIR相对提交目录解释，也可以写绝对Linux路径，路径有空格会正确处理。

从放有run.pbs的目录直接提交：

```bash
sbatch run.pbs
```

队列cpu，作业时限无限（Slurm --time=0，受分区及账号限制）。默认40节点×45MPI进程、每进程1线程，每节点申请170GiB。脚本会自动核对cr.bin、检查必要输入、只编译选中的程序、运行并将输出留在DATA_DIR。文件名run.pbs是历史名称，脚本由Slurm sbatch提交。

Work1 defaults: IncludeTMI=true, IncludeVector=false, IncludeDiffTMI=true, IncludeDiffVector=true, matching the current near-surface TMI and satellite endpoint-difference data. Enable IncludeVector only when all three single-point vector input files are supplied. The submission script checks every required file and then runs a single-process --validate-only pass before starting the full inversion. This pass reads and validates inputs without constructing a design matrix or calculating basis functions.

Work2读取Results_RSCHAModel.dat和Parameters_Coordinates.dat，输出Results_ForwardResults.dat；Work3读取Parameters_Modelling.dat和Parameters_SHModel.dat，外部SH及侧边配置文件可选，输出Results_RSCHAModel.dat及谱和元数据。Work1初值、权重和先验文件必须准备。两测区请分别用data/Work1及data/Work2存放资料和结果。

Intel或GCC编译器通过Makefile自动检测，MPI需加载匹配的MPI模块；模块名称按集群提供的环境。常规MPI启动使用mpirun，集群要求srun时可将脚本下方MPI_LAUNCHER默认值改为srun并配置集群PMI/PMIx。cr.bin在本目录，不依赖另一版本或根目录CR-Atlas/CRTab。

既有实验验收与 release 属于旧源码快照，保持原记录。当前 Work2 的本机合成回归见进度条说明；上传源码不包含 Windows 可执行文件。

## 根据桌面集群样例修正编译器初始化

先加载/opt/intel/bin/compilervars.sh intel64，再检测C++编译器并传给make。MPI版也会在存在时加载桌面Makefile指定的Intel MPI 2019.7.217的intel64/bin/mpivars.sh。桌面Makefile原来是Fortran配置，mpiifort不能用来编译本项目C++；对应MPI C++使用mpiicpc（Classic）或mpiicpx（oneAPI），OpenMP版使用icpc/icpx，另支持已有GNU环境。

无需在操作区添加加载命令，仍然选择WORK=1/2/3和DATA_DIR后sbatch run.pbs。新增日志会打印加载路径和实际C++编译器版本。make clean不再因为未加载编译器而在解析Makefile时退出。保留set -eo pipefail、本地查表、cpu无限时限和原并行资源。

桌面参考文件原样保留；本机没有编译或计算，集群实际编译仍须实际运行验证。

## 旧集群标准库兼容

构建采用C++14，不再使用std::clamp、std::filesystem和std::make_unique。新RSCHACompatibility.h提供等同比较逻辑的clamp，以及本项目所需的路径、目录和文件状态接口；Windows使用宽字符系统接口，Linux使用getcwd/stat/mkdir。数值容差、根求解策略、查表数据与计算公式未改。旧验收仅绑定旧源哈希，不在本次重绑定。集群实际编译需重新提交验证。


## Fortran-style MPI inversion architecture (2026-10-15)

The MPI sum helpers retain their public names but use MPI_Gatherv, explicit rank-order summation on rank 0, and MPI_Bcast. The numerical sum helpers do not call MPI_Allreduce or MPI_Reduce; Work2 uses an integer-only MPI_Allreduce for progress counters. Root sum scratch storage is bounded to 64 MiB using blocks; counts and displacements remain within the legacy MPI int range. The design matrix stays local to each rank. Final predictions use indexed MPI_Gatherv and MPI_Bcast, with validation of complete, unique global row ownership.

Work1 distributes each data family in contiguous blocks, including the three components of each differenced vector pair. Kernel completion messages follow all-rank barriers. Rank-0 local progress and elapsed-time stage messages are flushed immediately; each rank also emits a final kernel-ready marker.

For fixed data weights (RobustSwitch=0), the preconditioner and RHS are constructed once and BiCGSTAB carries its shadow residual, current residual, search direction, alpha, omega, and beta across iterations. A recursively converged residual is verified with the actual weighted normal operator. When robust weights change, the Krylov state is restarted and the weighted system is rebuilt. Sigma, Huber/Tukey choices, damping, prior model, inactive coefficients, input filenames, output layouts, and solver function signature are retained. Iteration trajectories can differ from the previous one-step restart solver; historical numerical acceptance records remain tied to their original source snapshots.

run.pbs still selects WORK=1/2/3 and DATA_DIR. Work1 now performs the existing one-rank input validation, followed by --mpi-check on the full allocation with a 240-second communication timeout, before any root calculation or design matrix construction. --mpi-check is an optional diagnostic argument requiring no input data or atlas. A failed or timed-out check aborts the submission script. WORK=2/3 retain their interfaces and use the new sum helpers.

Work2 gathers only each rank's owned entries for the ten VBT fields, instead of collecting sparse full-length vectors from every rank. Geographic coordinate filtering, basis evaluation, tensor rotation, and the 13-column output format remain unchanged.

### Work2 正演进度条

Work2 的 `forward local VBT calculation` 阶段现在由 rank 0 输出全局进度条。所有 MPI rank 按相同的本地工作轮次调用一次 `MPI_Allreduce` 汇总已完成坐标数，因此不会因不同 rank 的计算速度不同而错位；计算结果、坐标分配、MPI 汇集和输出文件格式不变。

输出示例：

```text
[Forward progress] [##############--------------------------] 35.00% ...
```

进度条约每阶段 25 次更新，并显示完成数、百分比、全局点速率和估计剩余时间。进度条只覆盖本地 V/B/T kernel 计算；进入 `forward owned-value gather and broadcast` 后仍使用原有阶段日志。若某个单点 kernel 长时间不返回，进度条会在该点完成后的下一次同步更新，不能作为单点内部的实时心跳。

Work2 精简（2026-10-17）：坐标筛选改为一次遍历原地压缩，保持原顺序和筛选阈值；合并重复的基函数初始化，移除多余临时变量及未使用的头文件，结果写入使用缓冲换行。MPI 的 10 个 LocalCart 结果数组只在 rank 0 分配，每个非主进程减少 80 × 有效坐标数 字节的数组存储；474496 个坐标时约为 36.2 MiB/非主进程。

已用本机 Intel oneAPI C++14/OpenMP 编译，并用小型合成数据完成修改前后回归：OpenMP 使用 2、4、8 线程，MPI 使用 1、3、6 进程，覆盖 3、41、96 个有效坐标、内外阶数大小的三种关系、无 Mehler 项、无效坐标筛选及空任务进程。共 7 组对比的 13 列结果文件逐字节一致，全部进度完成到 100%。集群 Intel 2019 和完整 1800 进程作业仍需实际提交验证。

Only small synthetic regression and communication tests are performed locally. No application dataset is inverted or forwarded locally. No executable is included in the upload patch. The installed Intel 2019 cluster runtime and the full 40-node/1800-rank job must still be checked by actual submission; local checks do not certify that cluster run.
