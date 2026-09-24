# LRFb 静态审阅与代码整理

日期：2026-09-24。

## 本次完成的工作

审阅覆盖描述子、局部参考系、MBS、FV 编码、IP / MSRC-12 驱动、基线适配、
回归备选分支、MEX 源码及依赖边界。统一整理了 **64 个 MATLAB 文件和 7 个 C++/头文件**。
第三方库检查了接入方式，并对其中 33 个 MATLAB 文件做语法检查；第三方源文件和二进制原样保留。

完成了缩进、运算符间距、多语句拆行、长行折行、LF 换行、尾部空格和冗余空行清理；
补充函数帮助、数据维度、MAT 文件读写和历史分支说明；修正了将平方误差写成 Rosenbrock /
logistic regression、将时间重采样写成弧长重采样、将 FV 池化写成稀疏编码等误导性注释。
删除了核心描述子中引用不存在变量的废弃注释片段，保留实验驱动中的历史备选菜单。

没有修改函数接口、计算表达式、实验参数、随机抽样顺序或分类器选择。
以下运行和数学问题仍然存在，不能把格式检查通过理解为仓库已经能够复现论文。

## 需要单独修复的问题

P1 表示可能导致崩溃、错误计算或关键实验无法执行；P2 表示边界条件、可移植性或维护问题。
这是静态审阅的已确认问题清单，不是对所有数值输入的完整正确性证明。

| 优先级 | 位置 / 触发条件 | 发现及影响 |
|---|---|---|
| P1 | `mbs/src/tricircumcenter3d.cpp` 的 `mexFunction` | 无条件写 `plhs[1]`、`plhs[2]`，而 `Estimate_Frenet` / `estimate_integral` 通常只请求一个输出。未按 `nlhs` 分配输出，有越界写风险；同时缺少 `nrhs`、类型、尺寸检查。`Determine_segment.cpp` 同样没有输入/输出数量防护。 |
| P1 | `mbs/src/tricircumcenter3d.cpp` 的 `orient2dadapt` | 最后一个误差界分支之后存在没有返回值的执行路径。接近共线且需要进一步精确运算时，非 void 函数可能落到底部，结果未定义。 |
| P1 | `encoding/GeneFisherCodeJointPyramid_whole.m` 与两个 `run.m` | `jointNum=2` 将已经按列合并的 T-by-14 描述子再沿行切成两半，输出维数为论文公式的两倍。与原始版本一致，但不能把它描述为与论文布局完全一致。 |
| P1 | 各 `GeneFisherCodeJointPyramid*` / 稀疏编码器，`pcaFlag=1` | GMM/字典学习使用经过特征值白化的投影样本，后续编码却只乘 `PcaM`，白化行被注释；训练码本与编码数据不在同一个尺度。默认主流程 `pcaFlag=0` 不经过该分支。 |
| P1 | 两套 `costFuncRegMultPartGp_v1_42.m` | 第二正则项是 `sqrt(sum(theta.^4))`，导数应含 `2*theta.^3/sqrt(sum(theta.^4))`；现有实现缺少系数 2（另有 epsilon 近似）。例如单元素 theta=2 时该项为 4，导数为 4，而现有代码约为 2。影响调用此目标的回归备选分支。 |
| P1 | `experiments/MSRC12/costFunc.m`，只传六个参数 | 缺省 `preTinitialTheta=0` 后立即 `reshape(..., D, C)`，在 `D*C>1` 时失败。`trainBinRegression_whole` 的调用正是六参形式；该训练器目前位于 MSRC-12 驱动的注释备选块。 |
| P1 | 两套回归目标函数，请求第三输出 | `ddf = ddgradient` 中 `ddgradient` 从未定义。当前 `minimize` 用的是目标值和梯度，但函数签名对第三输出的承诺不成立。 |
| P1 | `baselines/Temporal_SSM.m` 及 `GeneBaselineDB('SSM')` | 默认 SSM 分支调用 `distance_matrix_norm2`，该依赖未打包。其他模式的 `feature_dist_matching`、`distance_matrix_fd`、`distance_matrix_norm1`、`hist_cost_2` 也未提供。缺少依赖时该适配器不能直接运行。 |
| P1 | `mbs/splitting_curve_2D.m` / `Determine_xy.m` | `structure`、`testOctant`、`Determine_xy_sub` 未打包；部分路径还依赖只有 Windows 预编译版的 `MelkmanConvexHull`。3-D 主路径通常不经过它们，2-D 备选路径不具备完整依赖。 |
| P1 | `experiments/IP/scfv_pooling_ts.m` | `(feaSet-B*sc_codes)*sc_codes'` 得到 D-by-K 的字典原子统计，后续却按 1:T 的时间帧索引取列；T>K 时可能越界，T<=K 时索引含义也不一致。它是历史稀疏 FV 备选方案。 |
| P2 | 两套 `predictBinRegression.m`，测试样本少于一个 batch | 循环没有执行时，后续读取 `batchtimes+1` 使用未定义变量。部分批次编码器/训练器还假定始终存在尾批次文件；整除与不足一批都需单独测试。 |
| P2 | `encoding/rand_sampling_ts.m` 及调用它的采样循环 | 重复索引最多补一遍，单序列需要超过两倍帧数的样本时会越界；每序列取样数舍入为 0，或所有特征均为 0 时，上层 `while lastNsmp < nsmp` 可能不前进。 |
| P2 | `Myrotm2quat`、参考系/不变量、网格和池化归一化 | 180° 旋转的 `qw=0`、共线点的叉积范数为 0、常量轨迹的范围为 0、空轨迹，以及全零 FV 的 L2 范数为 0，都缺少一致的退化处理；可能产生 Inf/NaN 或索引错误。 |
| P2 | `mbs/transform_firstOctant.m` | 三处 `-pi < min_octant < -3*pi/4` 是链式比较，MATLAB 会先得到逻辑值再比较。MISS_HIT 报出三项 high-severity 检查结果。本次保留并记录，没有通过压制规则隐藏。 |
| P2 | `mbs/src/unit_reco.h` / `unit_reco.cpp` | `U`、`L` 数组容量为 1000，但从下标 1 开始存储且递增时没有容量检查；足够长且不断保留极值点的输入可能越界。网格宽度 1000 并不等于输入点数上限。 |
| P2 | `GeneSC`、旧驱动及备选训练器 | `GeneSC` 声明的 `SC_DB` 在返回前被 clear，正常使用依赖保存到磁盘的 RRV 文件；旧菜单仍使用硬编码服务器路径和 `eval`。IP 驱动会先运行回归训练，再由线性 SVM 覆盖预测结果。 |

原 README 的 MEX 构建示例同时编译 `Determine_segment.cpp` 和 `unit_reco.cpp`，
但前者已直接 include 后者，会形成重复定义。本次已将**文档示例**改为只编译入口文件；
没有改 C++ 的 include 结构，也没有声称 MEX 已成功编译。

对比基线的来源和缺失脚本仍以 [Fisher 检索报告](FISHER_CODE_SEARCH_2026-09-24.md) 为准。
`GeneBaselineDB` 中的单尺度 SSM 不是论文的 Multiscale SSM 原实现。

## 验证与保真

- 整理前建立完整归档，路径为
  `/Users/perryshao/Documents/ChatGPT/旧项目整理/backups/LRFb-before-polish-2026-09-24.tar.gz`。
- MATLAB：使用 MISS_HIT 0.9.44 解析全部 97 个 `.m` 文件（含第三方）；
  对自有 64 个文件执行项目格式规则检查。
- 自有 MATLAB 文件与备份逐个比对语法树、有效词法标记；仅忽略注释、布局换行及改成换行的语句末逗号。
  字符串、常量、函数名、转置运算和抑制输出的分号均参与比对。
- 7 个 C++/头文件通过 clang-format 23.1.1 检查，并用 Clang 原始词法标记比较前后代码。
  该比较不预处理、不链接，也不代替 MEX 编译测试。
- 第三方文件和原有 MEX 文件与备份逐字节比较。详情见
  `polish-verification.json` 和 `../tools/polish_manifest.tsv`。
- 最终结果：64 个 MATLAB 文件格式检查零问题；97 个 MATLAB 文件语法解析通过；
  191 个文件的保真比较全部通过（64 个 MATLAB、7 个 C++/头文件、120 个字节不变文件）。
  329 行行尾空格已清除。长语句拆行和帮助文档补充使源码总行数从 7067 增为 7363，
  因此总行数减少不作为本次整理的目标。
- `polish-style.json` 保存零格式问题结果；`polish-lint.json` 保存仍存在的三项链式比较警告。
- 本机未找到可调用的 MATLAB / Octave，本次没有执行真实轨迹数据、回归训练或数值复现实验。

完成本次代码 polish 时尚未初始化 Git。后续仓库发布使用整理后的文件作为首次提交；修改前归档保留在本地，前后哈希随仓库保存，可用于审阅。
