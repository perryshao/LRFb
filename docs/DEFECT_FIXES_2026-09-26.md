# LRFb 已确认缺陷修复与回归验证

日期：2026-09-26。根据[修复前联合复核](VALIDATION_RECHECK_2026-09-26.md)落实源码修复。
[完整结果 JSON](DEFECT_FIXES_2026-09-26.json)记录执行命令、源码哈希、数值误差和检查结果。

## 修复结果

| 已确认问题 | 本次修复 | 回归证据 |
| --- | --- | --- |
| 圆心 MEX 单输出越界 | 只分配请求的 1–3 个输出，未请求的 xi/eta 存于局部变量；提前检查参数数量和实 double 三元素输入 | 精确大小的 1/2/3 输出数组通过 ASan/UBSan；100 个随机三角形与独立线性方程解一致；Octave 实际调用通过 |
| GMM MEX 三输出越界 | 按输出数量分配均值、协方差、先验、似然和后验；支持 1–5 输出，拒绝不支持的数量 | 精确大小的 1–5 输出数组通过 ASan/UBSan；Octave 对 single/double、可选输出数值一致性、非法数量的检查通过 |
| 180° 四元数 Inf/NaN | 按迹或最大对角元素选择稳定计算分支，单位化并统一非负标量分量 | Python 对照 113 个旋转，Octave 用 Rodrigues 公式独立重建 30 个旋转，含三主轴、任意轴、精确及邻近半周旋转 |
| 空池化/零向量归一化 NaN | 仅当拼接向量范数大于零时归一化；空、全遮挡及零 Fisher 向量保持零 | Python/Octave 边界回归通过；删除了结果未使用、且会在遮挡过滤前调用 Fisher 的逐帧计算 |
| 采样量超过两倍帧数时越界 | 按需循环同一次随机排列，保持原有均分计数规则及旧有效范围内的顺序；非法数量和空训练序列明确报错 | 0–50000 请求、超过两倍帧数、排列循环及非法输入检查通过 |
| 回归正则项梯度错误 | IP 和 MSRC-12 均采用目标函数的准确导数 `2*theta.^3/sqrt(sum(theta.^4))`，零组导数设为零；行范数使用零次梯度 | 两个副本各覆盖标量、多类、多组和零参数；Python 18 例有限差分最大误差约 `1.70e-9`，Octave 最大误差约 `2.00e-9` |
| 小型 IP 输入残留 8000 预分配槽 | 保存前按实际训练/测试计数裁剪 | 实际加载器直接保存 9 条训练、6 条测试；标签及后续流水线无需测试端裁剪 |

回归函数同时修正了多组边界：按分组大小的累积和取连续块，MSRC-12 不再硬编码两个组。
这保持一组/两组情形的既有目标值，并纠正三组以上的重叠索引。
第一正则项原先在梯度中单独加 epsilon、但目标函数没有对应平滑项，现改为精确导数和零次梯度。
第三输出原本引用不存在的 Hessian，现明确报 `LRFb:hessianUnsupported`，不伪造 Hessian。
本次没有更改论文参数、`jointNum=2` 时间行切分或分类器选择。

GMM 是对随库网关的局部修复；没有升级 VLFeat，也没有改 GMM/Fisher 的底层公式。
两处 MEX 均需从新源码重建才会生效。

## 实际验证

- Python/C++ 共 **19 组检查全部通过**，没有保留“预期错误复现”作为成功条件。
- 正常原生夹具、圆心 1–3 输出、GMM 1–5 输出均无 ASan/UBSan 报错。
  GMM 的动态探针在真实 Octave 进程中执行，并采用严格的调用者输出数组大小。
- 六个 Octave MEX 重新编译，直接执行修复后的圆心和 GMM 源码；
  **51 项跨语言对照全部通过，最大绝对误差 `3.4589553443709065e-12`**。
- 真实 IP 16000 条记录的训练/测试划分核对通过，15 条抽样全部生成有限描述子。
- IP 真实子集和 MSRC-12 合成数据分别用 **9 条训练、6 条测试**，直接执行原加载器、
  GeneSC、50000 目标采样器、K=4 GMM/Fisher、LIBSVM 训练预测和模型保存重载。
  测试中没有重复训练槽位，没有裁剪加载器输出，没有给 MEX 补输出缓冲区。
- 9 个变更相关 MATLAB 文件的 MISS_HIT lint 检查通过。

Octave 仅保留 `save -v7.3` 到 v5 MAT 及 `.mat` 扩展名的格式适配。
IP 小样本命中 1/6，MSRC-12 合成样本命中 6/6；它们只用于检查计算与数据接口。
修复前使用重复训练槽位，修复后直接采样，因此不能据此比较识别性能。
Python 合成 IP/MSRC-12 的 K=64 流水线也通过，保留原 53760 维 FV。

首次运行发现新增梯度测试的局部变量覆盖了命令行参数，已改名；首轮日志保留在证据目录。
修正后的完整 Python 和 Octave 运行均退出为 0；两份项目在运行期间的已有文件均未被测试改写。
所有 sanitizer 日志均检查了错误文本。未进行内存泄漏或所有可能输入的认证。

## MATLAB 中使用修复后的 MEX

**仓库的旧平台 MEX 二进制没有覆盖。仅更新 `.cpp`/`.c` 不会更新已编译的 MEX。**
本次实际重建和验证的是 Octave MEX；不能直接把它作为 MATLAB MEX 使用。
在目标 MATLAB 平台准备好对应的编译器、完整 VLFeat 0.9.20 头文件和同架构 C 库后，
可将两个修复后的网关编译到独立目录，例如：

```matlab
root = pwd; % 当前目录为 LRFb
outputDir = fullfile(root, 'build', 'fixed-mex');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
fullVl = '/absolute/path/to/full/vlfeat-0.9.20';
vlLibrary = '/absolute/path/to/compatible/libvl.dylib'; % Linux .so / Windows .lib
mex('-R2017b', '-outdir', outputDir, fullfile(root, 'mbs/src/tricircumcenter3d.cpp'));
mex('-R2017b', '-outdir', outputDir, ['-I' fullVl], ['-I' fullfile(fullVl, 'toolbox')], ...
    fullfile(root, 'thirdparty/vlfeat-0.9.20/toolbox/gmm/vl_gmm.c'), vlLibrary);
setup_path;
addpath(outputDir, '-begin');
clear tricircumcenter3d vl_gmm;
which tricircumcenter3d -all
which vl_gmm -all
```

确认两个 `which` 的首项均为新构建产物，并确保动态库可加载。
这是 MATLAB 重建说明，本轮没有执行 MATLAB 编译或 MATLAB ABI 验证。
其他旧平台模块仍须满足对应平台要求，特别是没有现成 Apple-silicon MATLAB MEX。

## 重跑和边界

通用入口见 [validation/README.md](../validation/README.md)。完整本地证据在：
`/Users/perryshao/Documents/ChatGPT/旧项目整理/LRFb-fixes-20260926/`。

其中 `native/results.json`、`octave/results.json`、编译日志、sanitizer 日志、
`rerun.py` 与 `execution.json` 可直接复核；`before/` 保存修改前的源码和文档。
`working-before.json` / `github-before.json` 是修复前快照；
`before-files.json` / `other-before-files.json` 是最终验证开始时各副本的快照。

未验证完整 MATLAB 驱动、全量真实基准、MATLAB 工具箱、IP 完整回归训练、PCA 或可选基线。
原局部参考系的退化/零帧处理、历史 2-D 分支及截断几何谓词等其他审阅问题未在本次修复范围内。
本次通过不等于论文复现或所有历史分支均已修复。
