# LRFb 联合验证复核（2026-09-26）

> 历史记录：本页描述修复前状态。已确认的 MEX、四元数、池化、采样、梯度和 IP 空槽问题
> 后续修复及新的验证结果见 [2026-09-26 修复报告](DEFECT_FIXES_2026-09-26.md)；本页原始证据保留。

参照 IID 的 Python 数值对照、C++ ASan/UBSan 与 Octave 原代码执行方法，
本次重新编译并运行 LRFb 已有验证套件，未直接沿用 9 月 25 日的结果。
[机器可读结果](VALIDATION_RECHECK_2026-09-26.json)包含实际命令、退出码、源码哈希和全部检查项。

## 结论与结果

**适配后的主方法小型流水线可运行，核心数值对照通过；未经适配的原 MEX 仍有两处输出越界。**
这次是验证复核，没有修复历史算法，也没有替换随项目保存的 MEX 二进制。

| 层次 | 本次实际执行结果 |
| --- | --- |
| Python 与原 C/C++ | 19 组中，15 组正常检查通过，4 组按预期复现缺陷；没有意外断言失败 |
| 独立数值对照 | 100 个三角形圆心、30 条 MBS 轨迹、18 组 Fisher 向量，以及描述子、池化等检查通过 |
| C++ 内存检测 | 正常夹具通过 ASan/UBSan；圆心单输出复现 `stack-buffer-overflow`，GMM 三输出复现 `heap-buffer-overflow` |
| Octave 10.3.0 | 重建六个 MEX，51 项与 Python 的数值比较全部通过，最大绝对误差 `7.985723193826288e-12` |
| 真实 IP 数据 | 核对 16000 条记录，训练/测试各 8000；15 条固定抽样均产生有限描述子 |
| 原 `.m` 流水线 | IP 真实子集及 MSRC-12 合成输入的加载、GeneSC、GMM/Fisher、LIBSVM 训练预测和模型保存重载完成 |
| 源码保护 | 运行前后两份项目全部已有文件哈希一致（排除 Git 元数据及 Python 缓存）；Octave 另记录并核对 135 个研究/依赖源码文件 |

两个执行阶段均退出为 0。这里的成功退出只说明断言符合预期，
**4 组已知缺陷复现不是算法正确性通过项**；GMM 越界由 Octave 阶段另外确认。
正常夹具的 sanitizer 检查不覆盖所有输入，且没有进行内存泄漏认证。

Python 合成 IP/MSRC-12 流水线各使用 18 个训练和 6 个测试样本，K=64；
原 `jointNum=2` 行切分产生 53760 维 FV，而论文布局为 26880 维。本次保留并验证这一历史行为。

Octave 两组各使用 9 个独立训练和 6 个测试样本、K=4，FV 为 3360 维。
为运行原 50000 目标采样器，仅重复训练描述子槽位：IP 423 槽、MSRC-12 207 槽。
IP 预测命中 4/6，MSRC-12 命中 6/6；Python 微型真实 IP 分类为 1/6。
这些小样本结果仅说明接口连通，不能作为论文准确率或实现间性能差异的证据。

## 本次重新确认的缺陷

1. `mbs/src/tricircumcenter3d.cpp:459–462` 不检查输出数量便写入第 2、3 个输出槽。
   原局部参考系函数只请求一个输出，ASan 捕获越界写入。
2. `thirdparty/vlfeat-0.9.20/toolbox/gmm/vl_gmm.c:341` 无条件写入第 5 个输出槽。
   原编码器只请求三个输出，真实网关的 ASan 探针捕获堆缓冲区越界。
3. `Myrotm2quat` 遇到 180° 旋转、`fv_pooling_ts` 遇到空序列时产生非有限值。
4. `rand_sampling_ts` 从 6 帧序列请求 13 帧时索引越界。
5. IP 回归目标第二正则项的解析梯度约为 2，中心有限差分约为 4。
6. IP 加载器预分配 8000 个槽位，小型输入留下的空槽使 `getLabels` 报错；
   本次核对的完整数据恰好各 8000 条，不受该小样本问题影响。

这些都是对现有问题的再次确认，没有把缺陷修复混入验证工作。
建议先修复两个 MEX 输出契约并重新编译，再处理数值边界和梯度问题。

## 测试适配与覆盖边界

为继续验证流水线，临时网关分别为圆心和 GMM 提供 3/5 个输出槽，并释放多余输出。
Octave 不支持原 `save -v7.3`，测试路径上的适配器改用 v5 MAT 并补 `.mat` 扩展名；
小型 IP 夹具裁剪空槽，采样器使用上述训练槽位重复。适配器只存在于测试路径及输出目录。
因此，51 项通过验证的是明确适配后的执行，不证明原网关安全或原驱动可直接运行。

Python 的部分 MATLAB 函数使用受限 AST 执行器；其他流程使用显式 Python 转写。
Octave 执行实际 `.m`，但分段和部分几何比较共享 C++ 内核，不能视为完全独立的算法证明。
完整 `run.m`、IP 的 SVM 前置回归训练、PCA、其他基线、2-D 分支、完整真实数据评估、
MATLAB 工具箱/ABI 和旧二进制兼容性仍不在此次覆盖内。MSRC-12 使用合成输入。

## 重跑与证据

通用命令及环境要求见 [validation/README.md](../validation/README.md)。
本次使用 Python 3.11.7、NumPy 1.26.4、SciPy 1.11.4、MISS_HIT 0.9.44、Clang 与 Octave 10.3.0。
Python 复用了本机 Anaconda 数值包及已有虚拟环境中的纯 Python MISS_HIT；
Octave 使用 IID 的临时 Conda 环境及现有 runner 的显式重定位配置，没有安装或升级依赖。

完整本地证据目录：
`/Users/perryshao/Documents/ChatGPT/旧项目整理/LRFb-validation-20260926/`。

- `rerun.py`、`execution.json`：重跑入口、精确命令、时间和退出码。
- `native/results.json`、`octave/results.json`：各层详细结果。
- `native/one-output-asan.log`、`octave/gmm-asan.log`：两处越界的实际捕获日志。
- `native/build.log`、`octave/build.log`、`octave/octave.log`：编译及原 `.m` 执行日志。
- `before-files.json`、`other-before-files.json`：执行前两份项目文件哈希。

仓库只新增本复核报告和 JSON，并在 README 添加入口；实验数据和编译产物保存在上述独立目录。
