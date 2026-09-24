# Fisher vector 代码检索记录

检索日期：2026-09-24。对象：LRFb 及本地旧工程、用户指定的三个 SMB 路径。

本记录的逐字节相同比对发生在当天代码 polish 之前；后续格式与注释变更见
[代码审阅记录](CODE_REVIEW.md) 和 `../tools/polish_manifest.tsv`。
本文件及配套哈希表保留为检索时的历史证据。

## 结论

**Fisher vector 主实现存在，本地和远程均已实际读到。** 需要区分三件事：

1. LRFb 的 GMM / Fisher vector / 时间金字塔主流程已经在项目中。
2. VLFeat 的底层 C 实现也存在于本地 HBPL 和 9592… 远程共享；LRFb 只打包了接口、帮助文件和预编译二进制，没有打包完整 C 库源码。
3. 表 1 的 DI、II、SRVF、Multiscale SSM 描述子接入 FV 的**原始实验脚本**，本次检索尚未定位。不能据此断言其不存在或已丢失。

上次 Claude Code 的实际表述是“原来把这些对比方法接入 Fisher vector 流程的代码没找到”，并非主编码器不存在。其后检索因连接中断、共享掉线及会话额度耗尽而没有完成交接。

## 已确认的主流程

```text
experiments/{IP,MSRC12}/run.m
  -> GeneSC / RRV_DB.mat、RRV_SAMPLES.mat
  -> encoding/GeneFisherCodeJointPyramid_whole.m
       -> vl_gmm
       -> encoding/fv_pooling_ts.m
            -> vl_fisher(..., 'Improved')
  -> LIBSVM
```

| 层次 | 当前 LRFb 路径 | 核查结果 |
|---|---|---|
| GMM 学习、训练/测试集编码 | `encoding/GeneFisherCodeJointPyramid_whole.m` | 与 `_moved_HBPL_remote/MicrosoftGestureEvaluatingCode/` 原件逐字节一致 |
| 时间金字塔 FV | `encoding/fv_pooling_ts.m` | 与迁移保存的 IP、MSRC-12 两份原件逐字节一致 |
| MATLAB 到 C 的接口 | `thirdparty/vlfeat-0.9.20/toolbox/fisher/vl_fisher.c` | 调用 `vl_fisher_encode`；不是底层算法本身 |
| MATLAB 帮助文件 | `thirdparty/vlfeat-0.9.20/toolbox/fisher/vl_fisher.m` | 仅帮助文本，计算由同名 MEX 承担 |
| 预编译程序 | `thirdparty/vlfeat-0.9.20/toolbox/mex/*/vl_fisher.*` | 六种旧平台文件已在项目中；文件存在不代表在当前 MATLAB 上可加载 |

远程 IP 的 `GeneFisherCodeJointPyramid_whole.m` 与当前 LRFb 仅有一行 `fprintf` 的 `test data` / `testing data` 差异，没有计算差异。

## 完整 VLFeat 源码位置

本地根目录：

```text
/Users/perryshao/Documents/Projects/HBPL/third_party/vlfeat-0.9.20/
```

远程根目录（已重新挂载并读取）：

```text
/Volumes/9592df65-5b85-4531-8d5e-3ec20c74afc8/work/vlfeat-0.9.20/
```

两处的以下文件均经 SHA-256 比对一致：

- `vl/fisher.c`：真正的 Fisher vector 计算与归一化实现。
- `vl/fisher.h`：接口和标志定义。
- `vl/gmm.c`：GMM 实现。
- `toolbox/fisher/vl_fisher.c`：MATLAB MEX 包装。
- `toolbox/fisher/vl_fisher.m`：帮助文本。

例如，`vl/fisher.c` 的 SHA-256 为
`f7c3299934cfb84d37e1cf1631605f023b36cbc866df1882dcf3f92aa653b9b9`。
完整路径、字节数和哈希见同目录 `FISHER_CODE_HASHES_2026-09-24.tsv`。

本次没有搬动或替换这些原件。若后续需要从源码编译，应从该完整 VLFeat 树准备依赖，而不是仅复制 `fisher.c` 两个文件。

## 远程与迁移文件

| 位置 | 结果 |
|---|---|
| `/Volumes/9592df65-5b85-4531-8d5e-3ec20c74afc8/work/IPEvaluatingCode/` | 找到 `GeneFisherCodeJointPyramid_whole.m`、`GeneFisherCodeJointPyramid.m`、`GeneFisherCodeJointPyramid_2mod.m`、`GeneScFisherCodeJointPyramid.m` |
| `/Users/perryshao/Documents/Projects/_moved_HBPL_remote/MicrosoftGestureEvaluatingCode/` | 找到 MSRC-12 编码器及 `fv_pooling_ts.m`；它们已不在原共享目录的相同位置 |
| `/Users/perryshao/Documents/Projects/_moved_HBPL_remote/IPEvaluatingCode/` | 找到 `fv_pooling_ts.m`、`scfv_pooling_ts.m` |
| `/Users/perryshao/Documents/Projects/HBPL/extra/ip_dataset/` | 有多版本编码器，`__remote.m` 是保存的远程变体；不能只按文件名把它们当成相同版本 |
| `/Volumes/work` 与 `/Volumes/OneTouch/MatlabProjects/work` | 本次索引的 6486 个相对路径完全一致；顶层列表及抽检的三个源码文件也一致，与已有记录中的共享别名关系相符 |

IP 的 `GeneFisherCodeJointPyramid.m` 读取 `SC_DB` / `SC_SAMPLES` 并使用 FV 编码；这是仍然存在的 shape-context 接入路径。`GeneShapeContextJointPyramid.m` 也保留相关调用。它们的存在不能证明其他四项基线的原始接入脚本已找到。

另检查了远程 `IPEvaluatingCode/run.asv`：它保存了 `GeneSC -> GeneFisherCodeJointPyramid(..., pcaFlag=1)` 的旧流程，随后运行二元回归分类器，SVM 处于注释状态。它没有 DI / II / SRVF / Multiscale SSM 的选择分支，而且引用未定义的 `jointGroup`，不能作为可直接运行的论文表 1 驱动。

## 其他 Fisher 实现

- 本地 `Projects/Work/yael_v438/matlab/yael_fisher.c`，以及 TSSM、远程共享里的 YAEL 副本。
- `/Users/perryshao/Documents/Papers/IEEE T. CSVT-2018/Fisher Vectors/inria_fisher_v1/inria_fisher_v1/compute_fisher.m`：图像描述子示例，调用 `yael_fisher`，不是 LRFb 现有驱动调用的实现。
- Work 的 `NTU3DActionEvaluatingCode/recognition_ssm_NTUA_bat.m` 包含 SSM/RRV 接入 Fisher 编码的旧分支，但针对 NTU 数据，不能直接认定为 IP 表 1 的实验代码。

## 论文核对及未解决部分

附件第 4 页（印刷页 36118）给出 GMM K=64、D=14、15 个时间分段及 FV 长度 `2*K*D*15`。第 5 页（36119）明确写明表 1 各描述子共用编码流程和线性 SVM；本次同时检查了该页的渲染图。

因此，原 `baselines/README.md` 根据 Table 1 / Table 2 的 86.87% 相同而推测直接引用结果，证据不足，已经删除该推测。仍应寻找原始基线实验脚本。`GeneBaselineDB.m` 是整理时新写的适配代码；其中 SSM 是单尺度替代实现，不能当成已找回的 Multiscale SSM 原码，也没有复现准确率证据。

## 检索范围与限制

- 本地文件名：Documents、Downloads、Desktop 中 Fisher / FV 相关文件。
- 本地内容：Documents 下的 Projects、Papers、GitHub，以及 Desktop、Downloads 中的 `.m`、`.asv` 等代码；补查 `GeneFisher` 调用、`vl_fisher`、`yael_fisher`、`fv_pooling`、DI / II / SRVF 等描述子关键词。
- 三个 SMB 地址均已成功挂载。分别建立所列工作目录的文件名索引，6486 / 4999 / 6486 条；这是所选源码、Fisher 文件及归档文件名的计数，不是所有文件总数。
- 在两套不同目录树中递归搜索 MATLAB/C++/Python 源码与自动保存文件；另外检查 `.m~` / `.bak`。OneTouch 别名通过独立索引和抽样内容比对核验。
- 相关检索清单与远程命中行保存在本地 `docs/fisher-search-evidence-2026-09-24/`，可复核检索范围。该目录含工作站文件索引及会话摘录，不纳入 Git 仓库。
- 没有挂载或检索整个 Time Machine `.sparsebundle`，没有解包所有归档，也没有穷尽共享工作目录以外的历史磁盘备份。本次的“尚未定位”限于所述范围。
- 本次为源码定位、静态调用核查和哈希比对，没有运行 MATLAB 实验、重新编译 MEX，或验证论文准确率。没有修改算法、删除或移动原始文件。
