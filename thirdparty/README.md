# Third-party code

These toolboxes are **not** covered by this project's GPL-3.0 licence. Each remains under
its own authors' copyright and terms. They are bundled because the pipeline calls into
them directly; each was reduced to the code it needs. Origins and md5s are in
`../tools/provenance.tsv`.

| Toolbox | Needed by | Licence |
|---|---|---|
| `vlfeat-0.9.20` (Vedaldi & Fulkerson) — only `vl_gmm` / `vl_fisher`: help stubs, C sources and the prebuilt mex + `libvl` for six platforms | `GeneFisherCodeJointPyramid_whole`, `fv_pooling_ts` — the GMM and the improved Fisher vectors | `COPYING` (BSD-2-Clause) |
| `libsvm-3.17` (Chang & Lin) — C core + `matlab/` + `windows/` | `svmtrain` / `svmpredict`, the linear SVM of both papers' tables (`-t 0`) | `COPYRIGHT` (BSD-3-Clause) |
| `netlab/dist2.m` (Ian Nabney, netlab 3.3) | `GeneSC` and `GeneBaselineDB`: the mean-distance scale | netlab licence (BSD-style); no licence file bundled |
| `ndSparse` (Matt Jacobson) | `experiments/MSRC12/sc3d_compute.m`, the true 3-D shape context | © Xoran Technologies 2010 in the header; MATLAB File Exchange; no licence file bundled |
| `ScSPM/large_scale_svm`, `ScSPM/sparse_coding` (Jianchao Yang) | menu alternatives in `run.m`: `li2nsvm_multiclass_lbfgs` (2nd linear SVM), `reg_sparse_coding` (sparse-coding encoders) | no licence file; research code, © Jianchao Yang / NEC Labs America in headers |
| `Stochastic_Bosque` (Leo / Phillip M. Feldman) | menu alternative in `run.m`: random-forest classifier | `license.txt` (BSD-2-Clause) |

`vlfeat`, `libsvm`, the ScSPM subsets and the Stochastic_Bosque tree are copies of what
the TSSM / HBPL projects consolidated from the same workspace; `libsvm-3.17` and `ScSPM`
were taken from the already-trimmed copies in `Projects/TSSM/thirdparty`.

> ⚠️ Three of these ship without a licence file (`netlab/dist2.m`, `ndSparse`, ScSPM), so
> their redistribution terms are unclear. That is fine while the project is private;
> revisit before publishing it.

## vlfeat

Only the two functions the pipeline uses were kept. `vl_gmm.m` / `vl_fisher.m` are help
text; the work is done by `toolbox/mex/<platform>/vl_gmm.*` and `vl_fisher.*`, which load
`libvl` (`libvl.so`, `libvl.dylib`, `vl.dll`) from the same folder. `setup_path` adds the
folder for the running platform. There is no Apple-silicon (`maca64`) build: install
VLFeat from <https://www.vlfeat.org/> and call its `vl_setup` instead. Rebuilding the mex
from `vl_gmm.c` / `vl_fisher.c` also needs the full VLFeat source tree.

## Not bundled

| Toolbox | Why |
|---|---|
| full `vlfeat-0.9.20` (36 MB) | only `vl_gmm` / `vl_fisher` are called |
| `Kalman`, `C3D` | reached only from `preprocess_bat`, a stale commented line in the MSRC-12 driver from the IID era |
| `HMMall`, `bnt`, `pmtk3`, `yael` | no caller in either driver |
