# LRFb — Describing Local Reference Frames for 3-D Motion Trajectory Recognition

MATLAB/C++ implementation of the **local reference frame based (LRFb) descriptor** for
3-D point trajectories: local reference frames estimated with **3-D maximal blurred
segments (MBS)**, their **rotation (quaternion) and square-root velocity**, **Fisher-vector
encoding with a temporal pyramid** on a **GMM**, and a **linear SVM** — together with the
comparison descriptors of Table 1.

This directory is a consolidated, **code-only** reconstruction of a project that was
spread across two networked MATLAB workspaces, local folders and the HBPL clean-up.
Intermediate results (`.mat`), figures and raw datasets are deliberately **not**
included — see [Provenance](#provenance) and [Datasets](#datasets).

The [2026-09-24 Fisher-code search](docs/FISHER_CODE_SEARCH_2026-09-24.md)
verifies local and remote encoder copies, locates the full VLFeat C source, and
distinguishes the existing FV pipeline from baseline experiment scripts that have
not yet been located.

The [code review](docs/CODE_REVIEW.md) records the 2026-09-24 formatting and comment
cleanup, verification, and unresolved runtime issues. See [CONTRIBUTING.md](CONTRIBUTING.md)
for the formatting tools and conventions.

The [2026-09-26 defect fixes](docs/DEFECT_FIXES_2026-09-26.md) repair the two MEX
output overflows, half-turn quaternions, empty pooling, oversampling, regression
gradients and small-IP loader padding. All 19 Python/native checks and 51
Octave/Python comparisons pass; the Octave pipeline uses the repaired gateways
directly. See [validation/README.md](validation/README.md) to rerun the checks.
Historical evidence remains in the [numerical report](docs/NUMERICAL_VALIDATION.md),
[Octave report](docs/OCTAVE_VALIDATION.md) and
[pre-fix rerun](docs/VALIDATION_RECHECK_2026-09-26.md).
This is not full MATLAB certification or reproduction of the paper. The archived
MEX binaries are unchanged: rebuild the two corrected gateways before MATLAB use.

---

## Paper

Z. Shao, Y. Li, Y. Guo, and X. Zhou.
"Describing Local Reference Frames for 3-D Motion Trajectory Recognition."
*IEEE Access*, vol. 6, pp. 36115–36121, 2018.
[doi:10.1109/ACCESS.2018.2849690](https://doi.org/10.1109/ACCESS.2018.2849690)

Evaluated on **InteractPlay (IP)**, 16 hand gestures, stereo-tracked hands (88.95 %), and
**MSRC-12**, 12 Kinect gestures, left + right hand joints only (90.38 %).

## Licence

GPL-3.0 — see [`LICENSE`](LICENSE), the same licence as the companion IID and TSSM
projects. The toolboxes under `thirdparty/` keep their own licences; see
[`thirdparty/README.md`](thirdparty/README.md). The descriptor files and
`TrjResizeTime.m` carry a `by skaegy` header.

---

## Method at a glance

```
  3-D trajectory Γ(t), t = 1..N   (each hand; resampled to 64 frames on IP, 128 on MSRC-12)
        │
        ├─ lrf/, mbs/      Estimate_Frenet: digitise on a 1000-cell grid (curve_grid[_mean]),
        │                  split into maximal blurred segments of width v (splitting_curve_3D
        │                  → Determine_segment mex), take the left/right key points L(t), R(t),
        │                  circumcircle of Γ(L), Γ(t), Γ(R) (tricircumcenter3d mex)
        │                  →  F(t) = [t n b]  per point                          (§II-A, Eq. 3)
        │
        ├─ descriptor/     RRVdescriptor_BasedonFrames: rotation between adjacent frames as a
        │                  unit quaternion (Myrotm2quat) + square-root velocity ẋ/√‖ẋ‖
        │                  projected by that rotation  →  s_t ∈ R^7 per hand     (§II-B, Eq. 4–8)
        │                  GeneSC puts the two hands side by side  →  D = 14
        │
        ├─ encoding/       GeneFisherCodeJointPyramid_whole: GMM (vl_gmm, K = 64) on 50 000
        │                  sampled s_t; fv_pooling_ts: improved Fisher vector (vl_fisher) on
        │                  every segment of a temporal pyramid 2^0..2^3 (15 segments),
        │                  concatenated and L2-normalised                          (§III-A-1)
        │
        └─ run.m           linear SVM, LIBSVM  svmtrain(..., '-c 1 -t 0 -b 1')
```

## Paper → code

| Paper | Step | Code |
|---|---|---|
| §II-A, Eq. (3), Fig. 3 | MBS decomposition of the digitised trajectory | `mbs/splitting_curve_3D.m` → `Determine_segment` (mex, `mbs/src`) |
| §II-A | key points L(t), R(t); circumcircle; `n`, `t`, `b = t × n` | `lrf/Estimate_Frenet.m` (+ `tricircumcenter3d` mex) |
| §III-A-2 | "digitalize … with 1000 grids" | `lrf/curve_grid_mean.m` (MSRC-12), `lrf/curve_grid.m` (IP) |
| §II-B, Eq. (4)–(8) | quaternion of the frame-to-frame rotation + local SRV | `descriptor/RRVdescriptor_BasedonFrames.m`, `descriptor/Myrotm2quat.m` |
| §II-B, Eq. (5) | SRVF `ẋ/√‖ẋ‖` | line `VelNorm_End = …` in `RRVdescriptor_BasedonFrames.m` |
| §III-A-2 | two hands concatenated, D = 14 | `experiments/*/GeneSC.m` → `RRV_DB.mat`, `RRV_SAMPLES.mat` |
| §III-A-1 | GMM, K = 64 | `encoding/GeneFisherCodeJointPyramid_whole.m` (`vl_gmm`) |
| §III-A-1 | FV on a temporal pyramid, Z = 3, average pooling | `encoding/fv_pooling_ts.m` (`vl_fisher …'Improved'`) |
| §III-A-1 | linear SVM, default parameters | `svmtrain` / `svmpredict` block in `experiments/*/run.m` |
| Table 1 | Differential / Integral invariants, SRVF, multiscale SSM, 3-D shape context | `baselines/` — see [`baselines/README.md`](baselines/README.md) |
| Fig. 5, Fig. 7 | confusion matrices | `confusion_matrix` at the end of each `run.m` |

## Code vs. paper

The historical computations below are preserved. Formatting and comments were
polished on 2026-09-24; see `tools/polish_manifest.tsv` for before/after hashes.

| Item | Paper | Code as last saved |
|---|---|---|
| MBS width v | 8 (§III-A-2; the response to reviewers says 5–10 perform alike) | **5**: `width = 5` in `MSRC12/GeneSC.m`, `Estimate_Frenet(…, 5)` in `IP/GeneSC.m` |
| Frame rotation | Eq. (8) R(t) = F⁻¹(t−1) F(t) | `RotM = F(t+1) / F(t)` (forward neighbour, = F(t+1)F(t)⁻¹); the quaternion is taken of `RotM'`, the velocity is `v * RotM` |
| FV length | 2·K·D·15 with D = 14 | `run.m` passes `jointNum = length(bodyJoints) = 2`, so `GeneFisherCodeJointPyramid_whole` splits each gesture's **rows** (time) into two halves and encodes each half separately: 2 × 2·K·14·15. Pass `1` for the paper's layout |
| Resampling | not stated | `TrjResizeTime` to 64 frames (`IP/load_IPtxt_bat.m`) and 128 frames (`MSRC12/load_data_bat.m`) |
| Classifier | linear SVM | IP: `run.m` **also** trains the binary-regression classifier (`trainBinRegression_whole`) just before the SVM; the SVM runs last and overwrites `predict_label`, so the reported number is the SVM's |

Other things worth knowing:

- `Myrotm2quat` now handles 180° rotations with a stable component-based formula.
  `RotM` itself remains undefined if a frame row is all zero.
- `fv_pooling_ts` preserves empty/fully occluded/zero pooled vectors as zeros and
  no longer computes the unused per-frame Fisher pass.
- `-g 2.79e-4` in the `svmtrain` options has no effect with `-t 0` (linear kernel).
- MATLAB before R2018a has its own `svmtrain` in the Statistics Toolbox; `setup_path` puts
  LIBSVM on the path, but check `which svmtrain` if results look odd.

---

## Directory layout

```
LRFb/
├── lrf/          Estimate_Frenet (MSRC-12 version), curve_grid, curve_grid_mean,
│                 remove_stapoint, fillstapoint
├── mbs/          splitting_curve_3D/_2D, Determine_xy, transform_firstOctant
│   ├── src/        Determine_segment.cpp, tricircumcenter3d.cpp (+ unit_reco, memory pools)
│   └── bin/        prebuilt mex: a64, w32, w64, maci64 (+ MelkmanConvexHull w32/w64)
├── descriptor/   RRVdescriptor_BasedonFrames, Myrotm2quat   — the LRFb descriptor
├── encoding/     GeneFisherCodeJointPyramid_whole, fv_pooling_ts, rand_sampling_ts
├── utils/        TrjResizeTime
├── baselines/    Table 1 descriptors + GeneBaselineDB (new adapter, see its README)
├── experiments/
│   ├── IP/         InteractPlay driver and dataset-specific helpers
│   └── MSRC12/     MSRC-12 driver and dataset-specific helpers
├── thirdparty/   vlfeat (vl_gmm / vl_fisher only), libsvm-3.17, netlab dist2, ndSparse,
│                 ScSPM subsets, Stochastic_Bosque
├── tools/        provenance.tsv — origin and md5 of every file
└── setup_path.m
```

`experiments/*` keep their own copies of helpers whose IP and MSRC-12 versions differ:
`Estimate_Frenet` (IP: 2 arguments, `curve_grid`; MSRC-12: 3 arguments, `curve_grid_mean`
— the latter is the one in `lrf/`), `GeneSC`, `getLabels`, `sc3d_compute` (IP: pseudo 3-D
shape context, 135 bins; MSRC-12: true 3-D, 810 bins), and the binary-regression files.
Helpers whose two versions differed only cosmetically were merged into the shared
folders (the MSRC-12 copy was kept): `GeneFisherCodeJointPyramid_whole` (one `fprintf`
string) and `rand_sampling_ts` (`[a b]` vs `cat(2,a,b)`).

### Menu alternatives in the drivers

The drivers are *menus*: alternatives sit in commented blocks and one is uncommented. The
blocks as last saved run the paper's pipeline. The other blocks belong to the parallel
HBPL / RRV work and are **not** in this paper, but their code is included so every block
still resolves:

| Block in `run.m` | Code |
|---|---|
| batched FV, 3-D shape context (`SC_DB`), two modalities, sparse-coding FV | `GeneFisherCodeJointPyramid`, `…_2mod`, `GeneShapeContextJointPyramid`, `GeneScFisherCodeJointPyramid`, `GeneScCodeJointPyramid` (IP) |
| cost-sensitive binary regression | `trainBinRegression*`, `predictBinRegression*`, `costFunc`, `costFuncRegMultPartGp_v1_42`, `minimize` |
| second linear SVM | `li2nsvm_multiclass_lbfgs` / `_fwd` (thirdparty/ScSPM) |
| random forest | `Stochastic_Bosque` (thirdparty) |

The batched variants still carry the server paths they were written for
(`/home/data/IPdataset/…`, `C:/Users/perry/…`), and `GeneScCodeJointPyramid` loads a
learned dictionary `Results/reg_sc_b256_*.mat` that is not bundled. The paper's path
(`GeneFisherCodeJointPyramid_whole`) writes everything into the current folder.

## Running

```matlab
cd LRFb
setup_path                 % adds lrf/, mbs/, descriptor/, encoding/, baselines/ and thirdparty
cd experiments/MSRC12      % or experiments/IP
run                        % the driver; it starts with `clear all`
```

`experiments/` is intentionally **not** on the path (see `setup_path.m`). The drivers
write `DB.mat`, `SAMPLES.mat`, `RRV_DB.mat`, `RRV_SAMPLES.mat`, `traindata.mat`,
`testdata.mat` into the current folder, and `load modelForTest` expects a
`modelForTest.mat` written by the encoder in the same run.

Requirements: MATLAB with the Statistics Toolbox (`pca`, only if `pcaFlag = 1`). The
`SSM` baseline also needs the Image Processing Toolbox (`imfilter`).

## Datasets

Not bundled.

| Folder | Dataset | What the loader expects |
|---|---|---|
| `experiments/IP` | InteractPlay (Just, Bernier & Marcel, BMVC 2004) | `Hand_Data/Data.mat`: a cell array {trajectory, [gesture id, person id]} prepared by Y. Guo. `load_IPtxt_bat` ignores its `BAT_FOLDER` argument and reads this file; persons 1–10 train, 11, 13–15 and 17–22 test; gestures 17 and 18 are skipped. The reader for the raw `.txt` files (`readIPtxt`) is kept in a commented block |
| `experiments/MSRC12` | MSRC-12 Kinect Gesture (Fothergill et al., CHI 2012) | `MSRC12_Skeleton.mat` holding `SmthTrj{class, subject}{joint, instance}`: the sequences already segmented at the labelled action points. Odd subjects train, even subjects test; class 6 / subject 18 is skipped. The script that produced this file was not found |

Copies seen during consolidation (not touched): `Hand_Data/` in
`~/Documents/Projects/Work/IPEvaluatingCode/` and on the 9592… share, the raw IP data in
`HBPL/extra/ip_dataset/IPdataset/`, and `MSRC12_Skeleton.mat` (240 MB) in
`~/Documents/Papers/IEEE Access/MicrosoftGestureDataset-RC/`.

---

## Building the mex files

| Binary | Source | Prebuilt for |
|---|---|---|
| `Determine_segment`, `tricircumcenter3d` | `mbs/src` | a64, w32, w64, maci64 |
| `MelkmanConvexHull` | **none** — only reached from `splitting_curve_2D` (2-D input) | w32, w64 |
| `vl_gmm`, `vl_fisher` | `thirdparty/vlfeat-0.9.20/toolbox/{gmm,fisher}/*.c` + full VLFeat tree | a64, glx, maci, maci64, w32, w64 |
| `svmtrain`, `svmpredict` | `thirdparty/libsvm-3.17/matlab` (`make.m`) | w32 (`matlab/`), w64 (`windows/`) |
| `best_cut_node`, `mx_eval_cartree` | `thirdparty/Stochastic_Bosque/…/mx_files` | a64, w64 |

```matlab
cd mbs/src
mex Determine_segment.cpp  % already includes unit_reco.cpp
mex tricircumcenter3d.cpp
```

There are no Apple-silicon (`maca64`) builds of anything; rebuild from source there.
`Determine_segment.cpp` shares its origin with IID. The circumcenter gateways
now carry project-specific output-contract fixes; they are not byte-identical.
For the repaired MATLAB gateways and path precedence, see the
[rebuild instructions](docs/DEFECT_FIXES_2026-09-26.md#matlab-中使用修复后的-mex).

---

## Provenance

Every file's origin and original md5 is in [`tools/provenance.tsv`](tools/provenance.tsv):
185 verbatim copies, 2 modified, 1 new at initial consolidation. These are the
original-source hashes; `tools/polish_manifest.tsv` tracks the subsequent cosmetic edits.

The paper's code is the **Feb 2018 state of the second SMB share**
`smb://perry-System-Product-Name._smb._tcp.local/9592df65-…/work` (the Ubuntu server the
drivers' `/home/data/…` paths point to):
`IPEvaluatingCode/` and `MicrosoftGestureEvaluatingCode/`. On 2026-09-21 the HBPL clean-up
**moved** most of those files to `~/Documents/Projects/_moved_HBPL_remote/` (dates
preserved) and copied them into `HBPL/extra/{ip_dataset,msrc12_gesture}` — there the share
versions are the `*__remote.m` files and the plain names are the older local versions.
Files still on the share were read from the share; moved files from
`_moved_HBPL_remote`. The share versions were chosen because they are the ones that match
the paper: two hands only, `GeneFisherCodeJointPyramid_whole`, K = 64, Z = 3, linear SVM.
The local versions (in the Trash since the HBPL clean-up; the HBPL copies have lost their
dates) run the RRV / binary-regression experiments instead.

Other sources: `Determine_xy`, `transform_firstOctant`, `MelkmanConvexHull`, the w32 mex
and the Table 1 invariants from the top level of the first share (`MatlabProjects/work`,
also mounted as `work`); the maci64 mex from IID; `Temporal_SSM`, `Log_hogcalculator`,
`libsvm-3.17` and `ScSPM` from TSSM.

At initial consolidation, none of the paper's files were modified. The two `modified` rows are `Temporal_SSM.m` and
`Log_hogcalculator.m`, taken from TSSM, which had fixed them (argument defaults, a
`nargin` test). New: `baselines/GeneBaselineDB.m`, plus `setup_path.m` and the READMEs.

No MATLAB was available during consolidation, so **nothing here has been run**.
The later code review found unresolved dependencies in the legacy 2-D MBS and SSM
baseline branches, plus runtime defects in the MEX and optional regression paths.
See [docs/CODE_REVIEW.md](docs/CODE_REVIEW.md); the earlier static dependency check
was not sufficient to establish that every alternative is runnable.

Not included:

- `preprocess_bat` (Kalman smoothing of C3D mocap), a stale commented line in the MSRC-12
  driver from the IID era;
- the rest of the old `ToZP/` RRV toolbox (`Func_RRVdescriptor`, `RRV_AngTrj*`, `CalSIGN`,
  …) — the RRV descriptor of Guo et al. (IEEE TCyb 2018), which HBPL keeps;
- all other MSRC-12 files from the SSM / IID era (DTW, HMM, clustering, retrieval), which
  live in IID and TSSM;
- `.mat` caches (`RRV_DB`, `SC_DB`, `traindata`, `Results/reg_sc_*`, `train_test*/`),
  figures and datasets.

### Relation to the other consolidated projects

| Project | Shares with LRFb |
|---|---|
| IID (`Integral-Invariants`) | the MBS machinery (`splitting_curve_*`, `Determine_segment`, `tricircumcenter3d`) and `Estimate_Frenet`; the Table 1 integral / differential invariants |
| TSSM (`Temporal-Self-Similarity-Description`) | the SSM baseline |
| HBPL | the same IP / MSRC-12 folders under `extra/`, labelled there as outside HBPL's papers — they are this paper's experiments. Note that `HBPL/extra/msrc12_gesture` has no `ToZP/`, so `RRVdescriptor_BasedonFrames` does not resolve there |
