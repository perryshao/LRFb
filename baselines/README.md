# Table 1 baselines

Table 1 of the paper compares LRFb with five existing per-frame descriptors on the IP
dataset. All six go through the same Fisher-vector temporal-pyramid encoding and the same
linear SVM.

| Table 1 row | Acc. | Implementation | Origin |
|---|---|---|---|
| Differential Invariants [6] | 41.50 | `descriptor_comp.m` — curvature, torsion and their derivatives (4-D) | share `work/` top level (also in IID) |
| Integral Invariants [4] | 53.58 | `integral_invariant.m` (+ `estimate_integral`, `normalization`) — 2-D | share `work/` top level (also in IID) |
| SRVF [20] | 72.64 | one line in `GeneBaselineDB.m`: `ẋ / √‖ẋ‖` (paper Eq. 5, world frame) | new |
| Multiscale SSM [23] | 86.87 | `Temporal_SSM.m` + `Log_hogcalculator.m` — single-scale SSM of raw xyz, 150-D log-polar HOG | TSSM `ssm/`, `descriptor/` |
| 3D shape context [13] | 65.83 | `experiments/*/sc3d_compute.m` on the estimated LRFs | the experiment folders |
| **LRFb** | **88.95** | `descriptor/RRVdescriptor_BasedonFrames.m` | — |

## What is original and what is not

**The original wiring has not yet been located.** In both experiment folders only the 3-D shape context is
still connected: its `sc3d_compute` call sits, commented out, in `GeneSC.m`, which writes
`SC_DB.mat`, and `experiments/IP/GeneFisherCodeJointPyramid.m` /
`GeneShapeContextJointPyramid.m` encode `SC_DB`. For the other four rows no generator was
found in any source location (the two SMB shares, local `Projects/Work`, HBPL, IID, TSSM).
The available files do not establish whether the other blocks were overwritten or
survive in another backup. The bundled GMM/Fisher implementation does not establish
that the original baseline experiment wiring has been recovered.

**`GeneBaselineDB.m` is new code**, written during consolidation and **never run** (no
MATLAB was available). It reconnects the descriptor implementations above to the
pipeline by writing them into `RRV_DB.mat` / `RRV_SAMPLES.mat` in exactly the layout
`GeneSC` produces:

```matlab
cd experiments/IP
load_IPtxt_bat(BAT_FOLDER);  [trainGID,testGID] = getLabels();
GeneBaselineDB('DI');                              % instead of GeneSC(bodyJoints)
GeneFisherCodeJointPyramid_whole(jointNum, ntotalbh, numClusters(1), 0);
% ... the svmtrain / svmpredict block of run.m, unchanged
```

`method` is one of `'DI'`, `'II'`, `'SRVF'`, `'SSM'`, `'SC3D'`. For `'SSM'` (150-D) and
`'SC3D'` (135-D on IP, 810-D on MSRC-12) pass `pcaFlag = 1` to the encoder, as the
original `GeneFisherCodeJointPyramid` did for `SC_DB`.

Parameters the paper does not state were set as follows. They are guesses; change them in
`GeneBaselineDB.m`:

| Baseline | Setting | Why |
|---|---|---|
| II | `integral_invariant(xyz, 6, 0.005)` | the ASL setting in IID / TSSM |
| SSM | `Temporal_SSM(xyz, 5, 1, 0)`, Log-HOG defaults | the Euclidean SSM-raw setting of TSSM |
| SC3D | 18 × 9 × 5 bins, `r_inner = 0.1`, `r_outer = 2` | the values in `GeneSC.m` |
| LRF width (SC3D) | 5 | same as `GeneSC.m` |

## Caveats

- The SSM adapter is not currently self-contained: `Temporal_SSM(..., 5, 1, 0)`
  requires `distance_matrix_norm2`, which is not bundled. Other metric modes
  also reference missing helpers. Resolve these dependencies before using the
  reconstructed baseline adapter.
- **Multiscale SSM [23] is not this code.** [23] is Guo, Li & Shao, *IEEE TII* 2017. Its
  multiscale implementation has not yet been located. The paper explicitly states
  that Table 1 uses the same encoding procedure and linear SVM for the compared
  descriptors (p. 36119). Equal numbers in Tables 1 and 2 do not establish that the
  results were quoted. `'SSM'` here is a substitute: Shao's single-scale SSM
  + log-polar HOG from the TSSM project.
- `Temporal_SSM.m` and `Log_hogcalculator.m` are TSSM's copies, which carry two small
  fixes made during that consolidation (defaults for the short call forms, a `nargin`
  test). `Log_hogcalculator` needs the Image Processing Toolbox (`imfilter`).
- `descriptor_comp`, `integral_invariant`, `estimate_integral` and `normalization` are the
  untouched originals from the share. IID's copies differ only in trailing whitespace.
