# Python, native memory and Octave validation

This suite follows IID's layered validation approach. It executes current LRFb
source and checks independent numerical references, native ASan/UBSan contracts,
and original MATLAB functions in Octave. Tests require successful output contracts
and numerical behavior, including half-turn rotations, empty pooling, oversampling
and regression gradients.

## Run

Requires Python 3.11, Clang/Clang++, and the **full VLFeat 0.9.20 source tree**.
The full source recovered locally is under HBPL `third_party/vlfeat-0.9.20/`;
LRFb bundles the gateways and historical binaries, not the complete C library.

```sh
python3.11 -m venv /tmp/lrfb-validation-env
/tmp/lrfb-validation-env/bin/pip install -r validation/requirements.txt
/tmp/lrfb-validation-env/bin/python validation/run.py \
  --vlfeat-root /path/to/full/vlfeat-0.9.20 \
  --output /tmp/lrfb-validation \
  --ip-data /path/to/Hand_Data/Data.mat

/tmp/lrfb-validation-env/bin/python validation/octave_validate.py \
  --octave /path/to/octave-cli \
  --vlfeat-root /path/to/full/vlfeat-0.9.20 \
  --native-work /tmp/lrfb-validation \
  --output /tmp/lrfb-octave-validation \
  --ip-data /path/to/Hand_Data/Data.mat
```

`--ip-data` is optional; without it IP uses synthetic data. Output directories
must be outside the source repository. `--reuse-build` for the Python/native
runner and `--native-work` for the Octave runner verify native source fingerprints
before reusing a build. Without `--native-work`, the Octave runner builds it itself.
The macOS GMM ASan/UBSan probe injects Clang's matching ASan runtime into Octave.
It currently supports macOS only, and requires all 1-5 output calls to succeed.
Leak detection is not part of these checks.

On this machine, the existing IID Conda Octave 10.3.0 environment is at
`/tmp/iid-octave`. Supply `--octave /tmp/iid-octave/bin/octave-cli` and
`--conda-prefix /tmp/iid-octave` to enable the runner's child-process compiler
relocation workaround. Normal installations should omit `--conda-prefix`.
Temporary environments may disappear after system cleanup.

## Checks

- A limited MISS_HIT AST interpreter reads and executes the actual quaternion,
  RRV, temporal pooling, sampling, label and regression-objective source.
  Its indexing and core supported semantics are checked; unsupported syntax fails.
  It is not a general MATLAB interpreter. Preprocessing and pipeline orchestration
  are explicit Python bridges, additionally cross-checked in Octave.
- Actual C++ circumcenters are compared with independent linear-system solutions.
  Native MBS is compared with independent convex-hull widths, including the
  recognizer's monotone-axis constraint. Shared kernels in cross-language checks
  are not independent algorithm proofs.
- VLFeat C Fisher output is compared with independent diagonal-GMM score formulas.
  The arm64 library build disables x86 SIMD. The bundled LIBSVM source is used.
- Native memory checks allocate exact caller-sized output arrays: circumcenter
  1/2/3 outputs and GMM 1/2/3/4/5 outputs. The production single-output geometry
  bridge and the actual three-output GMM encoder call execute repaired gateways.
- Regression cases cover half-turn/near-half-turn rotations, empty/occluded/zero
  Fisher vectors, arbitrary cyclic oversampling, invalid sampler inputs, and
  finite-difference gradients for both IP and MSRC-12 (including zero parameters).
- Octave rebuilds six MEX modules directly from project source, with **no gateway
  padding adapters**. It checks single/double GMM output equivalence, optional
  outputs, parameter errors and independent quaternion reconstruction.
- Actual IP/MSRC-12 loaders, GeneSC, the 50000-target encoder, saved descriptors,
  linear SVM and model reload run on small fixtures. IP's loader trims its own
  unused allocation; the test does not trim cells. The fixed sampler receives
  the original training cells; the test does not repeat or reweight them.

The Octave path retains one format adapter: `octave_compat/save.m` substitutes
v5 MAT for unsupported MATLAB v7.3 and supplies implicit `.mat` filename extensions.
It is only added during validation. Original historical MEX files are never replaced.
MATLAB deployment must rebuild the changed gateways; see the instructions below.

## Scope

Python synthetic runs use 18 train/6 test gestures per dataset, K=64, `[1 2 4 8]`
and historical `jointNum=2`. All training descriptors feed GMM in that bridge.
The real-IP check validates the 16000-record split and probes 15 selected records;
its tiny K=4 classifier is an interface test. Octave uses 9 independent train/6 test
samples per dataset, K=4 and the original 50000-target sampling loop. Its 51
cross-language comparisons use the same fitted GMM where needed.

No reported tiny-sample accuracy is a benchmark result. Full `run.m`, complete
IP regression training, PCA, alternative baselines, 2-D MBS, full benchmark data,
MATLAB toolboxes/ABI, v7.3 serialization and shipped old binaries remain outside
this suite's certification. Logs, MAT data and compiled products stay outside
the repository.

## Rebuilding the repaired MATLAB gateways

Changing source does not update the archived MEX binaries. Octave MEX files cannot
be used as MATLAB MEX files. On the target MATLAB platform, prepare a supported
compiler, full VLFeat 0.9.20 headers and a matching C library, then build into a
separate directory:

```matlab
root = pwd; % Run from the LRFb project root.
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

Both functions must resolve to the new builds first, and the dynamic library must
be loadable. These are MATLAB build instructions, not a claim that the Octave suite
certifies MATLAB compilation or ABI compatibility. Other required MEX modules must
also match the target platform.
