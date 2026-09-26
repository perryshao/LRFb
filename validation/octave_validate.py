"""Execute real Octave MEX and original LRFb MATLAB functions, with only MAT-format compatibility."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import subprocess

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.io import loadmat, savemat

from native import ROOT, build, reuse
from pipeline import (
    describe,
    encode,
    fisher_reference,
    frames,
    grid,
    pool_reference,
    resample,
    split,
)
from run import make_engine


def json_value(value):
    """Serialize numeric arrays returned by scipy.io.loadmat."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported report value: {type(value).__name__}")


def quote(path):
    return str(path).replace("'", "''")


def octave_env(prefix):
    env = dict(os.environ)
    if prefix:
        if platform.system() != "Darwin":
            raise ValueError("The IID Conda relocation workaround is macOS-specific")
        env.update(
            OCTAVE_HOME=str(prefix),
            OCTAVE_EXEC_HOME=str(prefix),
            CPPFLAGS="-DLRFB_VALIDATION_BUILD",
            CFLAGS="-O2",
            CXXFLAGS="-O2 -std=c++17",
            XTRA_CFLAGS="-fPIC",
            XTRA_CXXFLAGS="-fPIC",
            LDFLAGS=f"-L{prefix}/lib -Wl,-rpath,{prefix}/lib",
            DL_LDFLAGS="-bundle -Wl,-undefined,dynamic_lookup",
        )
    return env


def run_octave(args):
    work = args.output.resolve()
    if work == ROOT or ROOT in work.parents:
        raise ValueError("Output must be outside the repository")
    work.mkdir(parents=True, exist_ok=True)
    native_work = args.native_work.resolve() if args.native_work else work / "native"
    native = (
        reuse(native_work, args.vlfeat_root.resolve())
        if args.native_work
        else build(native_work, args.vlfeat_root.resolve())
    )
    engine = make_engine(native)
    env = octave_env(args.conda_prefix)
    original_paths = [
        p
        for p in ROOT.rglob("*")
        if p.suffix in (".m", ".c", ".cpp", ".h")
        and ".git" not in p.parts
        and "validation" not in p.relative_to(ROOT).parts
    ]
    source_hashes = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in original_paths
    }
    test_hashes = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (ROOT / "validation").rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    }
    # Test exact caller-sized storage without padding or numerical substitutions.
    gmm_prefix = (
        "#define mexFunction tested_gmm_gateway\n#include "
        + json.dumps(str(ROOT / "thirdparty/vlfeat-0.9.20/toolbox/gmm/vl_gmm.c"))
        + "\n#undef mexFunction\n"
    )
    (work / "gmm_output_probe.c").write_text(
        gmm_prefix
        + """
#include <stdlib.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mxArray **exact = (mxArray**)calloc(nlhs, sizeof(mxArray*));
  tested_gmm_gateway(nlhs, exact, nrhs, prhs);
  for(int i=0;i<nlhs;i++) plhs[i]=exact[i];
  free(exact);
}
"""
    )

    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    expression = f"addpath('{quote(ROOT / 'validation')}');build_octave('{quote(ROOT)}','{quote(work)}','{quote(args.vlfeat_root.resolve())}','{quote(native_work / ('libvl' + suffix))}');"
    compiled = subprocess.run(
        [str(args.octave), "--no-gui", "--quiet", "--eval", expression],
        env=env,
        capture_output=True,
        text=True,
        timeout=240,
    )
    (work / "build.log").write_text(compiled.stdout + compiled.stderr)
    if compiled.returncode:
        raise RuntimeError("Octave MEX build failed; see " + str(work / "build.log"))
    print("Built six real Octave MEX modules", flush=True)
    asan_expression = (
        "mkoctfile('--mex','-g','-fsanitize=address,undefined',"
        f"'-I{quote(ROOT / 'validation/octave_include')}', '-I{quote(args.vlfeat_root.resolve())}',"
        f"'-I{quote(args.vlfeat_root.resolve() / 'toolbox')}', '{quote(work / 'gmm_output_probe.c')}',"
        f"'-L{quote(native_work)}', '-lvl', '-o', '{quote(work / 'vl_gmm_output_probe')}');"
    )
    probe_build = subprocess.run(
        [str(args.octave), "--no-gui", "--quiet", "--eval", asan_expression],
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    (work / "gmm-asan-build.log").write_text(probe_build.stdout + probe_build.stderr)
    if probe_build.returncode:
        raise RuntimeError("GMM sanitizer build failed")
    asan_env = dict(env, ASAN_OPTIONS="detect_leaks=0:halt_on_error=1", UBSAN_OPTIONS="halt_on_error=1")
    if platform.system() == "Darwin":
        runtime = subprocess.check_output(
            ["clang", "-print-file-name=libclang_rt.asan_osx_dynamic.dylib"], text=True
        ).strip()
        asan_env["DYLD_INSERT_LIBRARIES"] = runtime
    else:
        raise NotImplementedError(
            "This Octave ASan host-injection probe currently supports macOS only"
        )
    probe = subprocess.run(
        [
            str(args.octave),
            "--no-gui",
            "--quiet",
            "--eval",
            f"addpath('{quote(work)}'); for n=1:5; out=cell(1,n); [out{{:}}]=vl_gmm_output_probe(rand(2,20),2); assert(all(cellfun(@(v) all(isfinite(v(:))),out))); end; disp('PASS GMM outputs 1 through 5');",
        ],
        env=asan_env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    (work / "gmm-asan.log").write_text(probe.stdout + probe.stderr)
    assert probe.returncode == 0 and "PASS GMM outputs 1 through 5" in probe.stdout, probe.stderr[-2000:]
    print("Passed fixed GMM output contract (1-5 outputs) with ASan/UBSan", flush=True)

    rng = np.random.default_rng(240925)
    t = np.linspace(0, 5, 69)
    raw = np.c_[np.cos(t), np.sin(t), 0.2 * t]
    repeated = np.insert(raw, [12, 35], raw[[11, 34]], axis=0)
    real = loadmat(args.ip_data)["Data"] if args.ip_data else None
    indices = []
    ip = np.empty((15, 2), dtype=object)
    for label in range(1, 4):
        if real is not None:
            ids = np.vstack([x.ravel() for x in real[:, 1]])
            mask = ids[:, 0] == label
            train = np.flatnonzero(mask & np.isin(ids[:, 1], np.arange(1, 11)))[::100][:3]
            test = np.flatnonzero(mask & ~np.isin(ids[:, 1], np.arange(1, 11)))[:2]
            indices.extend([*map(int, train), *map(int, test)])
        else:
            for s in range(5):
                t = np.linspace(0, 4 + 0.7 * label, 69 + s)
                h = (1 + 0.01 * s) * np.c_[np.cos(t), np.sin(t), (0.15 + 0.12 * label) * t]
                ip[(label - 1) * 5 + s, 0] = np.c_[h, 0.8 * h + 0.1]
                ip[(label - 1) * 5 + s, 1] = np.array([[label, s + 1 if s < 3 else s + 8]])
    if real is not None:
        ip = real[indices].copy()
    msrc = np.empty((3, 2), dtype=object)
    msrc_raw = [[], []]
    for label in range(1, 4):
        for subject in range(2):
            trajectories = np.empty((12, 3 if subject == 0 else 2), dtype=object)
            for index in np.ndindex(trajectories.shape):
                trajectories[index] = np.empty((0, 0))
            for sample in range(trajectories.shape[1]):
                t = np.linspace(0, 4 + 0.7 * label, 69 + sample + subject)
                h = (1 + 0.008 * (sample + subject * 3)) * np.c_[
                    np.cos(t), np.sin(t), (0.15 + 0.12 * label) * t
                ]
                trajectories[7, sample] = h
                trajectories[11, sample] = 0.8 * h + np.array([0.1, 0.2, 0.3])
                msrc_raw[subject].append(np.c_[trajectories[7, sample], trajectories[11, sample]])
            msrc[label - 1, subject] = trajectories
    x = rng.normal(size=(7, 31))
    means = rng.normal(size=(7, 3))
    cov = rng.uniform(0.2, 2, size=(7, 3))
    priors = np.array([[0.3, 0.3, 0.4]])
    savemat(
        work / "fixtures.mat",
        dict(
            raw_curve=raw,
            repeated_curve=repeated,
            ip_data=ip,
            msrc_data=msrc,
            fisher_input=x,
            fisher_means=means,
            fisher_cov=cov,
            fisher_priors=priors,
        ),
    )
    expression = (
        f"addpath('{quote(ROOT / 'validation')}');run_octave('{quote(ROOT)}','{quote(work)}');"
    )
    with (work / "octave.log").open("w") as log:
        process = subprocess.run(
            [str(args.octave), "--no-gui", "--quiet", "--eval", expression],
            env=env,
            stdout=log,
            stderr=log,
            timeout=600,
        )
    if process.returncode:
        raise RuntimeError(
            "Original-source Octave execution failed; see " + str(work / "octave.log")
        )
    print("Original-source Octave stages completed; checking against Python", flush=True)
    result = loadmat(work / "octave_results.mat", simplify_cells=True)
    errors = {}

    def compare(name, got, expected, tolerance=1e-9):
        got, expected = np.asarray(got), np.asarray(expected)
        assert got.shape == expected.shape, (name, got.shape, expected.shape)
        assert_allclose(got, expected, atol=tolerance, rtol=tolerance, err_msg=name)
        errors[name] = float(np.max(abs(got - expected)))

    resized = resample(raw, 64)
    compare("resampling", result["resized"], resized)
    compare("grid", result["grid_points"], grid(resized, "IP"), 0)
    compare("segments", result["segments"], np.asarray(split(grid(resized, "IP"), native)), 0)
    f, _ = frames(resized, "IP", native)
    compare("local_frames", result["local_frames"], f)
    compare(
        "descriptor",
        result["core_descriptor"],
        engine.call("RRVdescriptor_BasedonFrames", [f[2:-2], resized[2:-2]]),
    )
    compare("fisher", result["fisher_actual"], fisher_reference(x, means, cov, priors).ravel())
    compare(
        "temporal_pool",
        result["pool_actual"],
        pool_reference(x, means, cov, priors, [1, 2, 4, 8]).ravel(),
    )
    dataset_results = []
    for output in result["outputs"]:
        dataset = output["name"]
        if dataset == "IP":
            train_raw = [a for a, b in ip if int(b.ravel()[1]) <= 10]
            test_raw = [a for a, b in ip if int(b.ravel()[1]) > 10]
            expected_train_labels = [int(b.ravel()[0]) for a, b in ip if int(b.ravel()[1]) <= 10]
            expected_test_labels = [int(b.ravel()[0]) for a, b in ip if int(b.ravel()[1]) > 10]
        else:
            train_raw, test_raw = msrc_raw
            expected_train_labels, expected_test_labels = np.repeat([1, 2, 3], 3), np.repeat(
                [1, 2, 3], 2
            )
        assert_array_equal(output["train_labels"], expected_train_labels)
        assert_array_equal(output["test_labels"], expected_test_labels)
        for which, curves in [("train", train_raw), ("test", test_raw)]:
            actual = output["original_descriptors" if which == "train" else "test_descriptors"]
            for i, curve in enumerate(curves):
                reference, _ = describe(curve, dataset, native, engine)
                compare(f"{dataset}_{which}_descriptor_{i}", actual[i], reference)
        gmm = output["gmm"]
        params = gmm["means"], gmm["covariances"], gmm["priors"]
        expected = encode(output["original_descriptors"][0], params, engine)
        compare(dataset + "_encoded_first_train", output["first_train_code"], expected)
        for i, desc in enumerate(output["test_descriptors"]):
            compare(
                f"{dataset}_encoded_test_{i}",
                output["test_codes"][:, i],
                encode(desc, params, engine),
            )
        dataset_results.append(
            dict(
                dataset=dataset,
                unique_train=len(train_raw),
                test=len(test_raw),
                training_cells=len(train_raw),
                clusters=4,
                original_sampling_target=50000,
                prediction=np.asarray(output["prediction"]).tolist(),
                accuracy_percent=float(np.asarray(output["accuracy"]).ravel()[0]),
                source="real IP subset" if dataset == "IP" and real is not None else "synthetic",
            )
        )
    safe = subprocess.run(
        [str(native_work / "native_checks")], capture_output=True, text=True, timeout=30
    )
    contract = subprocess.run(
        [str(native_work / "native_checks"), "output-contract"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    (work / "native-safe.log").write_text(safe.stdout + safe.stderr)
    (work / "native-output-contract.log").write_text(contract.stdout + contract.stderr)
    assert safe.returncode == 0
    assert contract.returncode == 0, contract.stderr
    for name, digest in source_hashes.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    report = dict(
        octave_executed=True,
        matlab_executed=False,
        paper_reproduced=False,
        octave_findings=result["report"],
        numerical_comparisons=len(errors),
        max_errors=errors,
        pipelines=dataset_results,
        real_IP_indices_zero_based=indices,
        native_memory_checks={
            "valid_fixture": "passed ASan/UBSan",
            "circumcenter_outputs_1_to_3": "passed ASan/UBSan with exact-sized output storage",
            "gmm_outputs_1_to_5": "passed ASan/UBSan with exact-sized output storage",
        },
        original_source_sha256=source_hashes,
        validation_source_sha256=test_hashes,
        adapters=[
            "v5 MAT save and MATLAB filename extension compatibility",
        ],
        excluded=[
            "MATLAB ABI",
            "full benchmark",
            "IP regression prepass",
            "PCA and alternative baselines",
        ],
    )
    if args.ip_data:
        report["IP_data_sha256"] = hashlib.sha256(args.ip_data.read_bytes()).hexdigest()
    (work / "results.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=json_value) + "\n"
    )
    print(
        json.dumps(
            {
                "comparisons": len(errors),
                "max_error": max(errors.values()),
                "pipelines": dataset_results,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--octave", type=Path, required=True)
    parser.add_argument(
        "--conda-prefix", type=Path, help="Opt-in workaround for the IID macOS Conda installation"
    )
    parser.add_argument("--vlfeat-root", type=Path, required=True)
    parser.add_argument(
        "--native-work",
        type=Path,
        help="Reuse a fingerprint-verified native build from validation/run.py",
    )
    parser.add_argument("--ip-data", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run_octave(parser.parse_args())
