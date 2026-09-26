"""Run source-derived checks and explicitly bridged LRFb pipelines without MATLAB."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.io import loadmat, savemat
from scipy.spatial.transform import Rotation
from scipy.spatial import ConvexHull

from matlab_subset import Engine, scalar
from native import build, reuse, ROOT
from pipeline import describe, encode, fisher_reference, matlab_round, pool_reference, resample


def make_engine(native):
    e = Engine()
    for relative in [
        "descriptor/Myrotm2quat.m",
        "descriptor/RRVdescriptor_BasedonFrames.m",
        "encoding/fv_pooling_ts.m",
        "encoding/rand_sampling_ts.m",
    ]:
        e.load(ROOT / relative)
    e.builtins.update(
        vl_fisher=native.fisher,
        round=matlab_round,
        max=lambda x: np.max(x),
        cat=lambda d, *x: np.concatenate(x, axis=int(d) - 1),
    )
    return e


def cells(items):
    result = np.empty((1, len(items)), dtype=object)
    for i, value in enumerate(items):
        result[0, i] = value
    return result


def digital_width(points):
    """Independent hull oracle with the recognizer's monotone-axis constraint."""
    delta = np.diff(points, axis=0)
    monotone = [
        (np.all(delta[:, axis] >= 0) or np.all(delta[:, axis] <= 0)) and np.ptp(points[:, axis]) > 0
        for axis in (0, 1)
    ]
    if not any(monotone):
        return float("inf")
    points = np.unique(points, axis=0)
    if len(points) < 3 or np.linalg.matrix_rank(points - points[0]) < 2:
        return 0.0
    hull = points[ConvexHull(points).vertices]
    edges = np.roll(hull, -1, axis=0) - hull
    widths = []
    for dx, dy in edges:
        if dx and monotone[0]:
            widths.append(np.ptp(points[:, 1] - (dy / dx) * points[:, 0]))
        if dy and monotone[1]:
            widths.append(np.ptp(points[:, 0] - (dx / dy) * points[:, 1]))
    return min(widths)


def execute(args):
    output = args.output.resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError("Validation output must be outside the source repository")
    output.mkdir(parents=True, exist_ok=True)
    args.vlfeat_root = args.vlfeat_root.resolve()
    native = (
        reuse(output, args.vlfeat_root) if args.reuse_build else build(output, args.vlfeat_root)
    )
    engine = make_engine(native)
    results = []
    rng = np.random.default_rng(240925)

    def record(name, status="passed", **data):
        results.append(dict(check=name, status=status, **data))
        print(name, status, json.dumps(data, ensure_ascii=False)[:500], flush=True)

    if args.real_probe is not None:
        data = loadmat(args.ip_data)["Data"]
        idx = args.real_probe
        desc, orthogonality = describe(data[idx, 0], "IP", native, engine)
        np.savez(output / ("real-" + str(idx) + ".npz"), descriptor=desc, ids=data[idx, 1])
        print(json.dumps(dict(index=idx, shape=list(desc.shape), orthogonality=orthogonality)))
        return 0

    # Fingerprint the actual sources, including external VLFeat; no stale result claims.
    sources = [
        p
        for p in ROOT.rglob("*")
        if ".git" not in p.parts and p.suffix in (".m", ".cpp", ".h", ".c")
    ]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    vlhash = {
        str(p.relative_to(args.vlfeat_root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((args.vlfeat_root / "vl").glob("*"))
        if p.suffix in (".c", ".h")
    }
    (output / "source-hashes.json").write_text(
        json.dumps({"repository": hashes, "external_vlfeat": vlhash}, indent=2)
    )

    with tempfile.TemporaryDirectory(dir=output) as tmp:
        probe = Path(tmp) / "probe.m"
        probe.write_text("function [a,b,c] = probe(x,y)\na=x(:);\nb=x/y;\nc=x(2:end,[1 3]);\nend\n")
        evaluator = Engine()
        evaluator.load(probe)
        x = np.arange(1.0, 10.0).reshape(3, 3)
        y = np.diag([2.0, 3.0, 4.0])
        a, b, c = evaluator.call("probe", [x, y], 3)
        assert_array_equal(a.ravel(), x.ravel(order="F"))
        assert_allclose(b @ y, x)
        assert_array_equal(c, x[1:, [0, 2]])
        probe.write_text("function a = probe(x)\na=unknown_function(x);\nend\n")
        evaluator.load(probe)
        try:
            evaluator.call("probe", [x])
        except KeyError:
            pass
        else:
            raise AssertionError("Unsupported function was silently accepted")
    record(
        "AST_evaluator_semantics",
        checks=[
            "column-major reshape",
            "one-based indexing",
            "end",
            "matrix right division",
            "unsupported function fails",
        ],
    )

    safe = subprocess.run(
        [str(output / "native_checks")], capture_output=True, text=True, timeout=30
    )
    (output / "native-safe.log").write_text(safe.stdout + safe.stderr)
    assert safe.returncode == 0, safe.stderr
    record("native_ASan_UBSan_valid_fixture")
    contract = subprocess.run(
        [str(output / "native_checks"), "output-contract"], capture_output=True, text=True, timeout=30
    )
    (output / "output-contract-asan.log").write_text(contract.stdout + contract.stderr)
    assert contract.returncode == 0, contract.stderr
    record("circumcenter_output_contract", outputs=[1, 2, 3], sanitizer="ASan/UBSan")

    errors = []
    for _ in range(100):
        tri = rng.normal(size=(3, 3))
        a, b, c = tri
        u, v = b - a, c - a
        normal = np.cross(u, v)
        expected = np.linalg.solve(np.vstack([u, v, normal]), [u @ u / 2, v @ v / 2, 0])
        errors.append(float(np.max(np.abs(native.circum(tri) - expected))))
    assert max(errors) < 1e-10
    record(
        "native_circumcenter_vs_linear_system",
        cases=100,
        max_absolute_error=max(errors),
        outputs_requested=1,
    )

    matches = 0
    for _ in range(30):
        p = np.cumsum(rng.integers(-6, 7, size=(30, 2)), axis=0).astype(float)
        expected = len(p)
        for i in range(2, len(p) + 1):
            if digital_width(p[:i]) >= 5:
                expected = i - 1
                break
        got = native.segment(p)
        assert got == expected, (got, expected, p.tolist())
        matches += 1
    record("native_MBS_vs_independent_convex_hull_width", cases=matches)

    # Test the MATLAB descriptor itself through parsed source, not a hand translation.
    maxerr = 0
    for n in (5, 19, 60):
        rotations = Rotation.from_rotvec(
            np.c_[np.zeros(n), np.zeros(n), np.linspace(0, 0.8, n)]
        ).as_matrix()
        f = rotations.transpose(0, 2, 1).reshape(n, 9)
        t = np.linspace(0, 1, n)
        tr = np.c_[t, t * t, np.sin(t)]
        got = engine.call("RRVdescriptor_BasedonFrames", [f, tr])
        velocity = np.gradient(tr, axis=0)
        velocity /= np.sqrt(np.linalg.norm(velocity, axis=1, keepdims=True) + np.finfo(float).eps)
        expected = []
        for i in range(n):
            relative = rotations[min(i + 1, n - 1)] @ rotations[i].T
            q = Rotation.from_matrix(relative.T).as_quat()[[3, 0, 1, 2]]
            if q[0] < 0:
                q *= -1
            expected.append(np.r_[q, velocity[i] @ relative])
        maxerr = max(maxerr, float(np.max(np.abs(got - expected))))
        assert_allclose(got, expected, atol=2e-13)
    record("MATLAB_AST_descriptor_vs_Scipy_rotation", cases=3, max_absolute_error=maxerr)

    errors = []
    for d in (1, 7, 14):
        for k in (1, 3):
            for n in (0, 1, 11):
                x = rng.normal(size=(d, n))
                m = rng.normal(size=(d, k))
                cov = rng.uniform(0.2, 2, size=(d, k))
                p = rng.uniform(0.1, 1, k)
                p /= p.sum()
                got = native.fisher(x, m, cov, p)
                expected = fisher_reference(x, m, cov, p)
                assert_allclose(got, expected, atol=2e-12, rtol=2e-10)
                errors.append(float(np.max(np.abs(got - expected))))
    record(
        "native_VLFeat_Fisher_vs_independent_score_formula",
        cases=len(errors),
        max_absolute_error=max(errors),
    )

    m = rng.normal(size=(7, 3))
    cov = rng.uniform(0.2, 2, size=m.shape)
    p = np.array([0.3, 0.3, 0.4])
    for n in (1, 5, 31, 60):
        x = rng.normal(size=(7, n))
        got = engine.call("fv_pooling_ts", [x, m, cov, p, "Improved", np.array([[1, 2, 4, 8]])])
        expected = pool_reference(x, m, cov, p, [1, 2, 4, 8])
        assert_allclose(got, expected, atol=2e-12)
    record("MATLAB_AST_temporal_pool_vs_independent_bins", cases=4, includes_empty_bins=True)

    rotation_cases = [np.eye(3), *[np.diag(row) for row in ([1, -1, -1], [-1, 1, -1], [-1, -1, 1])]]
    for axis in ([1, 2, 3], [-3, 1, 2], [1, -4, 2]):
        axis = np.asarray(axis) / np.linalg.norm(axis)
        rotation_cases.extend(Rotation.from_rotvec(axis * angle).as_matrix() for angle in
                              (np.pi - 1e-8, np.pi, np.pi + 1e-8))
    rotation_cases.extend(Rotation.random(100, random_state=rng).as_matrix())
    for rotation in rotation_cases:
        q = np.asarray(engine.call("Myrotm2quat", [rotation])).ravel()
        assert np.isfinite(q).all() and q[0] >= 0
        assert_allclose(np.linalg.norm(q), 1, atol=2e-14)
        assert_allclose(Rotation.from_quat(q[[1, 2, 3, 0]]).as_matrix(), rotation, atol=2e-13)
    for x in (np.zeros((2, 0)), np.full((2, 3), np.nan), np.zeros((2, 4))):
        z = engine.call("fv_pooling_ts", [x, np.zeros((2, 1)), np.ones((2, 1)),
                                          np.ones(1), "Improved", np.array([[1, 2]])])
        assert z.shape == (12, 1) and np.isfinite(z).all()
        if x.shape[1] == 0 or np.isnan(x).all():
            assert_array_equal(z, np.zeros((12, 1)))
    record("rotation_and_empty_pool_regressions", rotations=len(rotation_cases),
           pool_cases=["empty", "fully occluded", "zero descriptors"])

    labels_engine = Engine()
    labels_engine.load(ROOT / "experiments/IP/getLabels.m")
    train_cells = np.empty((2, 3), dtype=object)
    test_cells = np.empty((2, 2), dtype=object)
    for arr, classes in [(train_cells, [2, 1, 3]), (test_cells, [3, 2])]:
        for i, label in enumerate(classes):
            arr[0, i] = np.array([[label, i + 1]])
            arr[1, i] = np.zeros((64, 6))
    labels_engine.builtins["load_workspace"] = (
        lambda args: {"TRAJDB": train_cells} if args[0] == "DB" else {"TRAJSAMPLES": test_cells}
    )
    train_y, test_y = labels_engine.call("getLabels", [], 2)
    assert_array_equal(train_y.ravel(), [2, 1, 3])
    assert_array_equal(test_y.ravel(), [3, 2])
    record("MATLAB_AST_getLabels", training=3, testing=2, label_order="preserved")

    # Deterministic permutations isolate index/count semantics from RNG differences.
    samples = cells([np.arange(18.0).reshape(6, 3), np.arange(18.0, 36.0).reshape(6, 3)])
    engine.builtins["randperm"] = lambda n: (np.arange(int(n)) + 1).reshape(1, -1)
    for requested in (0, 2, 12, 16, 24, 26, 100, 50000):
        sample = engine.call("rand_sampling_ts", [samples, requested])
        per_sequence = requested // 2
        expected = np.hstack([np.asarray(seq)[np.arange(per_sequence) % 6].T for seq in samples.ravel()])
        assert_array_equal(sample, expected)
    for invalid, requested in [(cells([]), 1), (cells([np.empty((0, 3))]), 1),
                               (samples, -1), (samples, 1.5), (samples, float("nan"))]:
        try:
            engine.call("rand_sampling_ts", [invalid, requested])
        except ValueError as exc:
            assert "LRFb:" in str(exc)
        else:
            raise AssertionError("Invalid sampler input was accepted")
    record("MATLAB_AST_sampler_regression", requests=[0, 2, 12, 16, 24, 26, 100, 50000],
           invalid_cases=5, legacy_order_preserved=True)

    gradient_errors = []
    objective_errors = []
    for dataset in ("IP", "MSRC12"):
        cost = Engine()
        cost.load(ROOT / f"experiments/{dataset}/costFuncRegMultPartGp_v1_42.m")
        for groups, classes, lambdas in [([1], 1, [0, 1]), ([2, 3], 2, [0.4, 0.7]),
                                          ([1, 2, 3], 3, [0.3, 0.8])]:
            dimensions = sum(groups)
            theta = rng.normal(size=(dimensions, classes))
            x = rng.normal(size=(dimensions, 5))
            y = rng.normal(size=(5, classes))
            if dimensions == 1:
                theta[:] = 2; x[:] = 0; y[:] = 0
            cases = [theta, np.zeros_like(theta), theta.copy()]
            cases[-1][:groups[0], 0] = 0
            for values in cases:
                flat = values.reshape(-1, 1, order="F")
                cost_args = [flat, x, y, np.array([lambdas]), classes, np.array([groups])]
                f, g = cost.call("costFuncRegMultPartGp_v1_42", cost_args, 2)
                cuts = np.cumsum([0] + groups)
                objective = np.sum((x.T @ values - y)**2) + lambdas[0] * np.linalg.norm(values, axis=1).sum()
                objective += lambdas[1] * sum(np.sqrt(np.sum(values[a:b]**4, axis=0)).sum()
                                              for a, b in zip(cuts[:-1], cuts[1:]))
                assert_allclose(scalar(f), objective, atol=2e-12)
                objective_errors.append(abs(scalar(f) - objective))
                fd = np.zeros_like(flat)
                for index in range(flat.size):
                    plus, minus = flat.copy(), flat.copy()
                    plus[index] += 1e-5; minus[index] -= 1e-5
                    fd[index] = (scalar(cost.call("costFuncRegMultPartGp_v1_42", [plus] + cost_args[1:])) -
                                 scalar(cost.call("costFuncRegMultPartGp_v1_42", [minus] + cost_args[1:]))) / 2e-5
                assert_allclose(g, fd, atol=2e-7, rtol=2e-7)
                gradient_errors.append(float(np.max(np.abs(g - fd))))
    record("MATLAB_AST_regression_gradient_finite_difference", datasets=["IP", "MSRC12"],
           cases=len(gradient_errors), max_gradient_error=max(gradient_errors),
           max_objective_error=max(objective_errors))

    # Resampling checked on a polynomial whose cubic spline has an analytic value.
    t = np.arange(1.0, 8.0)
    curve = np.c_[t, t * t, t**3]
    got = resample(curve, 64)
    dense_len = 100 * (len(t) - 1) + 1
    ids = np.abs(np.arange(1, dense_len + 1)[:, None] - np.arange(1, 65) * dense_len / 64).argmin(
        axis=0
    )
    positions = 1 + ids / 100
    assert_allclose(got, np.c_[positions, positions**2, positions**3], atol=2e-12)
    assert_array_equal(matlab_round(np.array([-0.5, 0.5, -1.5, 1.5])), [-1, 1, -2, 2])
    record("resampling_bridge_analytic_cubic", length=64, rounding="half away from zero")

    # Complete mathematical route on nondegenerate synthetic two-hand trajectories.
    for dataset in ("IP", "MSRC12"):
        descriptions, labels, orthogonality = [], [], []
        for label in range(1, 4):
            for sample_id in range(8):
                t = np.linspace(0, 4 + 0.7 * label, 69 + sample_id)
                h = np.c_[np.cos(t), np.sin(t), (0.15 + 0.12 * label) * t]
                h *= 1 + 0.008 * sample_id
                trajectory = np.c_[h, h * (0.65 + 0.06 * label) + np.array([0.1, 0.2, 0.3])]
                desc, err = describe(trajectory, dataset, native, engine)
                descriptions.append(desc)
                labels.append(label)
                orthogonality.append(err)
        train_indices = [i for i in range(24) if i % 8 < 6]
        test_indices = [i for i in range(24) if i % 8 >= 6]
        x = np.concatenate([descriptions[i].T for i in train_indices], axis=1)
        means, cov, priors, ll = native.gmm(x, 64)
        assert np.isfinite(ll) and np.all(cov > 0)
        assert_allclose(priors.sum(), 1, atol=1e-12)
        encoded = np.array([encode(desc, (means, cov, priors), engine) for desc in descriptions])
        assert encoded.shape == (24, 53760)
        assert_allclose(np.linalg.norm(encoded, axis=1), np.sqrt(2), atol=2e-12)
        labels = np.array(labels)
        pred, sums = native.classify(
            encoded[train_indices], labels[train_indices], encoded[test_indices]
        )
        assert np.isin(pred, labels).all()
        assert_allclose(sums, 1, atol=2e-12)
        confusion = np.zeros((3, 3), dtype=int)
        for actual, prediction in zip(labels[test_indices], pred):
            confusion[actual - 1, int(prediction) - 1] += 1
        assert confusion.sum() == len(test_indices)
        record(
            dataset + "_synthetic_pipeline",
            train=18,
            test=6,
            clusters=64,
            descriptor_shape=list(descriptions[0].shape),
            fisher_dimension=53760,
            paper_dimension=26880,
            frame_orthogonality_error=max(orthogonality),
            confusion=confusion.tolist(),
            accuracy=float(np.mean(pred == labels[test_indices])),
            caveat="Synthetic smoke only; Python orchestration, all training descriptors instead of 50000 random samples; no MATLAB I/O or IP regression prepass",
        )

    if args.ip_data:
        data = loadmat(args.ip_data)["Data"]
        ids = np.vstack([v.ravel() for v in data[:, 1]])
        keep = ~np.isin(ids[:, 0], [17, 18])
        training = keep & np.isin(ids[:, 1], np.arange(1, 11))
        testing = keep & ~training
        record(
            "real_IP_loader_protocol",
            data_shape=list(data.shape),
            training=int(training.sum()),
            testing=int(testing.sum()),
            fixed_8000_allocation_safe_for_this_file=bool(
                training.sum() >= 8000 and testing.sum() >= 8000
            ),
            sha256=hashlib.sha256(args.ip_data.read_bytes()).hexdigest(),
        )
        selected = []
        for label in (1, 2, 3):
            selected.extend(np.flatnonzero(training & (ids[:, 0] == label))[::100][:3].tolist())
            selected.extend(np.flatnonzero(testing & (ids[:, 0] == label))[:2].tolist())
        valid, failures = [], []
        for idx in selected:
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--output",
                str(output),
                "--vlfeat-root",
                str(args.vlfeat_root),
                "--reuse-build",
                "--ip-data",
                str(args.ip_data),
                "--real-probe",
                str(idx),
            ]
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=60)
            except subprocess.TimeoutExpired:
                failures.append({"index": idx, "reason": "timeout"})
                continue
            (output / ("real-" + str(idx) + ".log")).write_text(result.stdout + result.stderr)
            if result.returncode:
                failures.append(
                    {"index": idx, "returncode": result.returncode, "tail": result.stderr[-500:]}
                )
            else:
                valid.append(idx)
        record(
            "real_IP_descriptor_probe",
            "observed_failures" if failures else "passed",
            tested=len(selected),
            finite=len(valid),
            failures=failures,
        )
        train_valid = [i for i in valid if training[i]]
        test_valid = [i for i in valid if testing[i]]
        if (
            len(train_valid) >= 6
            and len(test_valid) >= 3
            and len(np.unique(ids[train_valid, 0])) >= 2
        ):
            descriptors = {
                i: np.load(output / ("real-" + str(i) + ".npz"))["descriptor"] for i in valid
            }
            x = np.concatenate([descriptors[i].T for i in train_valid], axis=1)
            m, cov, p, _ = native.gmm(x, 4)
            train = np.array([encode(descriptors[i], (m, cov, p), engine) for i in train_valid])
            test = np.array([encode(descriptors[i], (m, cov, p), engine) for i in test_valid])
            pred, sums = native.classify(train, ids[train_valid, 0], test)
            assert_allclose(sums, 1, atol=2e-12)
            record(
                "real_IP_small_bridged_pipeline",
                train=len(train_valid),
                test=len(test_valid),
                clusters=4,
                accuracy=float(np.mean(pred == ids[test_valid, 0])),
                caveat="Selected finite subset only; not paper reproduction",
            )
        else:
            record(
                "real_IP_small_bridged_pipeline",
                "blocked",
                reason="Insufficient finite descriptors after original-source geometry checks",
            )

    # Verify v5 MAT shape/labels across the Python bridge only.
    fixture = output / "io-smoke.mat"
    savemat(
        fixture, {"descriptors": np.arange(28.0).reshape(2, 14), "labels": np.array([[1], [2]])}
    )
    recovered = loadmat(fixture)
    assert_array_equal(recovered["labels"], [[1], [2]])
    assert recovered["descriptors"].shape == (2, 14)
    fixture.unlink()
    record(
        "Python_MAT_v5_round_trip",
        caveat="Does not validate MATLAB v7.3 save/load or cell-array workspace side effects",
    )
    report = {
        "schema_version": 1,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "matlab_executed": False,
        "paper_reproduced": False,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "results": results,
        "AST_source_functions": engine.calls,
        "scope": "Original MATLAB subset execution + original native C/C++ + explicit Python orchestration; fixed source with regression assertions; original binaries not replaced",
    }
    (output / "results.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print("Results:", output / "results.json")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vlfeat-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ip-data", type=Path)
    parser.add_argument(
        "--reuse-build",
        action="store_true",
        help="Reuse binaries only when their original sources are unchanged",
    )
    parser.add_argument("--real-probe", type=int, help=argparse.SUPPRESS)
    raise SystemExit(execute(parser.parse_args()))
