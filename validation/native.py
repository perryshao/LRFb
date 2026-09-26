"""Build original numerical sources locally; no existing MEX binaries are replaced."""
from pathlib import Path
import ctypes as c
import json
import hashlib
import platform
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def build_fingerprint(vlfeat):
    paths = list((ROOT / "mbs/src").glob("*"))
    paths += list((ROOT / "thirdparty/libsvm-3.17").glob("*.cpp"))
    paths += list((ROOT / "thirdparty/libsvm-3.17").glob("*.h"))
    paths += [
        HERE / name
        for name in ("native.py", "native_bridge.cpp", "mex.h", "matrix.h", "vlfeat-arm64-cpuid.h")
    ]
    paths += list((vlfeat / "vl").glob("*.c")) + list((vlfeat / "vl").glob("*.h"))
    return {
        str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()
    }


def reuse(output, vlfeat):
    recorded = json.loads((output / "build-fingerprint.json").read_text())
    if recorded != build_fingerprint(vlfeat):
        raise ValueError("Native sources changed: rerun without --reuse-build")
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    return Native(output / ("liblrfb" + suffix), output / ("libvl" + suffix))


def build(output, vlfeat):
    output.mkdir(parents=True, exist_ok=True)
    wrappers = []
    for stem, symbol in [
        ("Determine_segment", "segment_gateway"),
        ("tricircumcenter3d", "circum_gateway"),
    ]:
        target = output / (stem + "_wrapper.cpp")
        target.write_text(
            "#define mexFunction "
            + symbol
            + "\n#include "
            + json.dumps(str(ROOT / "mbs/src" / (stem + ".cpp")))
            + "\n"
        )
        wrappers.append(str(target))
    sources = wrappers + [
        str(HERE / "native_bridge.cpp"),
        str(ROOT / "thirdparty/libsvm-3.17/svm.cpp"),
    ]
    common = ["clang++", "-std=c++11", "-I", str(HERE), "-I", str(ROOT / "thirdparty/libsvm-3.17")]
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    native = output / ("liblrfb" + suffix)
    vl = output / ("libvl" + suffix)
    commands = [
        common + ["-O1", "-shared", "-fPIC", *sources, "-o", str(native)],
        common
        + [
            "-O1",
            "-g",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            "-DVALIDATE_MAIN",
            *sources,
            "-o",
            str(output / "native_checks"),
        ],
    ]
    flags = ["-DVL_DISABLE_SSE2", "-DVL_DISABLE_AVX", "-DVL_DISABLE_THREADS"]
    if platform.machine().lower() in ("arm64", "aarch64"):
        flags += ["-include", str(HERE / "vlfeat-arm64-cpuid.h")]
    commands.append(
        [
            "clang",
            "-shared",
            "-fPIC",
            "-O2",
            *flags,
            "-I",
            str(vlfeat),
            *map(str, sorted((vlfeat / "vl").glob("*.c"))),
            "-lm",
            "-o",
            str(vl),
        ]
    )
    (output / "build-commands.json").write_text(json.dumps(commands, indent=2))
    with (output / "build.log").open("w") as log:
        for command in commands:
            subprocess.run(command, stdout=log, stderr=log, check=True, timeout=180)
    (output / "build-fingerprint.json").write_text(json.dumps(build_fingerprint(vlfeat), indent=2))
    return Native(native, vl)


class Native:
    def __init__(self, native, vl):
        self.core = c.CDLL(str(native))
        self.vl = c.CDLL(str(vl))
        self.core.segment.argtypes = [c.c_void_p, c.c_int, c.c_double]
        self.core.segment.restype = c.c_int
        self.core.circum.argtypes = [c.c_void_p, c.c_void_p]
        self.core.circum.restype = None
        self.core.classify.argtypes = [
            c.c_void_p,
            c.c_void_p,
            c.c_int,
            c.c_int,
            c.c_void_p,
            c.c_int,
            c.c_void_p,
            c.c_void_p,
        ]
        self.core.classify.restype = c.c_int
        self.vl.vl_fisher_encode.argtypes = [
            c.c_void_p,
            c.c_uint,
            c.c_void_p,
            c.c_ulonglong,
            c.c_ulonglong,
            c.c_void_p,
            c.c_void_p,
            c.c_void_p,
            c.c_ulonglong,
            c.c_int,
        ]
        self.vl.vl_fisher_encode.restype = c.c_ulonglong
        self.vl.vl_gmm_new.argtypes = [c.c_uint, c.c_ulonglong, c.c_ulonglong]
        self.vl.vl_gmm_new.restype = c.c_void_p
        self.vl.vl_gmm_delete.argtypes = [c.c_void_p]
        self.vl.vl_gmm_delete.restype = None
        self.vl.vl_gmm_cluster.argtypes = [c.c_void_p, c.c_void_p, c.c_ulonglong]
        self.vl.vl_gmm_cluster.restype = c.c_double
        self.vl.vl_gmm_set_max_num_iterations.argtypes = [c.c_void_p, c.c_ulonglong]
        self.vl.vl_gmm_set_max_num_iterations.restype = None
        self.vl.vl_get_rand.restype = c.c_void_p
        self.vl.vl_rand_seed.argtypes = [c.c_void_p, c.c_uint]
        self.vl.vl_rand_seed.restype = None
        for name in ("means", "covariances", "priors"):
            fn = getattr(self.vl, "vl_gmm_get_" + name)
            fn.argtypes = [c.c_void_p]
            fn.restype = c.POINTER(c.c_double)

    def segment(self, points, width=5):
        p = np.asfortranarray(points, dtype=float)
        if p.ndim != 2 or p.shape[1] != 2 or not 2 <= len(p) < 999 or not np.isfinite(p).all():
            raise ValueError("Validation bridge requires finite N-by-2 input, 2 <= N < 999")
        return self.core.segment(p.ctypes.data, len(p), width)

    def circum(self, points):
        p = np.ascontiguousarray(points, dtype=float)
        if p.shape != (3, 3) or not np.isfinite(p).all():
            raise ValueError("Expected three finite 3-D points")
        out = np.empty(3)
        self.core.circum(p.ctypes.data, out.ctypes.data)
        return out

    def fisher(self, x, means, cov, priors, mode="Improved"):
        if mode != "Improved":
            raise NotImplementedError(mode)
        x, means, cov, priors = [np.asfortranarray(a, dtype=float) for a in (x, means, cov, priors)]
        d, n = x.shape
        k = means.shape[1]
        assert means.shape == cov.shape == (d, k) and priors.size == k
        assert all(np.isfinite(a).all() for a in (x, means, cov, priors))
        assert np.all(cov > 0) and np.all(priors > 0)
        out = np.zeros(2 * d * k)
        self.vl.vl_fisher_encode(
            out.ctypes.data,
            2,
            means.ctypes.data,
            d,
            k,
            cov.ctypes.data,
            priors.ctypes.data,
            x.ctypes.data,
            n,
            3,
        )
        return out.reshape(-1, 1)

    def gmm(self, x, k):
        x = np.asfortranarray(x, dtype=float)
        assert np.isfinite(x).all() and x.shape[1] >= k
        self.vl.vl_rand_seed(self.vl.vl_get_rand(), 240925)
        model = self.vl.vl_gmm_new(2, x.shape[0], k)
        try:
            self.vl.vl_gmm_set_max_num_iterations(model, 100)
            likelihood = self.vl.vl_gmm_cluster(model, x.ctypes.data, x.shape[1])
            params = []
            for name, size in [
                ("means", x.shape[0] * k),
                ("covariances", x.shape[0] * k),
                ("priors", k),
            ]:
                arr = np.ctypeslib.as_array(
                    getattr(self.vl, "vl_gmm_get_" + name)(model), shape=(size,)
                ).copy()
                params.append(arr if name == "priors" else arr.reshape(x.shape[0], k, order="F"))
            return (*params, likelihood)
        finally:
            self.vl.vl_gmm_delete(model)

    def classify(self, train, labels, test):
        train, labels, test = [np.ascontiguousarray(a, dtype=float) for a in (train, labels, test)]
        assert train.ndim == test.ndim == 2 and train.shape[1] == test.shape[1]
        assert len(labels) == len(train)
        assert all(np.isfinite(a).all() for a in (train, labels, test))
        pred, sums = np.empty(len(test)), np.empty(len(test))
        classes = self.core.classify(
            train.ctypes.data,
            labels.ctypes.data,
            *train.shape,
            test.ctypes.data,
            len(test),
            pred.ctypes.data,
            sums.ctypes.data
        )
        assert classes == len(np.unique(labels))
        return pred, sums
