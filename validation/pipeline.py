"""Explicit Python bridge for MATLAB orchestration, not a replacement MATLAB runtime.

Preserves legacy grid scaling, overlapping MBS traversal, and row partitioning.
Guard failures stop the bridge rather than silently repairing a degenerate input.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial.distance import cdist
from scipy.special import logsumexp


def matlab_round(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def resample(curve, length):
    curve = np.asarray(curve, dtype=float)
    if curve.shape[0] < curve.shape[1]:
        curve = curve.T
    n = len(curve)
    dense = CubicSpline(np.arange(1, n + 1), curve)(np.arange(10 * (n - 1) + 1) / 10 + 1)
    dense2 = CubicSpline(np.arange(1, len(dense) + 1), dense)(
        np.arange(10 * (len(dense) - 1) + 1) / 10 + 1
    )
    targets = np.arange(1, length + 1) * len(dense2) / length
    indices = np.abs(np.arange(1, len(dense2) + 1)[:, None] - targets).argmin(axis=0)
    return dense2[indices]


def grid(curve, dataset, mean_distance=None):
    if dataset == "IP":
        return matlab_round((curve - curve.min(axis=0)) / (np.ptp(curve, axis=0).max() / 1000))
    return matlab_round((curve / mean_distance) / (mean_distance / 1000))


def split(curve, native, width=5):
    """Line-by-line control-flow bridge of splitting_curve_3D; indices are 1-based."""
    n = len(curve)
    if n < 2:
        raise ValueError("Original splitter requires at least two points")
    k = 1

    def end_of(points):
        xy = native.segment(points[:, [0, 1]], width)
        xz = native.segment(points[:, [0, 2]], width)
        return min(xy or 1000, xz or 1000)

    while k < n:
        k += 1
        end = end_of(curve[:k])
        if end < k:
            break
    begin = 1
    if not 1 <= end <= n:
        raise IndexError(f"Original first segment endpoint {end} outside 1:{n}")
    segments = [(begin, end)]
    steps = 0
    while k < n:
        while k < n:
            begin += 1
            steps += 1
            if begin >= k or steps > n * n:
                raise RuntimeError("Original MBS traversal cannot progress with a valid segment")
            if end_of(curve[begin - 1 : k]) == k - begin + 1:
                break
        while k < n:
            k += 1
            if end_of(curve[begin - 1 : k]) < k - begin + 1:
                break
        segments.append((begin, k - 1))
    return segments


def frames(curve, dataset, native, mean_distance=None):
    # MATLAB removes rows with NaN in the first coordinate only.
    clean = curve[~np.isnan(curve[:, 0])]
    repeats = np.flatnonzero(np.linalg.norm(np.diff(clean, axis=0), axis=1) == 0) + 1
    clean = np.delete(clean, repeats, axis=0)
    quantized = grid(clean, dataset, mean_distance)
    segments = split(quantized, native)
    n = len(clean)
    left, right = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    for i, (b, e) in enumerate(segments):
        first = 1 if i == 0 else segments[i - 1][1] + 1
        left[first - 1 : e] = b
        last = n if i == len(segments) - 1 else segments[i + 1][0] - 1
        right[b - 1 : last] = e
    output = np.zeros((n, 9))
    for i in range(2, n - 2):
        if not (1 <= left[i] <= n and 1 <= right[i] <= n):
            raise IndexError("Unassigned MBS key point")
        l, p, r = clean[left[i] - 1], clean[i], clean[right[i] - 1]
        t1, t2 = p - l, r - p
        t1, t2 = t1 / np.linalg.norm(t1), t2 / np.linalg.norm(t2)
        center = native.circum(np.array([l, p, r])) + l
        normal = center - p
        normal /= np.linalg.norm(normal)
        binormal = np.cross(t1, t2)
        binormal /= np.linalg.norm(binormal)
        output[i] = np.r_[np.cross(normal, binormal), normal, binormal]
    for i in repeats:
        output = np.insert(output, i, output[i - 1], axis=0)
    return output, segments


def describe(trajectory, dataset, native, engine):
    count = 64 if dataset == "IP" else 128
    hands = [resample(trajectory[:, i : i + 3], count) for i in (0, 3)]
    scale = max(cdist(hand, hand).mean() for hand in hands)
    descriptions = []
    diagnostics = []
    for hand in hands:
        hand = hand[~np.isnan(hand[:, 0])]
        # GeneSC explicitly emits seven zero columns for an absent hand.
        if not np.any(hand):
            descriptions.append(np.zeros((len(hand) - 4, 7)))
            diagnostics.append(0.0)
            continue
        f, segments = frames(hand, dataset, native, scale)
        f = f[2:-2]
        bases = f.reshape(-1, 3, 3).transpose(0, 2, 1)
        diagnostics.append(float(np.max(np.abs(bases.transpose(0, 2, 1) @ bases - np.eye(3)))))
        description = engine.call("RRVdescriptor_BasedonFrames", [f, hand[2:-2]])
        if not np.isfinite(description).all():
            raise ValueError("Nonfinite original-source descriptor")
        descriptions.append(description)
    return np.hstack(descriptions), max(diagnostics)


def fisher_reference(x, means, cov, priors):
    """Independent vectorized diagonal-GMM score formula with VLFeat thresholds."""
    d, n = x.shape
    if not n:
        return np.zeros((2 * means.size, 1))
    p = np.asarray(priors).ravel()
    z = (x.T[:, None, :] - means.T[None, :, :]) / np.sqrt(cov.T)[None, :, :]
    logq = np.log(p)[None, :] - 0.5 * (
        d * np.log(2 * np.pi) + np.log(cov).sum(axis=0)[None, :] + (z * z).sum(axis=2)
    )
    logq[:, p < 1e-6] = -np.inf
    q = np.exp(logq - logsumexp(logq, axis=1, keepdims=True))
    q[(q < 1e-6) | (p[None, :] < 1e-6)] = 0
    u = (q[:, :, None] * z).sum(axis=0) / (n * np.sqrt(p)[:, None])
    v = (q[:, :, None] * (z * z - 1)).sum(axis=0) / (n * np.sqrt(2 * p)[:, None])
    code = np.r_[u.ravel(), v.ravel()]
    code = np.sign(code) * np.sqrt(np.abs(code))
    code /= max(np.linalg.norm(code), 1e-12)
    return code[:, None]


def pool_reference(x, means, cov, priors, pyramid):
    blocks = []
    for bins in pyramid:
        ids = np.ceil(np.arange(1, x.shape[1] + 1) / (x.shape[1] / bins))
        for b in range(1, bins + 1):
            chunk = x[:, ids == b]
            chunk = chunk[:, ~np.isnan(chunk[-1])]
            blocks.append(fisher_reference(chunk, means, cov, priors).ravel())
    code = np.concatenate(blocks)
    norm = np.linalg.norm(code)
    return (code / norm if norm > 0 else code)[:, None]


def encode(descriptor, parameters, engine, partitions=2):
    means, covariance, priors = parameters
    frames_per_part = len(descriptor) // partitions
    codes = []
    for part in range(partitions):
        x = descriptor[part * frames_per_part : (part + 1) * frames_per_part].T
        codes.append(
            engine.call(
                "fv_pooling_ts",
                [x, means, covariance, priors, "Improved", np.array([[1, 2, 4, 8]])],
            )
        )
    return np.concatenate(codes).ravel()
