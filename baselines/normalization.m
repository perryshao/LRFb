function [Curve_xyz, max_limit, min_original] = normalization(Curve_xyz)
% NORMALIZATION  Translate and scale a trajectory by its largest coordinate range.
%   Returns normalized coordinates, the scale and the original minima.

n = size(Curve_xyz, 1);
m = size(Curve_xyz, 2);
Curve_limit = max(Curve_xyz, [], 1) - min(Curve_xyz, [], 1);
min_original = min(Curve_xyz, [], 1);
max_limit = max(Curve_limit);
Curve_xyz = (Curve_xyz - repmat(min(Curve_xyz, [], 1), n, 1)) ./ repmat(max_limit, n, m);
