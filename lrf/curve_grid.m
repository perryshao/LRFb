function [curvegrid, grid] = curve_grid(curve_xyz, width)
% CURVE_GRID  Quantize a trajectory using its largest coordinate range.
%   curve_xyz is N-by-D; width is the number of grid intervals.
%   grid is the scalar spacing; curvegrid contains rounded coordinates.

n = size(curve_xyz, 1);
curve_limit = max(curve_xyz, [], 1) - min(curve_xyz, [], 1);
max_limit = max(curve_limit);
grid = max_limit / width;
curvegrid = (curve_xyz - repmat(min(curve_xyz, [], 1), n, 1)) ./ grid;
curvegrid = round(curvegrid);
