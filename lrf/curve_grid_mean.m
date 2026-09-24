function [curvegrid, grid] = curve_grid_mean(curve_xyz, grids, mean_dist_global)
% CURVE_GRID_MEAN  Quantize coordinates with the historical mean-distance scaling.
%   Applies both curve_xyz / mean_dist_global and division by grid,
%   where grid = mean_dist_global / grids. This scaling is preserved.

n = size(curve_xyz, 1);
% curve_limit = max(curve_xyz,[],1)-min(curve_xyz,[],1);
% max_limit = max(curve_limit);
curve_xyz_n = curve_xyz / mean_dist_global;
grid = mean_dist_global / grids;
curvegrid = curve_xyz_n ./ grid;
curvegrid = round(curvegrid);
