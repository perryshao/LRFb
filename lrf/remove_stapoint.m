function [regular_xyz, stapoint_index] = remove_stapoint(curve_xyz)
% REMOVE_STAPOINT  Remove NaN rows and consecutive repeated trajectory points.
%   NaN detection uses the first coordinate, as in the original pipeline.
%   stapoint_index records repeated points for fillstapoint.

NaN_index = isnan(curve_xyz);
curve_xyz(NaN_index(:, 1) == 1, :) = [];
double(curve_xyz);
samples = size(curve_xyz, 1);
regular_xyz = [];
regular_xyz(end + 1, :) = curve_xyz(1, :);
stapoint_index = [];
for i = 2:samples
    velocity = norm((curve_xyz(i, :) - curve_xyz(i - 1, :)));
    if velocity > 0
        regular_xyz(end + 1, :) = curve_xyz(i, :); %#ok<AGROW>
    else
        stapoint_index(end + 1) = i; %#ok<AGROW>
    end
end
