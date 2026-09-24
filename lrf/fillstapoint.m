function smooth_curve = fillstapoint(smooth_curve, staindex)
% FILLSTAPOINT  Reinsert repeated points at the recorded trajectory indices.
%   Each inserted row copies its predecessor; indices refer to the input
%   trajectory before repeated points were removed.

n = length(staindex);
for i = 1:n
    smooth_curve = [smooth_curve(1:staindex(i) - 1, :); smooth_curve(staindex(i) - 1, ...
        :); smooth_curve(staindex(i):end, :)]; % insert operation
end
