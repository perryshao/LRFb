function ResizeTrj = TrjResizeTime(Trj, Length)
% TRJRESIZETIME  Resample an N-by-3 trajectory in the time-index domain.
%   Two spline interpolation passes precede selection of Length rows.
%   This does not resample by arc length. Original author: skaegy.

if nargin < 1
    error('Need more input');
end

if size(Trj, 1) < size(Trj, 2)
    Trj = Trj';
end
n = length(Trj);
Interp_Trj1(:, 1) = interp1(1:n, Trj(:, 1), 1:0.1:n, 'spline');
Interp_Trj1(:, 2) = interp1(1:n, Trj(:, 2), 1:0.1:n, 'spline');
Interp_Trj1(:, 3) = interp1(1:n, Trj(:, 3), 1:0.1:n, 'spline');
Interp_Trj2(:, 1) = interp1(1:length(Interp_Trj1), Interp_Trj1(:, 1), 1:0.1:length(Interp_Trj1), 'spline');
Interp_Trj2(:, 2) = interp1(1:length(Interp_Trj1), Interp_Trj1(:, 2), 1:0.1:length(Interp_Trj1), 'spline');
Interp_Trj2(:, 3) = interp1(1:length(Interp_Trj1), Interp_Trj1(:, 3), 1:0.1:length(Interp_Trj1), 'spline');

s_seg = 1:length(Interp_Trj2);
S_seg = repmat(s_seg, Length, 1);
Step = repmat([(length(Interp_Trj2) / Length):(length(Interp_Trj2) / Length):length(Interp_Trj2)], length(s_seg), 1);
clear s_seg
[~, min_idx] = min(abs(S_seg' - Step));

ResizeTrj = Interp_Trj2(min_idx, :);
end
