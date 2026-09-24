function quat = Myrotm2quat(rotm)
% MYROTM2QUAT  Convert a 3-by-3 rotation matrix to [w x y z].
%   Original author: skaegy (skaegy@gmail.com).
%   Uses the trace-based formula; the 180-degree singularity is preserved.

qw1 = cos(acos((trace(rotm) - 1) / 2) / 2);
qw = sqrt(1 + trace(rotm)) / 2;
qx = (rotm(3, 2) - rotm(2, 3)) / (4 * qw);
qy = (rotm(1, 3) - rotm(3, 1)) / (4 * qw);
qz = (rotm(2, 1) - rotm(1, 2)) / (4 * qw);

quat = [qw qx qy qz];
end
