function quat = Myrotm2quat(rotm)
% MYROTM2QUAT  Convert a 3-by-3 rotation matrix to unit [w x y z].
%   Original author: skaegy (skaegy@gmail.com).
%   Select a well-conditioned component to include rotations near 180 degrees.
%   The scalar component is nonnegative; at exactly 180 degrees either sign
%   represents the same rotation, with the dominant axis chosen positive.

tr = trace(rotm);
if tr > 0
    s = 2 * sqrt(1 + tr);
    quat = [s / 4, (rotm(3, 2) - rotm(2, 3)) / s, ...
            (rotm(1, 3) - rotm(3, 1)) / s, (rotm(2, 1) - rotm(1, 2)) / s];
elseif rotm(1, 1) >= rotm(2, 2) && rotm(1, 1) >= rotm(3, 3)
    s = 2 * sqrt(1 + rotm(1, 1) - rotm(2, 2) - rotm(3, 3));
    quat = [(rotm(3, 2) - rotm(2, 3)) / s, s / 4, ...
            (rotm(1, 2) + rotm(2, 1)) / s, (rotm(1, 3) + rotm(3, 1)) / s];
elseif rotm(2, 2) >= rotm(3, 3)
    s = 2 * sqrt(1 + rotm(2, 2) - rotm(1, 1) - rotm(3, 3));
    quat = [(rotm(1, 3) - rotm(3, 1)) / s, (rotm(1, 2) + rotm(2, 1)) / s, ...
            s / 4, (rotm(2, 3) + rotm(3, 2)) / s];
else
    s = 2 * sqrt(1 + rotm(3, 3) - rotm(1, 1) - rotm(2, 2));
    quat = [(rotm(2, 1) - rotm(1, 2)) / s, (rotm(1, 3) + rotm(3, 1)) / s, ...
            (rotm(2, 3) + rotm(3, 2)) / s, s / 4];
end
quat = quat / sqrt(sum(quat.^2));
if quat(1) < 0
    quat = -quat;
end
end
