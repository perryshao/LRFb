function descriptor = RRVdescriptor_BasedonFrames(frames, trajectory)
% RRVDESCRIPTOR_BASEDONFRAMES  Describe adjacent local reference frames.
%   frames is N-by-9, with tangent, normal and binormal in successive triples.
%   trajectory is N-by-3; output is N-by-7: quaternion then rotated SRV.
%   The final frame is paired with itself. Original author: skaegy.

if nargin < 3
    U = eye(3);
    SIGN = ones(1, 3);
end
T_joint_END = trajectory;
%% Compute square-root velocity and adjacent-frame rotations.
% Square-root velocity in the input coordinates.
Vel_joint_End = gradient((T_joint_END)')';
VelNorm_End = Vel_joint_End ./ repmat(sqrt((sqrt(sum(Vel_joint_End.^2, 2))) + eps), 1, 3);

current_frame = frames;
next_frame = [frames(2:end, :); frames(end, :)];
n = size(frames, 1);
descriptor = zeros(n, 7);
for i = 1:n
    % Historical tangent-axis intermediates.
    a = current_frame(i, 1:3);
    b = next_frame(i, 1:3); % use x axis
    % Map the current frame to its forward neighbor.
    RotM = reshape(next_frame(i, :), 3, 3) / reshape(current_frame(i, :), 3, 3);

    quat = Myrotm2quat(RotM');

    Velocity_axis = VelNorm_End(i, :) * RotM;
    descriptor(i, :) = [quat Velocity_axis];
end

end
