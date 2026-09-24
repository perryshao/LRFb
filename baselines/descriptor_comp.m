function approx_des = descriptor_comp(curve_xyz)
% DESCRIPTOR_COMP  Compute differential invariants for an N-by-3 trajectory.
%   Returns N-by-4 rows: curvature, transformed curvature derivative,
%   torsion and transformed torsion derivative. Boundary rows remain zero.

samples = size(curve_xyz, 1);
K = [];
T = [];
K_d = [];
T_d = []; % initial the value of the descriptors
K(1:2) = 0;
T(1:2) = 0;
K_d(1:2) = 0;
T_d(1:2) = 0; %% beginning and ending value of descriptors are set to 0
% Torsion derivative computation use the Curvature's the second order derivative approximation
T_d_2(1:2) = 0;
K_d_2(1:2) = 0;
K(samples - 1:samples) = 0;
T(samples - 1:samples) = 0;
K_d(samples - 1:samples) = 0;
T_d(samples - 1:samples) = 0;
T_d_2(samples - 1:samples) = 0;
K_d_2(samples - 1:samples) = 0;
%% descriptor approximation computation for curvature and torsion
for i = 3:samples - 2 % there will get samples-6 descriptors
    %% distance computing
    a = norm(curve_xyz(i, :) - curve_xyz(i - 1, :));
    b = norm(curve_xyz(i + 1, :) - curve_xyz(i, :));
    c = norm(curve_xyz(i + 1, :) - curve_xyz(i - 1, :));
    d = norm(curve_xyz(i + 2, :) - curve_xyz(i + 1, :));
    e = norm(curve_xyz(i + 2, :) - curve_xyz(i, :));
    f = norm(curve_xyz(i + 2, :) - curve_xyz(i - 1, :));
    g = norm(curve_xyz(i - 1, :) - curve_xyz(i - 2, :));
    n = norm(curve_xyz(i, :) - curve_xyz(i - 2, :));
    m = norm(curve_xyz(i + 1, :) - curve_xyz(i - 2, :));

    %% compute h value: the distance between a point and a line
    l = curve_xyz(i + 1, :) - curve_xyz(i - 1, :);
    pl = curve_xyz(i, :) - curve_xyz(i - 1, :);
    tem = cross(pl, l);
    h = norm(tem) / norm(l);
    %% compute area for abc to get curvature
    s = (a + b + c) / 2;
    abc_area = sqrt(s * (s - a) * (s - b) * (s - c));
    K(i) = 4 * abc_area / (a * b * c);
    %% compute two kinds of H respectively
    terahedron_def = [curve_xyz(i, :) 1;
                      curve_xyz(i - 1, :) 1;
                      curve_xyz(i + 1, :) 1;
                      curve_xyz(i + 2, :) 1;
      ];
    terahedron_def_det = det(terahedron_def);
    terahedron_gnm = [curve_xyz(i, :) 1;
                      curve_xyz(i - 1, :) 1;
                      curve_xyz(i + 1, :) 1;
                      curve_xyz(i - 2, :) 1;
      ];
    terahedron_gnm_det = det(terahedron_gnm);
    V_abcdef = terahedron_def_det / factorial(3);
    V_abcgnm = terahedron_gnm_det / factorial(3);
    H_def = 3 * V_abcdef / abc_area;
    H_gnm = 3 * V_abcgnm / abc_area;
    %% get torsion
    T(i) = 3 * H_def / (d * e * f * K(i)) + 3 * H_gnm / (g * n * m * K(i));
end
%% compute descriptors for their derivatives respectively
for i = 3:samples - 2
    a = norm(curve_xyz(i, :) - curve_xyz(i - 1, :));
    b = norm(curve_xyz(i + 1, :) - curve_xyz(i, :));
    c = norm(curve_xyz(i + 1, :) - curve_xyz(i - 1, :));
    d = norm(curve_xyz(i + 2, :) - curve_xyz(i + 1, :));
    e = norm(curve_xyz(i + 2, :) - curve_xyz(i, :));
    f = norm(curve_xyz(i + 2, :) - curve_xyz(i - 1, :));
    g = norm(curve_xyz(i - 1, :) - curve_xyz(i - 2, :));
    n = norm(curve_xyz(i, :) - curve_xyz(i - 2, :));
    m = norm(curve_xyz(i + 1, :) - curve_xyz(i - 2, :));
    %% compute h value: the distance between a point and a line
    l = curve_xyz(i + 1, :) - curve_xyz(i - 1, :);
    pl = curve_xyz(i, :) - curve_xyz(i - 1, :);
    tem = cross(pl, l);
    h = norm(tem) / norm(l);
    r = 2 * a + 2 * b - 2 * d - 3 * h + g;
    % get another version of curvature for computing the torsion derivative
    % according to Mireille Boutin paper (Numerically Invariant Signature Curves)
    K_d_2(i) = (K(i + 1) - K(i - 1)) / c;
    K_d(i) = 3 * (K(i + 1) - K(i - 1)) / (2 * a + 2 * b + d + g);
    T_d_2(i) = 4 * ((T(i + 1) - T(i - 1) + r * (T(i) * K_d_2(i) / (6 * K(i)))) / (2 * a + 2 * b * 2 * d + h + g));
    T_d(i) = 4 * ((T(i + 1) - T(i - 1) + r * (T(i) * K_d(i) / (6 * K(i)))) / (2 * a + 2 * b * 2 * d + h + g));
end
%% compute descriptors for two global features s and r
S(1) = 0;
R = []; %% initial the value of the descriptors
s = 0;
r = 0; %% cumulant initial value for computing cumulant value of arc_length along time series
arc_length = 0;
center_length = 0; %% initial the value of arc_length and geometric center
for i = 2:samples
    arc_length = arc_length + norm(curve_xyz(i, :) - curve_xyz(i - 1, :));
end
for i = 2:samples
    s = s + norm(curve_xyz(i, :) - curve_xyz(i - 1, :));
    S(end + 1) = s / arc_length;
end
c = sum(curve_xyz) / samples;
for i = 1:samples
    center_length = center_length + norm(curve_xyz(i, :) - c);
end
for i = 1:samples
    r = samples * norm(curve_xyz(i, :) - c);
    R(end + 1) = r / center_length;
end
%% normalize the K_d and T_d to eliminate the scaling effect
K_d(3:end - 2) = (abs(K_d(3:end - 2)) ./ K_d(3:end - 2)) .* sqrt(abs(K_d(3:end - 2)));
T_d_2(3:end - 2) = (abs(T_d_2(3:end - 2)) ./ T_d_2(3:end - 2)) .* sqrt(abs(T_d_2(3:end - 2)));
% approx_des=[K' K_d' T' T_d' S' R'];
% Torsion derivative computation use the Curvature's the second order derivative approximation
% approx_des=[K' K_d' T' T_d_2' S' R'];

% nanflag = isnan(curve_xyz);
% if all(~nanflag)
%     real_index = isfinite(approx_des);
%     approx_des(~real_index) = 0;
% else
%     approx_des(find(approx_des >= Inf)) = 0;
%     approx_des(find(approx_des <= -Inf)) = 0;
% end

% discard the global parameter in multiple trajectories recognition -- perry 28/02/13
approx_des = [K' K_d' T' T_d_2'];
