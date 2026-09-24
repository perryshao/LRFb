function area = estimate_integral(Curve_xyz, Ck, Cl, Cr, Radius, pos)
% ESTIMATE_INTEGRAL  Estimate two local integral invariants at a trajectory point.
%   Curve_xyz is the normalized trajectory; Ck, Cl and Cr define its
%   local frame. Radius sets the integration neighborhood; pos indexes Ck.

area = 0; %#ok<NASGU>
T1 = Ck - Cl;
T2 = Cr - Ck;
T1 = T1 / norm(T1);
T2 = T2 / norm(T2);

% B = cross(T1,T2); B= B/norm(B);
% T = (T1+T2)/norm(T1+T2);% modified by perry according to reviews of PR
% N = cross(B,T);

Cc = tricircumcenter3d(Cl, Ck, Cr) + Cl; % get the circumcenter
N = (Cc - Ck) / norm(Cc - Ck);
B = cross(T1, T2);
B = B / norm(B);
T = cross(N, B);

D = -(B * Ck');

% [~, grid] = curve_grid(Segment_curve,500);
% grid = 1/2000; %grid resolution
grid = 2 * Radius / 500;

%% get the projection on osculating plane
t_dnom = B * B';
t_nom = Curve_xyz * B' + D;
t = t_nom / t_dnom;
Po = Curve_xyz - t * B;
%% set to zero when the kernel reach the bound
% setzero_flag = 0;
%% tune the Radius when there are double features
% r_k = pos;r_r = 0;Po_length = size(Po,1);flag=0;
% while (r_r <= Radius && r_k < Po_length)
%     if (Po(r_k+1,:)-Po(r_k,:))*T'< 0
%         flag = flag +1;
%         if flag > 2
%             Radius = norm(Po(r_k+1,:)-Po(pos,:));
%             break;
%         end
%     end
%     r_k = r_k+1;
%     r_r = norm(Po(r_k,:)-Po(pos,:));
% end

%% get the area for osculating plane projection
r_l = 0;
l_k = pos;
Po_length = size(Po, 1);
area_o = [];

%% new method of approximation area
while (r_l <= Radius && l_k > 1)
    l_increment = Po(l_k, :) - Po(l_k - 1, :);
    l_increment_norm = norm((Po(l_k, :) - Po(l_k - 1, :)));
    costheta = ((Po(l_k, :) - Po(l_k - 1, :)) * T') / (l_increment_norm * norm(T));
    grid_l = grid / costheta;
    l_direction = l_increment / l_increment_norm;
    n_l = fix(l_increment_norm / grid_l);
    %     margin = abs(l_increment_norm/grid_l - n_l);
    m = abs(n_l);
    P_l = zeros(m, 3);

    for i = sign(n_l):sign(n_l):n_l
        P_l(abs(i), :) = Po(l_k, :) - i * grid_l * l_direction;
    end
    if isempty(P_l)
        l_k = l_k - 1;
        r_l = norm(Po(l_k, :) - Po(pos, :));
        continue;
    end
    P_l_vector = (P_l - repmat(Po(pos, :), m, 1));
    P_l_vector_norm = sqrt(sum(P_l_vector.^2, 2));
    P_l_vector((P_l_vector_norm >= Radius), :) = [];
    if isempty(P_l_vector)
        l_k = l_k - 1;
        r_l = norm(Po(l_k, :) - Po(pos, :));
        continue;
    end
    n = size(P_l_vector, 1);
    height_u = sqrt(abs(repmat(Radius^2, n, 1) - (P_l_vector * (-T)').^2));
    height_b = P_l_vector * N';
    height = height_u - height_b;
    %         area_o(end+1) = sum(grid*height)+height(end)*margin*grid;
    area_o(end + 1) = sum(grid * height) * sign(grid_l);
    l_k = l_k - 1;
    r_l = norm(Po(l_k, :) - Po(pos, :));
end

% set to zero when the kernel reach the bound
% if l_k <= 1 && r_l< Radius
%     setzero_flag = 1;
% end

r_k = pos;
r_r = 0;
while (r_r <= Radius && r_k < Po_length)
    r_increment = Po(r_k + 1, :) - Po(r_k, :);
    r_increment_norm = norm((Po(r_k + 1, :) - Po(r_k, :)));
    costheta = ((Po(r_k + 1, :) - Po(r_k, :)) * T') / (r_increment_norm * norm(T));
    grid_r = grid / costheta;
    r_direction = r_increment / r_increment_norm;
    n_r = fix(r_increment_norm / grid_r);
    %     margin = abs(r_increment_norm/grid_r - n_r);
    m = abs(n_r);
    P_r = zeros(m, 3);

    for i = sign(n_r):sign(n_r):n_r
        P_r(abs(i), :) = Po(r_k, :) + i * grid_r * r_direction;
    end

    m = size(P_r, 1);
    if isempty(P_r)
        r_k = r_k + 1;
        r_r = norm(Po(r_k, :) - Po(pos, :));
        continue;
    end
    P_r_vector = (P_r - repmat(Po(pos, :), m, 1));
    P_r_vector_norm = sqrt(sum(P_r_vector.^2, 2));
    P_r_vector((P_r_vector_norm >= Radius), :) = [];
    if isempty(P_r_vector)
        r_k = r_k + 1;
        r_r = norm(Po(r_k, :) - Po(pos, :));
        continue;
    end

    n = size(P_r_vector, 1);
    height_u = sqrt(abs(repmat(Radius^2, n, 1) - (P_r_vector * T').^2));
    height_b = P_r_vector * N';
    height = height_u - height_b;
    %         area_o(end+1) = sum(grid*height)+height(end)*margin*grid;
    area_o(end + 1) = sum(grid * height) * sign(grid_r);
    r_k = r_k + 1;
    r_r = norm(Po(r_k, :) - Po(pos, :));
end

% set to zero when the kernel reach the bound
% if r_k >= Po_length && r_r < Radius
%     setzero_flag = 1;
% end

area_o(isnan(area_o)) = [];
area_o = sum(area_o) / (pi * Radius^2);

% if setzero_flag==1
%     area_o=0.5;
% end

if area_o > 0.5
    area_o = 1 - area_o;
end

%% get the area for rectifying plane projection

D = -(N * Ck');
t_dnom = N * N';
t_nom = Curve_xyz * N' + D;
t = t_nom / t_dnom;
Pr = Curve_xyz - t * N;
%% set to zero when the kernel reach the bound
% setzero_flag = 0;
%% tune the Radius dynamically
% r_k = pos;r_r = 0;Pr_length = size(Pr,1);flag=0;
% while (r_r <= Radius && r_k < Pr_length)
%     if (Pr(r_k+1,:)-Pr(r_k,:))*T'< 0
%         flag = flag +1;
%         if flag > 2
%             Radius = norm(Pr(r_k+1,:)-Pr(pos,:));
%             break;
%         end
%     end
%     r_k = r_k+1;
%     r_r = norm(Pr(r_k,:)-Pr(pos,:));
% end
%% new method of approximation area - retifying plane projection
r_l = 0;
l_k = pos;
Pr_length = size(Pr, 1);
area_r = [];
while (r_l <= Radius && l_k > 1)
    l_increment = Pr(l_k, :) - Pr(l_k - 1, :);
    l_increment_norm = norm((Pr(l_k, :) - Pr(l_k - 1, :)));
    costheta = ((Pr(l_k, :) - Pr(l_k - 1, :)) * T') / (l_increment_norm * norm(T));
    grid_l = grid / costheta;
    l_direction = l_increment / l_increment_norm;
    n_l = fix(l_increment_norm / grid_l);
    %     margin = abs(l_increment_norm/grid_l - n_l);
    m = abs(n_l);
    P_l = zeros(m, 3);

    for i = sign(n_l):sign(n_l):n_l
        P_l(abs(i), :) = Pr(l_k, :) - i * grid_l * l_direction;
    end
    m = size(P_l, 1);
    if isempty(P_l)
        l_k = l_k - 1;
        r_l = norm(Pr(l_k, :) - Pr(pos, :));
        continue;
    end
    P_l_vector = (P_l - repmat(Pr(pos, :), m, 1));
    P_l_vector_norm = sqrt(sum(P_l_vector.^2, 2));
    P_l_vector((P_l_vector_norm >= Radius), :) = [];
    if isempty(P_l_vector)
        l_k = l_k - 1;
        r_l = norm(Pr(l_k, :) - Pr(pos, :));
        continue;
    end

    n = size(P_l_vector, 1);
    height_u = sqrt(abs(repmat(Radius^2, n, 1) - (P_l_vector * (-T)').^2));
    height_b = P_l_vector * B';
    height = height_u - height_b;
    %         area_r(end+1) = sum(grid*height)+height(end)*margin*grid;
    area_r(end + 1) = sum(grid * height) * sign(grid_l);
    l_k = l_k - 1;
    r_l = norm(Pr(l_k, :) - Pr(pos, :));
end

% set to zero when the kernel reach the bound
% if l_k <= 1 && r_l< Radius
%     setzero_flag = 1;
% end

r_k = pos;
r_r = 0;
while (r_r <= Radius && r_k < Pr_length)
    r_increment = Pr(r_k + 1, :) - Pr(r_k, :);
    r_increment_norm = norm((Pr(r_k + 1, :) - Pr(r_k, :)));
    costheta = ((Pr(r_k + 1, :) - Pr(r_k, :)) * T') / (r_increment_norm * norm(T));
    grid_r = grid / costheta;
    r_direction = r_increment / r_increment_norm;
    n_r = fix(r_increment_norm / grid_r);
    %     margin = abs(r_increment_norm/grid_r - n_r);
    m = abs(n_r);
    P_r = zeros(m, 3);

    for i = sign(n_r):sign(n_r):n_r
        P_r(abs(i), :) = Pr(r_k, :) + i * grid_r * r_direction;
    end

    m = size(P_r, 1);
    if isempty(P_r)
        r_k = r_k + 1;
        r_r = norm(Pr(r_k, :) - Pr(pos, :));
        continue;
    end
    P_r_vector = (P_r - repmat(Pr(pos, :), m, 1));
    P_r_vector_norm = sqrt(sum(P_r_vector.^2, 2));
    P_r_vector((P_r_vector_norm >= Radius), :) = [];
    if isempty(P_r_vector)
        r_k = r_k + 1;
        r_r = norm(Pr(r_k, :) - Pr(pos, :));
        continue;
    end

    n = size(P_r_vector, 1);
    height_u = sqrt(abs(repmat(Radius^2, n, 1) - (P_r_vector * T').^2));
    height_b = P_r_vector * B';
    height = height_u - height_b;
    %         area_r(end+1) = sum(grid*height)+height(end)*margin*grid;
    area_r(end + 1) = sum(grid * height) * sign(grid_r);
    r_k = r_k + 1;
    r_r = norm(Pr(r_k, :) - Pr(pos, :));
end
% set to zero when the kernel reach the bound
% if r_k >= Po_length && r_r < Radius
%     setzero_flag = 1;
% end

area_r(isnan(area_r)) = [];
area_r = abs(sum(area_r) / (pi * Radius^2));

% if setzero_flag==1
%     area_r=0.5;
% end
% if ~isempty(P_r) && (P_r(end,:) - Pr(pos,:))*T' < 0
%     area_r = 0.5-area_r;
% end

area = [area_o area_r];

%% get the area for normal plane projection
% D = -(T*Ck');
% t_dnom = T*T';
% t_nom = Curve_xyz*T'+D;
% t = t_nom/t_dnom;
% Pn = Curve_xyz-t*T;
% setzero_flag=0;

% %% tune the Radius dynamically
% % r_k = pos;r_r = 0;Pn_length = size(Pn,1);flag=0;
% % while (r_r <= 1 && r_k < Pn_length)
% %     if (Pn(r_k+1,:)-Pn(r_k,:))*B'< 0
% %         flag = flag +1;
% %         if flag > 1
% %             Radius = norm(Pn(r_k+1,:)-Pn(pos,:));
% %             break;
% %         end
% %     end
% %     r_k = r_k+1;
% %     r_r = norm(Pn(r_k,:)-Pn(pos,:));
% % end
% %% new method of approximation area - normal plane
% r_l = 0;l_k = pos;P_l = [];
% Pn_length = size(Pn,1);
% area_n = [];
% while (r_l <= Radius && l_k > 1)

%     l_increment = Pn(l_k,:)-Pn(l_k-1,:);
%     l_increment_norm = norm((Pn(l_k,:)-Pn(l_k-1,:)));
%     costheta = ((Pn(l_k,:)-Pn(l_k-1,:))*B')/(l_increment_norm*norm(B));
%     grid_l = grid/costheta;
%     l_direction = l_increment/l_increment_norm;
%     n_l = fix(l_increment_norm/grid_l);
%     margin = abs(l_increment_norm/grid_l - n_l);
%     m = abs(n_l);
%     P_l = zeros(m,3);

%     for i=sign(n_l):sign(n_l):n_l
%         P_l(abs(i),:) = Pn(l_k,:)-i*grid_l*l_direction;
%     end
%     m = size(P_l,1);
%     if isempty(P_l)
%         l_k = l_k-1;
%         r_l = norm(Pn(l_k,:)-Pn(pos,:));
%         continue;
%     end
%     P_l_vector = (P_l - repmat(Pn(pos,:),m,1));
%     P_l_vector_norm = sqrt(sum(P_l_vector.^2,2));
%     P_l_vector((P_l_vector_norm >= Radius),:) = [];
%     if isempty(P_l_vector)
%         l_k = l_k-1;
%         r_l = norm(Pn(l_k,:)-Pn(pos,:));
%         continue;
%     end
%     n = size(P_l_vector,1);
%     height_u = sqrt(abs(repmat(Radius^2,n,1) - (P_l_vector*(-B)').^2));
%     height_b = P_l_vector*N';
%     height = height_u-height_b;
%     area_n(end+1) = sum(grid*height)*sign(grid_l);

%     l_k = l_k-1;
%     r_l = norm(Pn(l_k,:)-Pn(pos,:));
% end
% % set to zero when the kernel reach the bound
% if l_k <= 1 && r_l< Radius
%     setzero_flag = 1;
% end

% r_r = 0;P_r = [];r_k = pos;
% while  (r_r <= Radius && r_k < Pn_length),

%     r_increment = Pn(r_k+1,:)-Pn(r_k,:);
%     r_increment_norm = norm((Pn(r_k+1,:)-Pn(r_k,:)));
%     costheta = ((Pn(r_k+1,:)-Pn(r_k,:))*B')/(r_increment_norm*norm(B));
%     grid_r = grid/costheta;
%     r_direction = r_increment/r_increment_norm;
%     n_r = fix(r_increment_norm/grid_r);
%     margin = abs(r_increment_norm/grid_r - n_r);
%     m = abs(n_r);
%     P_r = zeros(m,3);

%     for i=sign(n_r):sign(n_r):n_r
%         P_r(abs(i),:) = Pn(r_k,:)+i*grid_r*r_direction;
%     end

%     m = size(P_r,1);
%     if isempty(P_r)
%         r_k = r_k+1;
%         r_r = norm(Pn(r_k,:)-Pn(pos,:));
%         continue;
%     end
%     P_r_vector = (P_r - repmat(Pn(pos,:),m,1));
%     P_r_vector_norm = sqrt(sum(P_r_vector.^2,2));
%     P_r_vector((P_r_vector_norm >= Radius),:) = [];
%     if isempty(P_r_vector)
%         r_k = r_k+1;
%         r_r = norm(Pn(r_k,:)-Pn(pos,:));
%         continue;
%     end

%     n = size(P_r_vector,1);
%     height_u = sqrt(abs(repmat(Radius^2,n,1) - (P_r_vector*B').^2));
%     height_b = P_r_vector*N';
%     height = height_u-height_b;
%     area_n(end+1) = sum(grid*height)*sign(grid_r);

%     r_k = r_k+1;
%     r_r = norm(Pn(r_k,:)-Pn(pos,:));
% end
% % set to zero when the kernel reach the bound
% if r_k >= Po_length && r_r < Radius
%     setzero_flag = 1;
% end

% area_n(isnan(area_n)) = [];
% area_n = sum(area_n)/(pi*Radius^2);

% if setzero_flag == 1
%     area_n=0.5;
% end

% area_n=0.5-area_n;

% area = [area_o area_n];
