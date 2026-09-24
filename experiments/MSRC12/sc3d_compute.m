function [BH_theta, BH_alpha, mean_dist] = sc3d_compute(Bsamp, Tsamp, mean_dist, nbins_theta, nbins_alpha, ...
    nbins_r, r_inner, r_outer, out_vec)
% SC3D_COMPUTE  Build a joint azimuth/polar/radius histogram.
%   BH_theta contains the full histogram; BH_alpha is a zero placeholder.
%
% compute (r,theta,alpha) histograms for points along boundary
%
% Bsamp is 3 x nsamp (x and y ,z coords.)
% Tsamp is 9 x nsamp (Frenet_Frames)
% out_vec is 1 x nsamp (0 for inlier, 1 for outlier)
%
% mean_dist is the mean distance, used for length normalization
% Pass [] for mean_dist to estimate it from the inlier coordinates.
%
% outliers are not counted in the histograms, but they do get
% assigned a histogram
%

% compute r
r_array = real(sqrt(dist2(Bsamp', Bsamp'))); % real is needed to
% prevent bug in Unix version

nsamp = size(Bsamp, 2);
in_vec = out_vec == 0;
theta_array = zeros(nsamp, nsamp);
alpha_array = zeros(nsamp, nsamp);
% compute,theta arrays
for i = 1:nsamp
    Relative_xyz = bsxfun(@minus, Bsamp, Bsamp(:, i));
    %     Rot_Matrix = reshape(Tsamp(:,10),3,3);
    Rot_Matrix = reshape(Tsamp(:, i), 3, 3);
    Local_xyz = Rot_Matrix \ Relative_xyz;
    theta_array(i, :) = atan2(Local_xyz(2, :), Local_xyz(1, :))';
    r_array_xy = sqrt(sum(Local_xyz(1:2, :).^2, 1)); % get the sqrt(x2+y2) to compute alpha
    alpha_array(i, :) = atan2(r_array_xy, Local_xyz(3, :))';
end

% normalize distance by mean, ignoring outliers
if isempty(mean_dist)
    tmp = r_array(in_vec, :);
    tmp = tmp(:, in_vec);
    mean_dist = mean(tmp(:));
end
r_array_n = r_array / mean_dist;

% use a log. scale for binning the distances
r_bin_edges = logspace(log10(r_inner), log10(r_outer), nbins_r);
r_array_q = zeros(nsamp, nsamp);
for m = 1:nbins_r
    r_array_q = r_array_q + (r_array_n < r_bin_edges(m));
end
fz = r_array_q > 0; % flag all points inside outer boundary

% put all angles in [0,2pi) range
theta_array_2 = rem(rem(theta_array, 2 * pi) + 2 * pi, 2 * pi);

% quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
theta_array_q = 1 + floor(theta_array_2 / (2 * pi / nbins_theta));
alpha_array_q = 1 + floor(alpha_array / (pi / nbins_alpha));

%% 3d shape context are divided into 2 2d shape contexts
% nbins1=nbins_theta*nbins_r;
% BH_theta=zeros(nsamp,nbins1);
% for n=1:nsamp
%    fzn=fz(n,:)&in_vec;
%    Sn=sparse(theta_array_q(n,fzn),r_array_q(n,fzn),1,nbins_theta,nbins_r);
%    BH_theta(n,:)=Sn(:)';
% end

% nbins2=nbins_alpha*nbins_r;
% BH_alpha=zeros(nsamp,nbins2);
% for n=1:nsamp
%    fzn=fz(n,:)&in_vec;
%    Sn=sparse(alpha_array_q(n,fzn),r_array_q(n,fzn),1,nbins_alpha,nbins_r);
%    BH_alpha(n,:)=Sn(:)';
% end

%% Perry modified to real 3d shape context
nbins = nbins_theta * nbins_alpha * nbins_r;
BH_theta = zeros(nsamp, nbins);
for n = 1:nsamp
    fzn = fz(n, :) & in_vec;
    Sn = ndSparse.build([theta_array_q(n, fzn); alpha_array_q(n, fzn); r_array_q(n, fzn)]', 1, [nbins_theta, ...
        nbins_alpha, nbins_r]);
    BH_theta(n, :) = Sn(:)';
end
BH_alpha = 0;
