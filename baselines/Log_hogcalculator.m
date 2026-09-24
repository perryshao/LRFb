function [F_t, F] = Log_hogcalculator(img, radius, nbins_theta, nbins_r, ...
    nthet, issigned, normmethod)
% LOG_HOGCALCULATOR  Compute log-polar HOG descriptors from an SSM image.
%   [F_t, F] = Log_hogcalculator(img, radius, nbins_theta, nbins_r, ...
%       nthet, issigned, normmethod)
%   Optional parameters control the spatial bins and gradient histogram.
%   F_t stores per-time descriptors; F stores the underlying histograms.
%   Requires imfilter from the Image Processing Toolbox.
%   Author: Perry, December 2014. Based on Dalal and Triggs (CVPR 2005).

if nargin < 2
    % set default parameters value.
    radius = 30; % default
    %     radius = 15;
    nbins_theta = 8; % default
    %     nbins_theta = 5;
    nbins_r = 4; % default
    %     nbins_r = 2;
    nthet = 6; % default
    %     nthet = 4;
    overlap = 0.95; %#ok<NASGU>
    %     issigned = 'signed';%(0-2pi)
    issigned = 'unsigned'; % (0-pi) default
    normmethod = 'l2hys'; % default
else
    if nargin < 7
        error('Input parameters are not enough.');
    end
end
% the cell size counted by the number of pixels
celltheta = pi / nbins_theta;
cellro = log(radius) / nbins_r;
% check parameters's validity.
[M, N, K] = size(img);
% construct the indx matrix
indx_matrix = zeros(2, sum(1:M));
n = 0;
for i = 1:M
    for j = i:M
        n = n + 1;
        indx_matrix(:, n) = [i; j];
    end
end
% calculate gradient scale matrix.
hx = [-1, 0, 1];
hy = -hx';
gradscalx = imfilter(double(img), hx);
gradscaly = imfilter(double(img), hy);

if K > 1
    maxgrad = sqrt(double(gradscalx .* gradscalx + gradscaly .* gradscaly));
    [gradscal, gidx] = max(maxgrad, [], 3);
    gxtemp = zeros(M, N);
    gytemp = gxtemp;
    for kn = 1:K
        ttempx = gradscalx(:, :, kn);
        ttempy = gradscaly(:, :, kn);
        tmpindex = find(gidx == kn);
        gxtemp(tmpindex) = ttempx(tmpindex);
        gytemp(tmpindex) = ttempy(tmpindex);
    end
    gradscalx = gxtemp;
    gradscaly = gytemp;
else
    gradscal = sqrt(double(gradscalx .* gradscalx + gradscaly .* gradscaly));
end

% calculate gradient orientation matrix.
% plus small number for avoiding dividing zero.
gradscalxplus = gradscalx + ones(size(gradscalx)) * 0.0001;
gradorient = zeros(M, N);
% unsigned situation: orientation region is 0 to pi.
if strcmp(issigned, 'unsigned') == 1
    gradorient = ...
        atan(gradscaly ./ gradscalxplus);
    gradorient(gradorient < 0) = gradorient(gradorient < 0) + pi;
    or = 1;
else
    % signed situation: orientation region is 0 to 2*pi.
    if strcmp(issigned, 'signed') == 1
        idx = find(gradscalx >= 0 & gradscaly >= 0);
        gradorient(idx) = atan(gradscaly(idx) ./ gradscalxplus(idx));
        idx = find(gradscalx < 0);
        gradorient(idx) = atan(gradscaly(idx) ./ gradscalxplus(idx)) + pi;
        idx = find(gradscalx >= 0 & gradscaly < 0);
        gradorient(idx) = atan(gradscaly(idx) ./ gradscalxplus(idx)) + 2 * pi;
        or = 2;
    else
        error('Incorrect ISSIGNED parameter.');
    end
end
% calculate block slide step.
% xbstride = uint16(radius)*(1-overlap);
% ybstride = uint16(radius)*(1-overlap);
xbstride = 1;
xbstridend = M;
% ybstridend = M;
% calculate the total blocks number in the window detected, which is
ntotalbh = M / xbstride;
% generate the matrix hist3dbig for storing the 3-dimensions histogram. the
% matrix covers the whole image in the 'globalinterpolate' condition or
% covers the local block in the 'localinterpolate' condition. The matrix is
% bigger than the area where it covers by adding additional elements
% (corresponding to the cells) to the surround for calculation convenience.

hist3dbig = zeros(nbins_theta + 2, nbins_r + 2, nthet + 2);
F = zeros(1, ntotalbh * nbins_theta * nbins_r * nthet);
% generate the matrix for storing histogram of one block;
sF = zeros(1, nbins_theta * nbins_r * nthet);
% vote for histogram. there are two situations according to the interpolate
% condition('global' interpolate or local interpolate). The hist3d which is
% generated from the 'bigger' matrix hist3dbig is the final histogram.
xbstep = xbstride;
% rotate angle
rot_theta = -pi / 4;
% block slide loop
for btlx = 1:xbstep:xbstridend
    btly = btlx;
    % transform the indx_matrix to transformed matrix with respect to the
    % diagonal pixels
    t_matrix = indx_matrix - repmat(double([btlx; btly]), 1, size(indx_matrix, 2));
    % the pixel indexes of log-polar block
    btpixels = indx_matrix(:, 0 < sqrt(sum(t_matrix .* t_matrix)) & ...
                                     sqrt(sum(t_matrix .* t_matrix)) <= radius);
    % the pixel indexes of the trasformed log-polar block
    npixels = size(btpixels, 2);
    rt_matrix = [cos(rot_theta) -sin(rot_theta); ...
        sin(rot_theta) cos(rot_theta)] * (btpixels - repmat(double([btlx; btly]), 1, npixels));
    % adjust the negative value caused by accuracy of floating-point
    % operations.these value's scale is very small, usually at E-04 magnitude
    rt_matrix(abs(rt_matrix) < 10e-4) = 0;
    r_array = sqrt(sum(rt_matrix .* rt_matrix));
    ro = log(r_array);
    theta_array = atan2(rt_matrix(2, :), rt_matrix(1, :));
    for bi = 1:npixels
        i = btpixels(1, bi);
        j = btpixels(2, bi);
        gs = gradscal(i, j);
        go = gradorient(i, j);
        jorbj = ro(bi);
        iorbi = theta_array(bi);
        % calculate bin index of hist3dbig
        binx1 = floor((jorbj + cellro / 2) / cellro) + 1;
        biny1 = floor((iorbi + celltheta / 2) / celltheta) + 1;
        binz1 = floor((go + (or * pi / nthet) / 2) / (or * pi / nthet)) + 1;

        if gs < 1E-5
            continue;
        end

        binx2 = binx1 + 1;
        biny2 = biny1 + 1;
        binz2 = binz1 + 1;

        x1 = (binx1 - 1.5) * cellro; % don't need add 0.5 here
        %         x1 = (binx1-1.5)*cellro+0.5;
        y1 = (biny1 - 1.5) * celltheta;
        z1 = (binz1 - 1.5) * (or * pi / nthet);

        % trilinear interpolation.
        hist3dbig(biny1, binx1, binz1) = ...
            hist3dbig(biny1, binx1, binz1) + gs * ...
             (1 - (jorbj - x1) / cellro) * (1 - (iorbi - y1) / celltheta) * ...
            (1 - (go - z1) / (or * pi / nthet));
        hist3dbig(biny1, binx1, binz2) = ...
            hist3dbig(biny1, binx1, binz2) + gs * ...
             (1 - (jorbj - x1) / cellro) * (1 - (iorbi - y1) / celltheta) * ...
            ((go - z1) / (or * pi / nthet));
        hist3dbig(biny2, binx1, binz1) = ...
            hist3dbig(biny2, binx1, binz1) + gs * ...
            (1 - (jorbj - x1) / cellro) * ((iorbi - y1) / celltheta) * ...
            (1 - (go - z1) / (or * pi / nthet));
        hist3dbig(biny2, binx1, binz2) = ...
            hist3dbig(biny2, binx1, binz2) + gs * ...
            (1 - (jorbj - x1) / cellro) * ((iorbi - y1) / celltheta) * ...
            ((go - z1) / (or * pi / nthet));
        hist3dbig(biny1, binx2, binz1) = ...
            hist3dbig(biny1, binx2, binz1) + gs * ...
            ((jorbj - x1) / cellro) * (1 - (iorbi - y1) / celltheta) * ...
            (1 - (go - z1) / (or * pi / nthet));
        hist3dbig(biny1, binx2, binz2) = ...
            hist3dbig(biny1, binx2, binz2) + gs * ...
            ((jorbj - x1) / cellro) * (1 - (iorbi - y1) / celltheta) * ...
            ((go - z1) / (or * pi / nthet));
        hist3dbig(biny2, binx2, binz1) = ...
            hist3dbig(biny2, binx2, binz1) + gs * ...
            ((jorbj - x1) / cellro) * ((iorbi - y1) / celltheta) * ...
            (1 - (go - z1) / (or * pi / nthet));
        hist3dbig(biny2, binx2, binz2) = ...
            hist3dbig(biny2, binx2, binz2) + gs * ...
            ((jorbj - x1) / cellro) * ((iorbi - y1) / celltheta) * ...
            ((go - z1) / (or * pi / nthet));
    end

    % In the local interpolate condition. F is generated in this block
    % slide loop. hist3dbig should be cleared in each loop.
    if or == 2
        hist3dbig(:, :, 2) = hist3dbig(:, :, 2) + ...
            hist3dbig(:, :, nthet + 2);
        hist3dbig(:, :, (nthet + 1)) = ...
            hist3dbig(:, :, (nthet + 1)) + hist3dbig(:, :, 1);
    end
    hist3d = hist3dbig(2:(nbins_theta + 1), 2:(nbins_r + 1), 2:(nthet + 1));
    for ibin = 1:nbins_theta
        for jbin = 1:nbins_r
            idsF = nthet * ((ibin - 1) * nbins_r + jbin - 1) + 1;
            idsF = idsF:(idsF + nthet - 1);
            sF(idsF) = hist3d(ibin, jbin, :);
        end
    end
    iblock = ((btlx - 1) / xbstride) + 1;
    idF = (iblock - 1) * nbins_theta * nbins_r * nthet + 1;
    idF = idF:(idF + nbins_theta * nbins_r * nthet - 1);
    F(idF) = sF;
    hist3dbig(:, :, :) = 0;
end
% adjust the negative value caused by accuracy of floating-point
% operations.these value's scale is very small, usually at E-03 magnitude
% while others will be E+02 or E+03 before normalization.
F(F < 0) = 0;
% block normalization.
e = 0.001;
l2hysthreshold = 0.6;
% l2hysthreshold = 0.2;
fslidestep = nbins_r * nbins_theta * nthet;
switch normmethod
    case 'none'
    case 'l1'
        for fi = 1:fslidestep:size(F, 2)
            div = sum(F(fi:(fi + fslidestep - 1)));
            F(fi:(fi + fslidestep - 1)) = F(fi:(fi + fslidestep - 1)) / (div + e);
        end
    case 'l1sqrt'
        for fi = 1:fslidestep:size(F, 2)
            div = sum(F(fi:(fi + fslidestep - 1)));
            F(fi:(fi + fslidestep - 1)) = sqrt(F(fi:(fi + fslidestep - 1)) / (div + e));
        end
    case 'l2'
        for fi = 1:fslidestep:size(F, 2)
            sF = F(fi:(fi + fslidestep - 1)) .* F(fi:(fi + fslidestep - 1));
            div = sqrt(sum(sF) + e * e);
            F(fi:(fi + fslidestep - 1)) = F(fi:(fi + fslidestep - 1)) / div;
        end
    case 'l2hys'
        for fi = 1:fslidestep:size(F, 2)
            sF = F(fi:(fi + fslidestep - 1)) .* F(fi:(fi + fslidestep - 1));
            div = sqrt(sum(sF) + e * e);
            sF = F(fi:(fi + fslidestep - 1)) / div;
            sF(sF > l2hysthreshold) = l2hysthreshold;
            div = sqrt(sum(sF .* sF) + e * e);
            F(fi:(fi + fslidestep - 1)) = sF / div;
        end
    otherwise
        error('Incorrect NORMMETHOD parameter.');
end
F_t = reshape(F, fslidestep, ntotalbh)';
%% relocate the bins since that the center bins will be summed up to one bin
temp = 0;
nbins = nbins_r * nbins_theta - nbins_theta + 1; % the sum bins for single center
F_temp = zeros(ntotalbh, nbins * nthet);
for i = 1:nbins_r:nbins_r * nbins_theta
    temp = F_t(:, (i - 1) * nthet + 1:i * nthet) + temp;
    F_temp(:, 1:nthet) = temp;
end
for j = 0:nbins_theta - 1
    for i = 2 + j * nbins_r:2 + j * nbins_r + (nbins_r - 2)
        F_temp(:, (i - j - 1) * nthet + 1:(i - j) * nthet) = F_t(:, (i - 1) * nthet + 1:i * nthet);
    end
end
F_t = F_temp;
