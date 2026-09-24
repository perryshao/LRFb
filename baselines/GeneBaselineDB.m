function GeneBaselineDB(method)
% GENEBASELINEDB  Table 1 baselines, in the place of the LRFb descriptor.
%
%   NOT ORIGINAL CODE.  Added when the project was consolidated (2026-09-23)
%   and never run: no MATLAB was available.  The original wiring that fed
%   the Table 1 baselines into the Fisher-vector pipeline was not found in
%   any source location (only the 3-D shape-context lines survive, commented
%   out, in each GeneSC.m).  This adapter reconnects the original descriptor
%   implementations to that pipeline.  See baselines/README.md.
%
%   GeneBaselineDB(method) is a drop-in replacement for GeneSC(bodyJoints)
%   in experiments/*/run.m.  It reads DB.mat / SAMPLES.mat and writes
%   RRV_DB.mat / RRV_SAMPLES.mat in exactly GeneSC's layout: one T-by-(d*2)
%   matrix per gesture, the two hand trajectories side by side in columns,
%   the first and last two frames dropped.  GeneFisherCodeJointPyramid_whole
%   and the linear SVM in run.m then run unchanged.
%
%   method   Table 1 row                  per-frame descriptor              d
%   'DI'     Differential Invariants [6]  descriptor_comp                   4
%   'II'     Integral Invariants [4]      integral_invariant(xyz,6,0.005)   2
%   'SRVF'   SRVF [20]                    x'(t) / sqrt(|x'(t)|)             3
%   'SSM'    Multiscale SSM [23]          Temporal_SSM + Log_hogcalculator  150
%            (Shao's single-scale SSM + log-polar HOG from the TSSM project;
%             Guo et al.'s multiscale implementation [23] was not found)
%   'SC3D'   3D shape context [13]        sc3d_compute on the LRFs          135 (IP)
%            (the version in the current experiment folder is used)       810 (MSRC-12)
%
%   Run it from inside an experiment folder, after load_*_bat has written
%   DB.mat and SAMPLES.mat, e.g.
%       GeneBaselineDB('DI');
%       GeneFisherCodeJointPyramid_whole(jointNum, ntotalbh, numClusters(1), 0);
%   For 'SSM' and 'SC3D' pass pcaFlag = 1 to keep the Fisher vectors small.

load DB;
RRV_DB = compute_set(TRAJDB, method, 'training'); %#ok<NASGU>
save RRV_DB RRV_DB -v7.3;
clear RRV_DB TRAJDB;
load SAMPLES;
RRV_SAMPLES = compute_set(TRAJSAMPLES, method, 'testing'); %#ok<NASGU>
save RRV_SAMPLES RRV_SAMPLES -v7.3;
end

function OUT = compute_set(TRAJ, method, tag)
traj_no = [1; 2]; % the two hands, as in GeneSC
nsamp = size(TRAJ, 2);
OUT = cell(1, nsamp);
for i = 1:nsamp
    fprintf('%s descriptors for %s data %d/%d...\n', method, tag, i, nsamp);
    % scale used by GeneSC for the shape context (and by the MSRC-12 LRFs)
    mean_dist = zeros(1, length(traj_no));
    for n = traj_no'
        xyz = TRAJ{2, i}(:, (n - 1) * 3 + 1:n * 3);
        mean_dist(n) = mean(mean(real(sqrt(dist2(xyz, xyz)))));
    end
    mean_dist_global = max(mean_dist);
    m = 0;
    for n = traj_no'
        m = m + 1;
        xyz = TRAJ{2, i}(:, (n - 1) * 3 + 1:n * 3);
        xyz(isnan(xyz(:, 1)), :) = [];
        T = size(xyz, 1) - 4;
        if ~any(xyz(:))
            des = zeros(T, descriptor_length(method));
        else
            des = frame_descriptor(xyz, method, mean_dist_global);
            des = des(3:end - 2, :); % same trimming as GeneSC
        end
        d = size(des, 2);
        OUT{1, i}(:, (m - 1) * d + 1:m * d) = des;
    end
end
end

function des = frame_descriptor(xyz, method, mean_dist_global)
switch upper(method)
    case 'DI'
        des = descriptor_comp(xyz);
    case 'II'
        des = integral_invariant(xyz, 6, 0.005); % IID / TSSM default
    case 'SRVF'
        v = gradient(xyz')';
        des = v ./ repmat(sqrt(sqrt(sum(v.^2, 2)) + eps), 1, 3);
    case 'SSM'
        S = Temporal_SSM(xyz, 5, 1, 0); % Euclidean SSM of raw xyz
        S = floor((S / max(S(:))) * (2^16 - 1));
        des = Log_hogcalculator(S);
    case 'SC3D'
        if nargin('Estimate_Frenet') > 2 % MSRC-12 version
            F = Estimate_Frenet(xyz, 5, mean_dist_global);
        else % IP version
            F = Estimate_Frenet(xyz, 5);
        end
        F = F(3:end - 2, :);
        B = xyz(3:end - 2, :);
        [BH_theta, BH_alpha] = sc3d_compute(B', F', mean_dist_global, ...
            18, 9, 5, 0.1, 2, zeros(1, size(B, 1)));
        if isscalar(BH_alpha) % true 3-D shape context
            h = BH_theta;
        else % pseudo 3-D shape context
            h = cat(2, BH_theta, BH_alpha);
        end
        % pad back to full length so the common trimming applies
        des = zeros(size(xyz, 1), size(h, 2));
        des(3:end - 2, :) = full(h);
    otherwise
        error('GeneBaselineDB:method', 'Unknown method %s', method);
end
end

function d = descriptor_length(method)
switch upper(method)
    case 'DI'
        d = 4;
    case 'II'
        d = 2;
    case 'SRVF'
        d = 3;
    case 'SSM'
        d = 150;
    case 'SC3D'
        if nargin('Estimate_Frenet') > 2
            d = 810;
        else
            d = 135;
        end
    otherwise
        error('GeneBaselineDB:method', 'Unknown method %s', method);
end
end
