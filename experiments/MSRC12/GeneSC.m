function [SC_DB, SC_SAMPLES] = GeneSC(bodyJoints)
% GENESC  Generate the saved LRFb descriptors for the two hand trajectories.
%   Reads DB.mat / SAMPLES.mat; writes RRV_DB.mat / RRV_SAMPLES.mat.
%   Each gesture is a T-by-14 matrix with the hands concatenated in columns.
%   The legacy SC outputs and commented shape-context alternatives are retained.
%   Use the saved files: SC_DB is cleared before this function returns.

%% Load the training trajectories.
load DB;
samples_r = size(TRAJDB, 2);
SC_DB = cell(1, samples_r);
RRV_DB = cell(1, samples_r);
Joint_ID = bodyJoints;
% Parameters for computing shape context
mean_dist_global = []; % use [] to estimate scale from the data
nbins_theta = 18;
nbins_alpha = 9;
nbins_r = 5;
ndum1 = 0;
eps_dum = 0.15;
r_inner = 0.1;
r_outer = 2; % r_outer=2.5;
% total_bins=nbins_r*nbins_theta+nbins_r*nbins_alpha;  % pseudo 3d shape context
total_bins = nbins_r * nbins_theta * nbins_alpha; % true 3d shape context

traj_no = [1; 2]; % for data from Gy
% traj_no=[1;2;3;4;5];%
mean_dist = zeros(1, length(traj_no));
width = 5;
%% Compute training LRFb descriptors.
for i = 1:samples_r
    fprintf ('get the 3D shape context descriptors for training data %d/%d...\n', i, samples_r);

    for n = traj_no'
        marker_xyz = TRAJDB{2, i}(:, (n - 1) * 3 + 1:n * 3);
        r_array = real(sqrt(dist2(marker_xyz, marker_xyz))); % real is needed to  prevent bug in Unix version
        mean_dist(n) = mean(r_array(:));
    end
    mean_dist_global = max(mean_dist); % mean_dist_global is the mean distance, used for length normalization

    m = 0;
    for n = traj_no'
        m = m + 1;
        marker_xyz = TRAJDB{2, i}(:, (n - 1) * 3 + 1:n * 3);
        NaN_index = isnan(marker_xyz(:, 1));
        marker_xyz(NaN_index, :) = [];
        nsamp1 = size(marker_xyz, 1);
        % outliers on each iteration
        out_vec_1 = zeros(1, nsamp1 - 4); % beginning and ending elements are excluded
        T = nsamp1 - 4;
        if any(any(marker_xyz)) == 0
            SC_DB{1, i}(1:T, (m - 1) * total_bins + 1:m * total_bins) = 0;
            RRV_DB{1, i}(1:T, (m - 1) * 7 + 1:m * 7) = 0;
        else
            % Frenet Frames
            FrenetVector = Estimate_Frenet(marker_xyz, width, mean_dist_global);
            FVector = FrenetVector(3:end - 2, :);
            % compute shape contexts for (transformed) model
            Bsamp = marker_xyz(3:end - 2, :);
            %             [BH_theta,BH_alpha,mean_dist_1] = ...
            %             sc3d_compute(Bsamp', FVector', mean_dist_global, ...
            %                 nbins_theta, nbins_alpha, nbins_r, r_inner, r_outer, out_vec_1);
            %             SC_DB{1,i}((m-1)*T+1:m*T,:)  = cat(2,BH_theta, BH_alpha);
            %             SC_DB{1,i}(:,(m-1)*total_bins+1:m*total_bins)  = cat(2,BH_theta, BH_alpha);  % pseudo 3d shape
            %             context
            %             SC_DB{1,i}(:,(m-1)*total_bins+1:m*total_bins)  = BH_theta;% true 3d shape context
            RRV_DB{1, i}(:, (m - 1) * 7 + 1:m * 7) = RRVdescriptor_BasedonFrames(FVector, Bsamp);
        end
    end
end
save SC_DB SC_DB -v7.3;
save RRV_DB RRV_DB -v7.3;
clear SC_DB RRV_DB;
%% Compute test LRFb descriptors.
load SAMPLES;
samples_t = size(TRAJSAMPLES, 2);
SC_SAMPLES = cell(1, samples_t);
RRV_SAMPLES = cell(1, samples_t);
for i = 1:samples_t
    fprintf ('get the 3D shape context descriptors for testing data %d/%d...\n', i, samples_t);

    for n = traj_no'
        marker_xyz = TRAJSAMPLES{2, i}(:, (n - 1) * 3 + 1:n * 3);
        r_array = real(sqrt(dist2(marker_xyz, marker_xyz))); % real is needed to  prevent bug in Unix version
        mean_dist(n) = mean(r_array(:));
    end
    mean_dist_global = max(mean_dist); % mean_dist_global is the mean distance, used for length normalization

    m = 0;
    for n = traj_no'
        m = m + 1;
        marker_xyz = TRAJSAMPLES{2, i}(:, (n - 1) * 3 + 1:n * 3);
        NaN_index = isnan(marker_xyz(:, 1));
        marker_xyz(NaN_index, :) = [];
        nsamp1 = size(marker_xyz, 1);
        % outliers on each iteration
        out_vec_1 = zeros(1, nsamp1 - 4); % beginning and ending elements are excluded
        T = nsamp1 - 4;
        if any(any(marker_xyz)) == 0
            SC_SAMPLES{1, i}(1:nsamp1 - 4, (m - 1) * total_bins + 1:m * total_bins) = 0;
            SC_SAMPLES{1, i}(1:T, (m - 1) * total_bins + 1:m * total_bins) = 0;
            RRV_SAMPLES{1, i}(1:T, (m - 1) * 7 + 1:m * 7) = 0;
        else
            % Frenet Frames
            FrenetVector = Estimate_Frenet(marker_xyz, width, mean_dist_global);
            FVector = FrenetVector(3:end - 2, :);
            % compute shape contexts for (transformed) model
            Bsamp = marker_xyz(3:end - 2, :);
            %             [BH_theta,BH_alpha,mean_dist_1] = ...
            %             sc3d_compute(Bsamp', FVector', mean_dist_global, ...
            %                 nbins_theta, nbins_alpha, nbins_r, r_inner, r_outer, out_vec_1);
            %             SC_SAMPLES{1,i}((m-1)*T+1:m*T,:)  = cat(2,BH_theta, BH_alpha);
            %             SC_SAMPLES{1,i}(:,(m-1)*total_bins+1:m*total_bins)  = cat(2,BH_theta, BH_alpha);  % pseudo 3d
            %             shape context
            %             SC_SAMPLES{1,i}(:,(m-1)*total_bins+1:m*total_bins)  = BH_theta;% true 3d shape context
            RRV_SAMPLES{1, i}(:, (m - 1) * 7 + 1:m * 7) = RRVdescriptor_BasedonFrames(FVector, Bsamp);
        end
    end
end
save SC_SAMPLES SC_SAMPLES -v7.3;
save RRV_SAMPLES RRV_SAMPLES;
