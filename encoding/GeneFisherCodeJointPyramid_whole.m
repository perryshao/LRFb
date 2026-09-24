function GeneFisherCodeJointPyramid_whole(jointNum, ntotalbh, numClusters, pcaFlag)
% GENEFISHERCODEJOINTPYRAMID_WHOLE  Learn a GMM and encode complete LRFb datasets.
%   Reads RRV_DB.mat / RRV_SAMPLES.mat (cells of T-by-D descriptors).
%   Writes modelForTest.mat, traindata.mat and testdata.mat.
%   jointNum partitions descriptor rows; ntotalbh sets levels 0:ntotalbh.
%   Historical PCA whitening and row partitioning are preserved.

%% collect the visual words
load RRV_DB;
samples_r = size(RRV_DB, 2);

% parameters for GMM learning
nsmp = 50000;

% feature pooling parameters
pyramid = 2.^(0:ntotalbh); % spatial block number on each level of the pyramid

%% Learn the Gaussian mixture model.
currentTime = 0;
lastNsmp = 0;
% to avoid all(0) feature vector
while lastNsmp < nsmp
    currentTime = currentTime + 1;
    % randomly selecting local training features
    currentX{currentTime} = rand_sampling_ts(RRV_DB, nsmp);
    currentNsmp = size(currentX{currentTime}, 2); % recompute the nsmp after sampling
    emptyIndx = zeros(1, currentNsmp);
    for i = 1:currentNsmp
        if ~any(currentX{currentTime}(:, i))
            emptyIndx(i) = i;
        end
    end
    emptyIndx(emptyIndx == 0) = [];
    currentX{currentTime}(:, emptyIndx) = [];
    lastNsmp = lastNsmp + size(currentX{currentTime}, 2); % recompute the nsmp after sampling
end
X = [];
for i = 1:currentTime
    % X = [X currentX{currentTime}]; %big debug found by perry on 1/6/16
    X = [X currentX{i}];
end
clear currentX emptyIndx;

% X = rand_sampling_ts(DB_DESCRIPTORS, nsmp);

% do pca on X first;
if pcaFlag
    [coeff, ~, latent, ~, ~] = pca(X');
    PcaM = coeff(:, cumsum(latent) / sum(latent) < 0.98);
    X = PcaM' * X;
    X = X ./ sqrt(repmat(latent(1:size(PcaM, 2)), 1, lastNsmp)); % whiten the pca
else
    PcaM = zeros(3, 3);
end

[means, covariances, priors] = vl_gmm(X, numClusters);
clear X;

%% save key parameters for fisher vector encoding.
modelForTest.covariances = covariances;
modelForTest.means = means;
modelForTest.priors = priors;
modelForTest.PcaM = PcaM;
save modelForTest modelForTest;
clear modelForTest;
Length_fisherV = size(covariances, 1) * size(covariances, 2) * 2;
jointCodeLength = Length_fisherV * sum(pyramid);
% jointCodeLength = Length_fisherV;
%% Encode descriptors with the temporal Fisher pyramid.

disp('==================================================');
fprintf('Calculating the fisher vector...\n');
disp('==================================================');

traindata = zeros(jointCodeLength * jointNum, samples_r);
for iter1 = 1:samples_r
    fprintf ('computing and pooling Fisher Vectors for training data %d...\n', iter1);
    featsLength = size(RRV_DB{1, iter1}, 1);
    frameLength = floor(featsLength / jointNum);
    for m = 1:jointNum
        feats = RRV_DB{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        %         traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = vl_fisher(feats', means, covariances,
        %         priors,'Improved');
        traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, iter1) = fv_pooling_ts(feats', means, ...
            covariances, priors, 'Improved', pyramid);
    end
end
save traindata traindata -V7.3;
clear traindata;
clear RRV_DB;
load RRV_SAMPLES;
samples_t = size(RRV_SAMPLES, 2);
testdata = zeros(jointCodeLength * jointNum, samples_t);
tic;
for iter2 = 1:samples_t
    fprintf ('computing and pooling Fisher Vectors for test data %d...\n', iter2);
    featsLength = size(RRV_SAMPLES{1, iter2}, 1);
    frameLength = floor(featsLength / jointNum);
    for n = 1:jointNum
        feats = RRV_SAMPLES{1, iter2}((n - 1) * frameLength + 1:n * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        %         testdata((n-1)*jointCodeLength+1:n*jointCodeLength, iter2) = vl_fisher(feats', means, covariances,
        %         priors,'Improved');
        testdata((n - 1) * jointCodeLength + 1:n * jointCodeLength, iter2) = fv_pooling_ts(feats', means, ...
            covariances, priors, 'Improved', pyramid);
    end
end
sum_FvSPM_time = toc;
save testdata testdata -V7.3;
