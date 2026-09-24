function [trainGID, testGID] = GeneShapeContextJointPyramid(jointNum, ntotalbh, trainGID, testGID, pcaFlag)
% GENESHAPECONTEXTJOINTPYRAMID  Pool shape-context descriptors over a temporal pyramid.
%   Reads SC_DB.mat / SC_SAMPLES.mat and writes traindata.mat / testdata.mat.
%   Uses pv_pooling_ts, which sums descriptors within each temporal bin.

%% define the folder where the training and test data are saved
data_folder = '/home/data/IPdataset/train_test/'; % --for ubuntu
% data_folder = 'C:/Users/perry/Desktop/MyDocuments/work/IPEvaluatingCode/train_test/';

%% collect samples for PCA
load SC_DB;

% sample number of PCA
nsmp = 20000;
% feature pooling parameters
pyramid = 2.^(0:ntotalbh); % spatial block number on each level of the pyramid

currentTime = 0;
lastNsmp = 0;
% to avoid all(0) feature vector
while lastNsmp < nsmp
    currentTime = currentTime + 1;
    % randomly selecting local training features
    currentX{currentTime} = rand_sampling_ts(SC_DB, nsmp);
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

%% save key parameters for fisher vector encoding.
modelForTest.PcaM = PcaM;
save modelForTest modelForTest;
clear modelForTest;

Length_V = size(SC_DB{1}, 2);
samples_r = length(SC_DB);
jointCodeLength = Length_V * sum(pyramid);
%% calculate the pyramid vector

disp('==================================================');
fprintf('Calculating the pyramid vector...\n');
disp('==================================================');

%% for training data
% shuffle the training data before doing batching
% shuffleIndx = randperm(samples_r);
% trainGID = trainGID(shuffleIndx);
% SC_DB = SC_DB(1,shuffleIndx);

traindata = zeros(jointCodeLength * jointNum, samples_r);
for iter1 = 1:samples_r
    fprintf ('computing and pooling shape context Vectors for training data %d...\n', iter1);
    featsLength = size(SC_DB{1, iter1}, 1);
    frameLength = floor(featsLength / jointNum);
    for m = 1:jointNum
        feats = SC_DB{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, iter1) = pv_pooling_ts(feats', pyramid);
    end
end
save traindata traindata -v7.3;
clear traindata;
clear SC_DB;

load SC_SAMPLES;
samples_t = size(SC_SAMPLES, 2);
testdata = zeros(jointCodeLength * jointNum, samples_t);
tic;
for iter2 = 1:samples_t
    fprintf ('computing and pooling shape context Vectors for test data %d...\n', iter2);
    featsLength = size(SC_SAMPLES{1, iter2}, 1);
    frameLength = floor(featsLength / jointNum);
    for n = 1:jointNum
        feats = SC_SAMPLES{1, iter2}((n - 1) * frameLength + 1:n * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        %         testdata((n-1)*jointCodeLength+1:n*jointCodeLength, iter2) = vl_fisher(feats', means, covariances,
        %         priors,'Improved');
        testdata((n - 1) * jointCodeLength + 1:n * jointCodeLength, iter2) = pv_pooling_ts(feats', pyramid);
    end
end
sum_FvSPM_time = toc;
save testdata testdata -v7.3;
