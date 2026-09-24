function [trainGID, testGID] = GeneFisherCodeJointPyramid(jointNum, ntotalbh, numClusters, trainGID, testGID, ...
    batchsize, pcaFlag)
% GENEFISHERCODEJOINTPYRAMID  Encode descriptor datasets with a GMM and temporal Fisher pyramid.
%   Legacy batched variant; reads descriptor MAT files in the working folder.
%   Writes numbered training/test batches to the data_folder below.
%   Returns labels in the shuffled order used for the training batches.

%% define the folder where the training and test data are saved
data_folder = '/home/data/IPdataset/train_test_rrv/'; % --for ubuntu
data_folder = 'C:/Users/perry/Desktop/MyDocuments/work/MicrosoftGestureEvaluatingCode/train_test_rrv/';

%% collect the visual words
load RRV_DB;
samples_r = size(RRV_DB, 2);

% parameters for GMM learning
nsmp = 20000;

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

%% for training data
% shuffle the training data before doing batching
shuffleIndx_train = randperm(samples_r);
trainGID = trainGID(shuffleIndx_train);
RRV_DB = RRV_DB(1, shuffleIndx_train);

batchTimes = floor(samples_r / batchsize);
traindata = zeros(jointCodeLength * jointNum, batchsize);
for batchtimes = 1:batchTimes
    batch_num = 1;
    for iter1 = (batchtimes - 1) * batchsize + 1:batchtimes * batchsize
        fprintf ('computing and pooling Fisher Vectors for training data %d...\n', iter1);
        featsLength = size(RRV_DB{1, iter1}, 1);
        frameLength = floor(featsLength / jointNum);
        for m = 1:jointNum
            feats = RRV_DB{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
            if pcaFlag
                feats = feats * PcaM; % dimension reduction of pca;
                %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
            end
            %         traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = vl_fisher(feats', means,
            %         covariances, priors,'Improved');
            traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, batch_num) = fv_pooling_ts(feats', ...
                means, covariances, priors, 'Improved', pyramid);
        end
        batch_num = batch_num + 1;
    end
    eval([strcat('traindata', num2str(batchtimes)) '=traindata;']);
    eval(['save ' data_folder strcat('traindata', num2str(batchtimes)) ' ' strcat('traindata', ...
        num2str(batchtimes)) ' -v7.3;']);
    eval(['clear ' strcat('traindata', num2str(batchtimes))]);
end
clear traindata;

traindata = zeros(jointCodeLength * jointNum, samples_r - batchtimes * batchsize);
batch_num = 1;
for iter1 = batchtimes * batchsize + 1:samples_r
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
        traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, batch_num) = fv_pooling_ts(feats', means, ...
            covariances, priors, 'Improved', pyramid);
    end
    batch_num = batch_num + 1;
end
eval([strcat('traindata', num2str(batchtimes + 1)) '=traindata;']);
eval(['save ' data_folder strcat('traindata', num2str(batchtimes + 1)) ' ' strcat('traindata', ...
    num2str(batchtimes + 1)) ' -v7.3;']);
eval(['clear ' strcat('traindata', num2str(batchtimes + 1))]);
clear RRV_DB;
%% for test data
load RRV_SAMPLES;
samples_t = size(RRV_SAMPLES, 2);

% shuffle the testing data before doing batching
shuffleIndx_test = randperm(samples_t);
testGID = testGID(shuffleIndx_test);
RRV_SAMPLES = RRV_SAMPLES(1, shuffleIndx_test);

batchTimes = floor(samples_t / batchsize);
testdata = zeros(jointCodeLength * jointNum, batchsize);
for batchtimes = 1:batchTimes
    batch_num = 1;
    for iter1 = (batchtimes - 1) * batchsize + 1:batchtimes * batchsize
        fprintf ('computing and pooling Fisher Vectors for testing data %d...\n', iter1);
        featsLength = size(RRV_SAMPLES{1, iter1}, 1);
        frameLength = floor(featsLength / jointNum);
        for m = 1:jointNum
            feats = RRV_SAMPLES{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
            if pcaFlag
                feats = feats * PcaM; % dimension reduction of pca;
                %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
            end
            %         testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = vl_fisher(feats', means, covariances,
            %         priors,'Improved');
            testdata((m - 1) * jointCodeLength + 1:m * jointCodeLength, batch_num) = fv_pooling_ts(feats', ...
                means, covariances, priors, 'Improved', pyramid);
        end
        batch_num = batch_num + 1;
    end
    eval([strcat('testdata', num2str(batchtimes)) '=testdata;']);
    eval(['save ' data_folder strcat('testdata', num2str(batchtimes)) ' ' strcat('testdata', ...
        num2str(batchtimes)) ' -v7.3;']);
    eval(['clear ' strcat('testdata', num2str(batchtimes))]);
end
clear testdata;

testdata = zeros(jointCodeLength * jointNum, samples_t - batchtimes * batchsize);
batch_num = 1;
for iter1 = batchtimes * batchsize + 1:samples_t
    fprintf ('computing and pooling Fisher Vectors for testing data %d...\n', iter1);
    featsLength = size(RRV_SAMPLES{1, iter1}, 1);
    frameLength = floor(featsLength / jointNum);
    for m = 1:jointNum
        feats = RRV_SAMPLES{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        %         testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = vl_fisher(feats', means, covariances,
        %         priors,'Improved');
        testdata((m - 1) * jointCodeLength + 1:m * jointCodeLength, batch_num) = fv_pooling_ts(feats', means, ...
            covariances, priors, 'Improved', pyramid);
    end
    batch_num = batch_num + 1;
end
clear RRV_SAMPLES;
eval([strcat('testdata', num2str(batchtimes + 1)) '=testdata;']);
eval(['save ' data_folder strcat('testdata', num2str(batchtimes + 1)) ' ' strcat('testdata', ...
    num2str(batchtimes + 1)) ' -v7.3;']);
eval(['clear ' strcat('testdata', num2str(batchtimes + 1))]);
