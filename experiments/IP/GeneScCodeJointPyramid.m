function [trainGID, testGID] = GeneScCodeJointPyramid(jointNum, ntotalbh, nBases, trainGID, testGID, batchsize, pcaFlag)
% GENESCCODEJOINTPYRAMID  Pool sparse codes over temporal pyramid bins.
%   Legacy experiment using dictionaries in Results/reg_sc_b256_*.mat.
%   Keep batch paths consistent with the selected training/prediction helpers.

%% define the folder where the training and test data are saved
data_folder = '/home/data/IPdataset/train_test_Sc/'; % --for ubuntu
% data_folder = 'C:/Users/perry/Desktop/MyDocuments/work/IPEvaluatingCode/train_test/';
%% collect the visual words
%% define the folder where the training and test data are saved
load SC_DB;
samples_r = size(SC_DB, 2);

% dictionary training for sparse coding
% nBases =128;
nsmp = 20000;
beta = 1e-5; % a small regularization for stabilizing sparse coding
num_iters = 50;

% feature pooling parameters
pyramid = 2.^(0:ntotalbh); % spatial block number on each level of the pyramid
gamma = 0.15;
% knn = 200;
knn = 0; % find the k-nearest neighbors for approximate sparse coding
% if set 0, use the standard sparse coding
jointCodeLength = nBases * sum(pyramid);
%% learning sparse coding dictionary
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
    X = [X currentX{i}];
end
clear currentX emptyIndx;

% do pca on X first;
if pcaFlag
    [coeff, ~, latent, ~, ~] = pca(X');
    PcaM = coeff(:, cumsum(latent) / sum(latent) < 0.99);
    X = PcaM' * X;
    X = X ./ sqrt(repmat(latent(1:size(PcaM, 2)), 1, lastNsmp)); % whiten the pca
else
    PcaM = zeros(3, 3);
end

nsmp = size(X, 2); % recompute the nsmp after sampling
% batch_size = floor(nsmp/1); % batch size when learning sparse codes
batch_size = 512;

% check if there exist a mat file for bases
if exist('Results/reg_sc_b256_20171108T110255.mat', 'file')
    load Results/reg_sc_b256_20171108T110255.mat
    clear X;
else
    [B, S, stat] = reg_sparse_coding(X, nBases, eye(nBases), beta, gamma, num_iters, batch_size);
    clear X;
end

%  [B, S, stat] = reg_sparse_coding(X, nBases, eye(nBases), beta, gamma, num_iters,batch_size);
%  clear X;

%% save key parameters for fisher vector encoding.
modelForTest.B = B;
modelForTest.PcaM = PcaM;

save modelForTest modelForTest;
clear modelForTest;

%% calculate the sparse coding feature

disp('==================================================');
fprintf('Calculating the sparse coding feature...\n');
fprintf('Regularization parameter: %f\n', gamma);
disp('==================================================');

% %% for training data
% % shuffle the training data before doing batching
% %shuffleIndx = randperm(samples_r);
% %trainGID = trainGID(shuffleIndx);
% %SC_DB = SC_DB(1,shuffleIndx);

% batchTimes = floor(samples_r/batchsize);
% traindata=zeros(jointCodeLength*jointNum,batchsize);
% for batchtimes = 1:batchTimes
%     batch_num = 1;
%     for iter1 = (batchtimes-1)*batchsize+1:batchtimes*batchsize
%         fprintf ('computing and pooling Fisher Vectors for training data %d...\n',iter1);
%         featsLength = size(SC_DB{1,iter1},1);
%         frameLength = floor(featsLength/jointNum);
%         for m = 1:jointNum
%             feats = SC_DB{1,iter1}((m-1)*frameLength+1:m*frameLength,:);
%             if pcaFlag
%                 feats = feats*PcaM;% dimension reduction of pca;
%                 %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
%             end
%             if knn,
%               traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_approx_pooling_ts(feats', B, pyramid,
%               gamma, knn);
%             else
%               traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_pooling_ts(feats', B, pyramid, gamma);
%             end
%         end
%   batch_num = batch_num+1;
%     end
%     eval([strcat('traindata',num2str(batchtimes)) '=traindata;']);
%     eval(['save ' data_folder strcat('traindata',num2str(batchtimes)) ' ' strcat('traindata',num2str(batchtimes)) '
%     -v7.3;']);
%     eval(['clear ' strcat('traindata',num2str(batchtimes))]);
% end
% clear traindata;

% traindata=zeros(jointCodeLength*jointNum,samples_r-batchtimes*batchsize);
% batch_num = 1;
% for iter1 = batchtimes*batchsize+1:samples_r
%         fprintf ('computing and pooling Fisher Vectors for training data %d...\n',iter1);
%         featsLength = size(SC_DB{1,iter1},1);
%         frameLength = floor(featsLength/jointNum);
%         for m = 1:jointNum
%             feats = SC_DB{1,iter1}((m-1)*frameLength+1:m*frameLength,:);
%             if pcaFlag
%                 feats = feats*PcaM;% dimension reduction of pca;
%                 %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
%             end
%             if knn,
%               traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_approx_pooling_ts(feats', B, pyramid,
%               gamma, knn);
%             else
%               traindata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_pooling_ts(feats', B, pyramid, gamma);
%             end
%         end
%   batch_num = batch_num+1;
% end
% eval([strcat('traindata',num2str(batchtimes+1)) '=traindata;']);
% eval(['save ' data_folder  strcat('traindata',num2str(batchtimes+1)) ' ' strcat('traindata',num2str(batchtimes+1)) '
% -v7.3;']);
% eval(['clear ' strcat('traindata',num2str(batchtimes+1))]);
% clear SC_DB;

% %% for test data
% % shuffle the training data before doing batching
% %shuffleIndx = randperm(samples_t);
% %testGID = testGID(shuffleIndx);
% %SC_SAMPLES = SC_SAMPLES(1,shuffleIndx);

% load SC_SAMPLES;
% samples_t = size(SC_SAMPLES,2);
% batchTimes = floor(samples_t/batchsize);
% testdata=zeros(jointCodeLength*jointNum,batchsize);
% for batchtimes = 1:batchTimes
%     batch_num = 1;
%     for iter1 = (batchtimes-1)*batchsize+1:batchtimes*batchsize
%         fprintf ('computing and pooling Fisher Vectors for testing data %d...\n',iter1);
%         featsLength = size(SC_SAMPLES{1,iter1},1);
%         frameLength = floor(featsLength/jointNum);
%         for m = 1:jointNum
%             feats = SC_SAMPLES{1,iter1}((m-1)*frameLength+1:m*frameLength,:);
%             if pcaFlag
%                 feats = feats*PcaM;% dimension reduction of pca;
%                 %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
%             end
%             if knn,
%               testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_approx_pooling_ts(feats', B, pyramid,
%               gamma, knn);
%             else
%               testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_pooling_ts(feats', B, pyramid, gamma);
%             end
%         end
%   batch_num = batch_num+1;
%     end
%     eval([strcat('testdata',num2str(batchtimes)) '=testdata;']);
%     eval(['save ' data_folder strcat('testdata',num2str(batchtimes)) ' ' strcat('testdata',num2str(batchtimes)) '
%     -v7.3;']);
%     eval(['clear ' strcat('testdata',num2str(batchtimes))]);
% end
% clear testdata;

% testdata=zeros(jointCodeLength*jointNum,samples_t-batchtimes*batchsize);
% batch_num = 1;
% for iter1 = batchtimes*batchsize+1:samples_t
%         fprintf ('computing and pooling Fisher Vectors for testing data %d...\n',iter1);
%         featsLength = size(SC_SAMPLES{1,iter1},1);
%         frameLength = floor(featsLength/jointNum);
%         for m = 1:jointNum
%             feats = SC_SAMPLES{1,iter1}((m-1)*frameLength+1:m*frameLength,:);
%             if pcaFlag
%                 feats = feats*PcaM;% dimension reduction of pca;
%                 %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
%             end
%             if knn,
%               testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_approx_pooling_ts(feats', B, pyramid,
%               gamma, knn);
%             else
%               testdata((m-1)*jointCodeLength+1:m*jointCodeLength, iter1) = sc_pooling_ts(feats', B, pyramid, gamma);
%             end
%         end
%   batch_num = batch_num+1;
% end
% clear SC_SAMPLES;
% eval([strcat('testdata',num2str(batchtimes+1)) '=testdata;']);
% eval(['save ' data_folder strcat('testdata',num2str(batchtimes+1)) ' ' strcat('testdata',num2str(batchtimes+1)) '
% -v7.3;']);
% eval(['clear ' strcat('testdata',num2str(batchtimes+1))]);

%% sparse coding on whole dataset
traindata = zeros(jointCodeLength * jointNum, samples_r);
for iter1 = 1:samples_r
    fprintf ('computing and pooling sparce codes for training data %d...\n', iter1);
    featsLength = size(SC_DB{1, iter1}, 1);
    frameLength = floor(featsLength / jointNum);
    for m = 1:jointNum
        feats = SC_DB{1, iter1}((m - 1) * frameLength + 1:m * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        if knn
            traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, iter1) = sc_approx_pooling_ts(feats', ...
                B, pyramid, gamma, knn);
        else
            traindata((m - 1) * jointCodeLength + 1:m * jointCodeLength, iter1) = sc_pooling_ts(feats', B, ...
                pyramid, gamma);
        end
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
    fprintf ('computing and pooling sparce codes for test data %d...\n', iter2);
    featsLength = size(SC_SAMPLES{1, iter2}, 1);
    frameLength = floor(featsLength / jointNum);
    for n = 1:jointNum
        feats = SC_SAMPLES{1, iter2}((n - 1) * frameLength + 1:n * frameLength, :);
        if pcaFlag
            feats = feats * PcaM; % dimension reduction of pca;
            %             feats= feats./sqrt(repmat(latent(1:size(PcaM,2)),1,frameLength))'; % whiten the pca
        end
        if knn
            testdata((n - 1) * jointCodeLength + 1:n * jointCodeLength, iter2) = sc_approx_pooling_ts(feats', ...
                B, pyramid, gamma, knn);
        else
            testdata((n - 1) * jointCodeLength + 1:n * jointCodeLength, iter2) = sc_pooling_ts(feats', B, ...
                pyramid, gamma);
        end
    end
end
sum_ScSPM_time = toc;
save testdata testdata -v7.3;
