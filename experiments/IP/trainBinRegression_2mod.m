function theta = trainBinRegression_2mod(trainGID, lambda, jointNum, batchsize)
% TRAINBINREGRESSION_2MOD  Train regression on concatenated shape-context and RRV batches.
%   The two directories must contain aligned batches and training labels.

%% define the folder where the training and test data are saved
data_folder1 = '/home/data/IPdataset/train_test_sc/'; % --for ubuntu
% data_folder1 = 'C:/Users/perry/Desktop/MyDocuments/work/IPEvaluatingCode/train_test_cs/';
data_folder2 = '/home/data/IPdataset/train_test_rrv/'; % --for ubuntu
% data_folder2 = 'C:/Users/perry/Desktop/MyDocuments/work/IPEvaluatingCode/train_test_rrv/';
% Labels and training data
numTrain = length(trainGID);
classNum = length(unique(trainGID));
Y = zeros(numTrain, classNum);
for i = 1:classNum
    Y(trainGID == i, i) = 1;
end
% due to memory limitation, we train the regression model with a set of batch training.

batchTimes = floor(numTrain / batchsize);
eval(['load ' data_folder1 'traindata1']);
D1 = size(traindata1, 1);
eval(['load ' data_folder2 'traindata1']);
D2 = size(traindata1, 1);
modality = [D1, D2];
initialTheta = zeros((D1 + D2) * classNum, 1); % vec(W)
clear traindata1;
%% refine the training
% load preTinitialTheta;
% initialTheta = preTinitialTheta;
% load initialTheta
iterNum = 40;
for iter = 1:iterNum
    shuffleIndx = randperm(batchTimes);
    for batchtimes = shuffleIndx
        eval(['load ' data_folder1 strcat('traindata', num2str(batchtimes))]);
        eval(['Xbatch= ' strcat('traindata', num2str(batchtimes)) ';']);
        eval(['load ' data_folder2 strcat('traindata', num2str(batchtimes))]);
        eval(['Xbatch=[Xbatch;' strcat('traindata', num2str(batchtimes)) ']' ';']);
        Ybatch = Y((batchtimes - 1) * batchsize + 1:batchtimes * batchsize, :);
        fprintf('Iteration %d Batch %d\n', iter, batchtimes)
        % for multiple modality
        %         [initialTheta, J, c] = minimize(initialTheta, 'costFunctionReg', 50, Xbatch, Ybatch,
        %         lambda,classNum,jointNum, modalityNum,preTinitialTheta);
        % for single modality
        [initialTheta, ~, ~] = minimize(initialTheta, 'costFuncRegMultPartGp_v1_42', 1, Xbatch, Ybatch, ...
            lambda, classNum, modality);
        % [initialTheta, ~, ~] = minimize_fullbatch(initialTheta, 'costFuncRegMultPartGp_v1', 50, Xbatch, Ybatch, Y,
        % lambda,classNum,jointNum,batchTimes,preTinitialTheta);
        eval(['clear ' strcat('traindata', num2str(batchtimes))]);
    end
    eval(['load ' data_folder1 strcat('traindata', num2str(batchTimes + 1))]);
    eval(['Xbatch= ' strcat('traindata', num2str(batchTimes + 1)) ';']);
    eval(['load ' data_folder2 strcat('traindata', num2str(batchTimes + 1))]);
    eval(['Xbatch=[Xbatch;' strcat('traindata', num2str(batchTimes + 1)) ']' ';']);
    Ybatch = Y(batchTimes * batchsize + 1:end, :);
    fprintf('Iteration %d Batch %d \n', iter, batchTimes + 1)
    [initialTheta, ~, ~] = minimize(initialTheta, 'costFuncRegMultPartGp_v1_42', 1, Xbatch, Ybatch, lambda, ...
        classNum, modality);
    % [initialTheta, ~, ~] = minimize_fullbatch(initialTheta, 'costFuncRegMultPartGp_v2', 50, Xbatch, Ybatch,
    % lambda,classNum,jointNum,batchTimes,preTinitialTheta);
    eval(['clear ' strcat('traindata', num2str(batchTimes + 1))]);
end
save initialTheta initialTheta -v7.3;
theta = initialTheta;
theta = reshape(theta, [], classNum);
