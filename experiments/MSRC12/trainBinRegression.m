function theta = trainBinRegression(trainGID, lambda, jointNum, batchsize)
% TRAINBINREGRESSION  Train the legacy regression classifier from numbered MAT batches.
%   trainGID contains class labels; lambda sets the regularization weights.
%   The batch directory must match the selected feature encoder.

%% define the folder where the training and test data are saved
% data_folder = '/home/data/IPdataset/train_test_rrv/'; %--for ubuntu
data_folder = 'C:/Users/perry/Desktop/MyDocuments/work/MicrosoftGestureEvaluatingCode/train_test_rrv/';
% Labels and training data
numTrain = length(trainGID);
classNum = length(unique(trainGID));
Y = zeros(numTrain, classNum);
for i = 1:classNum
    Y(trainGID == i, i) = 1;
end
% due to memory limitation, we train the regression model with a set of batch training.

batchTimes = floor(numTrain / batchsize);
eval(['load ' data_folder 'traindata1']);
D = size(traindata1, 1);
initialTheta = zeros(D * classNum, 1); % vec(W)
preTinitialTheta = zeros(D * classNum, 1);
clear traindata1;
modality = D;
%% refine the training
% load preTinitialTheta;
% initialTheta = preTinitialTheta;
% load initialTheta

iterNum = 40;
for iter = 1:iterNum
    shuffleIndx = randperm(batchTimes);
    for batchtimes = shuffleIndx
        eval(['load ' data_folder strcat('traindata', num2str(batchtimes))]);
        eval(['Xbatch= ' strcat('traindata', num2str(batchtimes)) ';']);
        Ybatch = Y((batchtimes - 1) * batchsize + 1:batchtimes * batchsize, :);
        fprintf('Iteration %d Batch %d\n', iter, batchtimes)
        % for multiple modality
        %         [initialTheta, J, c] = minimize(initialTheta, 'costFunctionReg', 50, Xbatch, Ybatch,
        %         lambda,classNum,jointNum, modalityNum,preTinitialTheta);
        % for single modality
        [initialTheta, ~, ~] = minimize(initialTheta, 'costFunc', 1, Xbatch, Ybatch, lambda, classNum, ...
            jointNum, preTinitialTheta);
        %         [initialTheta, ~, ~] = minimize(initialTheta, 'costFuncRegMultPartGp_v1_42', 1, Xbatch, Ybatch,
        %         lambda,classNum,modality);
        % [initialTheta, ~, ~] = minimize_fullbatch(initialTheta, 'costFuncRegMultPartGp_v2', 50, Xbatch, Ybatch, Y,
        % lambda,classNum,jointNum,batchTimes,preTinitialTheta);
        eval(['clear ' strcat('traindata', num2str(batchtimes))]);
    end
    eval(['load ' data_folder strcat('traindata', num2str(batchTimes + 1))]);
    eval(['Xbatch= ' strcat('traindata', num2str(batchTimes + 1)) ';']);
    Ybatch = Y(batchTimes * batchsize + 1:end, :);
    fprintf('Iteration %d Batch %d \n', iter, batchTimes + 1)
    [initialTheta, ~, ~] = minimize(initialTheta, 'costFunc', 1, Xbatch, Ybatch, lambda, classNum, jointNum, ...
        preTinitialTheta);
    %     [initialTheta, ~, ~] = minimize(initialTheta, 'costFuncRegMultPartGp_v1_42', 1, Xbatch, Ybatch,
    %     lambda,classNum,modality);
    % [initialTheta, ~, ~] = minimize_fullbatch(initialTheta, 'costFuncRegMultPartGp_v2', 50, Xbatch, Ybatch,
    % lambda,classNum,jointNum,batchTimes,preTinitialTheta);
    eval(['clear ' strcat('traindata', num2str(batchTimes + 1))]);
end
save initialTheta initialTheta -v7.3;
theta = initialTheta;
theta = reshape(theta, D, classNum);
