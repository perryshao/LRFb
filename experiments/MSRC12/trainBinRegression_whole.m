function theta = trainBinRegression_whole(X, trainGID, lambda, jointNum)
% TRAINBINREGRESSION_WHOLE  Train squared-error regression on an in-memory feature matrix.
%   X is D-by-N and trainGID is N-by-1; returns a D-by-C weight matrix.
%   This legacy alternative is separate from the final linear SVM.

% Labels and training data
numTrain = length(trainGID);
classNum = length(unique(trainGID));
Y = zeros(numTrain, classNum);
for i = 1:classNum
    Y(trainGID == i, i) = 1;
end

%% Regularized squared-error regression.

% % Set Options
% options = optimset('GradObj', 'on', 'MaxIter', 400);
% preTheta = 0;
% D = size(X,1);
% initialTheta = zeros(D, classNum);
% % Optimize
% [theta, J, exit_flag] = ...
%     fminunc(@(t)(costFunctionReg(t, X, Y, lambda, classNum,jointNum, modulaNum,preTheta)), initialTheta, options);

% due to memory limitation, we train the regression model with a set of batch training.
D = size(X, 1);
N = size(X, 2);
J = D / jointNum;
% shuffle the training data before minibatch training
shuffleIndx = randperm(N);
Y = Y(shuffleIndx, :);
X = X(:, shuffleIndx);
batchTimes = 1;
batchsize = floor(N / batchTimes);
initialTheta = zeros(D * classNum, 1); % vec(W)
preTinitialTheta = zeros(D * classNum, 1);
%% refine the training
% load preTinitialTheta;
% initialTheta = preTinitialTheta;
% load initialTheta
iterNum = 10;
for iter = 1:iterNum
    for batchtimes = 1:batchTimes
        Xbatch = X(:, (batchtimes - 1) * batchsize + 1:batchtimes * batchsize);
        Ybatch = Y((batchtimes - 1) * batchsize + 1:batchtimes * batchsize, :);
        fprintf('Iteration %d Batch times %d\n', iter, batchtimes)
        % for multiple modality
        %         [initialTheta, J, c] = minimize(initialTheta, 'costFunctionReg', 50, Xbatch, Ybatch,
        %         lambda,classNum,jointNum, modalityNum,preTinitialTheta);
        % for single modality
        [initialTheta, ~, ~] = minimize(initialTheta, 'costFunc', 50, Xbatch, Ybatch, lambda, classNum, ...
            jointNum, preTinitialTheta);
        save initialTheta initialTheta;
        save iterNum iter;
    end
    %     Xbatch = X(:,batchtimes*batchsize+1:end);
    %     Ybatch = Y(batchtimes*batchsize+1:end,:);
    %     fprintf('Iteration %d Batch times %d \n',iter,batchtimes+1)
    %     [initialTheta, J, c] = minimize(initialTheta, 'costFuncRegMex', 50, Xbatch, Ybatch, lambda,classNum,jointNum,
    %     modalityNum, preTheta);
end

theta = initialTheta;
theta = reshape(theta, D, classNum);
