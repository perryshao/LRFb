function predictLabel = predictBinRegression_2mod(theta, testGID, batchsize)
% PREDICTBINREGRESSION_2MOD  Predict regression labels from two aligned feature modalities.
%   Concatenates shape-context and RRV features before computing class scores.

data_folder1 = '/home/data/IPdataset/train_test_sc/'; % --for ubuntu
data_folder1 = 'C:/Users/perry/Desktop/MyDocuments/work/MicrosoftGestureEvaluatingCode/train_test_sc/';
data_folder2 = '/home/data/IPdataset/train_test_rrv/'; % --for ubuntu
data_folder2 = 'C:/Users/perry/Desktop/MyDocuments/work/MicrosoftGestureEvaluatingCode/train_test_rrv/';
N = length(testGID);
batchTimes = floor(N / batchsize);
predictLabel = zeros(N, 1);
for batchtimes = 1:batchTimes
    eval(['load ' data_folder1 strcat('testdata', num2str(batchtimes))]);
    eval(['Xbatch= ' strcat('testdata', num2str(batchtimes)) ';']);
    eval(['load ' data_folder2 strcat('testdata', num2str(batchtimes))]);
    eval(['Xbatch=[Xbatch;' strcat('testdata', num2str(batchtimes)) ']' ';']);
    fprintf('Testing Batch times %d\n', batchtimes)
    [~, predictLabel((batchtimes - 1) * batchsize + 1:batchtimes * batchsize)] = max(Xbatch' * theta, [], 2);
    eval(['clear ' strcat('testdata', num2str(batchtimes))]);
end
eval(['load ' data_folder1 strcat('testdata', num2str(batchTimes + 1))]);
eval(['Xbatch= ' strcat('testdata', num2str(batchTimes + 1)) ';']);
eval(['load ' data_folder2 strcat('testdata', num2str(batchTimes + 1))]);
eval(['Xbatch=[Xbatch;' strcat('testdata', num2str(batchTimes + 1)) ']' ';']);
fprintf('Testing Batch times %d \n', batchTimes + 1)
[~, predictLabel(batchTimes * batchsize + 1:end)] = max(Xbatch' * theta, [], 2);
eval(['clear ' strcat('testdata', num2str(batchTimes + 1))]);
