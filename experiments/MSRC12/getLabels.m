function [trainGID, testGID] = getLabels()
% GETLABELS  Read gesture labels from the prepared DB.mat and SAMPLES.mat files.

%% for data from Gy
load DB;
load SAMPLES;
train_num = size(TRAJDB, 2);
test_num = size(TRAJSAMPLES, 2);
trainGID = zeros(train_num, 1);
testGID = zeros(test_num, 1);
for n = 1:size(TRAJDB, 2)
    trainGID(n) = TRAJDB{1, n}(1);
end
for n = 1:size(TRAJSAMPLES, 2)
    testGID(n) = TRAJSAMPLES{1, n}(1);
end
