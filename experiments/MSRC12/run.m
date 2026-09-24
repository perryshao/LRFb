% RUN  Execute the MSRC12 LRFb experiment from this directory.
%   Run setup_path from the project root first. This script overwrites
%   intermediate MAT files in the working directory.
%   Active path: LRFb -> GMM/Fisher temporal pyramid -> linear SVM.
%   Commented blocks retain historical alternative experiments.

% modified to support running in occlusion situation -- modified by perry
clear all
%% define the experiment times and class numbers
EXPERIMENT_TIMES = 1;
experiment_num = EXPERIMENT_TIMES;
CLASS_NUM = 12; % class numbers for classification task
%% define the directory of c3d data and corresponding numbers of directories
distance_matrix_diff = cell(1, EXPERIMENT_TIMES);
confusion_matrix_sc = cell(1, EXPERIMENT_TIMES);
compu_time_svm = zeros(1, EXPERIMENT_TIMES);
compu_time_diff = zeros(1, EXPERIMENT_TIMES);
recog_ratio_sc = zeros(1, EXPERIMENT_TIMES);

%% define the joint name
RANK = '19';
RKNE = '18';
LANK = '15';
LKNE = '14';
LFWT = '13';
RFWT = '17';
LSHO = '5';
RSHO = '9';
LELB = '6';
LWRA = '8';
RELB = '10';
RWRA = '12';
STRN = '1';
HEAD = '4';
C7 = '3';
T10 = '2';
ENSEMBLE = 'ENSEMBLE';
%% load targets of joints
% joints_no = {LWRA,RWRA,LANK,RANK,...
%              LELB,RELB,LKNE,RKNE,STRN};
% joints_no = {LWRA,RWRA};
joints_no = {RANK, RKNE, LANK, LKNE, LFWT, RFWT, LSHO, RSHO, LELB, LWRA, RELB, RWRA, STRN, HEAD, C7, T10};
ROOT = ENSEMBLE;
esemble_no = {ENSEMBLE};
bodyJoints = {LWRA; RWRA};
shuffle_sort = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29;
                 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
%% load data initially for first running
load_data_bat(bodyJoints, shuffle_sort);
[trainGID, testGID] = getLabels();
% preprocess_bat(joints_no,0,0);
GeneSC(bodyJoints); % generate SC_DB and SC_SAMPLES RRV_DB and RRV_SAMPLES
ntotalbh = 3; % l = 0,1,2,3 (L=3) blocks are 2^(l)
numClusters = [64, 64];
batchsize = 256;
jointNum = length(bodyJoints);
% temporal pyramid based on pooling fisher codes
% [trainGID,testGID]=GeneFisherCodeJointPyramid(1, ntotalbh, numClusters(1),trainGID,testGID,batchsize,0);% encoding
% descriptors
%  [trainGID,testGID] = GeneFisherCodeJointPyramid_2mod(1, ntotalbh,numClusters,trainGID,testGID,batchsize,0);  %
%  encoding 2 modality of descriptors by fisher vectors
GeneFisherCodeJointPyramid_whole(jointNum, ntotalbh, numClusters(1), 0);

% classify gesture by bin regression of 1 feature modality
% lambda = [0.0,0.3];
% theta = trainBinRegression(trainGID,lambda,jointNum,batchsize);
% predict_label = predictBinRegression(theta,testGID,batchsize);
% load modelForTest; modelForTest.theta = theta; save modelForTest modelForTest; clear modelForTest;

% classify gesture by bin regression without batch training
% lambda = [0.0,0.3];
% load traindata.mat;theta = trainBinRegression_whole(traindata,trainGID,lambda,jointNum);clear traindata;
% load testdata.mat;predict_label = predictBinRegression_whole(testdata,theta);clear testdata;

% classify gesture by bin regression of 2 modality feature
% lambda = [0.0,0.3];
% theta = trainBinRegression_2mod(trainGID,lambda,jointNum,batchsize);
% predict_label = predictBinRegression_2mod(theta,testGID,batchsize);
% load modelForTest; modelForTest.theta = theta; save modelForTest modelForTest; clear modelForTest;

% classify gesture by svm
% 1st type of svm linear classifier
load traindata.mat;
model = svmtrain(trainGID, traindata', '-c 1 -g 2.79e-4 -t 0 -b 1');
clear traindata;
load testdata.mat;
[predict_label, accuracy, dec_values] = svmpredict(testGID, testdata', model, '-b 1');
clear testdata;
load modelForTest;
modelForTest.model = model;
save modelForTest modelForTest;
clear modelForTest;

% 2nd type of svm linear classifier
% lambda = 0.1;     % regularization parameter for w
% load traindata.mat;[w, b, class_name] = li2nsvm_multiclass_lbfgs(traindata',trainGID, lambda);clear traindata;
% load testdata.mat;[predict_label, ~] = li2nsvm_multiclass_fwd(testdata', w, b, class_name);clear testdata;

% classify gestures by random forest
% load traindata.mat;Random_Forest = Stochastic_Bosque(traindata',trainGID,'ntrees',200);clear traindata;
% load testdata.mat;[f_output,~]= eval_Stochastic_Bosque(testdata',Random_Forest);clear testdata;
% predict_label = f_output;

confusion_matrix = zeros(CLASS_NUM, CLASS_NUM);
for i = 1:CLASS_NUM
    for j = 1:CLASS_NUM
        confusion_matrix(i, j) = length(find(testGID == i & predict_label == j));
    end
end
recog_ratio_sc(experiment_num) = trace(confusion_matrix) / sum(confusion_matrix(:));
recog_ratio_final_sc = mean(recog_ratio_sc) %#ok<NOPTS>
confusion_matrix_sc{1, experiment_num} = confusion_matrix;
