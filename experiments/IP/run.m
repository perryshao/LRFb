% RUN  Execute the IP LRFb experiment from this directory.
%   Run setup_path from the project root first. This script overwrites
%   intermediate MAT files in the working directory.
%   Active path: LRFb -> GMM/Fisher temporal pyramid -> linear SVM.
%   Commented blocks retain historical alternative experiments.

clear all
%% define the experiment times and class numbers
EXPERIMENT_TIMES = 1;
experiment_num = EXPERIMENT_TIMES;
CLASS_NUM = 16; % class numbers for classification task
% define the directory of c3d data and corresponding numbers of directories
% BAT_FOLDER = 'NTU3D_skeletons/';  %LOCATION OF SKELETON FILES --- for windows
BAT_FOLDER = '/home/data/IPdataset/database/'; % LOCATION OF SKELETON FILES -- for ubuntu
% BAT_FOLDER = 'IPdataset/interactplay/database/'; % for windows
file_ext = '.txt';

distance_matrix_3sc = cell(1, EXPERIMENT_TIMES);
confusion_matrix_3sc = cell(1, EXPERIMENT_TIMES);
distance_matrix_diff = cell(1, EXPERIMENT_TIMES);
confusion_matrix_diff = cell(1, EXPERIMENT_TIMES);
compu_time_svm = zeros(1, EXPERIMENT_TIMES);
compu_time_diff = zeros(1, EXPERIMENT_TIMES);
recog_ratio_diff = zeros(1, EXPERIMENT_TIMES);
recog_ratio_smml = zeros(1, EXPERIMENT_TIMES);
% define the joint name
HEAD = '1';
HAND_L = '2';
HAND_R = '3';
RIST = '4';
% joints number in the IP dataset
joints_no = {HEAD; HAND_L; HAND_R; RIST};
bodyJoints = {HAND_L; HAND_R};
%% load the data from database
load_IPtxt_bat(BAT_FOLDER);
[trainGID, testGID] = getLabels();
GeneSC(bodyJoints); % generate SC_DB and SC_SAMPLES
ntotalbh = 3; % l = 0,1,2,3 (L=3) blocks are 2^(l)
numClusters = [64, 64];
% nBases = 256;% sparse coding base number
batchsize = 256;
jointNum = length(bodyJoints);
% temporal pyramid based on pooling fisher codes
% [trainGID,testGID] = GeneFisherCodeJointPyramid(1, ntotalbh, numClusters,trainGID,testGID,batchsize,0);% encoding
% descriptors by fisher vectors
% [trainGID,testGID] = GeneShapeContextJointPyramid(1, ntotalbh,trainGID,testGID,0);% doing pyramid on shape context
% descriptors
%  [trainGID,testGID]=GeneScFisherCodeJointPyramid(1,ntotalbh,nBases,trainGID,testGID,batchsize,0);% encoding
%  descriptors by sparse coding based fisher vectors
% [trainGID,testGID] = GeneFisherCodeJointPyramid_2mod(1, ntotalbh,numClusters,trainGID,testGID,batchsize,0); % encoding
% 2 modality of descriptors by fisher vectors
GeneFisherCodeJointPyramid_whole(jointNum, ntotalbh, numClusters(1), 0);
% [trainGID,testGID] = GeneScCodeJointPyramid(jointNum,ntotalbh,nBases,trainGID,testGID,batchsize,0);
save trainGID trainGID;
save testGID testGID;
numTrain = length(trainGID);
numTest = length(testGID);

% classify gesture by a single modality feature
% load trainGID;load testGID;
% lambda = [0.0,0.3];
% theta = trainBinRegression(trainGID,lambda,1,batchsize);
% predict_label = predictBinRegression(theta,testGID,batchsize);
% load modelForTest; modelForTest.theta = theta; save modelForTest modelForTest; clear modelForTest;

% classify gesture by a single modality feature with non-batch training
lambda = [0.0, 0.5];
load traindata.mat;
theta = trainBinRegression_whole(traindata, trainGID, lambda, jointNum);
clear traindata;
load testdata.mat;
predict_label = predictBinRegression_whole(testdata, theta);
clear testdata;
load modelForTest;
modelForTest.theta = theta;
save modelForTest modelForTest;
clear modelForTest;

% classify gesture by bin regression of 2 modality feature
% load trainGID;load testGID;
% lambda = [0.0,0.3];
% theta = trainBinRegression_2mod(trainGID,lambda,1,batchsize);
% predict_label = predictBinRegression_2mod(theta,testGID,batchsize);
% load modelForTest; modelForTest.theta = theta; save modelForTest modelForTest; clear modelForTest;

% %1st type of svm linear classifier
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

% tic;
% load testdata.mat;[predict_label, ~] = li2nsvm_multiclass_fwd(testdata', w, b, class_name);clear testdata;
% sum_time_svm = toc;
% compu_time_ssm(experiment_num) =sum_time_svm/length(predict_label);

% Random Forest classifier
% load traindata.mat;Random_Forest = Stochastic_Bosque(traindata',trainGID,'ntrees',200);clear traindata;
% load testdata.mat;[f_output, f_votes]= eval_Stochastic_Bosque(testdata',Random_Forest);clear testdata;
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
