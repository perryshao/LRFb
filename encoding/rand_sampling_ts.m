function [X] = rand_sampling_ts(TRAJDB_DES, num_smp)
% RAND_SAMPLING_TS  Sample descriptor rows for codebook training.
%   TRAJDB_DES contains T-by-D matrices; X contains sampled D-by-N columns.
%   The requested count is rounded to an equal count per sequence.

num_training = length(TRAJDB_DES); % num of images
num_per_training = round(num_smp / num_training);
num_smp = num_per_training * num_training;
dimFea = size(TRAJDB_DES{1, 1}, 2);

X = zeros(dimFea, num_smp);
cnt = 0;

for ii = 1:num_training
    num_fea = size(TRAJDB_DES{1, ii}, 1);
    rndidx = randperm(num_fea);
    if num_per_training > max(rndidx)
        rndidx = cat(2, rndidx, rndidx(1:num_per_training - max(rndidx)));
    end
    X(:, cnt + 1:cnt + num_per_training) = TRAJDB_DES{1, ii}(rndidx(1:num_per_training), :)';
    cnt = cnt + num_per_training;
end;
