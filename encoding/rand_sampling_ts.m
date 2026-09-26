function [X] = rand_sampling_ts(TRAJDB_DES, num_smp)
% RAND_SAMPLING_TS  Sample descriptor rows for codebook training.
%   TRAJDB_DES contains nonempty T-by-D matrices in a row cell array.
%   The count is rounded to an equal count per sequence, as in the original.
%   If needed, cycle the same random permutation to fill any requested count.

num_training = length(TRAJDB_DES);
if num_training == 0
    error('LRFb:emptyTrainingSet', 'At least one training sequence is required.');
end
if ~isscalar(num_smp) || ~isfinite(num_smp) || num_smp < 0 || num_smp ~= floor(num_smp)
    error('LRFb:sampleCount', 'The sample count must be a finite nonnegative integer.');
end
num_per_training = round(num_smp / num_training);
num_smp = num_per_training * num_training;
dimFea = size(TRAJDB_DES{1, 1}, 2);
X = zeros(dimFea, num_smp);
cnt = 0;
if num_per_training == 0
    return;
end
for ii = 1:num_training
    num_fea = size(TRAJDB_DES{1, ii}, 1);
    if num_fea == 0 || size(TRAJDB_DES{1, ii}, 2) ~= dimFea
        error('LRFb:trainingSequence', 'Training sequences must be nonempty with equal feature dimensions.');
    end
    rndidx = randperm(num_fea);
    if num_per_training > num_fea
        rndidx = repmat(rndidx, 1, ceil(num_per_training / num_fea));
    end
    X(:, cnt + 1:cnt + num_per_training) = TRAJDB_DES{1, ii}(rndidx(1:num_per_training), :)';
    cnt = cnt + num_per_training;
end
end
