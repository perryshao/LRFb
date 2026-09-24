function [beta] = scfv_pooling_ts(feaSet, B, pyramid, gamma)
% SCFV_POOLING_TS  Compute the legacy sparse-coding Fisher pooling variant.
%   feaSet is D-by-T; B is the dictionary; gamma controls sparsity.
%   The atom/time indexing issue is documented in docs/CODE_REVIEW.md.
%
%   Adapted from pooling code by Jianchao Yang, NEC Research Lab America.
%   Mentor: Kai Yu. July 2008; revised May 2010.

dSize = size(B, 2);
nSmp = size(feaSet, 2);
sc_codes = zeros(dSize, nSmp);

% compute the local feature for each local feature
beta = 1e-4;
A = B' * B + 2 * beta * eye(dSize);
Q = -B' * feaSet;

for iter1 = 1:nSmp
    sc_codes(:, iter1) = L1QP_FeatureSign_yang(gamma, A, Q(:, iter1));
end

sc_codes = abs(sc_codes);
scfv_codes = (feaSet - B * sc_codes) * sc_codes'; % for sparse coding based fisher vectors
dSize = size(scfv_codes, 1);
% spatial levels
pLevels = length(pyramid);
% total spatial bins
tBins = sum(pyramid);

beta = zeros(dSize, tBins);
bId = 0;

for iter1 = 1:pLevels
    Unit = nSmp / pyramid(iter1);
    % find to which spatial bin each local descriptor belongs
    idxBin = ceil((1:nSmp) / Unit);

    for iter2 = 1:pyramid(iter1)
        bId = bId + 1;
        sidxBin = find(idxBin == iter2);
        if isempty(sidxBin)
            continue;
        end
        % beta(:, bId) = max(scfv_codes(:, sidxBin), [], 2);
        beta(:, bId) = sum(scfv_codes(:, sidxBin), 2);
    end
end

if bId ~= tBins
    error('Index number error!');
end

beta = beta(:);
beta = beta ./ sqrt(sum(beta.^2));
beta(isnan(beta)) = 0; % avoid NaN
