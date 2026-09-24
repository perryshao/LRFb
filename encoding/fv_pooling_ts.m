function [beta] = fv_pooling_ts(feaSet, means, covariances, priors, normalizeF, pyramid)
% FV_POOLING_TS  Compute and concatenate Fisher vectors over temporal bins.
%   feaSet is D-by-T; means/covariances are D-by-K; priors has K entries.
%   normalizeF is passed to vl_fisher (normally Improved).
%   pyramid lists bin counts, e.g. [1 2 4 8]. beta has 2*D*K*sum(pyramid) entries.
%   Each bin encodes its descriptors jointly; the concatenation is L2-normalized.
%   The unused per-frame FV pass is retained to preserve historical execution.
%
%   Adapted from pooling code by Jianchao Yang, NEC Research Lab America.
%   Mentor: Kai Yu. July 2008; revised May 2010.

dSize = size(covariances, 1) * size(covariances, 2) * 2;
nSmp = size(feaSet, 2);
fv_codes = zeros(dSize, nSmp);

% compute the local feature for each local feature
for iter1 = 1:nSmp
    fv_codes(:, iter1) = vl_fisher(feaSet(:, iter1), means, covariances, priors, normalizeF);
end

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

        % average pooling for occlusion

        RefeatSet = feaSet(:, sidxBin);
        RefeatSet(:, isnan(RefeatSet(end, :))) = []; % for occlusion
        if isempty(RefeatSet)
            beta(:, bId) = zeros(dSize, 1);
        else
            beta(:, bId) = vl_fisher(RefeatSet, means, covariances, priors, normalizeF);
        end

        % average pooling
        %         beta(:, bId) = vl_fisher(feaSet(:,sidxBin), means, covariances, priors,normalizeF);
        % max pooling
        %         beta(:, bId) = max(fv_codes(:,sidxBin),[],2);
    end
end

if bId ~= tBins
    error('Index number error!');
end

beta = beta(:);
beta = beta ./ sqrt(sum(beta.^2));
% beta(isnan(beta)) = 0;% avoid NaN
