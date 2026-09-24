function [beta] = pv_pooling_ts(feaSet, pyramid)
% PV_POOLING_TS  Sum descriptors within temporal bins and L2-normalize.
%   feaSet is D-by-T; pyramid lists the bin counts at each level.
%   beta is a column vector of length D*sum(pyramid).
%
%   Adapted from pooling code by Jianchao Yang, NEC Research Lab America.
%   Mentor: Kai Yu. July 2008; revised May 2010.

dSize = size(feaSet, 1);
nSmp = size(feaSet, 2);

% compute the local feature for each local feature
pv_codes = feaSet;

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
            %             beta(:, bId) = mean(RefeatSet,2);%average pooling
            %             beta(:, bId) = max(RefeatSet,[],2);
            beta(:, bId) = sum(RefeatSet, 2); % sum pooling
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
