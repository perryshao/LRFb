function [f, df, ddf] = costFuncRegMultPartGp_v1_42(initialTheta, X, Y, lambda, C, modality)
% COSTFUNCREGMULTPARTGP_V1_42  Evaluate the mixed-norm regression objective.
%   X is D-by-N, Y is N-by-C; modality gives contiguous feature-block sizes.
%   Regularizers are the row L2 sum and the sum of squared group L4 norms.
%   Zero rows/groups use zero (sub)gradients. A Hessian is not implemented.

if nargout > 2
    error('LRFb:hessianUnsupported', 'Only the objective and gradient are available.');
end
D = size(X, 1);
partGroup = length(modality);
modality = [0 cumsum(modality)];
theta = reshape(initialTheta, D, C);
residual = X' * theta - Y;
rowNorm = sqrt(sum(theta.^2, 2));
rowDenominator = rowNorm;
rowDenominator(rowDenominator == 0) = 1;
gradientRegularizationTerm1 = theta ./ repmat(rowDenominator, 1, C);
costRegularizationTerm2 = 0;
gradientRegularizationTerm2 = zeros(D, C);
for j = 1:partGroup
    block = theta(modality(j) + 1:modality(j + 1), :);
    groupNorm = sqrt(sum(block.^4, 1));
    costRegularizationTerm2 = costRegularizationTerm2 + sum(groupNorm);
    groupDenominator = groupNorm;
    groupDenominator(groupDenominator == 0) = 1;
    gradientRegularizationTerm2(modality(j) + 1:modality(j + 1), :) = ...
        2 * block.^3 ./ repmat(groupDenominator, modality(j + 1) - modality(j), 1);
end
f = sum(sum(residual.^2)) + lambda(1) * sum(rowNorm) + lambda(2) * costRegularizationTerm2;
gradient = 2 * X * residual + lambda(1) * gradientRegularizationTerm1 + ...
           lambda(2) * gradientRegularizationTerm2;
if nargout > 1
    df = reshape(gradient, D * C, 1);
end
end
