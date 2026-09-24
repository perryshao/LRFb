function [f, df, ddf] = costFuncRegMultPartGp_v1_42(initialTheta, X, Y, lambda, C, modality)
% COSTFUNCREGMULTPARTGP_V1_42  Evaluate the legacy mixed-norm regression objective.
%   X is D-by-N, Y is N-by-C, and modality contains feature-block sizes.
%   The gradient and third-output issues are documented in docs/CODE_REVIEW.md.
%   This historical objective is not the final linear-SVM classifier.

if nargin < 7
    preTinitialTheta = 0;
end

partGroup = length(modality);
groupCof = 1;
%% Compute the costJ of a particular choice of theta
% compute cost costJ
% X = DxN matrix, J and M is the partition parameters over joints and feature modalities
D = size(X, 1);
% J = D/jointNum; % the number of parts
modality = [0, modality];
% theta = DxC column vector
theta = reshape(initialTheta, D, C);
% costJ = single number
costJ = sum(sum((X' * theta - Y).^2));

costRegularizationTerm1 = sum(sqrt(sum(theta.^2, 2)));

%% part group1
tempJ = 0;
for j = 1:partGroup(1)
    tempM = sqrt(sum(theta(modality(j) + 1:modality(j) + modality(j + 1), :).^4, 1));
    tempJ = tempJ + tempM;
end
costRegularizationTerm2 = tempJ;
clear tempM tempJ;
%% compute the sum cost
costRegularizationTermRow = groupCof(1) * costRegularizationTerm2;
costRegularizationTermMain = sum(costRegularizationTermRow);

costJWithRegularization = costJ + lambda(1) * costRegularizationTerm1 + lambda(2) * costRegularizationTermMain;
% Compute the partial derivatives and set gradient to the partial
% derivatives of the cost w.r.t. each parameter in theta
%% compute the gradient
gradient = 2 * X * (X' * theta - Y);

clear X;

epsilon = 10e-8; % to avoid inf when divided by zero
gradientRegularizationTerm1 = repmat((sum(theta.^2, 2) + epsilon).^(-1 / 2), 1, C) .* theta;

%% part group1
gradientRegularizationTerm2 = zeros(D, C);
for j = 1:partGroup(1)
    gradientRegularizationTerm2(modality(j) + 1:modality(j) + modality(j + 1), :) = ...
        repmat(sum(theta(modality(j) + 1:modality(j) + modality(j + 1), :).^4, 1).^(1 / 4), modality(j + 1), 1) .* ...
                                                                              ((theta(modality(j) + ...
                                                                                  1:modality(j) + modality(j + ...
                                                                                      1), :).^3) ./ ...
                                                                              repmat(sum(theta(modality(j) + ...
                                                                                  1:modality(j) + modality(j + ...
                                                                                      1), :).^4, 1).^(3 / 4) + ...
                                                                                          epsilon, modality(j + 1), 1));
end

%% compute the sum gradient
gradient = gradient + lambda(1) * gradientRegularizationTerm1 + lambda(2) * gradientRegularizationTerm2;

f = costJWithRegularization;
gradient = reshape(gradient, D * C, 1); % vec(W)

if nargout > 1
    df = gradient;
end

if nargout > 2
    ddf = ddgradient;
end
