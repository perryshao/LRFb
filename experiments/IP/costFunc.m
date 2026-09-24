function [f, df, ddf] = costFunc(initialTheta, X, Y, lambda, C, jointNum)
% COSTFUNC  Evaluate the squared-error regression objective and its gradient.
%   X is D-by-N, Y is N-by-C, and initialTheta is the vectorized D-by-C model.
%   The active penalty is lambda(2) * sum(theta(:).^2).
%   Only the first two outputs are implemented; the legacy Hessian path is invalid.

D = size(X, 1);
J = D / jointNum; % the number of parts
% theta = DxC column vector
theta = reshape(initialTheta, D, C);
% costJ = single number
costJ = sum(sum((X' * theta - Y).^2));

costRegularizationTerm1 = sum(sum(theta.^2, 2));

%% compute the sum cost
costJWithRegularization = costJ + lambda(2) * costRegularizationTerm1;
% Compute the partial derivatives and set gradient to the partial
% derivatives of the cost w.r.t. each parameter in theta
%% compute the gradient
gradient = 2 * X * (X' * theta - Y);

clear X;

epsilon = 10e-8; % to avoid inf when divided by zero
gradientRegularizationTerm1 = 2 * theta;
%% compute the sum gradient
gradient = gradient + lambda(2) * gradientRegularizationTerm1;

f = costJWithRegularization;
gradient = reshape(gradient, D * C, 1); % vec(W)

if nargout > 1
    df = gradient;
end

if nargout > 2
    ddf = ddgradient;
end
