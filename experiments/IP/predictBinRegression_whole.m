function predictLabel = predictBinRegression_whole(testdata, theta)
% PREDICTBINREGRESSION_WHOLE  Predict labels by maximizing linear regression scores.
%   testdata is D-by-N; theta is D-by-C; predictLabel is N-by-1.

N = size(testdata, 2);
predictLabel = zeros(N, 1);
for i = 1:N
    [~, predictLabel(i)] = max(testdata(:, i)' * theta);
end
