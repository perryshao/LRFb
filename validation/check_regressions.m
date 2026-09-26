function report = check_regressions(root)
% CHECK_REGRESSIONS Exercise fixed contracts through the actual Octave MEX/.m APIs.
report = struct();
% Request every supported output count, using deterministic GMM initialization.
data = reshape(sin(1:80), 2, 40);
for kind = {'double', 'single'}
    x = cast(data, kind{1});
    options = {'Initialization', 'custom', 'InitMeans', x(:, [1 20]), ...
               'InitCovariances', ones(2, 2, kind{1}), 'InitPriors', cast([0.5 0.5], kind{1})};
    all_outputs = cell(1, 5);
    [all_outputs{:}] = vl_gmm(x, 2, options{:});
    assert(all(cellfun(@(v) all(isfinite(v(:))), all_outputs)));
    assert(isequal(size(all_outputs{5}), [2 40]));
    assert(max(abs(sum(all_outputs{5}, 1) - 1)) < 1e-5);
    for count = 1:5
        outputs = cell(1, count);
        [outputs{:}] = vl_gmm(x, 2, options{:});
        for k = 1:count
            assert(isequal(outputs{k}, all_outputs{k}));
        end
    end
end
try
    vl_gmm(data, 2);
    error('test:missingError', 'No-output GMM call was accepted.');
catch err
    assert(~isempty(strfind(err.message, 'Expected one to five output arguments')));
end
try
    outputs = cell(1, 6);
    [outputs{:}] = vl_gmm(data, 2);
    error('test:missingError', 'Six-output GMM call was accepted.');
catch err
    assert(~isempty(strfind(err.message, 'Expected one to five output arguments')));
end
report.gmm = '1-5 outputs; single/double; deterministic optional-output equivalence; invalid arity';
for count = 1:3
    outputs = cell(1, count);
    [outputs{:}] = tricircumcenter3d([3 4 5], [5 4 5], [3 6 5]);
    assert(isequal(outputs{1}, [1 1 0]));
    if count >= 2
        assert(outputs{2} == 0.5);
    end
    if count >= 3
        assert(outputs{3} == 0.5);
    end
end
assert_error(@() tricircumcenter3d([0 0], [1 0 0], [0 1 0]), 'LRFb:circumcenterInput');
assert_error(@() tricircumcenter3d(single([0 0 0]), [1 0 0], [0 1 0]), 'LRFb:circumcenterInput');
assert_error(@() tricircumcenter3d([0 0 0], [1 0 0]), 'LRFb:circumcenterArity');
report.circumcenter = '1-3 outputs and input contract';
% Rodrigues matrices provide an independent reconstruction oracle for half-turns.
axes = [eye(3); 1 2 3; -3 1 2; 1 -4 2];
rotation_cases = 0;
for axis = axes'
    axis = axis / norm(axis);
    skew = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    for angle = [0, 0.1, pi - 1e-8, pi, pi + 1e-8]
        rotation = cos(angle) * eye(3) + (1 - cos(angle)) * (axis * axis') + sin(angle) * skew;
        q = Myrotm2quat(rotation);
        assert(all(isfinite(q)) && abs(norm(q) - 1) < 1e-12 && q(1) >= 0);
        w = q(1); v = q(2:4)';
        crossv = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        reconstructed = (w^2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * crossv;
        assert(norm(reconstructed - rotation, 'fro') < 1e-12);
        rotation_cases = rotation_cases + 1;
    end
end
report.rotation_cases = rotation_cases;
empty = fv_pooling_ts(zeros(2, 0), zeros(2, 1), ones(2, 1), 1, 'Improved', [1 2]);
occluded = fv_pooling_ts(nan(2, 3), zeros(2, 1), ones(2, 1), 1, 'Improved', [1 2]);
zero = fv_pooling_ts([-1 1], 0, 1, 1, 'Improved', 1);
assert(isequal(empty, zeros(12, 1)) && isequal(occluded, empty));
assert(isequal(zero, zeros(2, 1)));
report.pooling = 'empty, fully occluded and exactly zero Fisher vectors remain zero';
sample_data = [(1:6)' (11:16)' (21:26)'];
for count = [0 1 6 8 12 13 101]
    sampled = rand_sampling_ts({sample_data}, count);
    assert(isequal(size(sampled), [3 count]));
    assert(all(ismember(sampled', sample_data, 'rows')));
    if count >= 6
        assert(isequal(sortrows(sampled(:, 1:6)'), sample_data));
        assert(isequal(sampled, sampled(:, mod(0:count - 1, 6) + 1)));
    end
end
assert_error(@() rand_sampling_ts({}, 1), 'LRFb:emptyTrainingSet');
assert_error(@() rand_sampling_ts({zeros(0, 3)}, 1), 'LRFb:trainingSequence');
assert_error(@() rand_sampling_ts({sample_data}, -1), 'LRFb:sampleCount');
assert_error(@() rand_sampling_ts({sample_data}, NaN), 'LRFb:sampleCount');
report.sampler = '0,1,6,8,12,13,101 rows; cyclic permutation; explicit invalid-input errors';
% Both experiment-specific objectives must differentiate the same stated cost.
max_gradient_error = 0;
for dataset = {'IP', 'MSRC12'}
    folder = fullfile(root, 'experiments', dataset{1});
    addpath(folder, '-begin');
    clear costFuncRegMultPartGp_v1_42;
    for groups = {[1], [2 3], [1 2 3]}
        dims = sum(groups{1});
        classes = 2;
        X = reshape(sin(1:dims * 4), dims, 4);
        Y = reshape(cos(1:4 * classes), 4, classes);
        theta = reshape(sin(1:dims * classes), [], 1);
        for values = {theta, zeros(size(theta))}
            values = values{1};
            [f, gradient] = costFuncRegMultPartGp_v1_42(values, X, Y, [0.3 0.7], classes, groups{1});
            assert(isfinite(f) && all(isfinite(gradient)));
            fd = zeros(size(values));
            for i = 1:numel(values)
                delta = zeros(size(values)); delta(i) = 1e-5;
                plus = costFuncRegMultPartGp_v1_42(values + delta, X, Y, [0.3 0.7], classes, groups{1});
                minus = costFuncRegMultPartGp_v1_42(values - delta, X, Y, [0.3 0.7], classes, groups{1});
                fd(i) = (plus - minus) / 2e-5;
            end
            max_gradient_error = max(max_gradient_error, max(abs(gradient - fd)));
            assert(max(abs(gradient - fd)) < 2e-7);
        end
    end
    [f, gradient] = costFuncRegMultPartGp_v1_42(2, 0, 0, [0 1], 1, 1);
    assert(f == 4 && gradient == 4);
    rmpath(folder);
end
report.max_gradient_error = max_gradient_error;
fprintf('PASS fixed gateways, half-turns, pooling, sampler and both regression objectives\n');
end

function assert_error(operation, identifier)
try
    unused = operation();
catch err
    assert(strcmp(err.identifier, identifier));
    return;
end
error('test:missingError', 'Expected %s.', identifier);
end
