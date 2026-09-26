function run_octave(root, work)
% RUN_OCTAVE Execute LRFb MATLAB functions with real Octave MEX modules.
oldpath = path;
oldpwd = pwd;
cleanup = onCleanup(@() restore_session(oldpath, oldpwd));
for folder = {'lrf', 'mbs', 'descriptor', 'encoding', 'utils', 'thirdparty/netlab'}
    addpath(fullfile(root, folder{1}));
end
addpath(work, '-begin');
cd(work);
load(fullfile(work, 'fixtures.mat'));
report = struct();
report.octave_version = version;
% Record Octave's real save limitation before adding the compatibility bridge.
try
    save(fullfile(work, 'unsupported.mat'), 'raw_curve', '-v7.3');
    error('test:unexpected', 'Octave unexpectedly accepted MATLAB v7.3 save.');
catch err
    assert(~strcmp(err.identifier, 'test:unexpected'));
    assert(~isempty(strfind(err.message, 'v7.3')));
    report.original_v73_error = err.message;
end
addpath(fullfile(root, 'validation', 'octave_compat'), '-begin');
% Execute repaired source gateways directly; no output-padding wrapper is used.
report.regressions = check_regressions(root);
resized = TrjResizeTime(raw_curve, 64);
grid_points = curve_grid(resized, 1000);
segments_struct = splitting_curve_3D(grid_points, 5);
segments = segments_struct.B_E_seg;
[regular, repeated_indices] = remove_stapoint(repeated_curve);
assert(isequal(fillstapoint(regular, repeated_indices), repeated_curve));
cc = tricircumcenter3d([0 0 0], [2 0 0], [0 2 0]);
[cc2, xi] = tricircumcenter3d([0 0 0], [2 0 0], [0 2 0]);
[cc3, xi3, eta] = tricircumcenter3d([0 0 0], [2 0 0], [0 2 0]);
assert(isequal(cc, [1 1 0]) && isequal(cc, cc2) && isequal(cc, cc3));
assert(xi == .5 && xi3 == .5 && eta == .5);
fisher_actual = vl_fisher(fisher_input, fisher_means, fisher_cov, fisher_priors, 'Improved');
pool_actual = fv_pooling_ts(fisher_input, fisher_means, fisher_cov, fisher_priors, 'Improved', [1 2 4 8]);
assert(all(isfinite([fisher_actual(:); pool_actual(:)])));
report.geometry_and_pool = 'executed original .m and rebuilt Octave MEX';
datasets = {'IP', 'MSRC12'};
outputs = cell(1, 2);
for ds = 1:2
    name = datasets{ds};
    experiment = fullfile(root, 'experiments', name);
    addpath(experiment, '-begin');
    clear Estimate_Frenet GeneSC getLabels load_data_bat load_IPtxt_bat costFuncRegMultPartGp_v1_42;
    case_dir = fullfile(work, name);
    mkdir(case_dir);
    cd(case_dir);
    if ds == 1
        mkdir('Hand_Data');
        Data = ip_data;
        save('Hand_Data/Data.mat', 'Data');
        load_IPtxt_bat('unused');
        original_db = load('DB.mat');
        original_test = load('SAMPLES.mat');
        assert(size(original_db.TRAJDB, 2) == 9);
        assert(size(original_test.TRAJSAMPLES, 2) == 6);
        report.small_IP_loader = '9 training and 6 test slots; no test-side trimming';
        assert(strcmp(which('Estimate_Frenet'), fullfile(experiment, 'Estimate_Frenet.m')));
        local_frames = Estimate_Frenet(resized, 5);
        core_descriptor = RRVdescriptor_BasedonFrames(local_frames(3:end - 2, :), resized(3:end - 2, :));
    else
        SmthTrj = msrc_data;
        save('MSRC12_Skeleton.mat', 'SmthTrj');
        load_data_bat({'8'; '12'}, [1; 2]);
        assert(strcmp(which('Estimate_Frenet'), fullfile(root, 'lrf', 'Estimate_Frenet.m')));
    end
    [train_labels, test_labels] = getLabels();
    % Actual train and test generators, including their saved intermediate files.
    assert(strcmp(which('GeneSC'), fullfile(experiment, 'GeneSC.m')));
    GeneSC({'8'; '12'});
    train_desc = load('RRV_DB.mat');
    test_desc = load('RRV_SAMPLES.mat');
    originals = train_desc.RRV_DB;
    assert(numel(originals) == numel(train_labels));
    assert(numel(test_desc.RRV_SAMPLES) == numel(test_labels));
    for i = 1:numel(originals)
        assert(all(isfinite(originals{i}(:))) && size(originals{i}, 2) == 14);
    end
    for i = 1:numel(test_desc.RRV_SAMPLES)
        assert(all(isfinite(test_desc.RRV_SAMPLES{i}(:))));
    end
    % The fixed sampler accepts the original small training set directly.
    % No descriptor-cell repetition or reweighting is applied by this test.
    rand('seed', 240925);
    GeneFisherCodeJointPyramid_whole(2, 3, 4, 0);
    train_codes = load('traindata.mat');
    test_codes = load('testdata.mat');
    model_file = load('modelForTest.mat');
    assert(size(train_codes.traindata, 1) == 3360);
    assert(size(train_codes.traindata, 2) == numel(train_labels));
    assert(all(isfinite(train_codes.traindata(:))) && all(isfinite(test_codes.testdata(:))));
    assert(max(abs(sqrt(sum(train_codes.traindata.^2, 1)) - sqrt(2))) < 1e-9);
    model = svmtrain(train_labels, train_codes.traindata', '-c 1 -g 2.79e-4 -t 0 -b 1');
    [prediction, accuracy, probability] = svmpredict(test_labels, test_codes.testdata', model, '-b 1');
    assert(all(isfinite(probability(:))) && max(abs(sum(probability, 2) - 1)) < 1e-9);
    modelForTest = model_file.modelForTest;
    modelForTest.model = model;
    save('modelForTest.mat', 'modelForTest');
    restored = load('modelForTest.mat');
    prediction2 = svmpredict(test_labels, test_codes.testdata', restored.modelForTest.model, '-b 1');
    assert(isequal(prediction, prediction2));
    outputs{ds} = struct('name', name, 'original_descriptors', {originals}, ...
                         'test_descriptors', {test_desc.RRV_SAMPLES}, 'train_labels', train_labels, ...
                         'test_labels', test_labels, 'gmm', model_file.modelForTest, ...
                         'first_train_code', train_codes.traindata(:, 1), ...
                         'test_codes', test_codes.testdata, 'prediction', prediction, 'accuracy', accuracy);
    fprintf('PASS %s loader, GeneSC, original 50000-sample encoder, SVM, model round trip\n', name);
    rmpath(experiment);
end
cd(work);
builtin('save', '-mat7-binary', fullfile(work, 'octave_results.mat'), 'report', 'outputs', ...
        'resized', 'grid_points', 'segments', 'local_frames', 'core_descriptor', ...
        'fisher_actual', 'pool_actual');
fprintf('PASS original .m execution with MAT-format compatibility only\n');
end

function restore_session(oldpath, oldpwd)
path(oldpath);
cd(oldpwd);
end
