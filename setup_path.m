function setup_path()
% SETUP_PATH  Add the LRFb project to the MATLAB path.
%
%   Run this once per MATLAB session from the LRFb root:
%       >> setup_path
%   then change into one experiment folder and run its driver, e.g.
%       >> cd experiments/MSRC12
%       >> run
%
%   experiments/ is deliberately NOT added to the path.  The IP and MSRC-12
%   folders each carry their own version of helpers that share a name
%   (Estimate_Frenet, GeneSC, getLabels, sc3d_compute, trainBinRegression*,
%   ...) and the versions differ.  MATLAB resolves the current folder first,
%   which is how the original experiments were run.  In particular
%   experiments/IP/Estimate_Frenet.m (2 arguments, curve_grid) shadows
%   lrf/Estimate_Frenet.m (3 arguments, curve_grid_mean) while you are in
%   the IP folder.
%
%   See README.md for the layout and for which mex targets exist.

root = fileparts(mfilename('fullpath'));

% Project code.
folders = {'lrf', 'mbs', 'descriptor', 'encoding', 'utils', 'baselines'};
for k = 1:numel(folders)
    add_or_warn(fullfile(root, folders{k}), false);
end

% Third-party code: only the folders the pipeline calls into.
tp = fullfile(root, 'thirdparty');
thirdparty = {
    'vlfeat-0.9.20/toolbox/gmm' % vl_gmm     (GMM for the Fisher vectors)
    'vlfeat-0.9.20/toolbox/fisher' % vl_fisher  (improved Fisher vectors)
    'libsvm-3.17/matlab' % svmtrain / svmpredict (linear SVM)
    'netlab' % dist2
    'ndSparse' % n-d sparse arrays, MSRC-12 sc3d_compute
    'ScSPM/large_scale_svm' % li2nsvm_multiclass_lbfgs  (menu alternative)
    'ScSPM/sparse_coding' % reg_sparse_coding         (menu alternative)
    'Stochastic_Bosque/Stochastic_Bosque' % random forest (menu alternative)
    'Stochastic_Bosque/Stochastic_Bosque/cartree'
    'Stochastic_Bosque/Stochastic_Bosque/cartree/mx_files'
    'Stochastic_Bosque/Stochastic_Bosque/ancillary'
    };
for k = 1:numel(thirdparty)
    add_or_warn(fullfile(tp, thirdparty{k}), false);
end

arch = computer('arch');

% vlfeat ships one mex folder per platform, each next to its libvl.
vlmex = struct('glnxa64', 'mexa64', 'glnx86', 'mexglx', 'maci', 'mexmaci', ...
               'maci64', 'mexmaci64', 'win32', 'mexw32', 'win64', 'mexw64');
if isfield(vlmex, arch)
    addpath(fullfile(tp, 'vlfeat-0.9.20', 'toolbox', 'mex', vlmex.(arch)), '-begin');
end

% libsvm: 32-bit mex live in matlab/, 64-bit Windows mex in windows/.
if strcmp(arch, 'win64')
    addpath(fullfile(tp, 'libsvm-3.17', 'windows'), '-begin');
end

% Compiled MBS helpers go in front of everything else.
addpath(fullfile(root, 'mbs', 'bin'), '-begin');

fprintf('LRFb path configured (root: %s)\n', root);

missing = check_mex();
if isempty(missing)
    fprintf('All core mex functions are available for %s.\n', arch);
else
    fprintf(2, 'Missing mex for %s: %s\n', arch, strjoin(missing, ', '));
    fprintf(2, 'See the "Building the mex files" section of README.md.\n');
end
end

function add_or_warn(p, recursive)
if exist(p, 'dir')
    if recursive
        addpath(genpath(p));
    else
        addpath(p);
    end
else
    warning('setup_path:missing', 'Folder not found, skipped: %s', p);
end
end

function missing = check_mex()
% CHECK_MEX  Report which compiled helpers cannot be resolved on this platform.
%   vl_gmm.m / vl_fisher.m are help stubs, so test for the mex file itself.
required = {'Determine_segment', 'tricircumcenter3d', 'vl_gmm', 'vl_fisher', 'svmtrain'};
missing = {};
for k = 1:numel(required)
    if exist(required{k}, 'file') ~= 3
        missing{end + 1} = required{k}; %#ok<AGROW>
    end
end
end
