function build_octave(root, work, vlroot, vllib)
% BUILD_OCTAVE Compile real Octave MEX files outside the research checkout.
inc = ['-I' fullfile(root, 'validation', 'octave_include')];
mkoctfile('--mex', inc, fullfile(root, 'mbs', 'src', 'Determine_segment.cpp'), ...
          '-o', fullfile(work, 'Determine_segment'));
mkoctfile('--mex', inc, fullfile(root, 'mbs', 'src', 'tricircumcenter3d.cpp'), ...
          '-o', fullfile(work, 'tricircumcenter3d'));
for name = {'fisher', 'gmm'}
    source = fullfile(root, 'thirdparty', 'vlfeat-0.9.20', 'toolbox', name{1}, ['vl_' name{1} '.c']);
    mkoctfile('--mex', inc, ['-I' vlroot], ['-I' fullfile(vlroot, 'toolbox')], ...
              source, ['-L' fileparts(vllib)], '-lvl', '-o', fullfile(work, ['vl_' name{1}]));
end
svm = fullfile(root, 'thirdparty', 'libsvm-3.17');
for name = {'svmtrain', 'svmpredict'}
    mkoctfile('--mex', fullfile(svm, 'matlab', [name{1} '.c']), ...
              fullfile(svm, 'svm.cpp'), fullfile(svm, 'matlab', 'svm_model_matlab.c'), ...
              '-o', fullfile(work, name{1}));
end
fprintf('PASS six Octave MEX modules built\n');
end
