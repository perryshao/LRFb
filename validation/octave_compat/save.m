function save(varargin)
% SAVE Test-only bridge for MATLAB MAT filenames and unsupported v7.3 saves.
%   Values are read from the original caller workspace. Only named variables
%   and character options used by the LRFb scripts are supported.
args = varargin;
filename_seen = false;
for i = 1:numel(args)
    assert(ischar(args{i}), 'Validation save bridge only accepts character arguments.');
    if strcmpi(args{i}, '-v7.3')
        args{i} = '-mat7-binary';
    elseif ~isempty(args{i}) && args{i}(1) ~= '-' && ~filename_seen
        [~, ~, ext] = fileparts(args{i});
        if isempty(ext)
            args{i} = [args{i} '.mat'];
        end
        filename_seen = true;
    end
end
quoted = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], args, 'UniformOutput', false);
expression = ['builtin(''save'', ''-mat7-binary'', ' strjoin(quoted, ', ') ');'];
evalin('caller', expression);
end
