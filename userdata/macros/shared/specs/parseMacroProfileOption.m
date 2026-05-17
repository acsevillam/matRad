function [profile,args] = parseMacroProfileOption(defaultProfile,varargin)
% parseMacroProfileOption Extract the macro execution profile from options.

if nargin < 1 || isempty(defaultProfile)
    defaultProfile = 'prod';
end
profile = normalizeProfile(defaultProfile);
args = varargin;
i = 1;
while i <= numel(args)
    if isstruct(args{i}) && isscalar(args{i}) && isfield(args{i},'profile')
        profile = normalizeProfile(args{i}.profile);
        args{i} = rmfield(args{i},'profile');
        if isempty(fieldnames(args{i}))
            args(i) = [];
        else
            i = i + 1;
        end
    elseif i < numel(args) && isTextScalar(args{i}) && ...
            strcmp(char(string(args{i})),'profile')
        profile = normalizeProfile(args{i + 1});
        args(i:i + 1) = [];
    else
        i = i + 1;
    end
end

end

function profile = normalizeProfile(profile)

if ~(ischar(profile) || (isstring(profile) && isscalar(profile)))
    error('planWorkflow:macros:InvalidMacroProfile', ...
        'Macro profile must be "prod" or "testing".');
end
profile = lower(char(string(profile)));
if ~any(strcmp(profile,{'prod','testing'}))
    error('planWorkflow:macros:InvalidMacroProfile', ...
        'Macro profile must be "prod" or "testing".');
end

end

function tf = isTextScalar(value)

tf = ischar(value) || (isstring(value) && isscalar(value));

end
