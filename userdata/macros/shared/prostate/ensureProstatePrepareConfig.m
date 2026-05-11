function prepareConfig = ensureProstatePrepareConfig(prepareConfig)
% ensureProstatePrepareConfig Validate or initialize prostate prepare settings.

if nargin < 1 || isempty(prepareConfig)
    prepareConfig = struct();
elseif ~isstruct(prepareConfig)
    error('planWorkflow:macros:InvalidPrepareConfig', ...
        'A workflowConfig.prepare struct is required.');
end

end
