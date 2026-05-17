function result = runWorkflowMacroSpec(specId,varargin)
% runWorkflowMacroSpec Run a local macro spec through MacroRunner.

macroRoot = fileparts(fileparts(mfilename('fullpath')));
specFolder = fullfile(macroRoot,'shared','specs');
if exist(specFolder,'dir') ~= 7
    error('planWorkflow:macros:MissingSpecCatalog', ...
        'Macro spec catalog folder not found: %s.',specFolder);
end
if ~any(strcmp(strsplit(path,pathsep),specFolder))
    addpath(specFolder);
end
ensurePlanWorkflowAvailable('runWorkflowMacroSpec');

[profile,args] = parseMacroProfileOption('prod',varargin{:});
spec = macroSpecCatalog(specId,'profile',profile);
result = planWorkflow.macros.MacroRunner.run(spec,args{:});

end
