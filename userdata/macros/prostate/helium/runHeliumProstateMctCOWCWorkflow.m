function runHeliumProstateMctCOWCWorkflow(varargin)
% runHeliumProstateMctCOWCWorkflow Run the helium 1_mct prostate COWC workflow.
% Usage: runHeliumProstateMctCOWCWorkflow('rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = findMacroRoot(macroFolder);
sharedFolder = fullfile(macroRoot,'shared','prostate');
if exist(sharedFolder,'dir') ~= 7
    error('planWorkflow:macros:MissingSharedMacro', ...
        'Shared prostate macro folder not found: %s.',sharedFolder);
end

if ~any(strcmp(strsplit(path,pathsep),sharedFolder))
    addpath(sharedFolder);
end

prepareConfig = heliumProstatePrepareConfig( ...
    prostateMctPrepareConfig('COWC_001'));

runProstateCOWCWorkflowCore(prepareConfig,varargin{:});

end

function macroRoot = findMacroRoot(startFolder)

macroRoot = startFolder;
while true
    if exist(fullfile(macroRoot,'helpers'),'dir') == 7
        return;
    end
    parentFolder = fileparts(macroRoot);
    if strcmp(parentFolder,macroRoot)
        error('planWorkflow:macros:MacroRootNotFound', ...
            'Could not locate userdata/macros root from %s.',startFolder);
    end
    macroRoot = parentFolder;
end

end
