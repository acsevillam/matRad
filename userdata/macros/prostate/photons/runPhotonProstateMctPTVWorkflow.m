function runPhotonProstateMctPTVWorkflow(varargin)
% runPhotonProstateMctPTVWorkflow Run the photon 1_mct prostate PTV workflow.
% Usage: runPhotonProstateMctPTVWorkflow('rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'))

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

prepareConfig = photonProstatePrepareConfig( ...
    prostateMctPrepareConfig('PTV_001'));

runProstatePTVWorkflowCore(prepareConfig,varargin{:});

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
