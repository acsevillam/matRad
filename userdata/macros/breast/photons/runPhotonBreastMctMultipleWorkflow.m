function runPhotonBreastMctMultipleWorkflow(varargin)
% runPhotonBreastMctMultipleWorkflow Run the photon breast 4136_mct multiple-plan workflow.
% Usage: runPhotonBreastMctMultipleWorkflow('rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'))

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = findMacroRoot(macroFolder);
sharedFolder = fullfile(macroRoot,'shared','breast');
if exist(sharedFolder,'dir') ~= 7
    error('planWorkflow:macros:MissingSharedMacro', ...
        'Shared breast macro folder not found: %s.',sharedFolder);
end

if ~any(strcmp(strsplit(path,pathsep),sharedFolder))
    addpath(sharedFolder);
end

prepareConfig = photonBreastPrepareConfig('4136_mct','comparison_001');

runBreastMultipleWorkflowCore(prepareConfig,varargin{:});

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
