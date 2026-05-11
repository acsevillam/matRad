function runPhotonProstateMultipleWorkflow(varargin)
% runPhotonProstateMultipleWorkflow Run the photon prostate multiple-plan workflow.
% Usage: runPhotonProstateMultipleWorkflow('caseID','3482','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

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

prepareConfig = struct();
prepareConfig.caseID = '3482';
prepareConfig.AcquisitionType = 'dicom';
prepareConfig.hlutFileName = 'matRad_default.hlut';
prepareConfig.description = 'prostate';
prepareConfig.plan_template = 'comparison_001';
prepareConfig.radiationMode = 'photons';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'none';
prepareConfig.plan_beams = '9F';
prepareConfig.dicomMetadata = struct();
prepareConfig.resolution = [3 3 3];
prepareConfig.n_cores = 8;

runProstateMultipleWorkflowCore(prepareConfig,varargin{:});

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
