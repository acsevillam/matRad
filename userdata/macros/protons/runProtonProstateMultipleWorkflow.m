function runProtonProstateMultipleWorkflow(varargin)
% runProtonProstateMultipleWorkflow Run the proton RBExD prostate multiple-plan workflow.
% Usage: runProtonProstateMultipleWorkflow('caseID','3482','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = fileparts(macroFolder);
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
prepareConfig.radiationMode = 'protons';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'constRBE';
prepareConfig.quantityOpt = 'RBExD';
prepareConfig.plan_beams = '2F';
prepareConfig.dicomMetadata = struct();
prepareConfig.resolution = [3 3 3];
prepareConfig.n_cores = 8;

runProstateMultipleWorkflowCore(prepareConfig,varargin{:});

end
