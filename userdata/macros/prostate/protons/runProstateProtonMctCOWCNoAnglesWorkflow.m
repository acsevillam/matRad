function result = runProstateProtonMctCOWCNoAnglesWorkflow(varargin)
% runProstateProtonMctCOWCNoAnglesWorkflow Configure and run a local planWorkflow MacroSpec.

% Path setup.
macroRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(macroRoot,'helpers'));

% Execution defaults.
% profile: prod uses production defaults; testing applies explicit overrides.
% openGui: true opens the planWorkflow GUI before running stages.
profile = 'prod';
openGui = true;

% MacroSpec selectors.
% site: anatomy/catalog key (breast, prostate, head_and_neck).
% particleType: radiation mode used by prepare.radiationMode.
% caseID: patient/case identifier resolved under userdata.
% robustness: canonical robust plan key or multiple-plan selector.
% samplingProfile: sampling profile (default, noAngles).
site = 'prostate';
particleType = 'protons';
caseID = '1_mct';
robustness = 'COWC';
samplingProfile = 'noAngles';

% Optimization scenario toggles.
% Controls robust scenarios used for precompute and optimization plans.
optimizationScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', true, ...
    'gantryActive', false, ...
    'couchActive', false);

% Sampling scenario toggles.
% Controls active uncertainty dimensions used by the sampling stage.
samplingScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', true, ...
    'gantryActive', false, ...
    'couchActive', false);

% Resolve the MacroSpec and run the complete workflow.
specId = macroSpecId( ...
    site,particleType,caseID,robustness,samplingProfile);
result = runWorkflowMacroSpec( ...
    specId,'profile',profile,'openGui',openGui, ...
    'optimizationScenario',optimizationScenario, ...
    'samplingScenario',samplingScenario,varargin{:});
end
