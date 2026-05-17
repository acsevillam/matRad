function jobResult = runBreastPhotonMctRobustnessJob(varargin)
% runBreastPhotonMctRobustnessJob Run breast macro parameter sets in series.

% Path setup.
macroRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(macroRoot,'helpers'));

% Job execution defaults.
% profile: prod uses production defaults; testing applies explicit overrides.
% openGui: serial jobs usually run non-interactively.
profile = 'testing';
openGui = false;
stopOnError = true;

% Base MacroSpec selectors.
% Parameter sets below override robustness while sharing this base context.
base = struct();
base.label = 'breast_4136mct';
base.site = 'breast';
base.particleType = 'photons';
base.caseID = '4136_mct';
base.samplingProfile = 'default';

% Base optimization scenario toggles.
base.optimizationScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);

% Base sampling scenario toggles.
base.samplingScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);

% Job parameter sets.
% Each parameter set is merged over the base before execution.
parameterSets = repmat(struct('label','','robustness',''),1,3);
parameterSets(1).label = 'PTV';
parameterSets(1).robustness = 'PTV';
parameterSets(2).label = 'COWC';
parameterSets(2).robustness = 'COWC';
parameterSets(3).label = 'STOCH';
parameterSets(3).robustness = 'STOCH';
%parameterSets(4).label = 'INTERVAL2';
%parameterSets(4).robustness = 'INTERVAL2';
%parameterSets(5).label = 'INTERVAL3';
%parameterSets(5).robustness = 'INTERVAL3';

% Resolve all parameter sets and run them in series.
job = struct();
job.id = 'breast.photons.4136_mct.robustness';
job.description = 'Breast photons 4136_mct robustness study.';
job.profile = profile;
job.openGui = openGui;
job.stopOnError = stopOnError;
job.base = base;
job.parameterSets = parameterSets;

jobResult = runWorkflowMacroJob(job,varargin{:});

end
