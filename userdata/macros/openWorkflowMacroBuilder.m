function openWorkflowMacroBuilder(varargin)
% openWorkflowMacroBuilder Open the planWorkflow GUI for macro creation.

macroRoot = fileparts(mfilename('fullpath'));
helperFolder = fullfile(macroRoot,'helpers');
specFolder = fullfile(macroRoot,'shared','specs');
if ~any(strcmp(strsplit(path,pathsep),helperFolder))
    addpath(helperFolder);
end
if ~any(strcmp(strsplit(path,pathsep),specFolder))
    addpath(specFolder);
end
ensurePlanWorkflowAvailable('openWorkflowMacroBuilder');

[specId,profile,args] = parseBuilderArgs(varargin{:});
spec = macroSpecCatalog(specId,'profile',profile);
[workflowConfig,~] = planWorkflow.macros.MacroRunner.build( ...
    spec,args{:});

workflow = planWorkflow.Workflow(workflowConfig);
cleanupObj = onCleanup(@() workflow.releaseMemory()); %#ok<NASGU>

try
    workflow.gui();
catch ME
    if strcmp(ME.identifier,'planWorkflow:gui:PlanEditor:Cancelled')
        return;
    end
    rethrow(ME);
end

if ~isempty(workflow.guiProgressReporter)
    workflow.guiProgressReporter.log( ...
        'Macro builder finished. No workflow stages were executed.');
    workflow.guiProgressReporter.setProgress(0, ...
        'Macro builder finished.');
end

end

function [specId,profile,args] = parseBuilderArgs(varargin)

specId = 'prostate.photons.3482.INTERVAL2';
profile = 'prod';
args = varargin;

if ~isempty(args) && isTextScalar(args{1}) && ...
        contains(char(string(args{1})),'.')
    specId = char(string(varargin{1}));
    args = args(2:end);
end

if numel(args) >= 2 && isTextScalar(args{1}) && ...
        strcmp(char(string(args{1})),'specId')
    specId = char(string(args{2}));
    args(1:2) = [];
end

[profile,args] = parseMacroProfileOption('prod',args{:});

end

function tf = isTextScalar(value)

tf = ischar(value) || (isstring(value) && isscalar(value));

end
