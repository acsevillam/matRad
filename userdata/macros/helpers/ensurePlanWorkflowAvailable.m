function ensurePlanWorkflowAvailable(callerName)
% ensurePlanWorkflowAvailable Ensure the planWorkflow submodule is on the MATLAB path.

if nargin < 1 || isempty(callerName)
    callerName = 'planWorkflow workflow';
end

if isPlanWorkflowAvailable()
    return;
end

candidateMatRadRoots = {};

if exist('MatRad_Config','class') == 8
    matRad_cfg = MatRad_Config.instance();
    candidateMatRadRoots{end + 1} = matRad_cfg.matRadRoot;
end

helperFolder = fileparts(mfilename('fullpath'));
macroRoot = fileparts(helperFolder);
userDataRoot = fileparts(macroRoot);
candidateMatRadRoots{end + 1} = fileparts(userDataRoot);

for i = 1:numel(candidateMatRadRoots)
    planWorkflowRoot = fullfile(candidateMatRadRoots{i},'submodules','planWorkflow');
    if exist(fullfile(planWorkflowRoot,'+planWorkflow'),'dir') == 7
        addpath(planWorkflowRoot);
    end
end

if ~isPlanWorkflowAvailable()
    error('planWorkflow:SubmoduleNotLoaded', ...
        ['planWorkflow.Workflow is not on the MATLAB path. Initialize ' ...
         'the planWorkflow submodule under <matRadRoot>/submodules/planWorkflow and ' ...
         'run matRad_rc from that checkout before calling %s.'],callerName);
end

end

function isAvailable = isPlanWorkflowAvailable()

isAvailable = ~isempty(which('planWorkflow.Workflow'));

end
