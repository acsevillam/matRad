function compatible = matRad_hasCompatibleCstMetadata(cstScenarios)
% matRad checks whether multi-scenario cst metadata can be merged
%
% call
%   compatible = matRad_hasCompatibleCstMetadata(cstScenarios)
%
% input
%   cstScenarios:   cell array with one matRad cst per CT scenario
%
% output
%   compatible:     logical scalar indicating whether structure IDs, names,
%                   types, properties, and objectives match across scenarios
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

if ~iscell(cstScenarios) || isempty(cstScenarios)
    matRad_cfg.dispError('cstScenarios must be a non-empty cell array.');
end

referenceCst = cstScenarios{1};
compatible = iscell(referenceCst) && size(referenceCst,2) >= 6;

for scenarioIx = 2:numel(cstScenarios)
    currentCst = cstScenarios{scenarioIx};
    compatible = compatible && iscell(currentCst) && ...
        size(currentCst,1) == size(referenceCst,1) && ...
        size(currentCst,2) >= 6 && ...
        isequal(currentCst(:,1),referenceCst(:,1)) && ...
        isequal(currentCst(:,2),referenceCst(:,2)) && ...
        isequal(currentCst(:,3),referenceCst(:,3)) && ...
        isequal(currentCst(:,5),referenceCst(:,5)) && ...
        isequal(currentCst(:,6),referenceCst(:,6));

    if ~compatible
        return;
    end
end
end
