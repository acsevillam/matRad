function model = matRad_createScenarioModel(ct, modelMetadata)
% matRad_createScenarioModel creates a matRad scenario model for a CT
%
% call
%   model = matRad_createScenarioModel(ct)
%   model = matRad_createScenarioModel(ct,modelMetadata)
%
% input
%   ct:             matRad CT struct
%   modelMetadata:  scenario model short name, metadata struct, or
%                   matRad_ScenarioModel instance
%
% output
%   model:          matRad_ScenarioModel instance
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

if nargin < 2 || isempty(modelMetadata)
    modelMetadata = 'nomScen';
end

model = matRad_ScenarioModel.create(modelMetadata, ct);

end
