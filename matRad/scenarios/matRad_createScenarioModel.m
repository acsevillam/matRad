function multScen = matRad_createScenarioModel(ct,scenarioModel,varargin)
% matRad_createScenarioModel creates scenario model instances.
%
% call
%   multScen = matRad_createScenarioModel(ct,scenarioModel)
%   multScen = matRad_createScenarioModel(ct,scenarioModel,'property',value,...)
%
% input
%   ct:                 matRad ct struct, used to initialize available CT
%                       scenario ids; can be empty for a default model
%   scenarioModel:      scenario creation method. Supported names are
%                       'nomScen', 'wcScen', 'impScen',
%                       'truncatedImpScen', 'rndScen', and 'random'
%   varargin:           optional property-value pairs assigned to the
%                       created scenario model instance. Unknown property
%                       names raise an error.
%
% output
%   multScen:           scenario model instance
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

if isstring(scenarioModel) && isscalar(scenarioModel)
    scenarioModel = char(scenarioModel);
end

switch scenarioModel
    case 'nomScen'
        multScen = matRad_NominalScenario(ct);
    case 'wcScen'
        multScen = matRad_WorstCaseScenarios(ct);
    case 'impScen'
        multScen = matRad_ImportanceScenarios(ct);
    case 'truncatedImpScen'
        multScen = matRad_TruncatedImportanceScenarios(ct);
    case {'rndScen','random'}
        multScen = matRad_RandomScenarios(ct);
    otherwise
        matRad_cfg.dispError('''%s'' not known as scenario type!',scenarioModel);
end

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        matRad_cfg.dispError('Scenario model options need to be property/value pairs.');
    end
    for i = 1:2:numel(varargin)
        propName = varargin{i};
        propValue = varargin{i+1};
        if isstring(propName) && isscalar(propName)
            propName = char(propName);
        end
        if matRad_ispropCompat(multScen,propName)
            multScen.(propName) = propValue;
        else
            matRad_cfg.dispError('Scenario model "%s" does not have property "%s".',scenarioModel,propName);
        end
    end
end

end
