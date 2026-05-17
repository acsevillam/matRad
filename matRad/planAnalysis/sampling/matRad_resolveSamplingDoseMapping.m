function doseMapping = matRad_resolveSamplingDoseMapping(ct, multScen, refScen)
% matRad_resolveSamplingDoseMapping validates multi-CT dose mapping
%
% call
%   doseMapping = matRad_resolveSamplingDoseMapping(ct,multScen,refScen)
%
% input
%   ct:          matRad ct struct
%   multScen:    matRad scenario model used for sampling
%   refScen:     reference CT scenario id
%
% output
%   doseMapping: struct describing whether sampled doses must be mapped to
%                the reference CT scenario before analysis
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

matRadCfg = MatRad_Config.instance();

doseMapping = struct();
doseMapping.enabled = false;
doseMapping.refScen = refScen;
doseMapping.method = 'none';

if ~(isfield(ct, 'numOfCtScen') && ct.numOfCtScen > 1)
    return
end

if ~isa(multScen, 'matRad_ScenarioModel')
    matRadCfg.dispError('Multi-CT sampling analysis requires a matRad_ScenarioModel instance.');
end

scenarioIds = multScen.scenarioIds();
ctScenIds = unique(arrayfun(@(id) multScen.getCtScenario(id), scenarioIds));
if all(ctScenIds == refScen)
    return
end

if ~isfield(ct, 'dvf') || isempty(ct.dvf)
    matRadCfg.dispError('Multi-CT sampling analysis requires pull DVFs to map sampled doses to the reference CT scenario.');
end

if ~matRad_hasPullDvf(ct)
    matRadCfg.dispError('Multi-CT sampling analysis requires pull DVFs to map sampled doses to the reference CT scenario.');
end

matRad_validateDvfReferenceScenario(ct, refScen, matRadCfg);
matRad_validateDvfCoverage(ct, ctScenIds, refScen, matRadCfg);

doseMapping.enabled = true;
doseMapping.method = 'pullDirectDoseMapping';

end

function matRad_validateDvfReferenceScenario(ct, refScen, matRadCfg)
dvfRefScen = [];
if isfield(ct, 'dvfMetadata')
    if isfield(ct.dvfMetadata, 'refScen') && ~isempty(ct.dvfMetadata.refScen)
        dvfRefScen = ct.dvfMetadata.refScen;
    elseif isfield(ct.dvfMetadata, 'referenceCtScen') && ~isempty(ct.dvfMetadata.referenceCtScen)
        dvfRefScen = ct.dvfMetadata.referenceCtScen;
    elseif isfield(ct.dvfMetadata, 'referenceScenario') && ~isempty(ct.dvfMetadata.referenceScenario)
        dvfRefScen = ct.dvfMetadata.referenceScenario;
    end
end

if ~isempty(dvfRefScen) && dvfRefScen ~= refScen
    matRadCfg.dispError('Pull DVFs were generated for reference CT scenario %d, but sampling uses reference CT scenario %d.', ...
                        dvfRefScen, refScen);
end
end

function matRad_validateDvfCoverage(ct, ctScenIds, refScen, matRadCfg)
for i = 1:numel(ctScenIds)
    ctScenId = ctScenIds(i);
    if ctScenId == refScen
        continue
    end

    if numel(ct.dvf) < ctScenId || isempty(ct.dvf{ctScenId})
        matRadCfg.dispError('ct.dvf must contain a pull deformation vector field for CT scenario %d.', ctScenId);
    end
end
end

function hasPull = matRad_hasPullDvf(ct)
dvfType = '';
if isfield(ct, 'dvfMetadata') && isfield(ct.dvfMetadata, 'dvfType') && ~isempty(ct.dvfMetadata.dvfType)
    dvfType = char(ct.dvfMetadata.dvfType);
elseif isfield(ct, 'dvfType') && ~isempty(ct.dvfType)
    dvfType = char(ct.dvfType);
end
hasPull = strcmp(dvfType, 'pull');
end
