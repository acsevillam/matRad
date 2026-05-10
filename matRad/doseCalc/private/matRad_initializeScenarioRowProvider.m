function [provider,dijForResolve] = matRad_initializeScenarioRowProvider( ...
    ct,cst,stf,pln,cfg,scenarioInfo,matRad_cfg,calculationName)
% matRad_initializeScenarioRowProvider initializes streaming row access
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

if nargin < 8 || isempty(calculationName)
    calculationName = 'streaming scenario dose';
end

provider = struct();
provider.secondPassStrategy = cfg.SecondPassStrategy;
provider.keepCache = cfg.KeepCache;
provider.cacheRoot = cfg.CacheRoot;
provider.cacheContext = [];

if isfield(cfg,'PrecomputedDij') && ~isempty(cfg.PrecomputedDij)
    provider.type = 'inMemory';
    provider.dij = cfg.PrecomputedDij;
    dijForResolve = cfg.PrecomputedDij;
    return;
end

provider.type = 'streaming';
provider.ct = ct;
provider.cst = cst;
provider.stf = stf;
provider.pln = pln;
provider.scenarioIds = scenarioInfo.scenarioIds(:);

firstScenarioId = scenarioInfo.scenarioIds(1);
pln_s = matRad_selectSingleScenarioPlan(pln,firstScenarioId);
dijTemplate = matRad_calcDoseInfluence(ct,cst,stf,pln_s);
provider.preloadedScenarioId = firstScenarioId;
provider.preloadedDij = dijTemplate;

dijForResolve = matRad_buildSyntheticDijForScenarioResolution( ...
    dijTemplate,pln,scenarioInfo.scenarioDijIx,matRad_cfg,calculationName);
end
