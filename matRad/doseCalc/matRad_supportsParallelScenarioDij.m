function [isSupported,reason] = matRad_supportsParallelScenarioDij(plnOrEngine)
% matRad_supportsParallelScenarioDij checks multi-scenario dij support
%
% call
%   isSupported = matRad_supportsParallelScenarioDij(pln)
%   isSupported = matRad_supportsParallelScenarioDij(engine)
%   [isSupported,reason] = matRad_supportsParallelScenarioDij(...)
%
% input
%   plnOrEngine:  matRad plan struct or configured matRad dose engine
%
% output
%   isSupported: true when the engine can calculate independent scenario
%                dose-influence matrices in parallel
%   reason:      empty when supported, otherwise a stable reason identifier
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

if nargin < 1 || isempty(plnOrEngine)
    matRad_cfg.dispError('A plan or dose engine is required.');
end

if isa(plnOrEngine,'DoseEngines.matRad_DoseEngineBase')
    engine = plnOrEngine;
elseif isstruct(plnOrEngine) && isfield(plnOrEngine,'propDoseCalc') && ...
        isa(plnOrEngine.propDoseCalc,'DoseEngines.matRad_DoseEngineBase')
    engine = plnOrEngine.propDoseCalc;
else
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(plnOrEngine);
end

isSupported = false;
reason = '';

if ~isa(engine,'DoseEngines.matRad_PencilBeamEngineAbstract')
    reason = 'notAnalyticalPencilBeam';
    return;
end

if matRad_ispropCompat(engine,'enableDijSampling') && ...
        logical(engine.enableDijSampling)
    reason = 'stochasticDijSampling';
    return;
end

isSupported = true;
end
