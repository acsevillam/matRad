function [isSupported, reason] = matRad_supportsParallelScenarioDij(engine)
% matRad_supportsParallelScenarioDij checks multi-scenario dij support
%
% call
%   isSupported = ScenarioBatch.Parallel.matRad_supportsParallelScenarioDij(engine)
%   [isSupported,reason] = ScenarioBatch.Parallel.matRad_supportsParallelScenarioDij(engine)
%
% input
%   engine:     configured matRad dose engine
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

if nargin < 1 || isempty(engine) || ...
        ~isa(engine, 'DoseEngines.matRad_DoseEngineBase')
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('A configured dose engine is required.');
end

isSupported = false;
reason = '';

if ~isa(engine, 'DoseEngines.matRad_PencilBeamEngineAbstract')
    reason = 'notAnalyticalPencilBeam';
    return
end

if matRad_ispropCompat(engine, 'enableDijSampling') && ...
        logical(engine.enableDijSampling)
    reason = 'stochasticDijSampling';
    return
end

isSupported = true;
end
