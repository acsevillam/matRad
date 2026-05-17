function dijIntervalContext = matRad_buildDoseIntervalDijContext(ct, cst, stf, ...
                                                                 pln, dijInterval, quantity, cfg)
% matRad_buildDoseIntervalDijContext builds nominal INTERVAL context
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

numBixels = size(dijInterval.center, 2);
robustDij = [];
if isfield(cfg, 'PrecomputedDij')
    robustDij = cfg.PrecomputedDij;
end
dijIntervalContext = ScenarioBatch.Dij.matRad_buildNominalDijContext( ...
                                                                     ct, cst, stf, pln, cfg, numBixels, robustDij);
dijIntervalContext.intervalQuantity = quantity.name;
dijIntervalContext.intervalQuantityField = quantity.field;
end
