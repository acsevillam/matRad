function dijInterval = matRad_buildDoseIntervalDij(ctx, quantity, cfg, intervalMode)
% matRad_buildDoseIntervalDij builds the compact INTERVAL payload
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

dijInterval = struct();
dijInterval.center = sparse(ctx.numVoxels, ctx.numBixels);
dijInterval.radius = sparse(ctx.numBixels, ctx.numBixels);
dijInterval.targetSubIx = ctx.targetRows(:);
dijInterval.OARSubIx = ctx.oarRows(:);
dijInterval.quantity = quantity.name;
dijInterval.quantityField = quantity.field;
dijInterval.quantityScale = quantity.scale;
dijInterval.optimizationQuantity = quantity.optimizationQuantity;
dijInterval.refScen = cfg.refScen;
dijInterval.scenarioDijIx = ctx.scenarioDijIx(:);
dijInterval.scenarioCtScenIds = ctx.scenarioCtScenIds(:);
dijInterval.scenarioWeights = ctx.scenarioWeights(:);
dijInterval.intervalMode = intervalMode;
dijInterval.radiusMode = cfg.RadiusMode;
dijInterval.precomputeMode = 'scenario-batch';
dijInterval.secondPassStrategy = cfg.SecondPassStrategy;

if cfg.CollectTiming
    dijInterval.timing = DoseInterval.matRad_buildDoseIntervalTiming(intervalMode, cfg, ctx);
end

if strcmp(intervalMode, 'INTERVAL3')
    dijInterval.OARRadiusRank = zeros(numel(ctx.oarRows), 1);
    dijInterval.OARRadiusFactor = cell(numel(ctx.oarRows), 1);
end
end
