function [stats,timing] = matRad_computeScenarioDoseBatchStats(quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels,matRad_cfg,options)
% matRad_computeScenarioDoseBatchStats computes weighted scenario dose rows
%
% call
%   [stats,timing] = matRad_computeScenarioDoseBatchStats(quantity, ...
%       scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
%       matRad_cfg,options)
%
% input
%   quantity:         struct describing the linear dij quantity, including
%                     matrixCell and scale
%   scenarioDijIx:    DIJ linear scenario indices used to select matrices
%                     from quantity.matrixCell
%   scenarioWeights:  scenario weights used for weighted statistics
%   scenarioMaps:     cell array with scenario-to-reference CT dose mappings
%   rows:             dose-grid linear row indices in the reference CT scenario
%   numBixels:        number of bixels in the influence matrices
%   matRad_cfg:       MatRad_Config instance for diagnostics
%   options:          optional struct controlling stored scenario rows,
%                     second moments, centered rows, extreme deltas,
%                     timing, and progress callbacks
%
% output
%   stats:            struct with requested batch statistics
%                     (scenarioRows, scenarioWeights, meanRows,
%                     secondMoment, centeredRows, extremeDeltaRows)
%   timing:           struct with extractMapSeconds, centerAccumSeconds,
%                     secondMomentSeconds, centeredRowsSeconds, and
%                     extremeDeltaSeconds
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

if nargin < 8 || isempty(options)
    options = struct();
end

options = setDefault(options,'StoreScenarioRows',true);
options = setDefault(options,'ComputeSecondMoment',false);
options = setDefault(options,'ComputeCenteredRows',false);
options = setDefault(options,'ComputeExtremeDelta',false);
options = setDefault(options,'CollectTiming',false);
options = setDefault(options,'ScenarioProgressFcn',[]);

if options.ComputeCenteredRows || options.ComputeExtremeDelta
    options.StoreScenarioRows = true;
end

scenarioWeights = scenarioWeights(:);
numScenarios = numel(scenarioDijIx);
numRows = numel(rows);

stats = struct();
if options.StoreScenarioRows
    stats.scenarioRows = cell(numScenarios,1);
else
    stats.scenarioRows = {};
end
stats.scenarioWeights = scenarioWeights;
stats.meanRows = sparse(numRows,numBixels);
stats.secondMoment = [];
stats.centeredRows = {};
stats.extremeDeltaRows = [];

if options.ComputeSecondMoment
    stats.secondMoment = sparse(numBixels,numBixels);
end

timing = initializeTiming();
for s = 1:numScenarios
    if ~isempty(options.ScenarioProgressFcn)
        options.ScenarioProgressFcn(s);
    end

    if options.CollectTiming
        sectionTimer = tic;
    end
    scenarioRows = matRad_getScenarioDoseRows(quantity,scenarioDijIx(s), ...
        scenarioMaps{s},rows,matRad_cfg);
    if options.CollectTiming
        timing.extractMapSeconds = timing.extractMapSeconds + toc(sectionTimer);
        sectionTimer = tic;
    end

    if options.StoreScenarioRows
        stats.scenarioRows{s} = scenarioRows;
    end

    stats.meanRows = stats.meanRows + scenarioWeights(s).*scenarioRows;
    if options.CollectTiming
        timing.centerAccumSeconds = timing.centerAccumSeconds + toc(sectionTimer);
    end

    if options.ComputeSecondMoment
        if options.CollectTiming
            sectionTimer = tic;
        end
        stats.secondMoment = stats.secondMoment + scenarioRows' * ...
            (scenarioWeights(s).*scenarioRows);
        if options.CollectTiming
            timing.secondMomentSeconds = timing.secondMomentSeconds + ...
                toc(sectionTimer);
        end
    end
end

if options.ComputeCenteredRows || options.ComputeExtremeDelta
    if options.CollectTiming
        sectionTimer = tic;
    end
    stats.centeredRows = cell(numScenarios,1);

    for s = 1:numScenarios
        stats.centeredRows{s} = stats.scenarioRows{s} - stats.meanRows;
    end

    if options.CollectTiming
        timing.centeredRowsSeconds = timing.centeredRowsSeconds + ...
            toc(sectionTimer);
    end
end

if options.ComputeExtremeDelta
    if options.CollectTiming
        sectionTimer = tic;
    end
    stats.extremeDeltaRows = sparse(numRows,numBixels);
    for s = 1:numScenarios
        stats.extremeDeltaRows = max(stats.extremeDeltaRows, ...
            abs(stats.centeredRows{s}));
    end
    if options.CollectTiming
        timing.extremeDeltaSeconds = timing.extremeDeltaSeconds + ...
            toc(sectionTimer);
    end
end
end

function timing = initializeTiming()
timing = struct();
timing.extractMapSeconds = 0;
timing.centerAccumSeconds = 0;
timing.secondMomentSeconds = 0;
timing.centeredRowsSeconds = 0;
timing.extremeDeltaSeconds = 0;
end

function s = setDefault(s,fieldName,value)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = value;
end
end
