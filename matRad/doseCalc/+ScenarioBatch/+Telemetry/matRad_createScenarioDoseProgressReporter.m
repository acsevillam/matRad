function reporter = matRad_createScenarioDoseProgressReporter(matRadCfg, cfg, ...
                                                              stageName, totalWorkItems, useDataQueue)
% matRad_createScenarioDoseProgressReporter creates scenario-batch progress output
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

if nargin < 5 || isempty(useDataQueue)
    useDataQueue = false;
end

reporter = struct('queue', [], 'enabled', false, 'update', @matRad_noopUpdate, ...
                  'cleanup', @matRad_noopCleanup);

if nargin < 4 || isempty(totalWorkItems) || totalWorkItems <= 0
    return
end

totalWorkItems = double(totalWorkItems);
summaryStep = max(1, ceil(totalWorkItems / 10));
progressListener = [];
progressState = containers.Map({'completedWorkItems'}, {0});
progressUpdate = @(detailText) matRad_logProgress(progressState, matRadCfg, cfg, ...
                                                  stageName, totalWorkItems, summaryStep, detailText);

if useDataQueue
    if exist('parallel.pool.DataQueue', 'class') ~= 8
        matRad_warnDataQueueUnavailable(matRadCfg, stageName);
        return
    end

    progressQueue = parallel.pool.DataQueue;
    progressListener = afterEach(progressQueue, progressUpdate);
    reporter.queue = progressQueue;
    reporter.cleanup = @() matRad_cleanupProgressListener(progressListener);
else
    reporter.update = progressUpdate;
end

reporter.enabled = true;
end

function matRad_logProgress(progressState, matRadCfg, cfg, stageName, totalWorkItems, ...
                            summaryStep, detailText)
completedWorkItems = progressState('completedWorkItems') + 1;
progressState('completedWorkItems') = completedWorkItems;
percentDone = 100 * completedWorkItems / totalWorkItems;

if matRad_isDetailedProgress(cfg)
    detailText = matRad_normalizeDetailText(detailText);
    if isempty(detailText)
        detailText = 'completed work item';
    end
    matRadCfg.dispInfo('matRad: %s progress %d/%d (%.0f%%): %s.\n', ...
                       stageName, completedWorkItems, totalWorkItems, percentDone, ...
                       detailText);
elseif completedWorkItems == 1 || completedWorkItems == totalWorkItems || ...
        mod(completedWorkItems, summaryStep) == 0
    matRadCfg.dispInfo('matRad: %s progress %d/%d (%.0f%%).\n', ...
                       stageName, completedWorkItems, totalWorkItems, percentDone);
end

drawnow('limitrate');
end

function matRad_cleanupProgressListener(progressListener)
drawnow;
if ~isempty(progressListener)
    delete(progressListener);
end
end

function detailText = matRad_normalizeDetailText(detailText)
if nargin < 1 || isempty(detailText)
    detailText = '';
elseif isstring(detailText) && isscalar(detailText)
    detailText = char(detailText);
elseif ~ischar(detailText)
    detailText = '';
end
end

function detailed = matRad_isDetailedProgress(cfg)
detailed = isfield(cfg, 'ProgressLevel') && strcmp(cfg.ProgressLevel, 'detailed');
end

function matRad_warnDataQueueUnavailable(matRadCfg, stageName)
persistent dataQueueWarningShown

if isempty(dataQueueWarningShown) || ~dataQueueWarningShown
    matRadCfg.dispInfo(['matRad: Parallel progress for %s is unavailable ', ...
                        'because parallel.pool.DataQueue is unavailable.\n'], stageName);
    dataQueueWarningShown = true;
end
end

function matRad_noopUpdate(varargin)
end

function matRad_noopCleanup()
end
