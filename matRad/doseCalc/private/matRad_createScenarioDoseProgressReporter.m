function reporter = matRad_createScenarioDoseProgressReporter(matRad_cfg,cfg, ...
    stageName,totalWorkItems,useDataQueue)
% matRad_createScenarioDoseProgressReporter creates streaming progress output
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

reporter = struct('queue',[],'enabled',false,'update',@noopUpdate, ...
    'cleanup',@noopCleanup);

if nargin < 4 || isempty(totalWorkItems) || totalWorkItems <= 0
    return;
end

totalWorkItems = double(totalWorkItems);
completedWorkItems = 0;
summaryStep = max(1,ceil(totalWorkItems/10));
progressListener = [];

if useDataQueue
    if exist('parallel.pool.DataQueue','class') ~= 8
        warnDataQueueUnavailable(matRad_cfg,stageName);
        return;
    end

    progressQueue = parallel.pool.DataQueue;
    progressListener = afterEach(progressQueue,@logProgress);
    reporter.queue = progressQueue;
    reporter.cleanup = @cleanupProgressListener;
else
    reporter.update = @logProgress;
end

reporter.enabled = true;

    function logProgress(detailText)
        completedWorkItems = completedWorkItems + 1;
        percentDone = 100*completedWorkItems/totalWorkItems;

        if isDetailedProgress(cfg)
            detailText = normalizeDetailText(detailText);
            if isempty(detailText)
                detailText = 'completed work item';
            end
            matRad_cfg.dispInfo('matRad: %s progress %d/%d (%.0f%%): %s.\n', ...
                stageName,completedWorkItems,totalWorkItems,percentDone, ...
                detailText);
        elseif completedWorkItems == 1 || completedWorkItems == totalWorkItems || ...
                mod(completedWorkItems,summaryStep) == 0
            matRad_cfg.dispInfo('matRad: %s progress %d/%d (%.0f%%).\n', ...
                stageName,completedWorkItems,totalWorkItems,percentDone);
        end

        drawnow('limitrate');
    end

    function cleanupProgressListener()
        drawnow;
        if ~isempty(progressListener)
            delete(progressListener);
        end
    end
end

function detailText = normalizeDetailText(detailText)
if nargin < 1 || isempty(detailText)
    detailText = '';
elseif isstring(detailText) && isscalar(detailText)
    detailText = char(detailText);
elseif ~ischar(detailText)
    detailText = '';
end
end

function detailed = isDetailedProgress(cfg)
detailed = isfield(cfg,'ProgressLevel') && strcmp(cfg.ProgressLevel,'detailed');
end

function warnDataQueueUnavailable(matRad_cfg,stageName)
persistent dataQueueWarningShown

if isempty(dataQueueWarningShown) || ~dataQueueWarningShown
    matRad_cfg.dispInfo(['matRad: Parallel progress for %s is unavailable ', ...
        'because parallel.pool.DataQueue is unavailable.\n'],stageName);
    dataQueueWarningShown = true;
end
end

function noopUpdate(varargin)
end

function noopCleanup()
end
