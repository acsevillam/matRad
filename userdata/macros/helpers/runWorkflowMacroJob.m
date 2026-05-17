function jobResult = runWorkflowMacroJob(job,varargin)
% runWorkflowMacroJob Run workflow macro parameter sets in series.

jobPlan = resolveWorkflowMacroJob(job,varargin{:});

jobResult = struct();
jobResult.id = jobPlan.id;
jobResult.description = jobPlan.description;
jobResult.stopOnError = jobPlan.stopOnError;
jobResult.plan = jobPlan;
jobResult.runs = repmat(emptyJobRunResult(),1,numel(jobPlan.runs));

for runIx = 1:numel(jobPlan.runs)
    runPlan = jobPlan.runs(runIx);
    fprintf('[%d/%d] Running macro job item "%s": %s\n', ...
        runIx,numel(jobPlan.runs),runPlan.label,runPlan.specId);

    runResult = emptyJobRunResult();
    runResult.index = runPlan.index;
    runResult.label = runPlan.label;
    runResult.specId = runPlan.specId;
    runResult.status = 'running';
    runResult.startedAt = char(datetime('now'));

    try
        runResult.result = runWorkflowMacroSpec( ...
            runPlan.specId,runPlan.args{:});
        runResult.status = 'completed';
    catch ME
        runResult.status = 'failed';
        runResult.error = struct( ...
            'identifier',ME.identifier, ...
            'message',ME.message);
        jobResult.runs(runIx) = runResult;
        if logical(jobPlan.stopOnError)
            rethrow(ME);
        end
    end

    runResult.finishedAt = char(datetime('now'));
    jobResult.runs(runIx) = runResult;
end

end

function runResult = emptyJobRunResult()

runResult = struct();
runResult.index = [];
runResult.label = '';
runResult.specId = '';
runResult.status = 'pending';
runResult.startedAt = '';
runResult.finishedAt = '';
runResult.error = [];
runResult.result = [];

end
