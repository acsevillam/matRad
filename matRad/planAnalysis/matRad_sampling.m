function [caSampRes, mSampDose, pln, resultGUInomScen] = matRad_sampling(ct, stf, cst, pln, w, structSel, multScen, varargin)
% matRad_randomSampling enables sampling multiple treatment scenarios
%
% call:
%   [caSampRes,mSampDose,pln,resultGUInomScen] =
%       matRad_sampling(ct,stf,cst,pln,w,structSel,multScen)
%   [...] = matRad_sampling(...,'dvhDoseWindow',dvhDoseWindow)
%   [...] = matRad_sampling(...,'dvhDoseGrid',dvhDoseGrid)
%
% input:
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
%   dvhDoseWindow: (optional) per-fraction dose window for the common DVH grid
%   dvhDoseGrid:   (optional) explicit per-fraction DVH dose grid
%   autoLimitWorkers: (optional) reduce the parallel pool size based on
%               available system memory and estimated memory per worker
%   workerMemorySafetyFactor: (optional) safety factor for worker memory
%   memoryReserveFraction: (optional) system memory fraction kept in reserve
%   minWorkerMemoryBytes: (optional) lower bound for per-worker memory
%   workerUpperBound: (optional) explicit upper bound for workers
%   calibrateWorkerMemory: (optional) calculate the first sampled scenario
%               serially to calibrate the worker memory estimate when
%               observed process memory is available
%   allowCalibrationToReduceWorkerMemory: (optional) allow a reliable
%               calibration to lower the static worker estimate
%   calibratedMinForwardDoseWorkerMemoryBytes: (optional) lower bound for a
%               calibrated transient forward dose worker memory estimate
%   minForwardDoseWorkerMemoryBytes: (optional) lower bound for transient
%               forward dose worker memory during sampling
%
% output:
%   caSampRes:         cell array of sampling results depicting plan parameter
%   mSampDose:         matrix holding the sampled doses, each row corresponds to
%                      one dose sample
%   pln:               matRad pln struct containing sampling information
%   resultGUInomScen:  resultGUI struct of the nominal scenario
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
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

if nargin < 7
    multScen = [];
elseif ischar(multScen) || (isstring(multScen) && isscalar(multScen))
    varargin = [{char(multScen)} varargin];
    multScen = [];
end

parser = inputParser;
parser.addParameter('dvhDoseWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x) >= 2));
parser.addParameter('dvhDoseGrid', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
parser.addParameter('autoLimitWorkers', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parser.addParameter('workerMemorySafetyFactor', 1.2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
parser.addParameter('memoryReserveFraction', 0.10, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
parser.addParameter('minWorkerMemoryBytes', 4 * 1024^3, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parser.addParameter('workerUpperBound', [], @(x) isempty(x) || ...
                    (isnumeric(x) && isscalar(x) && isfinite(x) && round(x) == x && x >= 1));
parser.addParameter('calibrateWorkerMemory', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parser.addParameter('allowCalibrationToReduceWorkerMemory', true, ...
                    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parser.addParameter('calibratedMinForwardDoseWorkerMemoryBytes', 4 * 1024^3, ...
                    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parser.addParameter('minForwardDoseWorkerMemoryBytes', 16 * 1024^3, ...
                    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parser.parse(varargin{:});

refScen = matRad_resolveSamplingReferenceScenario(ct);
cstEval = matRad_buildSamplingReferenceCst(cst, refScen);

% save nonSampling pln for nominal scenario calculation and add dummy fields
plnNominal = pln;
% create nominal scenario
plnNominal.multScen = matRad_ScenarioModel.create('nomScen', ct);
plnNominal.multScen.ctScenProb = [refScen 1];

if ct.numOfCtScen > 1 && isempty(multScen)
    matRad_cfg.dispWarning(['No explicit multi-CT sampling scenario model was provided; ', ...
                            'random sampling will use the reference CT scenario only.\n']);
end

% either use existing multScen struct or create new one
if exist('multScen', 'var') && ~isempty(multScen)
    pln.multScen = multScen;
else
    % create random scenarios for sampling
    pln.multScen = matRad_RandomScenarios(ct);
    pln.multScen.ctScenProb = [refScen 1];
    pln.multScen.nSamples = matRad_cfg.defaults.samplingScenarios;
end

scenarioIds = pln.multScen.scenarioIds();
numScenarios = numel(scenarioIds);

matRad_cfg.dispInfo('Using %d samples in total \n', numScenarios);

doseMapping = matRad_resolveSamplingDoseMapping(ct, pln.multScen, refScen);
if doseMapping.enabled
    matRad_cfg.dispInfo('matRad: Mapping sampled multi-CT dose cubes to reference CT scenario %d before analysis.\n', ...
                        refScen);
end

V = [];
% define voxels for sampling
if ~exist('structSel', 'var') || sum(size(structSel)) == 0
    for i = 1:size(cstEval, 1)
        if ~isempty(cstEval{i, 4}{1})
            V = [V cstEval{i, 4}(1)];
        end
    end
else
    for i = 1:size(cstEval, 1)
        for j = 1:numel(structSel)
            if strcmp(structSel{j}, cstEval{i, 2})
                V = [V cstEval{i, 4}{1}];
            end
        end
    end
end

if isempty(V)
    matRad_cfg.dispError('No voxels selected for sampling.');
end

% final voxel subset for sampling
subIx = unique(vertcat(V{:}));

% disable structures for DVH plotting which are not completely in subIx
for i = 1:size(cstEval, 1)
    if ~all(ismember(cstEval{i, 4}{1}, subIx))
        cstEval{i, 5}.Visible = false;
    end
end

% define variable for storing scenario doses
mSampDose   = single(zeros(numel(subIx), numScenarios, 1));
StorageInfo = whos('mSampDose');
matRad_cfg.dispInfo('matRad: Realizations variable will need: %f GB \n', StorageInfo.bytes / 1e9);

%% calculate nominal scenario
nomScenTimer     = tic;
resultGUInomScen = matRad_calcDoseForward(ct, cst, stf, plnNominal, w);
nomScenTime      = toc(nomScenTimer);
matRad_cfg.dispInfo('Finished nominal Scenario Calculation. Computation time: %f h \n', round(nomScenTime / 3600));

refVol = [2 5 50 95 98];
quantityVis = matRad_resolveDoseAnalysisQuantity(resultGUInomScen, pln, '');

refGy = linspace(0, max(resultGUInomScen.(quantityVis)(:)), 6);

dvhPoints = matRad_resolveSamplingDvhDoseGrid(resultGUInomScen.(quantityVis), ...
                                              parser.Results.dvhDoseWindow, ...
                                              parser.Results.dvhDoseGrid);

resultGUInomScen.dvh = matRad_calcDVH(cstEval, resultGUInomScen.(quantityVis), 'cum', dvhPoints);
dvhPoints            = resultGUInomScen.dvh(1).doseGrid;
nomQi                = matRad_calcQualityIndicators(cstEval, pln, resultGUInomScen.(quantityVis), refGy, refVol);

resultGUInomScen.qi  = nomQi;
resultGUInomScen.cst = cstEval;
resultGUInomScen.analysisQuantity = quantityVis;
resultGUInomScen.refScen = refScen;
resultGUInomScen.evaluationModeBase = 'perFraction';

samplingContext = matRad_buildSamplingContext(ct, stf, cst, pln, w, cstEval, subIx, ...
                                              dvhPoints, refGy, refVol, doseMapping, quantityVis);
samplingResourceConfig = matRad_samplingResourceConfig(parser.Results);
samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingContext, numScenarios, ...
                                                       StorageInfo.bytes, ...
                                                       samplingResourceConfig);

%% perform sampling
samplingOutput = matRad_runSampling(samplingContext, scenarioIds, mSampDose, ...
                                    nomScenTime, samplingMemoryEstimate, ...
                                    samplingResourceConfig, matRad_cfg);
mSampDose = samplingOutput.mSampDose;
caSampRes = samplingOutput.caSampRes;
samplingMemoryEstimate = samplingOutput.samplingMemoryEstimate;

%% add subindices
pln.subIx        = subIx;
pln.samplingReferenceCtScen = refScen;
pln.samplingDoseMapping = doseMapping;
pln.samplingMemoryEstimate = samplingMemoryEstimate;

end

function samplingOutput = matRad_runSampling(samplingContext, scenarioIds, mSampDose, ...
                                             nomScenTime, samplingMemoryEstimate, ...
                                             samplingResourceConfig, matRadCfg)

numScenarios = numel(scenarioIds);
sampleResults = cell(1, numScenarios);

[mSampDose, sampleResults, samplingMemoryEstimate] = ...
    matRad_calibrateSamplingWorkerMemoryIfNeeded(samplingContext, scenarioIds, ...
                                                 mSampDose, sampleResults, ...
                                                 samplingMemoryEstimate, ...
                                                 samplingResourceConfig, matRadCfg);
remainingScenarioIx = matRad_remainingSamplingScenarioIndices(sampleResults);

poolContext = matRad_prepareSamplingPool(matRadCfg, numel(remainingScenarioIx), ...
                                         samplingMemoryEstimate, samplingResourceConfig);
samplingMemoryEstimate = poolContext.samplingMemoryEstimate;

if isempty(remainingScenarioIx)
    caSampRes = [sampleResults{:}];
elseif poolContext.useParallel
    [mSampDose, sampleResults] = matRad_executeParallelSamplingScenarios( ...
        samplingContext, scenarioIds, remainingScenarioIx, mSampDose, sampleResults, ...
        nomScenTime, poolContext.pool, poolContext.parallelPlan, matRadCfg);
    caSampRes = [sampleResults{:}];
else
    [mSampDose, sampleResults] = matRad_executeSerialSamplingScenarios( ...
        samplingContext, scenarioIds, remainingScenarioIx, mSampDose, sampleResults, ...
        nomScenTime, matRadCfg);
    caSampRes = [sampleResults{:}];
end

if numel(caSampRes) ~= numScenarios
    matRadCfg.dispError('Sampling did not produce the expected number of scenario results.');
end

samplingOutput = struct('mSampDose', mSampDose, ...
                        'caSampRes', caSampRes, ...
                        'samplingMemoryEstimate', samplingMemoryEstimate);

end

function poolContext = matRad_prepareSamplingPool(matRadCfg, numScenarios, samplingMemoryEstimate, samplingResourceConfig)
poolContext = struct('parallelToolboxLicensed', false, ...
                     'useParallel', false, ...
                     'pool', [], ...
                     'parallelPlan', [], ...
                     'samplingMemoryEstimate', samplingMemoryEstimate);

if numScenarios < 2
    samplingMemoryEstimate.parallelPlan = matRad_planSamplingParallelTasks( ...
        numScenarios, samplingMemoryEstimate, samplingResourceConfig, matRadCfg);
    samplingMemoryEstimate.parallelPlan.useParallel = false;
    samplingMemoryEstimate.parallelPlan.fallbackReason = 'taskCount';
    poolContext.samplingMemoryEstimate = samplingMemoryEstimate;
    matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, ...
                                     samplingMemoryEstimate.parallelPlan, matRadCfg);
    return
end

try
    [poolContext.parallelToolboxLicensed, ~] = license('checkout', 'Distrib_Computing_Toolbox');
    if ~poolContext.parallelToolboxLicensed
        matRadCfg.dispWarning('Could not check out parallel computing toolbox. \n');
        matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, [], matRadCfg);
        return
    end

    if samplingResourceConfig.autoLimitWorkers
        parallelPlan = matRad_planSamplingParallelTasks(numScenarios, ...
                                                        samplingMemoryEstimate, ...
                                                        samplingResourceConfig, ...
                                                        matRadCfg);
        samplingMemoryEstimate.parallelPlan = parallelPlan;
        samplingMemoryEstimate.workerLimit = parallelPlan;
        matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, parallelPlan, matRadCfg);
        if ~parallelPlan.useParallel
            poolContext.samplingMemoryEstimate = samplingMemoryEstimate;
            return
        end
        poolContext.pool = matRad_configureParallelPoolSize( ...
            parallelPlan.workerUpperBound, 'sampling', matRadCfg);
        poolContext.useParallel = true;
        poolContext.parallelPlan = parallelPlan;
    else
        matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, [], matRadCfg);
        poolContext.pool = gcp(); % If no pool exists, create one.
        poolContext.useParallel = true;
    end

    poolContext.samplingMemoryEstimate = samplingMemoryEstimate;

    if poolContext.useParallel && isempty(poolContext.pool)
        matRadCfg.dispError('matRad: Could not start valid parallel pool. Please check your parallel computing toolbox installation. \n');
    end
catch
    matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, [], matRadCfg);
    poolContext.parallelToolboxLicensed = false;
    poolContext.useParallel = false;
    poolContext.pool = [];
    poolContext.samplingMemoryEstimate = samplingMemoryEstimate;
end
end

function [mSampDose, caSampRes] = matRad_executeParallelSamplingScenarios(samplingContext, scenarioIds, ...
                                                                          scenarioIndices, mSampDose, sampleResults, ...
                                                                          nomScenTime, pool, parallelPlan, matRadCfg)

if isempty(pool)
    poolSize = 1;
else
    poolSize = pool.NumWorkers;
end

logLevel = matRadCfg.logLevel;
numScenarios = numel(scenarioIndices);
totCompTime = ceil(numScenarios / poolSize) * nomScenTime * 1.35;
matRad_logSamplingTimeEstimate(totCompTime, matRadCfg);
matRadCfg.dispInfo(['matRad: Sampling parallel plan uses chunks of %d ', ...
                    'scenario(s) with %d worker(s).\n'], ...
                   parallelPlan.chunkSize, poolSize);

progressEnabled = matRad_startSamplingProgress(numScenarios, logLevel, matRadCfg);
[progressQueue, progressListener] = ...
    matRad_createSamplingProgressQueue(numScenarios, logLevel, matRadCfg); %#ok<ASGLU>

chunks = matRad_buildSamplingChunks(numScenarios, parallelPlan.chunkSize);
for chunkIx = 1:numel(chunks)
    chunkLocalIx = chunks{chunkIx};
    chunkScenarioIx = scenarioIndices(chunkLocalIx);
    chunkDoseColumns = single(zeros(size(mSampDose, 1), numel(chunkScenarioIx)));
    chunkResults = cell(1, numel(chunkScenarioIx));
    chunkScenarioIds = scenarioIds(chunkScenarioIx);
    matRadCfg.dispInfo('matRad: Sampling chunk %d/%d with %d scenario(s).\n', ...
                       chunkIx, numel(chunks), numel(chunkScenarioIx));
    parfor localIx = 1:numel(chunkScenarioIx)
        matRadCfgWorker = MatRad_Config.instance();
        matRadCfgWorker.logLevel = logLevel;

        [chunkDoseColumns(:, localIx), chunkResults{localIx}] = ...
            matRad_calculateSamplingScenario(samplingContext, chunkScenarioIds(localIx));

        if progressEnabled && logLevel > 2
            parfor_progress;
        end
        if ~isempty(progressQueue)
            send(progressQueue, 1);
        end
    end
    mSampDose(:, chunkScenarioIx) = chunkDoseColumns;
    sampleResults(chunkScenarioIx) = chunkResults;
    clear chunkDoseColumns chunkResults;
end

if progressEnabled && logLevel > 2
    parfor_progress(0);
end
caSampRes = sampleResults;

end

function [mSampDose, caSampRes] = matRad_executeSerialSamplingScenarios(samplingContext, scenarioIds, ...
                                                                        scenarioIndices, mSampDose, sampleResults, ...
                                                                        nomScenTime, matRadCfg)

numScenarios = numel(scenarioIndices);
totCompTime = numScenarios * nomScenTime * 1.1;
matRad_logSamplingTimeEstimate(totCompTime, matRadCfg);

for localIx = 1:numScenarios
    scenarioIx = scenarioIndices(localIx);
    [mSampDose(:, scenarioIx), sampleResults{scenarioIx}] = ...
        matRad_calculateSamplingScenario(samplingContext, scenarioIds(scenarioIx));

    if matRadCfg.logLevel > 2
        matRad_logSamplingProgress(localIx, numScenarios, matRadCfg);
        matRad_progress(localIx, numScenarios);
    end
end
caSampRes = sampleResults;

end

function [mSampDose, sampleResults, samplingMemoryEstimate] = ...
    matRad_calibrateSamplingWorkerMemoryIfNeeded(samplingContext, scenarioIds, ...
                                                 mSampDose, sampleResults, ...
                                                 samplingMemoryEstimate, ...
                                                 samplingResourceConfig, matRadCfg)
samplingMemoryEstimate.calibration = matRad_emptySamplingCalibration();
if ~samplingResourceConfig.calibrateWorkerMemory || isempty(scenarioIds)
    return
end

calibrationScenarioIx = 1;
beforeBytes = matRad_currentProcessMemoryBytes();
timer = tic;
[mSampDose(:, calibrationScenarioIx), sampleResults{calibrationScenarioIx}] = ...
    matRad_calculateSamplingScenario(samplingContext, scenarioIds(calibrationScenarioIx));
elapsedSeconds = toc(timer);
afterBytes = matRad_currentProcessMemoryBytes();

calibration = matRad_emptySamplingCalibration();
calibration.enabled = true;
calibration.scenarioIndex = calibrationScenarioIx;
calibration.scenarioId = scenarioIds(calibrationScenarioIx);
calibration.elapsedSeconds = elapsedSeconds;
calibration.beforeProcessMemoryBytes = beforeBytes;
calibration.afterProcessMemoryBytes = afterBytes;
calibration.reusedScenario = true;
calibration.staticWorkerBytes = samplingMemoryEstimate.rawWorkerBytes;
calibration.allowReduction = samplingResourceConfig.allowCalibrationToReduceWorkerMemory;
calibration.calibratedMinWorkerBytes = ...
    samplingResourceConfig.calibratedMinForwardDoseWorkerMemoryBytes;

measurementReliable = ~isempty(beforeBytes) && ~isempty(afterBytes) && ...
    isfinite(beforeBytes) && isfinite(afterBytes) && afterBytes >= beforeBytes;
calibration.measurementReliable = measurementReliable;

if measurementReliable
    measuredBytes = max(afterBytes - beforeBytes, 0) + ...
        samplingMemoryEstimate.resultBytesPerTask;
    calibration.measuredWorkerBytes = measuredBytes;
    calibratedWorkerBytes = max(measuredBytes, ...
        samplingResourceConfig.calibratedMinForwardDoseWorkerMemoryBytes);
    calibration.calibratedWorkerBytes = calibratedWorkerBytes;
    calibration.action = 'kept';

    if calibratedWorkerBytes > samplingMemoryEstimate.rawWorkerBytes
        samplingMemoryEstimate.rawWorkerBytes = calibratedWorkerBytes;
        samplingMemoryEstimate.estimateBasis = 'samplingForwardDoseCalibrated';
        calibration.usedForPlanning = true;
        calibration.action = 'raised';
    elseif samplingResourceConfig.allowCalibrationToReduceWorkerMemory && ...
            calibratedWorkerBytes < samplingMemoryEstimate.rawWorkerBytes
        samplingMemoryEstimate.rawWorkerBytes = calibratedWorkerBytes;
        samplingMemoryEstimate.estimateBasis = 'samplingForwardDoseCalibrated';
        calibration.usedForPlanning = true;
        calibration.action = 'lowered';
    end
else
    calibration.action = 'unreliable';
end
samplingMemoryEstimate.calibration = calibration;

matRad_logSamplingCalibration(calibration, samplingMemoryEstimate, matRadCfg);
matRad_logSamplingProgress(1, numel(scenarioIds), matRadCfg);
end

function calibration = matRad_emptySamplingCalibration()
calibration = struct();
calibration.enabled = false;
calibration.scenarioIndex = [];
calibration.scenarioId = [];
calibration.elapsedSeconds = [];
calibration.beforeProcessMemoryBytes = [];
calibration.afterProcessMemoryBytes = [];
calibration.measuredWorkerBytes = [];
calibration.staticWorkerBytes = [];
calibration.calibratedWorkerBytes = [];
calibration.calibratedMinWorkerBytes = [];
calibration.measurementReliable = false;
calibration.allowReduction = false;
calibration.usedForPlanning = false;
calibration.reusedScenario = false;
calibration.action = 'disabled';
end

function matRad_logSamplingCalibration(calibration, samplingMemoryEstimate, matRadCfg)
scenarioText = matRad_formatScenarioId(calibration.scenarioId);
switch char(calibration.action)
    case 'raised'
        matRadCfg.dispInfo(['matRad: Sampling worker memory calibration raised ', ...
                            'the worker estimate from %s to %s using scenario %s.\n'], ...
                           matRad_formatSamplingBytes(calibration.staticWorkerBytes), ...
                           matRad_formatSamplingBytes(samplingMemoryEstimate.rawWorkerBytes), ...
                           scenarioText);
    case 'lowered'
        matRadCfg.dispInfo(['matRad: Sampling worker memory calibration lowered ', ...
                            'the worker estimate from %s to %s using scenario %s ', ...
                            '(measured %s, calibrated floor %s).\n'], ...
                           matRad_formatSamplingBytes(calibration.staticWorkerBytes), ...
                           matRad_formatSamplingBytes(samplingMemoryEstimate.rawWorkerBytes), ...
                           scenarioText, ...
                           matRad_formatSamplingBytes(calibration.measuredWorkerBytes), ...
                           matRad_formatSamplingBytes(calibration.calibratedMinWorkerBytes));
    case 'unreliable'
        matRadCfg.dispInfo(['matRad: Sampling worker memory calibration reused ', ...
                            'scenario %s but kept the static worker estimate because ', ...
                            'process memory measurement was unavailable or non-monotonic.\n'], ...
                           scenarioText);
    otherwise
        matRadCfg.dispInfo(['matRad: Sampling worker memory calibration reused ', ...
                            'scenario %s and kept the static worker estimate %s.\n'], ...
                           scenarioText, ...
                           matRad_formatSamplingBytes(samplingMemoryEstimate.rawWorkerBytes));
end
end

function scenarioIndices = matRad_remainingSamplingScenarioIndices(sampleResults)
scenarioIndices = find(cellfun(@isempty, sampleResults));
end

function parallelPlan = matRad_planSamplingParallelTasks(numScenarios, ...
                                                         samplingMemoryEstimate, ...
                                                         samplingResourceConfig, ...
                                                         matRadCfg)
if numScenarios < 1
    numScenarios = 1;
end
parallelPlan = matRad_planMemoryLimitedParallelTasks( ...
    numScenarios, samplingMemoryEstimate.rawWorkerBytes, ...
    'stageName', 'sampling forward dose', ...
    'resultBytesPerTask', samplingMemoryEstimate.resultBytesPerTask, ...
    'accumulatorBytes', samplingMemoryEstimate.mainProcessOutputBytes, ...
    'reserveFraction', samplingResourceConfig.memoryReserveFraction, ...
    'safetyFactor', samplingResourceConfig.workerMemorySafetyFactor, ...
    'minWorkerMemoryBytes', samplingResourceConfig.minWorkerMemoryBytes, ...
    'workerUpperBound', samplingResourceConfig.workerUpperBound, ...
    'matRadCfg', matRadCfg);
end

function chunks = matRad_buildSamplingChunks(numScenarios, chunkSize)
chunkSize = max(1, min(numScenarios, floor(chunkSize)));
starts = 1:chunkSize:numScenarios;
chunks = cell(1, numel(starts));
for chunkIx = 1:numel(starts)
    chunkStart = starts(chunkIx);
    chunkEnd = min(numScenarios, chunkStart + chunkSize - 1);
    chunks{chunkIx} = chunkStart:chunkEnd;
end
end

function bytes = matRad_currentProcessMemoryBytes()
bytes = [];
try
    pid = feature('getpid');
    [status, output] = system(sprintf('ps -o rss= -p %d', pid));
    rssKbText = regexp(output, '\d+', 'match', 'once');
    if status == 0 && ~isempty(rssKbText)
        bytes = str2double(rssKbText) * 1024;
    end
catch
    bytes = [];
end
end

function text = matRad_formatScenarioId(scenarioId)
if isnumeric(scenarioId)
    text = mat2str(scenarioId);
else
    text = char(string(scenarioId));
end
end

function matRad_logSamplingTimeEstimate(totCompTime, matRadCfg)
try
    msg = ['Approximate Total calculation time: ', num2str(round(totCompTime / 3600)), ...
           'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n'];
catch
    msg = ['Approximate Total calculation time: ', num2str(round(totCompTime / 3600)), '\n'];
end
matRadCfg.dispInfo(msg);
end

function progressEnabled = matRad_startSamplingProgress(numScenarios, logLevel, matRadCfg)
if exist('parfor_progress', 'file') == 2 && logLevel > 2
    progressEnabled = true;
    parfor_progress(numScenarios);
else
    msg = ['matRad: Consider downloading parfor_progress function from the matlab central fileexchange ', ...
           'to get feedback from parfor loop.\n'];
    matRadCfg.dispInfo(msg);
    progressEnabled = false;
end
end

function [progressQueue, progressListener] = matRad_createSamplingProgressQueue(numScenarios, logLevel, matRadCfg)
progressQueue = [];
progressListener = [];
if logLevel <= 2 || exist('parallel.pool.DataQueue', 'class') ~= 8
    return
end

progressState = containers.Map({'finishedScenarios'}, {0});
progressQueue = parallel.pool.DataQueue;
progressListener = afterEach(progressQueue, ...
                             @(~) matRad_updateSamplingProgress(progressState, numScenarios, matRadCfg));
end

function matRad_updateSamplingProgress(progressState, numScenarios, matRadCfg)
finishedScenarios = progressState('finishedScenarios') + 1;
progressState('finishedScenarios') = finishedScenarios;
matRad_logSamplingProgress(finishedScenarios, numScenarios, matRadCfg);
end

function matRad_logSamplingProgress(finishedScenarios, totalScenarios, matRadCfg)
matRadCfg.dispInfo('matRad: Sampling progress: %d/%d scenarios.\n', ...
                   finishedScenarios, totalScenarios);
end

function samplingContext = matRad_buildSamplingContext(ct, stf, cst, pln, w, cstEval, subIx, ...
                                                       dvhPoints, refGy, refVol, doseMapping, quantity)

samplingContext = struct();
samplingContext.ct = ct;
samplingContext.stf = stf;
samplingContext.cst = cst;
samplingContext.pln = pln;
samplingContext.w = w;
samplingContext.cstEval = cstEval;
samplingContext.subIx = subIx;
samplingContext.dvhPoints = dvhPoints;
samplingContext.refGy = refGy;
samplingContext.refVol = refVol;
samplingContext.doseMapping = doseMapping;
samplingContext.quantity = quantity;

end

function samplingResourceConfig = matRad_samplingResourceConfig(parserResults)
samplingResourceConfig = struct();
samplingResourceConfig.autoLimitWorkers = logical(parserResults.autoLimitWorkers);
samplingResourceConfig.workerMemorySafetyFactor = parserResults.workerMemorySafetyFactor;
samplingResourceConfig.memoryReserveFraction = parserResults.memoryReserveFraction;
samplingResourceConfig.minWorkerMemoryBytes = parserResults.minWorkerMemoryBytes;
samplingResourceConfig.workerUpperBound = parserResults.workerUpperBound;
samplingResourceConfig.calibrateWorkerMemory = logical(parserResults.calibrateWorkerMemory);
samplingResourceConfig.allowCalibrationToReduceWorkerMemory = ...
    logical(parserResults.allowCalibrationToReduceWorkerMemory);
samplingResourceConfig.calibratedMinForwardDoseWorkerMemoryBytes = ...
    parserResults.calibratedMinForwardDoseWorkerMemoryBytes;
samplingResourceConfig.minForwardDoseWorkerMemoryBytes = ...
    parserResults.minForwardDoseWorkerMemoryBytes;
end

function samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingContext, numScenarios, sampleDoseStorageBytes, samplingResourceConfig)
inputBytes = matRad_variableBytes(samplingContext.ct) + ...
             matRad_variableBytes(samplingContext.stf) + ...
             matRad_variableBytes(samplingContext.cst) + ...
             matRad_variableBytes(samplingContext.pln) + ...
             matRad_variableBytes(samplingContext.w) + ...
             matRad_variableBytes(samplingContext.cstEval) + ...
             matRad_variableBytes(samplingContext.subIx) + ...
             matRad_variableBytes(samplingContext.dvhPoints) + ...
             matRad_variableBytes(samplingContext.refGy) + ...
             matRad_variableBytes(samplingContext.refVol) + ...
             matRad_variableBytes(samplingContext.doseMapping);
sampleDoseBytes = numel(samplingContext.subIx) * 4;
sampleResultBytes = 64 * 1024;
doseCubeBytes = prod(double(samplingContext.ct.cubeDim)) * 8;
doseMappingWorkspaceBytes = doseCubeBytes * double(samplingContext.doseMapping.enabled);
[workerBytes, workerDetails] = matRad_estimateSamplingForwardDoseWorkerBytes( ...
    samplingContext, inputBytes, doseCubeBytes, doseMappingWorkspaceBytes, ...
    samplingResourceConfig.minForwardDoseWorkerMemoryBytes);

samplingMemoryEstimate = struct();
samplingMemoryEstimate.estimateBasis = workerDetails.estimateBasis;
samplingMemoryEstimate.numSamples = numScenarios;
samplingMemoryEstimate.numVoxels = numel(samplingContext.subIx);
samplingMemoryEstimate.inputBytes = inputBytes;
samplingMemoryEstimate.doseResultProxyBytes = doseCubeBytes;
samplingMemoryEstimate.sampleDoseBytes = sampleDoseBytes;
samplingMemoryEstimate.sampleResultBytes = sampleResultBytes;
samplingMemoryEstimate.resultBytesPerTask = sampleDoseBytes + sampleResultBytes;
samplingMemoryEstimate.doseMappingWorkspaceBytes = doseMappingWorkspaceBytes;
samplingMemoryEstimate.forwardDoseWorkerDetails = workerDetails;
samplingMemoryEstimate.rawWorkerBytes = workerBytes;
samplingMemoryEstimate.sampleDoseStorageBytes = sampleDoseStorageBytes;
samplingMemoryEstimate.sampleResultStorageBytes = sampleResultBytes * numScenarios;
samplingMemoryEstimate.mainProcessOutputBytes = sampleDoseStorageBytes + ...
    samplingMemoryEstimate.sampleResultStorageBytes;
samplingMemoryEstimate.workerLimit = [];
samplingMemoryEstimate.parallelPlan = [];
samplingMemoryEstimate.calibration = matRad_emptySamplingCalibration();
end

function bytes = matRad_variableBytes(value)
info = whos('value');
bytes = info.bytes;
end

function matRad_logSamplingMemoryEstimate(samplingMemoryEstimate, parallelPlan, matRadCfg)
msg = ['matRad: Estimated sampling memory: output ', ...
       matRad_formatSamplingBytes(samplingMemoryEstimate.mainProcessOutputBytes), ...
       ', raw worker ', matRad_formatSamplingBytes(samplingMemoryEstimate.rawWorkerBytes), ...
       ', result ', matRad_formatSamplingBytes(samplingMemoryEstimate.resultBytesPerTask), ...
       '/scenario'];
if ~isempty(parallelPlan)
    msg = [msg, ', memory-limited workers ', num2str(parallelPlan.workerUpperBound), ...
           ' (usable ', matRad_formatSamplingBytes(parallelPlan.memoryBudgetBytes), ...
           ', effective worker ', matRad_formatSamplingBytes(parallelPlan.workerBytes), ...
           ', chunk ', num2str(parallelPlan.chunkSize)];
    if isfield(parallelPlan, 'memoryInfo') && isstruct(parallelPlan.memoryInfo) && ...
            isfield(parallelPlan.memoryInfo, 'source') && ~isempty(parallelPlan.memoryInfo.source)
        msg = [msg, ', source ', parallelPlan.memoryInfo.source];
    end
    if isfield(parallelPlan, 'allocatedCpuCount') && ~isempty(parallelPlan.allocatedCpuCount)
        msg = [msg, ', allocated CPUs ', num2str(parallelPlan.allocatedCpuCount)];
    end
    if isfield(parallelPlan, 'fallbackReason') && ~isempty(parallelPlan.fallbackReason)
        msg = [msg, ', fallback ', parallelPlan.fallbackReason];
    end
    msg = [msg, ')'];
end
matRadCfg.dispInfo([msg, '.\n']);
end

function text = matRad_formatSamplingBytes(bytes)
if isempty(bytes) || ~isnumeric(bytes) || ~isscalar(bytes) || ~isfinite(bytes)
    text = 'unavailable';
    return
end
units = {'B', 'KiB', 'MiB', 'GiB', 'TiB'};
value = double(bytes);
unitIx = 1;
while value >= 1024 && unitIx < numel(units)
    value = value / 1024;
    unitIx = unitIx + 1;
end
text = sprintf('%.2f %s', value, units{unitIx});
end

function refScen = matRad_resolveSamplingReferenceScenario(ct)

matRad_cfg = MatRad_Config.instance();

if isfield(ct, 'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end

validRefScen = isnumeric(refScen) && isscalar(refScen) && isfinite(refScen) && ...
    round(refScen) == refScen && refScen >= 1;
if ~validRefScen
    matRad_cfg.dispError('ct.refScen must be a positive integer scalar.');
end

if isfield(ct, 'numOfCtScen') && refScen > ct.numOfCtScen
    matRad_cfg.dispError('ct.refScen (%d) exceeds ct.numOfCtScen (%d).', refScen, ct.numOfCtScen);
end

refScen = double(refScen);

end

function cstEval = matRad_buildSamplingReferenceCst(cst, refScen)

matRad_cfg = MatRad_Config.instance();
cstEval = cst;

for i = 1:size(cstEval, 1)
    if numel(cstEval{i, 4}) < refScen
        matRad_cfg.dispError('Structure %s does not contain contours for reference CT scenario %d.', ...
                             cstEval{i, 2}, refScen);
    end
    cstEval{i, 4}{1} = cstEval{i, 4}{refScen};
end

end

function dvhDoseGrid = matRad_resolveSamplingDvhDoseGrid(doseCube, dvhDoseWindow, dvhDoseGrid)

if ~isempty(dvhDoseGrid)
    dvhDoseGrid = dvhDoseGrid(:)';
    if numel(dvhDoseGrid) < 2 || any(~isfinite(dvhDoseGrid)) || any(diff(dvhDoseGrid) <= 0)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('dvhDoseGrid must be a finite increasing vector.');
    end
    return
end

if ~isempty(dvhDoseWindow)
    doseWindow = dvhDoseWindow(:)';
    if numel(doseWindow) >= 2 && all(isfinite(doseWindow(1:2))) && doseWindow(2) > doseWindow(1)
        dvhDoseGrid = linspace(min(0, doseWindow(1)), doseWindow(2), 1000);
        return
    end

    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('dvhDoseWindow must contain finite increasing dose limits.');
end

maxDose = max(doseCube(:));
if ~isfinite(maxDose) || maxDose <= 0
    maxDose = 1;
end
dvhDoseGrid = linspace(0, 1.05 * maxDose, 1000);

end
