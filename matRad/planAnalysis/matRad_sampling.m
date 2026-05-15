function [caSampRes,mSampDose,pln,resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,w,structSel,multScen,varargin)
% matRad_sampling enables sampling multiple treatment scenarios
%
% call
%   [caSampRes,mSampDose,pln,resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,w,structSel)
%   [caSampRes,mSampDose,pln,resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,w,structSel,multScen)
%   [caSampRes,mSampDose,pln,resultGUInomScen] = matRad_sampling(...,'dvhDoseWindow',dvhDoseWindow)
%   [caSampRes,mSampDose,pln,resultGUInomScen] = matRad_sampling(...,'dvhDoseGrid',dvhDoseGrid)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
%   structSel:  (optional) cell array of structure names used to define the
%               sampled voxel subset. If empty, all structures are used.
%   multScen:   (optional) matRad scenario model used for sampling. If
%               empty, a random scenario model is created.
%
% input (optional Name-Value pairs)
%   varargin:       optional Name-Value pairs
%   dvhDoseWindow:  1x2 dose window used to create the common DVH dose grid
%   dvhDoseGrid:    explicit finite increasing dose grid used for all DVHs
%   autoLimitWorkers: automatically reduce the parallel pool size based on
%               available system memory and estimated memory per worker
%   workerMemorySafetyFactor: safety factor applied to the worker memory
%               estimate
%   memoryReserveFraction: fraction of total system memory kept in reserve
%               when autoLimitWorkers is enabled
%   minWorkerMemoryBytes: lower bound for the estimated per-worker memory
%               footprint before applying workerMemorySafetyFactor
%   workerUpperBound: explicit upper bound for the memory-limited worker
%               count. Empty means no explicit bound.
%
% note
%   Multi-CT sampled doses are always mapped to the reference CT scenario
%   before computing sampling statistics.
%
% output
%   caSampRes:         cell array of sampling results depicting plan parameter
%   mSampDose:         matrix holding the sampled doses, each row corresponds
%                      to one sampled voxel in pln.subIx and each column to
%                      one sampled scenario
%   pln:               matRad pln struct containing sampling information,
%                      including samplingMemoryEstimate diagnostics
%   resultGUInomScen:  resultGUI struct of the nominal scenario
%
% References
%   -
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
elseif ischar(multScen) || (isa(multScen,'string') && isscalar(multScen))
    varargin = [{char(multScen)} varargin];
    multScen = [];
end

% Dose calculation and parfor serialization allocate memory not captured by whos.
defaultMinWorkerMemoryBytes = 1024^3;

p = inputParser;
p.addParameter('dvhDoseWindow',[],@(x) isempty(x) || (isnumeric(x) && numel(x) >= 2));
p.addParameter('dvhDoseGrid',[],@(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('autoLimitWorkers',true,@(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.addParameter('workerMemorySafetyFactor',1.2,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addParameter('memoryReserveFraction',0.10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
p.addParameter('minWorkerMemoryBytes',defaultMinWorkerMemoryBytes, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('workerUpperBound',[],@(x) true);
p.parse(varargin{:});
workerUpperBound = validateOptionalWorkerUpperBound( ...
    p.Results.workerUpperBound,matRad_cfg);

if isfield(ct,'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end
if ~(isnumeric(refScen) && isscalar(refScen) && isfinite(refScen) && ...
        round(refScen) == refScen && refScen >= 1)
    matRad_cfg.dispError('ct.refScen must be a positive integer scalar.');
end
if isfield(ct,'numOfCtScen') && refScen > ct.numOfCtScen
    matRad_cfg.dispError('ct.refScen (%d) exceeds ct.numOfCtScen (%d).', ...
        refScen,ct.numOfCtScen);
end
refScen = double(refScen);

cstEval = cst;
for i = 1:size(cstEval,1)
    if numel(cstEval{i,4}) < refScen
        matRad_cfg.dispError('Structure %s does not contain contours for reference CT scenario %d.', ...
            cstEval{i,2},refScen);
    end
    cstEval{i,4}{1} = cstEval{i,4}{refScen};
end

% save nonSampling pln for nominal scenario calculation and add dummy fields
plnNominal = pln;
% create nominal scenario
plnNominal.multScen = matRad_NominalScenario(ct);
plnNominal.multScen.ctScenProb = [refScen 1];

% check for different ct scenarios
if ct.numOfCtScen > 1 && isempty(multScen)
    matRad_cfg.dispWarning(['No explicit multi-CT sampling scenario model was provided; ', ...
        'random sampling will use the reference CT scenario only.\n']);
end

% either use existing multScen struct or create new one
if ~isempty(multScen)
    pln.multScen = multScen;
else
    % create random scenarios for sampling
    pln.multScen = matRad_RandomScenarios(ct);
    pln.multScen.ctScenProb = [refScen 1];
    pln.multScen.nSamples = matRad_cfg.defaults.samplingScenarios;
end

numSamples = pln.multScen.totNumScen;
matRad_cfg.dispInfo('Using %d samples in total \n',numSamples);

doseMapping = matRad_resolveSamplingDoseMapping(ct,pln.multScen,refScen);
if doseMapping.enabled
    matRad_cfg.dispInfo('matRad: Mapping sampled multi-CT dose cubes to reference CT scenario %d before analysis.\n',refScen);
end

scenarioIds = pln.multScen.scenarioIds();
samplingCtScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id),scenarioIds);
if ~exist('structSel','var')
    structSel = {};
end

voxelSets = cell(1,size(cstEval,1));
numVoxelSets = 0;
for i = 1:size(cstEval,1)
    if isempty(structSel) || any(strcmp(structSel,cstEval{i,2}))
        if ~isempty(cstEval{i,4}{1})
            numVoxelSets = numVoxelSets + 1;
            voxelSets{numVoxelSets} = cstEval{i,4}{1};
        end
    end
end
if numVoxelSets == 0
    matRad_cfg.dispError('No voxels selected for sampling.');
end
subIx = unique(vertcat(voxelSets{1:numVoxelSets}));
for i = 1:size(cstEval,1)
    if ~all(ismember(cstEval{i,4}{1},subIx))
        cstEval{i,5}.Visible = false;
    end
end

mSampDose = single(zeros(numel(subIx),numSamples,1));
caSampRes(1,numSamples) = createEmptySamplingResult();
storageInfo = whos('mSampDose');
samplingDoseStorageBytes = storageInfo.bytes;

% check if parallel toolbox is installed and license can be checked out
try
    [parallelToolboxLicensed,~] = license('checkout','Distrib_Computing_Toolbox');
    if ~parallelToolboxLicensed
        matRad_cfg.dispWarning('Could not check out parallel computing toolbox. \n');
    end

catch
    parallelToolboxLicensed = false;
end

%% calculate nominal scenario
nomScenTimer     = tic;
resultGUInomScen = matRad_calcDoseForward(ct,cst,stf,plnNominal,w);
nomScenTime      = toc(nomScenTimer);
matRad_cfg.dispInfo('Finished nominal Scenario Calculation. Computation time: %f h \n',round(nomScenTime / 3600));

refVol = [2 5 50 95 98];
doseCubeNominal = resultGUInomScen.(pln.bioParam.quantityVis);
refGy = linspace(0,max(doseCubeNominal(:)),6);
dvhPoints = resolveDvhDoseGrid(doseCubeNominal,p.Results.dvhDoseWindow, ...
    p.Results.dvhDoseGrid);

resultGUInomScen.dvh = matRad_calcDVH(cstEval,doseCubeNominal,'cum',dvhPoints);
nomQi                = matRad_calcQualityIndicators(cstEval,pln,doseCubeNominal,refGy,refVol);

resultGUInomScen.qi  = nomQi;
resultGUInomScen.cst = cstEval;
resultGUInomScen.evaluationModeBase = 'perFraction';

samplingMemoryContext = struct();
samplingMemoryContext.ct = ct;
samplingMemoryContext.stf = stf;
samplingMemoryContext.cst = cst;
samplingMemoryContext.cstEval = cstEval;
samplingMemoryContext.pln = pln;
samplingMemoryContext.w = w;
samplingMemoryContext.subIx = subIx;
samplingMemoryContext.samplingCtScenIds = samplingCtScenIds;
samplingMemoryContext.dvhPoints = dvhPoints;
samplingMemoryContext.refGy = refGy;
samplingMemoryContext.refVol = refVol;
samplingMemoryContext.resultGUInomScen = resultGUInomScen;
samplingMemoryContext.doseMapping = doseMapping;
samplingMemoryContext.refScen = refScen;
samplingMemoryContext.samplingDoseStorageBytes = samplingDoseStorageBytes;
samplingMemoryContext.numSamples = numSamples;
samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingMemoryContext);

%% perform parallel sampling
if parallelToolboxLicensed
    % MatRad_Config is process-local; mirror relevant settings on workers.
    logLevel = matRad_cfg.logLevel;

    poolSizeLimit = [];
    memoryEstimate = [];
    if p.Results.autoLimitWorkers
        [poolSizeLimit,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(samplingMemoryEstimate.rawWorkerBytes, ...
            'numTasks',numSamples, ...
            'safetyFactor',p.Results.workerMemorySafetyFactor, ...
            'reserveFraction',p.Results.memoryReserveFraction, ...
            'minWorkerMemoryBytes',p.Results.minWorkerMemoryBytes, ...
            'workerUpperBound',workerUpperBound);
        samplingMemoryEstimate.workerLimit = memoryEstimate;
    end
    matRad_logEstimatedSamplingMemory(samplingMemoryEstimate,memoryEstimate,matRad_cfg);

    % Create or resize parallel pool on cluster
    pPool = configureSamplingPool(poolSizeLimit,matRad_cfg);

    if isempty(pPool)
        poolSize = 1;
    else
        poolSize = pPool.NumWorkers;
    end

    % rough estimate of total computation time
    totCompTime = ceil(numSamples / poolSize) * nomScenTime * 1.35;
    logEstimatedSamplingTime(totCompTime,matRad_cfg);

    progressQueue = [];
    progressListener = [];
    nFinishedScenarios = 0;
    if exist('parallel.pool.DataQueue','class') == 8 && logLevel > 2
        parforProgressEnabled = true;
        progressQueue = parallel.pool.DataQueue;
        matRad_cfg.dispInfo('Sampling progress: 0 scenarios of %d (0%%).\n',numSamples);
        progressListener = afterEach(progressQueue,@(~) updateSamplingProgress());
    elseif exist('parfor_progress', 'file') == 2 && logLevel > 2
        parforProgressEnabled = true;
        parfor_progress(numSamples);
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        parforProgressEnabled = false;
    end

    parfor i = 1:numSamples
        matRad_cfg_worker = MatRad_Config.instance();
        matRad_cfg_worker.logLevel = logLevel;

        [mSampDose(:,i),caSampRes(i)] = matRad_calculateSamplingScenario( ...
            ct,stf,cst,pln,w,cstEval,subIx,dvhPoints,refGy,refVol, ...
            samplingCtScenIds,doseMapping,refScen,i);

        if parforProgressEnabled && logLevel > 2
            if isempty(progressQueue)
                parfor_progress;
            else
                send(progressQueue,i);
            end
        end
    end

    if parforProgressEnabled && logLevel > 2 && isempty(progressQueue)
        parfor_progress(0);
    end
    if ~isempty(progressListener)
        delete(progressListener);
    end

else
    %% perform serial sampling
    matRad_logEstimatedSamplingMemory(samplingMemoryEstimate,[],matRad_cfg);

    % rough estimate of total computation time
    totCompTime = numSamples * nomScenTime * 1.1;
    logEstimatedSamplingTime(totCompTime,matRad_cfg);

    for i = 1:numSamples
        [mSampDose(:,i),caSampRes(i)] = matRad_calculateSamplingScenario( ...
            ct,stf,cst,pln,w,cstEval,subIx,dvhPoints,refGy,refVol, ...
            samplingCtScenIds,doseMapping,refScen,i);

        % Show progress
        if matRad_cfg.logLevel > 2
            matRad_progress(i,numSamples);
        end
    end
end

%% add subindices
pln.subIx = subIx;
pln.samplingMemoryEstimate = samplingMemoryEstimate;

    function updateSamplingProgress()
        nFinishedScenarios = nFinishedScenarios + 1;
        matRad_cfg.dispInfo('Sampling progress: %d scenarios of %d (%.0f%%).\n', ...
            nFinishedScenarios,numSamples, ...
            100 * nFinishedScenarios / numSamples);
        drawnow('limitrate');
    end

end

function sampleResult = createEmptySamplingResult()
sampleResult = struct( ...
    'bioParam',[], ...
    'ctScenId',[], ...
    'refScen',[], ...
    'doseMapping',[], ...
    'relRangeShift',[], ...
    'absRangeShift',[], ...
    'isoShift',[], ...
    'evaluationModeBase','perFraction', ...
    'dvh',[], ...
    'qi',[]);
end

function logEstimatedSamplingTime(totCompTime,matRad_cfg)
estimatedFinishText = formatEstimatedFinish(totCompTime);
if isempty(estimatedFinishText)
    matRad_cfg.dispInfo(['Approximate Total calculation time: ', ...
        num2str(round(totCompTime / 3600)), 'h.\n']);
else
    matRad_cfg.dispInfo(['Approximate Total calculation time: ', ...
        num2str(round(totCompTime / 3600)), ...
        'h. Estimated finish: ', estimatedFinishText, '\n']);
end
end

function estimatedFinishText = formatEstimatedFinish(totCompTime)
try
    estimatedFinishText = char(datetime('now') + seconds(totCompTime));
catch
    estimatedFinishText = '';
end
end

function dvhDoseGrid = resolveDvhDoseGrid(doseCube,dvhDoseWindow,dvhDoseGrid)
if ~isempty(dvhDoseGrid)
    dvhDoseGrid = dvhDoseGrid(:)';
    if numel(dvhDoseGrid) < 2 || any(~isfinite(dvhDoseGrid)) || ...
            any(diff(dvhDoseGrid) <= 0)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('dvhDoseGrid must be a finite increasing vector.');
    end
    return;
end

if ~isempty(dvhDoseWindow)
    doseWindow = dvhDoseWindow(:)';
    if numel(doseWindow) >= 2 && all(isfinite(doseWindow(1:2))) && ...
            doseWindow(2) > doseWindow(1)
        dvhDoseGrid = linspace(min(0,doseWindow(1)),doseWindow(2),1000);
        return;
    end
end

maxDose = max(doseCube(:));
if ~isfinite(maxDose) || maxDose <= 0
    maxDose = 1;
end
dvhDoseGrid = linspace(0,maxDose*1.05,1000);
end

function workerUpperBound = validateOptionalWorkerUpperBound(workerUpperBound,matRad_cfg)
if isempty(workerUpperBound)
    return;
end

if ~(isnumeric(workerUpperBound) && isscalar(workerUpperBound) && ...
        isfinite(workerUpperBound) && round(workerUpperBound) == workerUpperBound && ...
        workerUpperBound >= 1)
    matRad_cfg.dispError('workerUpperBound must be a positive integer scalar or empty.');
end
end

function p = configureSamplingPool(poolSizeLimit,matRad_cfg)
p = gcp('nocreate');

if isempty(p)
    if isempty(poolSizeLimit)
        p = gcp();
    else
        p = parpool(poolSizeLimit);
    end
elseif ~isempty(poolSizeLimit) && p.NumWorkers > poolSizeLimit
    matRad_cfg.dispWarning(['Reducing parallel pool from ', num2str(p.NumWorkers), ...
        ' to ', num2str(poolSizeLimit), ...
        ' worker(s) to keep sampling within the estimated available memory.\n']);
    delete(p);
    p = parpool(poolSizeLimit);
end
end
