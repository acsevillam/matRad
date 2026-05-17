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

%% perform sampling
[mSampDose, caSampRes] = matRad_executeSamplingScenarios(samplingContext, scenarioIds, mSampDose, ...
                                                         nomScenTime, matRad_cfg);

%% add subindices
pln.subIx        = subIx;
pln.samplingReferenceCtScen = refScen;
pln.samplingDoseMapping = doseMapping;

end

function [mSampDose, caSampRes] = matRad_executeSamplingScenarios(samplingContext, scenarioIds, mSampDose, ...
                                                                  nomScenTime, matRadCfg)

[parallelToolboxLicensed, pool] = matRad_prepareSamplingPool(matRadCfg);
numScenarios = numel(scenarioIds);

if parallelToolboxLicensed
    [mSampDose, caSampRes] = matRad_executeParallelSamplingScenarios(samplingContext, scenarioIds, ...
                                                                     mSampDose, nomScenTime, pool, matRadCfg);
else
    [mSampDose, caSampRes] = matRad_executeSerialSamplingScenarios(samplingContext, scenarioIds, ...
                                                                   mSampDose, nomScenTime, matRadCfg);
end

if numel(caSampRes) ~= numScenarios
    matRadCfg.dispError('Sampling did not produce the expected number of scenario results.');
end

end

function [parallelToolboxLicensed, pool] = matRad_prepareSamplingPool(matRadCfg)
pool = [];
try
    [parallelToolboxLicensed, ~] = license('checkout', 'Distrib_Computing_Toolbox');
    if ~parallelToolboxLicensed
        matRadCfg.dispWarning('Could not check out parallel computing toolbox. \n');
        return
    end

    pool = gcp(); % If no pool exists, create one.
    if isempty(pool)
        matRadCfg.dispError('matRad: Could not start valid parallel pool. Please check your parallel computing toolbox installation. \n');
    end
catch
    parallelToolboxLicensed = false;
    pool = [];
end
end

function [mSampDose, caSampRes] = matRad_executeParallelSamplingScenarios(samplingContext, scenarioIds, ...
                                                                          mSampDose, nomScenTime, pool, matRadCfg)

if isempty(pool)
    poolSize = 1;
else
    poolSize = pool.NumWorkers;
end

logLevel = matRadCfg.logLevel;
numScenarios = numel(scenarioIds);
totCompTime = ceil(numScenarios / poolSize) * nomScenTime * 1.35;
matRad_logSamplingTimeEstimate(totCompTime, matRadCfg);

progressEnabled = matRad_startSamplingProgress(numScenarios, logLevel, matRadCfg);
sampleResults = cell(1, numScenarios);

parfor i = 1:numScenarios
    matRadCfgWorker = MatRad_Config.instance();
    matRadCfgWorker.logLevel = logLevel;

    [mSampDose(:, i), sampleResults{i}] = matRad_calculateSamplingScenario(samplingContext, scenarioIds(i));

    if progressEnabled && logLevel > 2
        parfor_progress;
    end
end

if progressEnabled && logLevel > 2
    parfor_progress(0);
end
caSampRes = [sampleResults{:}];

end

function [mSampDose, caSampRes] = matRad_executeSerialSamplingScenarios(samplingContext, scenarioIds, ...
                                                                        mSampDose, nomScenTime, matRadCfg)

numScenarios = numel(scenarioIds);
totCompTime = numScenarios * nomScenTime * 1.1;
matRad_logSamplingTimeEstimate(totCompTime, matRadCfg);
sampleResults = cell(1, numScenarios);

for i = 1:numScenarios
    [mSampDose(:, i), sampleResults{i}] = matRad_calculateSamplingScenario(samplingContext, scenarioIds(i));

    if matRadCfg.logLevel > 2
        matRad_progress(i, numScenarios);
    end
end
caSampRes = [sampleResults{:}];

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
