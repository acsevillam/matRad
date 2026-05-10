function dijContext = matRad_buildNominalDijContext( ...
    ct,cst,stf,pln,cfg,compactNumBixels,robustDij)
% matRad_buildNominalDijContext builds a compatible nominal dij
%
% The returned dij is used as optimization context for scenario-free robust
% data (INTERVAL/PROB2). The compact robust payload remains attached to
% pln.propOpt; nominal objectives must evaluate a real nominal dose matrix.
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

if nargin < 7
    robustDij = [];
end

matRad_cfg = MatRad_Config.instance();
refScen = resolveReferenceScenario(cfg);
compactNumBixels = validateCompactBixelCount(compactNumBixels,matRad_cfg);

dijContext = extractCompatibleNominalDij(pln,robustDij,refScen, ...
    compactNumBixels,matRad_cfg);
if isempty(dijContext)
    dijContext = calculateCompatibleNominalDij(ct,cst,stf,pln,refScen, ...
        matRad_cfg);
end

dijContext = finalizeNominalDijContext(dijContext,refScen, ...
    compactNumBixels,matRad_cfg);
end

function refScen = resolveReferenceScenario(cfg)
if isstruct(cfg) && isfield(cfg,'refScen') && ~isempty(cfg.refScen)
    refScen = cfg.refScen;
else
    refScen = 1;
end
end

function compactNumBixels = validateCompactBixelCount(compactNumBixels,matRad_cfg)
if ~isnumeric(compactNumBixels) || ~isscalar(compactNumBixels) || ...
        ~isfinite(compactNumBixels) || compactNumBixels < 1 || ...
        fix(compactNumBixels) ~= compactNumBixels
    matRad_cfg.dispError(['Scenario-free nominal dij context requires a ', ...
        'positive integer compact bixel count.']);
end
compactNumBixels = double(compactNumBixels);
end

function dijContext = extractCompatibleNominalDij(pln,robustDij,refScen, ...
    compactNumBixels,matRad_cfg)
dijContext = [];
if isempty(robustDij) || ~isstruct(robustDij) || ...
        ~isfield(pln,'multScen') || isempty(pln.multScen) || ...
        ~isa(pln.multScen,'matRad_ScenarioModel')
    return;
end

nominalScenarioIds = pln.multScen.getNominalScenarioIds();
for scenarioId = nominalScenarioIds(:)'
    if pln.multScen.getCtScenario(scenarioId) ~= refScen
        continue;
    end

    scenarioDijIx = pln.multScen.getDijScenarioIndex(scenarioId);
    if ~isCompatibleRobustDijScenario(robustDij,scenarioDijIx, ...
            compactNumBixels)
        continue;
    end

    dijContext = extractRobustDijScenario(robustDij,scenarioDijIx);
    dijContext.nominalScenarioDijIx = scenarioDijIx;
    dijContext.nominalScenarioId = scenarioId;
    return;
end

if ~isempty(nominalScenarioIds)
    matRad_cfg.dispWarning(['No compatible nominal robust dij scenario was ', ...
        'found for refScen %d. Recalculating nominal dij with the robust ', ...
        'steering geometry.'],refScen);
end
end

function tf = isCompatibleRobustDijScenario(dij,scenarioDijIx,compactNumBixels)
tf = false;
if ~isfield(dij,'physicalDose') || ~iscell(dij.physicalDose) || ...
        numel(dij.physicalDose) < scenarioDijIx || ...
        isempty(dij.physicalDose{scenarioDijIx})
    return;
end

matrix = dij.physicalDose{scenarioDijIx};
if ~isnumeric(matrix) || ~ismatrix(matrix) || ...
        size(matrix,2) ~= compactNumBixels
    return;
end

if isfield(dij,'totalNumOfBixels') && ~isempty(dij.totalNumOfBixels) && ...
        dij.totalNumOfBixels ~= compactNumBixels
    return;
end

tf = true;
end

function dijOut = extractRobustDijScenario(dijIn,scenarioDijIx)
dijOut = dijIn;
fieldNames = fieldnames(dijOut);
for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    value = dijOut.(fieldName);
    if ~iscell(value) || numel(value) < scenarioDijIx
        continue;
    end

    scenarioValue = value{scenarioDijIx};
    if ~(isempty(scenarioValue) || isnumeric(scenarioValue) || ...
            islogical(scenarioValue))
        continue;
    end

    scenarioCell = cell(1,1,1);
    scenarioCell{1} = scenarioValue;
    dijOut.(fieldName) = scenarioCell;
end
end

function dijContext = calculateCompatibleNominalDij(ct,cst,stf,pln, ...
    refScen,matRad_cfg)
if isempty(stf)
    matRad_cfg.dispError(['Cannot build a compatible nominal dij context ', ...
        'without either an in-memory robust dij nominal scenario or the ', ...
        'robust steering information stf.']);
end

plnNominal = pln;
plnNominal.multScen = buildNominalScenarioModel(refScen);
if isfield(plnNominal,'propOpt') && isstruct(plnNominal.propOpt) && ...
        isfield(plnNominal.propOpt,'scen4D')
    plnNominal.propOpt = rmfield(plnNominal.propOpt,'scen4D');
end

dijContext = matRad_calcDoseInfluence(ct,cst,stf,plnNominal);
end

function dijContext = finalizeNominalDijContext(dijContext,refScen, ...
    compactNumBixels,matRad_cfg)
if ~isstruct(dijContext) || ~isfield(dijContext,'physicalDose') || ...
        ~iscell(dijContext.physicalDose) || isempty(dijContext.physicalDose) || ...
        isempty(dijContext.physicalDose{1})
    matRad_cfg.dispError(['Scenario-free nominal dij context must contain ', ...
        'physicalDose{1}.']);
end

if size(dijContext.physicalDose{1},2) ~= compactNumBixels
    matRad_cfg.dispError(['Scenario-free nominal dij context bixel count ', ...
        '(%d) is inconsistent with compact robust data (%d).'], ...
        size(dijContext.physicalDose{1},2),compactNumBixels);
end

if isfield(dijContext,'totalNumOfBixels') && ...
        ~isempty(dijContext.totalNumOfBixels) && ...
        dijContext.totalNumOfBixels ~= compactNumBixels
    matRad_cfg.dispError(['Scenario-free nominal dij metadata ', ...
        'totalNumOfBixels (%d) is inconsistent with compact robust data ', ...
        '(%d).'],dijContext.totalNumOfBixels,compactNumBixels);
end

dijContext.totalNumOfBixels = compactNumBixels;
if ~isfield(dijContext,'beamNum') || ...
        numel(dijContext.beamNum) ~= compactNumBixels
    dijContext.beamNum = ones(compactNumBixels,1);
else
    dijContext.beamNum = dijContext.beamNum(:);
end
if ~isfield(dijContext,'numOfBeams') || isempty(dijContext.numOfBeams)
    dijContext.numOfBeams = max(dijContext.beamNum);
end

dijContext.numOfScenarios = 1;
dijContext.scenarioModel = buildNominalScenarioModel(refScen);
end

function scenarioModel = buildNominalScenarioModel(refScen)
components = matRad_createScenarioComponents([0 0 0],0,0,{'ct'});
scenarioValues = zeros(1,numel(components));
ctScenIds = refScen;
scenProb = 1;
scenWeight = 1;
scenForProb = [ctScenIds scenarioValues];
linearMask = [1 1 1];
scenMask = true(1,1,1);

scenarioModel = matRad_NominalScenario();
scenarioModel.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
    scenProb,scenWeight,scenForProb,linearMask,scenMask);
end
