function ctx = matRad_resolveDoseIntervalInputs(ct,cst,pln,dij,cfg,intervalMode,matRad_cfg)
% matRad_resolveDoseIntervalInputs validates inputs for interval dose calculation
%
% call
%   ctx = matRad_resolveDoseIntervalInputs(ct,cst,pln,dij,cfg,intervalMode,matRad_cfg)
%
% input
%   ct:           matRad ct struct
%   cst:          matRad cst cell array
%   pln:          matRad pln struct with pln.multScen as matRad_ScenarioModel
%   dij:          robust dose influence struct containing the requested
%                 linear quantity as scenario cell matrices
%   cfg:          interval configuration struct; dose quantities are in Gy
%                 or Gy(RBE) according to the selected linear dij field
%   intervalMode: interval method identifier, either 'INTERVAL2' or 'INTERVAL3'
%   matRad_cfg:   MatRad_Config instance for diagnostics
%
% output
%   ctx:          validated context struct with quantity metadata, DIJ
%                 scenario indices, CT scenario ids, normalized scenario
%                 weights, selected target/OAR rows and CT mapping metadata
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

cfg = normalizeConfig(cfg,ct,pln,intervalMode,matRad_cfg);
quantity = resolveQuantity(dij,pln,cfg,matRad_cfg);
[scenarioDijIx,scenarioCtScenIds,scenarioWeights] = resolveScenarios(pln,dij,quantity.field,matRad_cfg);

cstDoseGrid = resizeCstToDoseGrid(cst,dij);
targetRows = resolveStructureRows(cstDoseGrid,cfg.targetStructSel,'TARGET',cfg.refScen);
oarRows = resolveStructureRows(cstDoseGrid,cfg.OARStructSel,'OAR',cfg.refScen);

ctx.cfg = cfg;
ctx.quantity = quantity;
ctx.scenarioDijIx = scenarioDijIx;
ctx.scenarioCtScenIds = scenarioCtScenIds;
ctx.scenarioWeights = scenarioWeights;
ctx.targetRows = targetRows;
ctx.oarRows = oarRows;
ctx.numVoxels = getNumDoseVoxels(dij,quantity.matrixCell,matRad_cfg);
ctx.numBixels = getNumBixels(dij,quantity.matrixCell,scenarioDijIx,matRad_cfg);
ctx.scenarioMaps = matRad_buildDoseIntervalScenarioMappings(ct,dij, ...
    scenarioCtScenIds,cfg.refScen,matRad_cfg);
end

function cfg = normalizeConfig(cfg,ct,pln,intervalMode,matRad_cfg)
if nargin < 1 || isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRad_cfg.dispError('Dose interval configuration must be a struct.');
end

cfg = setDefault(cfg,'Quantity',getDefaultQuantity(pln));
cfg = setDefault(cfg,'QuantityField',[]);
cfg = setDefault(cfg,'refScen',getDefaultReferenceScenario(ct));
cfg = setDefault(cfg,'targetStructSel',[]);
cfg = setDefault(cfg,'OARStructSel',[]);
cfg = setDefault(cfg,'UseParallel',false);
cfg = setDefault(cfg,'MemoryLimitMB',[]);
cfg = setDefault(cfg,'BatchSize',[]);
cfg = setDefault(cfg,'ProgressLevel','summary');
cfg = setDefault(cfg,'KMode','dynamic');
cfg = setDefault(cfg,'KMax',10);
cfg = setDefault(cfg,'RetentionThreshold',1.0);

cfg.refScen = validatePositiveInteger(cfg.refScen,'refScen',matRad_cfg);
if isfield(ct,'numOfCtScen') && cfg.refScen > ct.numOfCtScen
    matRad_cfg.dispError('refScen (%d) exceeds ct.numOfCtScen (%d).', ...
        cfg.refScen,ct.numOfCtScen);
end

if isfield(ct,'refScen') && ~isempty(ct.refScen) && ct.refScen ~= cfg.refScen
    matRad_cfg.dispError('Requested refScen %d does not match ct.refScen (%d).', ...
        cfg.refScen,ct.refScen);
end

cfg.UseParallel = logicalScalar(cfg.UseParallel,'UseParallel',matRad_cfg);

if ~isempty(cfg.MemoryLimitMB) && ...
   (~isnumeric(cfg.MemoryLimitMB) || ~isscalar(cfg.MemoryLimitMB) || ...
    ~isfinite(cfg.MemoryLimitMB) || cfg.MemoryLimitMB <= 0)
    matRad_cfg.dispError('MemoryLimitMB must be a positive finite scalar.');
end

if ~isempty(cfg.BatchSize)
    cfg.BatchSize = validatePositiveInteger(cfg.BatchSize,'BatchSize',matRad_cfg);
end

cfg.ProgressLevel = validateProgressLevel(cfg.ProgressLevel,matRad_cfg);

if strcmp(intervalMode,'INTERVAL3')
    cfg.KMax = validatePositiveInteger(cfg.KMax,'KMax',matRad_cfg);
    if isstring(cfg.KMode) && isscalar(cfg.KMode)
        cfg.KMode = char(cfg.KMode);
    end
    if ~ischar(cfg.KMode) || ~any(strcmpi(cfg.KMode,{'dynamic','static'}))
        matRad_cfg.dispError('KMode must be ''dynamic'' or ''static''.');
    end
    cfg.KMode = lower(cfg.KMode);
    if ~isnumeric(cfg.RetentionThreshold) || ~isscalar(cfg.RetentionThreshold) || ...
       ~isfinite(cfg.RetentionThreshold) || cfg.RetentionThreshold <= 0 || ...
       cfg.RetentionThreshold > 1
        matRad_cfg.dispError('RetentionThreshold must be in the interval (0,1].');
    end
end
end

function progressLevel = validateProgressLevel(progressLevel,matRad_cfg)
if isstring(progressLevel) && isscalar(progressLevel)
    progressLevel = char(progressLevel);
end

if ~ischar(progressLevel) || ...
   ~any(strcmpi(progressLevel,{'summary','detailed'}))
    matRad_cfg.dispError('ProgressLevel must be ''summary'' or ''detailed''.');
end

progressLevel = lower(progressLevel);
end

function value = getDefaultQuantity(pln)
value = 'physicalDose';
if isfield(pln,'bioParam')
    if isobject(pln.bioParam) && isprop(pln.bioParam,'quantityOpt') && ...
       ~isempty(pln.bioParam.quantityOpt)
        value = pln.bioParam.quantityOpt;
    elseif isstruct(pln.bioParam) && isfield(pln.bioParam,'quantityOpt') && ...
           ~isempty(pln.bioParam.quantityOpt)
        value = pln.bioParam.quantityOpt;
    end
end
end

function refScen = getDefaultReferenceScenario(ct)
if isfield(ct,'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function cfg = setDefault(cfg,fieldName,value)
if ~isfield(cfg,fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = value;
end
end

function value = validatePositiveInteger(value,name,matRad_cfg)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
   value < 1 || fix(value) ~= value
    matRad_cfg.dispError('%s must be a positive integer scalar.',name);
end
value = double(value);
end

function value = logicalScalar(value,name,matRad_cfg)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
    matRad_cfg.dispError('%s must be a logical scalar.',name);
end
value = logical(value);
end

function quantity = resolveQuantity(dij,pln,cfg,matRad_cfg)
quantityName = normalizeText(cfg.Quantity,'Quantity',matRad_cfg);
quantityField = normalizeOptionalText(cfg.QuantityField,'QuantityField',matRad_cfg);
scale = 1;

if ~isempty(quantityField)
    field = quantityField;
elseif hasLinearDijField(dij,quantityName)
    field = quantityName;
elseif strcmpi(quantityName,'physicalDose')
    field = 'physicalDose';
elseif strcmpi(quantityName,'RBExD')
    if isfield(dij,'RBE') && isnumeric(dij.RBE) && isscalar(dij.RBE) && isfinite(dij.RBE)
        field = 'physicalDose';
        scale = dij.RBE;
    else
        matRad_cfg.dispError(['RBExD interval calculation requires a scalar dij.RBE ' ...
            'or an explicit linear QuantityField.']);
    end
elseif any(strcmpi(quantityName,{'effect','BED'}))
    matRad_cfg.dispError('%s interval calculation requires an explicit linear QuantityField.', ...
        quantityName);
else
    matRad_cfg.dispError('Could not resolve linear dose interval quantity ''%s''.',quantityName);
end

if ~hasLinearDijField(dij,field)
    matRad_cfg.dispError('dij.%s must be a cell array of dose influence matrices.',field);
end

quantity.name = quantityName;
quantity.field = field;
quantity.scale = scale;
quantity.matrixCell = dij.(field);
quantity.optimizationQuantity = getDefaultQuantity(pln);
end

function textValue = normalizeText(value,name,matRad_cfg)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || isempty(value)
    matRad_cfg.dispError('%s must be a non-empty character vector.',name);
end
textValue = value;
end

function textValue = normalizeOptionalText(value,name,matRad_cfg)
if isempty(value)
    textValue = '';
    return;
end
textValue = normalizeText(value,name,matRad_cfg);
end

function tf = hasLinearDijField(dij,fieldName)
tf = isfield(dij,fieldName) && iscell(dij.(fieldName));
end

function [scenarioDijIx,scenarioCtScenIds,scenarioWeights] = resolveScenarios(pln,dij,quantityField,matRad_cfg)
if ~isfield(pln,'multScen') || isempty(pln.multScen)
    matRad_cfg.dispError('Dose interval calculation requires pln.multScen.');
end

if ~isa(pln.multScen,'matRad_ScenarioModel')
    matRad_cfg.dispError('Dose interval calculation requires a matRad_ScenarioModel instance.');
end

scenarioIds = pln.multScen.scenarioIds();
scenarioDijIx = arrayfun(@(id) pln.multScen.getDijScenarioIndex(id),scenarioIds);

if isempty(scenarioDijIx)
    matRad_cfg.dispError('No active scenarios found for dose interval calculation.');
end

quantityCells = dij.(quantityField);
if max(scenarioDijIx) > numel(quantityCells)
    matRad_cfg.dispError('Scenario indices exceed dij.%s cell array dimensions.',quantityField);
end

emptyScenario = cellfun(@isempty,quantityCells(scenarioDijIx));
if any(emptyScenario)
    matRad_cfg.dispError('dij.%s contains empty active scenarios.',quantityField);
end

scenarioCtScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id),scenarioIds);

scenarioWeights = resolveScenarioWeights(pln.multScen,scenarioDijIx,matRad_cfg);
scenarioDijIx = scenarioDijIx(:);
scenarioCtScenIds = scenarioCtScenIds(:);
end

function scenarioWeights = resolveScenarioWeights(multScen,scenarioDijIx,matRad_cfg)
if hasFieldOrProp(multScen,'scenWeight') && ~isempty(multScen.scenWeight)
    rawWeights = multScen.scenWeight(:);
elseif hasFieldOrProp(multScen,'scenProb') && ~isempty(multScen.scenProb)
    rawWeights = multScen.scenProb(:);
else
    matRad_cfg.dispError('Scenario model must provide scenWeight or scenProb.');
end

if numel(rawWeights) ~= numel(scenarioDijIx)
    matRad_cfg.dispError('Number of scenario weights is inconsistent with active scenarios.');
end
scenarioWeights = rawWeights;

scenarioWeights = scenarioWeights(:);
scenarioWeights(~isfinite(scenarioWeights) | scenarioWeights < 0) = 0;
if sum(scenarioWeights) <= 0
    matRad_cfg.dispError('Scenario weights must contain at least one positive finite value.');
end
scenarioWeights = scenarioWeights./sum(scenarioWeights);
end

function tf = hasFieldOrProp(value,fieldName)
tf = (isobject(value) && isprop(value,fieldName)) || ...
    (isstruct(value) && isfield(value,fieldName));
end

function cst = resizeCstToDoseGrid(cst,dij)
requiredGridFields = {'x','y','z'};
if isfield(dij,'ctGrid') && isfield(dij,'doseGrid') && ...
   all(isfield(dij.ctGrid,requiredGridFields)) && ...
   all(isfield(dij.doseGrid,requiredGridFields))
    cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z, ...
        dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
end
end

function rows = resolveStructureRows(cst,structureSelection,structureType,refScen)
if isempty(structureSelection)
    cstRows = matRad_resolveStructureSelection(cst,'all',[],structureType);
else
    cstRows = matRad_resolveStructureSelection(cst,'include',structureSelection,structureType);
end

rows = [];
for i = cstRows
    rows = [rows; getStructureVoxelIndices(cst,i,refScen)]; %#ok<AGROW>
end
rows = unique(rows(:),'stable');
end

function voxels = getStructureVoxelIndices(cst,rowIx,refScen)
voxels = [];
if size(cst,2) < 4 || isempty(cst{rowIx,4})
    return;
end

contours = cst{rowIx,4};
if iscell(contours) && numel(contours) >= refScen && ~isempty(contours{refScen})
    voxels = contours{refScen}(:);
elseif iscell(contours) && ~isempty(contours{1})
    voxels = contours{1}(:);
elseif isnumeric(contours)
    voxels = contours(:);
end
end

function numVoxels = getNumDoseVoxels(dij,matrixCell,matRad_cfg)
if isfield(dij,'doseGrid') && isfield(dij.doseGrid,'numOfVoxels')
    numVoxels = dij.doseGrid.numOfVoxels;
elseif isfield(dij,'doseGrid') && isfield(dij.doseGrid,'dimensions')
    numVoxels = prod(dij.doseGrid.dimensions);
else
    firstNonEmpty = find(~cellfun(@isempty,matrixCell(:)),1,'first');
    if isempty(firstNonEmpty)
        matRad_cfg.dispError('Could not determine number of dose voxels.');
    end
    numVoxels = size(matrixCell{firstNonEmpty},1);
end
end

function numBixels = getNumBixels(dij,matrixCell,scenarioDijIx,matRad_cfg)
if isfield(dij,'totalNumOfBixels') && ~isempty(dij.totalNumOfBixels)
    numBixels = dij.totalNumOfBixels;
else
    numBixels = size(matrixCell{scenarioDijIx(1)},2);
end

if any(cellfun(@(m) size(m,2) ~= numBixels,matrixCell(scenarioDijIx)))
    matRad_cfg.dispError('All active scenario matrices must have the same number of bixels.');
end
end
