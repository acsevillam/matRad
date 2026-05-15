function ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dij,cfg,calculationMode,matRad_cfg)
% matRad_resolveScenarioDoseInputs validates shared scenario-dose inputs
%
% call
%   ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dij,cfg,calculationMode,matRad_cfg)
%
% input
%   ct:           matRad ct struct
%   cst:          matRad cst cell array
%   pln:          matRad pln struct with pln.multScen as matRad_ScenarioModel
%   dij:          robust dose influence struct containing the requested
%                 linear quantity as scenario cell matrices
%   cfg:          robust scenario-dose configuration struct; dose quantities
%                 are in Gy or Gy(RBE) according to the selected linear dij
%                 field
%   calculationMode: method identifier, e.g. 'PROB2', 'INTERVAL2', or
%                    'INTERVAL3'
%   matRad_cfg:   MatRad_Config instance for diagnostics
%
% output
%   ctx:          validated context struct with quantity metadata, DIJ
%                 scenario indices, CT scenario ids selected by
%                 pln.propOpt.scen4D (default: 1), normalized scenario
%                 weights, selected target/OAR rows, selected structure
%                 voxels in the reference scenario, and CT mapping metadata
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

cfg = normalizeConfig(cfg,ct,pln,calculationMode,matRad_cfg);
quantity = resolveQuantity(dij,pln,cfg,matRad_cfg);
[scenarioDijIx,scenarioCtScenIds,scenarioWeights] = resolveScenarios(pln,dij,quantity.field,matRad_cfg);

cstDoseGrid = prepareCstForScenarioDoseRows(cst,dij);
[targetRows,targetStructures] = resolveStructureRows(cstDoseGrid, ...
    cfg.targetStructSel,'TARGET',cfg.refScen);
[oarRows,oarStructures] = resolveStructureRows(cstDoseGrid, ...
    cfg.OARStructSel,'OAR',cfg.refScen);

ctx.cfg = cfg;
ctx.quantity = quantity;
ctx.scenarioDijIx = scenarioDijIx;
ctx.scenarioCtScenIds = scenarioCtScenIds;
ctx.scenarioWeights = scenarioWeights;
ctx.targetRows = targetRows;
ctx.oarRows = oarRows;
ctx.targetStructures = targetStructures(:);
ctx.oarStructures = oarStructures(:);
ctx.targetStructRows = getSelectedStructureRows(targetStructures);
ctx.oarStructRows = getSelectedStructureRows(oarStructures);
ctx.cstDoseGrid = cstDoseGrid;
ctx.numVoxels = getNumDoseVoxels(dij,quantity.matrixCell,matRad_cfg);
ctx.numBixels = getNumBixels(dij,quantity.matrixCell,scenarioDijIx,matRad_cfg);
ctx.scenarioMaps = matRad_buildScenarioDoseMappings(ct,dij, ...
    scenarioCtScenIds,cfg.refScen,matRad_cfg);
end

function cfg = normalizeConfig(cfg,ct,pln,calculationMode,matRad_cfg)
if nargin < 1 || isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRad_cfg.dispError('Scenario dose configuration must be a struct.');
end

calculationMode = upper(normalizeText(calculationMode,'calculationMode',matRad_cfg));

cfg = setDefault(cfg,'Quantity',getDefaultQuantity(pln));
cfg = setDefault(cfg,'QuantityField',[]);
cfg = setDefault(cfg,'refScen',getDefaultReferenceScenario(ct));
cfg = setDefault(cfg,'targetStructSel',[]);
cfg = setDefault(cfg,'OARStructSel',[]);
cfg = setDefault(cfg,'UseParallel',false);
cfg = setDefault(cfg,'parallelOptions',[]);
cfg = setDefault(cfg,'CollectTiming',false);
cfg = setDefault(cfg,'MemoryLimitMB',[]);
cfg = setDefault(cfg,'BatchSize',[]);
cfg = setDefault(cfg,'ProgressLevel','summary');
cfg = setDefault(cfg,'RadiusMode','std');
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
cfg.CollectTiming = logicalScalar(cfg.CollectTiming,'CollectTiming',matRad_cfg);
matRad_doseParallelPoolOptions(cfg,matRad_cfg,'parallelOptions');

if ~isempty(cfg.MemoryLimitMB) && ...
   (~isnumeric(cfg.MemoryLimitMB) || ~isscalar(cfg.MemoryLimitMB) || ...
    ~isfinite(cfg.MemoryLimitMB) || cfg.MemoryLimitMB <= 0)
    matRad_cfg.dispError('MemoryLimitMB must be a positive finite scalar.');
end

if ~isempty(cfg.BatchSize)
    cfg.BatchSize = validatePositiveInteger(cfg.BatchSize,'BatchSize',matRad_cfg);
end

cfg.ProgressLevel = validateProgressLevel(cfg.ProgressLevel,matRad_cfg);
if any(strcmp(calculationMode,{'INTERVAL2','INTERVAL3'}))
    cfg.RadiusMode = validateRadiusMode(cfg.RadiusMode,matRad_cfg);
end

if strcmp(calculationMode,'INTERVAL3')
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

function radiusMode = validateRadiusMode(radiusMode,matRad_cfg)
if isstring(radiusMode) && isscalar(radiusMode)
    radiusMode = char(radiusMode);
end

if ~ischar(radiusMode) || ...
   ~any(strcmpi(radiusMode,{'std','extreme'}))
    matRad_cfg.dispError('RadiusMode must be ''std'' or ''extreme''.');
end

radiusMode = lower(radiusMode);
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
        matRad_cfg.dispError(['RBExD scenario dose calculation requires a scalar dij.RBE ' ...
            'or an explicit linear QuantityField.']);
    end
elseif any(strcmpi(quantityName,{'effect','BED'}))
    matRad_cfg.dispError('%s scenario dose calculation requires an explicit linear QuantityField.', ...
        quantityName);
else
    matRad_cfg.dispError('Could not resolve linear scenario dose quantity ''%s''.',quantityName);
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
    matRad_cfg.dispError('Scenario dose calculation requires pln.multScen.');
end

if ~isa(pln.multScen,'matRad_ScenarioModel')
    matRad_cfg.dispError('Scenario dose calculation requires a matRad_ScenarioModel instance.');
end

scenarioIds = pln.multScen.scenarioIds();
scenarioDijIx = arrayfun(@(id) pln.multScen.getDijScenarioIndex(id),scenarioIds);
scenarioCtScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id),scenarioIds);

selectedCtScenIds = resolveScen4DSelection(pln,scenarioCtScenIds,matRad_cfg);
scenarioMask = ismember(scenarioCtScenIds,selectedCtScenIds);

if ~any(scenarioMask)
    matRad_cfg.dispError('pln.propOpt.scen4D does not select any active CT scenarios.');
end

scenarioDijIx = scenarioDijIx(scenarioMask);
scenarioCtScenIds = scenarioCtScenIds(scenarioMask);

if isempty(scenarioDijIx)
    matRad_cfg.dispError('No active scenarios found for scenario dose calculation.');
end

quantityCells = dij.(quantityField);
if max(scenarioDijIx) > numel(quantityCells)
    matRad_cfg.dispError('Scenario indices exceed dij.%s cell array dimensions.',quantityField);
end

emptyScenario = cellfun(@isempty,quantityCells(scenarioDijIx));
if any(emptyScenario)
    matRad_cfg.dispError('dij.%s contains empty active scenarios.',quantityField);
end

scenarioWeights = resolveScenarioWeights(pln.multScen,scenarioMask,matRad_cfg);
scenarioDijIx = scenarioDijIx(:);
scenarioCtScenIds = scenarioCtScenIds(:);
end

function selectedCtScenIds = resolveScen4DSelection(pln,scenarioCtScenIds,matRad_cfg)
activeCtScenIds = unique(scenarioCtScenIds(:),'stable');
selectedCtScenIds = 1;

if ~isfield(pln,'propOpt') || ~isstruct(pln.propOpt) || ...
   ~isfield(pln.propOpt,'scen4D') || isempty(pln.propOpt.scen4D)
    validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg);
    return;
end

scen4D = pln.propOpt.scen4D;
if isstring(scen4D) && isscalar(scen4D)
    scen4D = char(scen4D);
end

if ischar(scen4D)
    scen4D = strtrim(scen4D);
    if strcmpi(scen4D,'all')
        selectedCtScenIds = activeCtScenIds;
        return;
    end
    matRad_cfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
        'vector of positive integer CT scenario ids.']);
end

if ~isnumeric(scen4D) || isempty(scen4D) || any(~isfinite(scen4D(:))) || ...
   any(scen4D(:) < 1) || any(fix(scen4D(:)) ~= scen4D(:))
    matRad_cfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
        'vector of positive integer CT scenario ids.']);
end

selectedCtScenIds = unique(double(scen4D(:)),'stable');
validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg);
end

function validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg)
missingCtScenIds = selectedCtScenIds(~ismember(selectedCtScenIds,activeCtScenIds));
if ~isempty(missingCtScenIds)
    matRad_cfg.dispError('Scenario dose scen4D selection includes inactive CT scenario ids: %s.', ...
        mat2str(missingCtScenIds(:)'));
end
end

function scenarioWeights = resolveScenarioWeights(multScen,scenarioMask,matRad_cfg)
if hasFieldOrProp(multScen,'scenWeight') && ~isempty(multScen.scenWeight)
    rawWeights = multScen.scenWeight(:);
elseif hasFieldOrProp(multScen,'scenProb') && ~isempty(multScen.scenProb)
    rawWeights = multScen.scenProb(:);
else
    matRad_cfg.dispError('Scenario model must provide scenWeight or scenProb.');
end

if numel(rawWeights) ~= numel(scenarioMask)
    matRad_cfg.dispError('Number of scenario weights is inconsistent with active scenarios.');
end
scenarioWeights = rawWeights(scenarioMask);

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

function cst = prepareCstForScenarioDoseRows(cst,dij)
cst = applyOverlapPrioritiesIfAvailable(cst);
cst = resizeCstToDoseGrid(cst,dij);
end

function cst = applyOverlapPrioritiesIfAvailable(cst)
hasPriority = cellfun(@(voiProperties) isstruct(voiProperties) && ...
    isfield(voiProperties,'Priority'),cst(:,5));
if all(hasPriority)
    cst = matRad_setOverlapPriorities(cst);
end
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

function [rows,structures] = resolveStructureRows(cst,structureSelection,structureType,refScen)
if isempty(structureSelection)
    cstRows = matRad_resolveStructureSelection(cst,'all',[],structureType);
else
    cstRows = matRad_resolveStructureSelection(cst,'include',structureSelection,structureType);
end

rows = [];
structures = struct('cstRow',cell(numel(cstRows),1), ...
    'type',cell(numel(cstRows),1),'voxelIx',cell(numel(cstRows),1));
for k = 1:numel(cstRows)
    rowIx = cstRows(k);
    voxels = unique(getStructureVoxelIndices(cst,rowIx,refScen),'stable');
    structures(k).cstRow = rowIx;
    structures(k).type = structureType;
    structures(k).voxelIx = voxels(:);
    rows = [rows; voxels(:)]; %#ok<AGROW>
end
rows = unique(rows(:),'stable');
end

function cstRows = getSelectedStructureRows(structures)
if isempty(structures)
    cstRows = zeros(0,1);
else
    cstRows = [structures.cstRow]';
end
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
