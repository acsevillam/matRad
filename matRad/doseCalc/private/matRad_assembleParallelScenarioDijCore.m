function dij = matRad_assembleParallelScenarioDijCore(scenarioDijs,scenarioIds,scenarioModel)
% matRad_assembleParallelScenarioDijCore assembles scenario-local dijs
%
% call
%   dij = matRad_assembleParallelScenarioDijCore(scenarioDijs,scenarioIds,scenarioModel)
%
% input
%   scenarioDijs:    cell array with one single-scenario dij per scenario
%   scenarioIds:     public scenario ids from the original scenario model
%   scenarioModel:   original multi-scenario model
%
% output
%   dij:             robust multi-scenario dij with matrix fields inserted
%                    at the original DIJ scenario indices and biological
%                    metadata inserted at the original CT scenario indices
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

matRad_cfg = MatRad_Config.instance();

validateAssemblerInputs(scenarioDijs,scenarioIds,scenarioModel,matRad_cfg);

scenarioIds = scenarioIds(:);
scenarioDijs = scenarioDijs(:);
numScenarios = numel(scenarioIds);
dij = scenarioDijs{1};

validateCommonDijMetadata(scenarioDijs,matRad_cfg);
bioCtFields = findBiologicalCtFields(scenarioDijs,matRad_cfg);
matrixFields = findScenarioMatrixFields(scenarioDijs,bioCtFields,matRad_cfg);
if ~any(strcmp(matrixFields,'physicalDose'))
    matRad_cfg.dispError('Parallel scenario dij assembly requires a physicalDose matrix field.');
end
containerSize = scenarioModel.getDijContainerSize();
ctContainerSize = getCtContainerSize(scenarioModel);

fieldNames = fieldnames(dij);
for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if any(strcmp(fieldName,matrixFields))
        dij.(fieldName) = cell(containerSize);
        for scenarioIx = 1:numScenarios
            targetScenIx = scenarioModel.getDijScenarioIndex(scenarioIds(scenarioIx));
            matrix = extractSingleScenarioMatrix(scenarioDijs{scenarioIx}, ...
                fieldName,matRad_cfg);
            validateScenarioMatrixDimensions(matrix,fieldName,dij,matRad_cfg);
            dij.(fieldName){targetScenIx} = matrix;
        end
    elseif any(strcmp(fieldName,bioCtFields))
        dij.(fieldName) = cell(ctContainerSize);
        for scenarioIx = 1:numScenarios
            ctScenId = scenarioModel.getCtScenario(scenarioIds(scenarioIx));
            value = extractSingleScenarioCtValue(scenarioDijs{scenarioIx}, ...
                fieldName,ctScenId,matRad_cfg);
            dij = insertCtScenarioValue(dij,fieldName,ctScenId,value, ...
                matRad_cfg);
        end
    else
        validateCompatibleNonScenarioField(scenarioDijs,fieldName,matRad_cfg);
    end
end

dij.scenarioModel = scenarioModel;
dij.scenarioIds = scenarioModel.scenarioIds();
dij.numOfScenarios = scenarioModel.numScenarios();

end

function validateAssemblerInputs(scenarioDijs,scenarioIds,scenarioModel,matRad_cfg)
if ~iscell(scenarioDijs) || isempty(scenarioDijs)
    matRad_cfg.dispError('Parallel scenario dij assembly requires a non-empty cell array of dijs.');
end

if numel(scenarioDijs) ~= numel(scenarioIds)
    matRad_cfg.dispError('Number of scenario dijs does not match number of scenario ids.');
end

if ~isa(scenarioModel,'matRad_ScenarioModel')
    matRad_cfg.dispError('Parallel scenario dij assembly requires a matRad_ScenarioModel.');
end

for i = 1:numel(scenarioDijs)
    if ~isstruct(scenarioDijs{i})
        matRad_cfg.dispError('Scenario dij %d is not a struct.',i);
    end
end
end

function validateCommonDijMetadata(scenarioDijs,matRad_cfg)
commonFields = {'ctGrid','doseGrid','numOfBeams','numOfRaysPerBeam', ...
    'totalNumOfBixels','totalNumOfRays','beamNum','rayNum','bixelNum', ...
    'minMU','maxMU','numOfParticlesPerMU'};

for i = 1:numel(commonFields)
    fieldName = commonFields{i};
    validateCompatibleNonScenarioField(scenarioDijs,fieldName,matRad_cfg);
end
end

function bioCtFields = findBiologicalCtFields(scenarioDijs,matRad_cfg)
fieldNames = fieldnames(scenarioDijs{1});
knownFields = {'ax','bx','abx','gamma','ixDose','vTissueIndex'};
bioCtFields = {};

for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if ~any(strcmp(fieldName,knownFields))
        continue;
    end

    if ~iscell(scenarioDijs{1}.(fieldName))
        continue;
    end

    for scenarioIx = 1:numel(scenarioDijs)
        if ~isfield(scenarioDijs{scenarioIx},fieldName)
            matRad_cfg.dispError('Scenario dij %d is missing field "%s".', ...
                scenarioIx,fieldName);
        end
        if ~iscell(scenarioDijs{scenarioIx}.(fieldName))
            matRad_cfg.dispError(['Biological CT field "%s" must be a ', ...
                'cell array in all scenario-local dijs.'],fieldName);
        end
    end

    bioCtFields{end+1} = fieldName; %#ok<AGROW>
end
end

function matrixFields = findScenarioMatrixFields(scenarioDijs,bioCtFields,matRad_cfg)
fieldNames = fieldnames(scenarioDijs{1});
matrixFields = {};

for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if ~isfield(scenarioDijs{1},fieldName) || ...
            ~iscell(scenarioDijs{1}.(fieldName)) || ...
            any(strcmp(fieldName,bioCtFields))
        continue;
    end

    isMatrixField = true;
    for scenarioIx = 1:numel(scenarioDijs)
        if ~isfield(scenarioDijs{scenarioIx},fieldName)
            matRad_cfg.dispError('Scenario dij %d is missing field "%s".', ...
                scenarioIx,fieldName);
        end

        try
            matrix = extractSingleScenarioMatrix(scenarioDijs{scenarioIx}, ...
                fieldName,matRad_cfg);
        catch
            isMatrixField = false;
            break;
        end

        expectedSize = [scenarioDijs{scenarioIx}.doseGrid.numOfVoxels ...
            scenarioDijs{scenarioIx}.totalNumOfBixels];
        if ~isequal(size(matrix),expectedSize)
            isMatrixField = false;
            break;
        end
    end

    if isMatrixField
        matrixFields{end+1} = fieldName; %#ok<AGROW>
    end
end
end

function ctContainerSize = getCtContainerSize(scenarioModel)
numCtScenarios = [];
if matRad_ispropCompat(scenarioModel,'numOfAvailableCtScen') && ...
        ~isempty(scenarioModel.numOfAvailableCtScen)
    numCtScenarios = scenarioModel.numOfAvailableCtScen;
elseif matRad_ispropCompat(scenarioModel,'numOfCtScen') && ...
        ~isempty(scenarioModel.numOfCtScen)
    numCtScenarios = scenarioModel.numOfCtScen;
elseif matRad_ispropCompat(scenarioModel,'ctScenProb') && ...
        ~isempty(scenarioModel.ctScenProb)
    numCtScenarios = max(scenarioModel.ctScenProb(:,1));
end

if isempty(numCtScenarios)
    scenarioIds = scenarioModel.scenarioIds();
    numCtScenarios = max(arrayfun(@(id) scenarioModel.getCtScenario(id), ...
        scenarioIds));
end

ctContainerSize = [numCtScenarios 1];
end

function matrix = extractSingleScenarioMatrix(dij,fieldName,matRad_cfg)
if ~isfield(dij,fieldName) || ~iscell(dij.(fieldName))
    matRad_cfg.dispError('Field "%s" is not a scenario matrix cell array.',fieldName);
end

populatedIx = find(~cellfun(@isempty,dij.(fieldName)(:)));
if numel(populatedIx) ~= 1
    matRad_cfg.dispError(['Field "%s" must contain exactly one populated ', ...
        'single-scenario matrix, but %d populated entries were found.'], ...
        fieldName,numel(populatedIx));
end

matrix = dij.(fieldName){populatedIx};
if ~(isnumeric(matrix) || islogical(matrix)) || ndims(matrix) > 2
    matRad_cfg.dispError('Field "%s" does not contain a numeric matrix.',fieldName);
end
end

function value = extractSingleScenarioCtValue(dij,fieldName,ctScenId,matRad_cfg)
if ~isfield(dij,fieldName) || ~iscell(dij.(fieldName))
    matRad_cfg.dispError('Field "%s" is not a biological CT cell array.', ...
        fieldName);
end

fieldValues = dij.(fieldName);
if ctScenId <= numel(fieldValues) && ~isempty(fieldValues{ctScenId})
    value = fieldValues{ctScenId};
    return;
end

populatedIx = find(~cellfun(@isempty,fieldValues(:)));
if numel(populatedIx) == 1
    value = fieldValues{populatedIx};
    return;
end

matRad_cfg.dispError(['Biological CT field "%s" does not contain a ', ...
    'unique value for CT scenario %d.'],fieldName,ctScenId);
end

function dij = insertCtScenarioValue(dij,fieldName,ctScenId,value,matRad_cfg)
if ctScenId > numel(dij.(fieldName))
    matRad_cfg.dispError(['Biological CT field "%s" references CT ', ...
        'scenario %d outside the assembled container.'],fieldName,ctScenId);
end

validateBiologicalCtValue(value,fieldName,dij.doseGrid.numOfVoxels,matRad_cfg);

if isempty(dij.(fieldName){ctScenId})
    dij.(fieldName){ctScenId} = value;
elseif ~isequaln(dij.(fieldName){ctScenId},value)
    matRad_cfg.dispError(['Biological CT field "%s" is not compatible ', ...
        'for CT scenario %d across scenario-local results.'], ...
        fieldName,ctScenId);
end
end

function validateBiologicalCtValue(value,fieldName,numDoseVoxels,matRad_cfg)
if isempty(value)
    matRad_cfg.dispError('Biological CT field "%s" contains an empty value.', ...
        fieldName);
end

if ~(isnumeric(value) || islogical(value)) || ~isvector(value)
    matRad_cfg.dispError(['Biological CT field "%s" must contain a ', ...
        'numeric or logical vector.'],fieldName);
end

switch fieldName
    case {'ax','bx','abx','gamma'}
        validateFullDoseGridVector(value,fieldName,numDoseVoxels,matRad_cfg);
    case 'ixDose'
        validateDoseIndexVector(value,fieldName,numDoseVoxels,matRad_cfg);
    case 'vTissueIndex'
        validateTissueIndexVector(value,fieldName,numDoseVoxels,matRad_cfg);
end
end

function validateFullDoseGridVector(value,fieldName,numDoseVoxels,matRad_cfg)
if numel(value) ~= numDoseVoxels
    matRad_cfg.dispError(['Biological CT field "%s" has %d entries, ', ...
        'but doseGrid.numOfVoxels is %d.'],fieldName,numel(value), ...
        numDoseVoxels);
end
end

function validateDoseIndexVector(value,fieldName,numDoseVoxels,matRad_cfg)
if islogical(value)
    if numel(value) ~= numDoseVoxels
        matRad_cfg.dispError(['Logical biological CT field "%s" has %d ', ...
            'entries, but doseGrid.numOfVoxels is %d.'],fieldName, ...
            numel(value),numDoseVoxels);
    end
    return;
end

if any(~isfinite(value(:))) || any(round(value(:)) ~= value(:)) || ...
        any(value(:) < 1) || any(value(:) > numDoseVoxels)
    matRad_cfg.dispError(['Biological CT field "%s" contains invalid ', ...
        'dose-grid indices.'],fieldName);
end
end

function validateTissueIndexVector(value,fieldName,numDoseVoxels,matRad_cfg)
if numel(value) > numDoseVoxels
    matRad_cfg.dispError(['Biological CT field "%s" has %d entries, ', ...
        'which exceeds doseGrid.numOfVoxels %d.'],fieldName,numel(value), ...
        numDoseVoxels);
end

if isnumeric(value) && (any(~isfinite(value(:))) || ...
        any(round(value(:)) ~= value(:)) || any(value(:) < 0))
    matRad_cfg.dispError(['Biological CT field "%s" contains invalid ', ...
        'tissue indices.'],fieldName);
end
end

function validateScenarioMatrixDimensions(matrix,fieldName,dij,matRad_cfg)
expectedSize = [dij.doseGrid.numOfVoxels dij.totalNumOfBixels];
if ~isequal(size(matrix),expectedSize)
    matRad_cfg.dispError(['Scenario matrix field "%s" has size [%d %d], ', ...
        'but expected [%d %d].'],fieldName,size(matrix,1),size(matrix,2), ...
        expectedSize(1),expectedSize(2));
end
end

function validateCompatibleNonScenarioField(scenarioDijs,fieldName,matRad_cfg)
if any(strcmp(fieldName,{'scenarioModel','scenarioIds','numOfScenarios'}))
    return;
end

if ~isfield(scenarioDijs{1},fieldName)
    return;
end

referenceValue = scenarioDijs{1}.(fieldName);
for scenarioIx = 2:numel(scenarioDijs)
    if ~isfield(scenarioDijs{scenarioIx},fieldName)
        matRad_cfg.dispError('Scenario dij %d is missing field "%s".', ...
            scenarioIx,fieldName);
    end

    if ~isequaln(referenceValue,scenarioDijs{scenarioIx}.(fieldName))
        matRad_cfg.dispError(['Scenario dij field "%s" is not compatible ', ...
            'across scenario-local results.'],fieldName);
    end
end
end
