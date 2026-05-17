function dij = matRad_assembleParallelScenarioDij(scenarioDijs, scenarioIds, scenarioModel)
% matRad_assembleParallelScenarioDij assembles scenario-local dijs
%
% call
%   dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij(scenarioDijs,scenarioIds,scenarioModel)
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

matRadCfg = MatRad_Config.instance();

matRad_validateAssemblerInputs(scenarioDijs, scenarioIds, scenarioModel, matRadCfg);

scenarioIds = scenarioIds(:);
scenarioDijs = scenarioDijs(:);
numScenarios = numel(scenarioIds);
dij = scenarioDijs{1};

matRad_validateCommonDijMetadata(scenarioDijs, matRadCfg);
bioCtFields = matRad_findBiologicalCtFields(scenarioDijs, matRadCfg);
matrixFields = matRad_findScenarioMatrixFields(scenarioDijs, bioCtFields, matRadCfg);
if ~any(strcmp(matrixFields, 'physicalDose'))
    matRadCfg.dispError('Parallel scenario dij assembly requires a physicalDose matrix field.');
end
containerSize = scenarioModel.getDijContainerSize();
ctContainerSize = matRad_getCtContainerSize(scenarioModel);

fieldNames = fieldnames(dij);
for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if any(strcmp(fieldName, matrixFields))
        dij.(fieldName) = cell(containerSize);
        for scenarioIx = 1:numScenarios
            targetScenIx = scenarioModel.getDijScenarioIndex(scenarioIds(scenarioIx));
            matrix = matRad_extractSingleScenarioMatrix(scenarioDijs{scenarioIx}, ...
                                                        fieldName, matRadCfg);
            matRad_validateScenarioMatrixDimensions(matrix, fieldName, dij, matRadCfg);
            dij.(fieldName){targetScenIx} = matrix;
        end
    elseif any(strcmp(fieldName, bioCtFields))
        dij.(fieldName) = cell(ctContainerSize);
        for scenarioIx = 1:numScenarios
            ctScenId = scenarioModel.getCtScenario(scenarioIds(scenarioIx));
            value = matRad_extractSingleScenarioCtValue(scenarioDijs{scenarioIx}, ...
                                                        fieldName, ctScenId, matRadCfg);
            dij = matRad_insertCtScenarioValue(dij, fieldName, ctScenId, value, ...
                                               matRadCfg);
        end
    else
        matRad_validateCompatibleNonScenarioField(scenarioDijs, fieldName, matRadCfg);
    end
end

dij.scenarioModel = scenarioModel;
dij.scenarioIds = scenarioModel.scenarioIds();
dij.numOfScenarios = scenarioModel.numScenarios();

end

function matRad_validateAssemblerInputs(scenarioDijs, scenarioIds, scenarioModel, matRadCfg)
if ~iscell(scenarioDijs) || isempty(scenarioDijs)
    matRadCfg.dispError('Parallel scenario dij assembly requires a non-empty cell array of dijs.');
end

if ~isnumeric(scenarioIds) || isempty(scenarioIds) || any(~isfinite(scenarioIds(:))) || ...
        any(scenarioIds(:) < 1) || any(fix(scenarioIds(:)) ~= scenarioIds(:))
    matRadCfg.dispError('Parallel scenario dij assembly requires positive integer scenario ids.');
end

if numel(scenarioDijs) ~= numel(scenarioIds)
    matRadCfg.dispError('Number of scenario dijs does not match number of scenario ids.');
end

if ~isa(scenarioModel, 'matRad_ScenarioModel')
    matRadCfg.dispError('Parallel scenario dij assembly requires a matRad_ScenarioModel.');
end

scenarioIds = scenarioIds(:);
if numel(unique(scenarioIds)) ~= numel(scenarioIds)
    matRadCfg.dispError('Parallel scenario dij assembly received duplicate scenario ids.');
end

modelScenarioIds = scenarioModel.scenarioIds();
modelScenarioIds = modelScenarioIds(:);
missingScenarioIds = setdiff(modelScenarioIds, scenarioIds);
extraScenarioIds = setdiff(scenarioIds, modelScenarioIds);
if ~isempty(missingScenarioIds) || ~isempty(extraScenarioIds)
    matRadCfg.dispError(['Parallel scenario dij assembly requires exactly the active ', ...
                         'scenario ids of the scenario model. Missing: %s. Extra: %s.'], ...
                        mat2str(missingScenarioIds'), mat2str(extraScenarioIds'));
end

for i = 1:numel(scenarioDijs)
    if ~isstruct(scenarioDijs{i})
        matRadCfg.dispError('Scenario dij %d is not a struct.', i);
    end
end
end

function matRad_validateCommonDijMetadata(scenarioDijs, matRadCfg)
commonFields = {'ctGrid', 'doseGrid', 'numOfBeams', 'numOfRaysPerBeam', ...
                'totalNumOfBixels', 'totalNumOfRays', 'beamNum', 'rayNum', 'bixelNum', ...
                'minMU', 'maxMU', 'numOfParticlesPerMU'};

for i = 1:numel(commonFields)
    fieldName = commonFields{i};
    matRad_validateCompatibleNonScenarioField(scenarioDijs, fieldName, matRadCfg);
end
end

function bioCtFields = matRad_findBiologicalCtFields(scenarioDijs, matRadCfg)
fieldNames = fieldnames(scenarioDijs{1});
knownFields = {'ax', 'bx', 'abx', 'gamma', 'ixDose', 'vTissueIndex'};
bioCtFields = {};

for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if ~any(strcmp(fieldName, knownFields))
        continue
    end

    if ~iscell(scenarioDijs{1}.(fieldName))
        continue
    end

    for scenarioIx = 1:numel(scenarioDijs)
        if ~isfield(scenarioDijs{scenarioIx}, fieldName)
            matRadCfg.dispError('Scenario dij %d is missing field "%s".', ...
                                scenarioIx, fieldName);
        end
        if ~iscell(scenarioDijs{scenarioIx}.(fieldName))
            matRadCfg.dispError(['Biological CT field "%s" must be a ', ...
                                 'cell array in all scenario-local dijs.'], fieldName);
        end
    end

    bioCtFields{end + 1} = fieldName; %#ok<AGROW>
end
end

function matrixFields = matRad_findScenarioMatrixFields(scenarioDijs, bioCtFields, matRadCfg)
fieldNames = fieldnames(scenarioDijs{1});
matrixFields = {};

for fieldIx = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIx};
    if ~isfield(scenarioDijs{1}, fieldName) || ...
            ~iscell(scenarioDijs{1}.(fieldName)) || ...
            any(strcmp(fieldName, bioCtFields))
        continue
    end

    isMatrixField = true;
    for scenarioIx = 1:numel(scenarioDijs)
        if ~isfield(scenarioDijs{scenarioIx}, fieldName)
            matRadCfg.dispError('Scenario dij %d is missing field "%s".', ...
                                scenarioIx, fieldName);
        end

        try
            matrix = matRad_extractSingleScenarioMatrix(scenarioDijs{scenarioIx}, ...
                                                        fieldName, matRadCfg);
        catch
            isMatrixField = false;
            break
        end

        expectedSize = [scenarioDijs{scenarioIx}.doseGrid.numOfVoxels ...
                        scenarioDijs{scenarioIx}.totalNumOfBixels];
        if ~isequal(size(matrix), expectedSize)
            isMatrixField = false;
            break
        end
    end

    if isMatrixField
        matrixFields{end + 1} = fieldName; %#ok<AGROW>
    end
end
end

function ctContainerSize = matRad_getCtContainerSize(scenarioModel)
scenarioIds = scenarioModel.scenarioIds();
numCtScenarios = max(arrayfun(@(id) scenarioModel.getCtScenario(id), ...
                              scenarioIds(:)));

ctContainerSize = [numCtScenarios 1];
end

function matrix = matRad_extractSingleScenarioMatrix(dij, fieldName, matRadCfg)
if ~isfield(dij, fieldName) || ~iscell(dij.(fieldName))
    matRadCfg.dispError('Field "%s" is not a scenario matrix cell array.', fieldName);
end

populatedIx = find(~cellfun(@isempty, dij.(fieldName)(:)));
if numel(populatedIx) ~= 1
    matRadCfg.dispError(['Field "%s" must contain exactly one populated ', ...
                         'single-scenario matrix, but %d populated entries were found.'], ...
                        fieldName, numel(populatedIx));
end

matrix = dij.(fieldName){populatedIx};
if ~(isnumeric(matrix) || islogical(matrix)) || ndims(matrix) > 2
    matRadCfg.dispError('Field "%s" does not contain a numeric matrix.', fieldName);
end
end

function value = matRad_extractSingleScenarioCtValue(dij, fieldName, ctScenId, matRadCfg)
if ~isfield(dij, fieldName) || ~iscell(dij.(fieldName))
    matRadCfg.dispError('Field "%s" is not a biological CT cell array.', ...
                        fieldName);
end

fieldValues = dij.(fieldName);
if ctScenId <= numel(fieldValues) && ~isempty(fieldValues{ctScenId})
    value = fieldValues{ctScenId};
    return
end

populatedIx = find(~cellfun(@isempty, fieldValues(:)));
if numel(populatedIx) == 1
    value = fieldValues{populatedIx};
    return
end

matRadCfg.dispError(['Biological CT field "%s" does not contain a ', ...
                     'unique value for CT scenario %d.'], fieldName, ctScenId);
end

function dij = matRad_insertCtScenarioValue(dij, fieldName, ctScenId, value, matRadCfg)
if ctScenId > numel(dij.(fieldName))
    matRadCfg.dispError(['Biological CT field "%s" references CT ', ...
                         'scenario %d outside the assembled container.'], fieldName, ctScenId);
end

matRad_validateBiologicalCtValue(value, fieldName, dij.doseGrid.numOfVoxels, matRadCfg);

if isempty(dij.(fieldName){ctScenId})
    dij.(fieldName){ctScenId} = value;
elseif ~isequaln(dij.(fieldName){ctScenId}, value)
    matRadCfg.dispError(['Biological CT field "%s" is not compatible ', ...
                         'for CT scenario %d across scenario-local results.'], ...
                        fieldName, ctScenId);
end
end

function matRad_validateBiologicalCtValue(value, fieldName, numDoseVoxels, matRadCfg)
if isempty(value)
    matRadCfg.dispError('Biological CT field "%s" contains an empty value.', ...
                        fieldName);
end

if ~(isnumeric(value) || islogical(value)) || ~isvector(value)
    matRadCfg.dispError(['Biological CT field "%s" must contain a ', ...
                         'numeric or logical vector.'], fieldName);
end

switch fieldName
    case {'ax', 'bx', 'abx', 'gamma'}
        matRad_validateFullDoseGridVector(value, fieldName, numDoseVoxels, matRadCfg);
    case 'ixDose'
        matRad_validateDoseIndexVector(value, fieldName, numDoseVoxels, matRadCfg);
    case 'vTissueIndex'
        matRad_validateTissueIndexVector(value, fieldName, numDoseVoxels, matRadCfg);
end
end

function matRad_validateFullDoseGridVector(value, fieldName, numDoseVoxels, matRadCfg)
if numel(value) ~= numDoseVoxels
    matRadCfg.dispError(['Biological CT field "%s" has %d entries, ', ...
                         'but doseGrid.numOfVoxels is %d.'], fieldName, numel(value), ...
                        numDoseVoxels);
end
end

function matRad_validateDoseIndexVector(value, fieldName, numDoseVoxels, matRadCfg)
if islogical(value)
    if numel(value) ~= numDoseVoxels
        matRadCfg.dispError(['Logical biological CT field "%s" has %d ', ...
                             'entries, but doseGrid.numOfVoxels is %d.'], fieldName, ...
                            numel(value), numDoseVoxels);
    end
    return
end

if any(~isfinite(value(:))) || any(round(value(:)) ~= value(:)) || ...
        any(value(:) < 1) || any(value(:) > numDoseVoxels)
    matRadCfg.dispError(['Biological CT field "%s" contains invalid ', ...
                         'dose-grid indices.'], fieldName);
end
end

function matRad_validateTissueIndexVector(value, fieldName, numDoseVoxels, matRadCfg)
if numel(value) > numDoseVoxels
    matRadCfg.dispError(['Biological CT field "%s" has %d entries, ', ...
                         'which exceeds doseGrid.numOfVoxels %d.'], fieldName, numel(value), ...
                        numDoseVoxels);
end

if isnumeric(value) && (any(~isfinite(value(:))) || ...
                        any(round(value(:)) ~= value(:)) || any(value(:) < 0))
    matRadCfg.dispError(['Biological CT field "%s" contains invalid ', ...
                         'tissue indices.'], fieldName);
end
end

function matRad_validateScenarioMatrixDimensions(matrix, fieldName, dij, matRadCfg)
expectedSize = [dij.doseGrid.numOfVoxels dij.totalNumOfBixels];
if ~isequal(size(matrix), expectedSize)
    matRadCfg.dispError(['Scenario matrix field "%s" has size [%d %d], ', ...
                         'but expected [%d %d].'], fieldName, size(matrix, 1), size(matrix, 2), ...
                        expectedSize(1), expectedSize(2));
end
end

function matRad_validateCompatibleNonScenarioField(scenarioDijs, fieldName, matRadCfg)
if any(strcmp(fieldName, {'scenarioModel', 'scenarioIds', 'numOfScenarios'}))
    return
end

if ~isfield(scenarioDijs{1}, fieldName)
    return
end

referenceValue = scenarioDijs{1}.(fieldName);
for scenarioIx = 2:numel(scenarioDijs)
    if ~isfield(scenarioDijs{scenarioIx}, fieldName)
        matRadCfg.dispError('Scenario dij %d is missing field "%s".', ...
                            scenarioIx, fieldName);
    end

    if ~isequaln(referenceValue, scenarioDijs{scenarioIx}.(fieldName))
        matRadCfg.dispError(['Scenario dij field "%s" is not compatible ', ...
                             'across scenario-local results.'], fieldName);
    end
end
end
