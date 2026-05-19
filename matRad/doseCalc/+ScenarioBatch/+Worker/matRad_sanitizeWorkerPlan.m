function pln = matRad_sanitizeWorkerPlan(pln, engine)
% matRad_sanitizeWorkerPlan removes runtime state from a worker plan
%
% call
%   pln = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(pln)
%   pln = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(pln,engine)
%
% input
%   pln:        matRad plan struct
%   engine:     optional dose engine whose public configuration
%               should be serialized into pln.propDoseCalc
%
% output
%   pln:        copy of the plan with nested dose parallelism disabled,
%               propDoseCalc converted from engine handle state to a struct,
%               and biological model handles converted to serializable structs
%
% note
%   This helper is intentionally conservative and only sanitizes worker
%   boundary state that is known to be unsafe or unnecessary to serialize.
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

if nargin < 2
    engine = [];
end

pln.propDoseCalc = matRad_buildSerialPropDoseCalc(pln, engine);
pln = matRad_sanitizeBiologicalModelFields(pln);
end

function propDoseCalc = matRad_buildSerialPropDoseCalc(pln, engine)
if isa(engine, 'DoseEngines.matRad_DoseEngineBase')
    propDoseCalc = matRad_doseEngineToPropDoseCalcStruct(engine);
elseif isfield(pln, 'propDoseCalc') && ...
        isa(pln.propDoseCalc, 'DoseEngines.matRad_DoseEngineBase')
    propDoseCalc = matRad_doseEngineToPropDoseCalcStruct(pln.propDoseCalc);
elseif isfield(pln, 'propDoseCalc') && isstruct(pln.propDoseCalc)
    propDoseCalc = pln.propDoseCalc;
else
    propDoseCalc = struct();
end

propDoseCalc.UseParallel = false;
end

function propDoseCalc = matRad_doseEngineToPropDoseCalcStruct(engine)
propDoseCalc = struct();
propDoseCalc.engine = char(engine.shortName);

propNames = matRad_publicSettableDoseEngineProperties(engine);
excludedProps = {'UseParallel', 'parallelOptions', 'multScen'};

for i = 1:numel(propNames)
    propName = propNames{i};
    if any(strcmp(propName, excludedProps))
        continue
    end
    try
        propDoseCalc.(propName) = engine.(propName);
    catch
        % Keep the worker plan conservative if a custom property cannot be
        % read outside the engine class.
    end
end
end

function pln = matRad_sanitizeBiologicalModelFields(pln)
[pln, ~] = matRad_sanitizeBioModelField(pln, 'bioModel');
[pln, bioParamMetadata] = matRad_sanitizeBioModelField(pln, 'bioParam');

if ~isfield(pln, 'bioModel') && ~isempty(bioParamMetadata)
    pln.bioModel = bioParamMetadata;
end

if isfield(pln, 'propDoseCalc') && isstruct(pln.propDoseCalc)
    pln.propDoseCalc = ...
        matRad_sanitizePropDoseCalcBiologicalModelFields(pln.propDoseCalc);
end
end

function [valueStruct, metadata] = matRad_sanitizeBioModelField(valueStruct, fieldName)
metadata = [];
if ~isfield(valueStruct, fieldName)
    return
end

[valueStruct.(fieldName), metadata] = ...
    matRad_sanitizeBioModelValue(valueStruct.(fieldName));
end

function propDoseCalc = matRad_sanitizePropDoseCalcBiologicalModelFields(propDoseCalc)
for elementIx = 1:numel(propDoseCalc)
    [propDoseCalc(elementIx), ~] = ...
        matRad_sanitizeBioModelField(propDoseCalc(elementIx), 'bioModel');
    [propDoseCalc(elementIx), bioParamMetadata] = ...
        matRad_sanitizeBioModelField(propDoseCalc(elementIx), 'bioParam');
    if ~isfield(propDoseCalc(elementIx), 'bioModel') && ...
            ~isempty(bioParamMetadata)
        propDoseCalc(elementIx).bioModel = bioParamMetadata;
    end
end
end

function [value, metadata] = matRad_sanitizeBioModelValue(value)
metadata = [];
if isa(value, 'matRad_BiologicalModel')
    metadata = matRad_bioModelToSerializableStruct(value);
    value = metadata;
end
end

function metadata = matRad_bioModelToSerializableStruct(bioModel)
metadata = struct();
metadata.model = matRad_scalarText(bioModel.model);

propNames = matRad_publicSerializableBioModelProperties(bioModel);
excludedProps = {'model', 'quantityOpt', 'quantityVis'};
for i = 1:numel(propNames)
    propName = propNames{i};
    if any(strcmp(propName, excludedProps))
        continue
    end
    try
        propValue = bioModel.(propName);
    catch
        continue
    end
    if matRad_isSerializableBioModelValue(propValue)
        metadata.(propName) = propValue;
    end
end
end

function propNames = matRad_publicSerializableBioModelProperties(bioModel)
propNames = {};
try
    metaClass = metaclass(bioModel);
    propertyList = metaClass.PropertyList;
catch
    return
end

for i = 1:numel(propertyList)
    metaProp = propertyList(i);
    if iscell(propertyList)
        metaProp = propertyList{i};
    end
    if ~matRad_isPublicAccess(metaProp.GetAccess) || ...
            ~matRad_isPublicAccess(metaProp.SetAccess) || ...
            matRad_metaFlag(metaProp, 'Constant') || ...
            matRad_metaFlag(metaProp, 'Dependent') || ...
            matRad_metaFlag(metaProp, 'Hidden')
        continue
    end
    propNames{end + 1} = metaProp.Name; %#ok<AGROW>
end

propNames = unique(propNames, 'stable');
end

function tf = matRad_isSerializableBioModelValue(value)
tf = ~isa(value, 'handle') && ~isa(value, 'function_handle');
end

function value = matRad_scalarText(value)
if isstring(value) && isscalar(value)
    value = char(value);
elseif ~ischar(value)
    value = '';
end
end

function propNames = matRad_publicSettableDoseEngineProperties(engine)
try
    propNames = matRad_publicSettableDoseEnginePropertiesFromMetaClass(engine);
catch
    propNames = matRad_fallbackPublicSettableDoseEngineProperties(engine);
end

propNames = unique(propNames, 'stable');
end

function propNames = matRad_publicSettableDoseEnginePropertiesFromMetaClass(engine)
propNames = {};
metaClass = metaclass(engine);
propertyList = metaClass.PropertyList;

for i = 1:numel(propertyList)
    metaProp = propertyList(i);
    if iscell(propertyList)
        metaProp = propertyList{i};
    end
    if ~matRad_isPublicAccess(metaProp.GetAccess) || ...
            ~matRad_isPublicAccess(metaProp.SetAccess) || ...
            matRad_metaFlag(metaProp, 'Constant') || ...
            matRad_metaFlag(metaProp, 'Dependent')
        continue
    end
    propNames{end + 1} = metaProp.Name; %#ok<AGROW>
end
end

function tf = matRad_isPublicAccess(accessValue)
if ischar(accessValue)
    tf = strcmp(accessValue, 'public');
else
    tf = strcmp(char(accessValue), 'public');
end
end

function tf = matRad_metaFlag(metaProp, fieldName)
try
    tf = logical(metaProp.(fieldName));
catch
    tf = false;
end
end

function propNames = matRad_fallbackPublicSettableDoseEngineProperties(engine)
candidateNames = { ...
                  'doseGrid', 'voxelSubIx', 'selectVoxelsInScenarios', 'bioParam', ...
                  'precision', 'enableGPU', ...
                  'keepRadDepthCubes', 'geometricLateralCutOff', 'dosimetricLateralCutOff', ...
                  'ssdDensityThreshold', 'useGivenEqDensityCube', 'ignoreOutsideDensities', ...
                  'numOfDijFillSteps', ...
                  'useCustomPrimaryPhotonFluence', 'kernelCutOff', 'randomSeed', ...
                  'intConvResolution', 'enableDijSampling', 'dijSampling', ...
                  'calcLET', 'calcBioDose', 'airOffsetCorrection', 'lateralModel', ...
                  'visBoolLateralCutOff', 'fineSampling', ...
                  'massDensity', 'p', 'alpha', 'beta', 'gammaNuc', 'Z', 'MM', 'radLength', ...
                  'phi0', 'epsilonTail', 'sigmaEnergy', 'modeWidth'};

propNames = {};
for i = 1:numel(candidateNames)
    if matRad_ispropCompat(engine, candidateNames{i})
        propNames{end + 1} = candidateNames{i}; %#ok<AGROW>
    end
end
end
