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
%   pln:        copy of the plan with nested dose parallelism disabled and
%               propDoseCalc converted from engine handle state to a struct
%
% note
%   This helper currently sanitizes propDoseCalc only. It is intentionally
%   conservative and does not claim to remove all possible runtime state
%   from arbitrary plan fields.
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
