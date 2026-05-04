function reportParameters = matRad_getSamplingReportScenarioParameters(multScen)
% matRad_getSamplingReportScenarioParameters derives report values from a scenario model
%
% call
%   reportParameters = matRad_getSamplingReportScenarioParameters(multScen)
%
% input
%   multScen: matRad scenario model used for sampling analysis
%
% output
%   reportParameters: struct with uncertainty values consumed by
%                     matRad_latexReport
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

if ~isa(multScen,'matRad_ScenarioModel')
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Sampling report requires a matRad_ScenarioModel instance.');
end

[shiftSize,maxAbsRangeShift,maxRelRangeShift] = getScenarioUncertaintyExtrema(multScen);

reportParameters.numOfShiftScen = multScen.totNumShiftScen;
reportParameters.shiftSize = shiftSize;
reportParameters.shiftGenType = getDimensionGenerationLabel(multScen,'setup');
reportParameters.shiftCombType = getDimensionCombinationLabel(multScen,'setup');

reportParameters.numOfRangeShiftScen = multScen.totNumRangeScen;
reportParameters.maxAbsRangeShift = maxAbsRangeShift;
reportParameters.maxRelRangeShift = maxRelRangeShift;
reportParameters.rangeCombType = getDimensionCombinationLabel(multScen,'range');
reportParameters.rangeGenType = getDimensionGenerationLabel(multScen,'range');
reportParameters.scenCombType = getScenarioCombinationLabel(multScen);

reportParameters.ctScen = logicalToText(multScen.numOfCtScen > 1);
reportParameters.rangeScen = logicalToText(multScen.totNumRangeScen > 1);
reportParameters.shiftScen = logicalToText(multScen.totNumShiftScen > 1);

reportParameters.shiftSD = multScen.shiftSD;
reportParameters.rangeAbsSD = multScen.rangeAbsSD;
reportParameters.rangeRelSD = multScen.rangeRelSD;

end

function [shiftSize,maxAbsRangeShift,maxRelRangeShift] = getScenarioUncertaintyExtrema(multScen)

scenarioIds = multScen.scenarioIds();
setupShifts = zeros(numel(scenarioIds),3);
rangeShifts = zeros(numel(scenarioIds),2);

for i = 1:numel(scenarioIds)
    setupShifts(i,:) = multScen.getSetupShift(scenarioIds(i));
    rangeShifts(i,:) = multScen.getRangeShift(scenarioIds(i));
end

shiftSize = maxFiniteValue(abs(setupShifts(:)),0);
maxAbsRangeShift = maxFiniteValue(abs(rangeShifts(:,1)),0);
maxRelRangeShift = maxFiniteValue(abs(rangeShifts(:,2)),0);

end

function label = getDimensionGenerationLabel(multScen,dimensionName)

if ~multScen.hasActiveScenarioDimension(dimensionName)
    label = 'inactive';
else
    label = char(multScen.name);
end

end

function label = getDimensionCombinationLabel(multScen,dimensionName)

if ~multScen.hasActiveScenarioDimension(dimensionName)
    label = 'inactive';
    return;
end

switch dimensionName
    case 'setup'
        if hasFieldOrProp(multScen,'combinations') && ~isempty(multScen.combinations)
            label = char(multScen.combinations);
        else
            label = char(multScen.name);
        end
    case 'range'
        if hasFieldOrProp(multScen,'combineRange')
            if multScen.combineRange
                label = 'combined';
            else
                label = 'separate';
            end
        else
            label = char(multScen.name);
        end
    otherwise
        label = char(multScen.name);
end

end

function label = getScenarioCombinationLabel(multScen)

activeDimensions = multScen.scenarioDimensionActive;
if ischar(activeDimensions)
    activeDimensions = {activeDimensions};
elseif isstring(activeDimensions)
    activeDimensions = cellstr(activeDimensions);
end
if isempty(activeDimensions)
    activeText = 'none';
else
    activeText = strjoin(activeDimensions,'+');
end

label = [char(multScen.name) ':' activeText];

end

function text = logicalToText(value)

if value
    text = 'true';
else
    text = 'false';
end

end

function maxValue = maxFiniteValue(values,defaultValue)

finiteValues = values(isfinite(values));
if isempty(finiteValues)
    maxValue = defaultValue;
else
    maxValue = max(finiteValues(:));
end

end

function tf = hasFieldOrProp(value,fieldName)

tf = (isobject(value) && isprop(value,fieldName)) || ...
    (isstruct(value) && isfield(value,fieldName));

end
