function [parameterLines, modelSummaryLines, componentSummaryLines] = matRad_buildScenarioReportParameters(multScen)
% matRad_buildScenarioReportParameters builds scenario model LaTeX summaries.
%
% call
%   [parameterLines,modelSummaryLines,componentSummaryLines] = matRad_buildScenarioReportParameters(multScen)
%
% input
%   multScen:             matRad scenario model
%
% output
%   parameterLines:       cell array with LaTeX macros for scenario metadata
%   modelSummaryLines:    cell array with a LaTeX table summarizing the model
%   componentSummaryLines: cell array with a LaTeX table summarizing active
%                         uncertainty components
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

if ~isa(multScen, 'matRad_ScenarioModel')
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Scenario report parameters require a matRad_ScenarioModel instance.');
end

modelName = helper_propertyText(multScen, 'name', 'Scenario Model');
storagePolicy = helper_propertyText(multScen, 'scenarioStoragePolicy', 'unknown');
activeDimensions = helper_joinText(multScen.scenarioDimensionActive);
ctScenarioIds = unique(multScen.scenarioCtScenIds(:))';
probabilities = multScen.scenProb(:);

parameterLines = helper_buildParameterLines(modelName, storagePolicy, ...
                                            activeDimensions, multScen.numScenarios());
modelSummaryLines = helper_buildModelSummaryLines(modelName, storagePolicy, ...
                                                  activeDimensions, ctScenarioIds, probabilities, multScen.numScenarios());
componentSummaryLines = helper_buildComponentSummaryLines(multScen);

end

function lines = helper_buildParameterLines(modelName, storagePolicy, activeDimensions, numScenarios)

lines = { ...
         ['\newcommand{\scenarioModelName}{', helper_latexText(modelName), '}']; ...
         ['\newcommand{\scenarioStoragePolicy}{', helper_latexText(storagePolicy), '}']; ...
         ['\newcommand{\scenarioActiveDimensions}{', helper_latexText(activeDimensions), '}']; ...
         ['\newcommand{\numOfTotalScen}{', num2str(numScenarios), '}']};

end

function lines = helper_buildModelSummaryLines(modelName, storagePolicy, ...
                                               activeDimensions, ctScenarioIds, probabilities, numScenarios)

if isempty(probabilities)
    probabilitySum = 'N.A.';
    probabilityMin = 'N.A.';
    probabilityMax = 'N.A.';
else
    probabilitySum = helper_numberText(sum(probabilities));
    probabilityMin = helper_numberText(min(probabilities));
    probabilityMax = helper_numberText(max(probabilities));
end

lines = { ...
         '\begin{table}[!h]'; ...
         '  \centering'; ...
         '  \label{table:scenarioModelSummary}'; ...
         '  \begin{tabular}{rl}'; ...
         '    \toprule'; ...
         ['    Scenario model & ', helper_latexText(modelName), ' \\']; ...
         ['    Storage policy & ', helper_latexText(storagePolicy), ' \\']; ...
         ['    Total scenarios & ', num2str(numScenarios), ' \\']; ...
         ['    Active dimensions & ', helper_latexText(activeDimensions), ' \\']; ...
         ['    Active CT scenario IDs & ', helper_latexText(helper_joinNumbers(ctScenarioIds)), ' \\']; ...
         ['    Probability sum & ', probabilitySum, ' \\']; ...
         ['    Minimum scenario probability & ', probabilityMin, ' \\']; ...
         ['    Maximum scenario probability & ', probabilityMax, ' \\']; ...
         '    \bottomrule'; ...
         '  \end{tabular}'; ...
         '\end{table}'};

end

function lines = helper_buildComponentSummaryLines(multScen)

lines = { ...
         '\begin{table}[!h]'; ...
         '  \centering'; ...
         '  \small'; ...
         '  \label{table:scenarioComponentSummary}'; ...
         '  \begin{tabular}{llllrrrr}'; ...
         '    \toprule'; ...
         '    Dimension & Component & Unit & Scale & Min. & Max. & Mean & Std. dev. \\'; ...
         '    \midrule'};

componentLines = helper_componentRows(multScen);
if isempty(componentLines)
    componentLines = {'    \multicolumn{8}{l}{No active numeric scenario components.} \\'};
end

lines = [lines; componentLines(:); { ...
                                    '    \bottomrule'; ...
                                    '  \end{tabular}'; ...
                                    '\end{table}'}];

end

function lines = helper_componentRows(multScen)

components = multScen.scenarioComponents;
scenarioValues = multScen.scenarioValues;
lines = cell(0);

for i = 1:numel(components)
    component = components(i);
    if isfield(component, 'active') && ~component.active
        continue
    end

    values = scenarioValues(:, i);
    lines{end + 1, 1} = sprintf('    %s & %s & %s & %s & %s & %s & %s & %s \\\\', ...
                                helper_latexText(helper_componentField(component, 'applicator', 'unknown')), ...
                                helper_latexText(helper_componentField(component, 'name', sprintf('component%d', i))), ...
                                helper_latexText(helper_componentField(component, 'unit', '')), ...
                                helper_numberText(helper_componentField(component, 'scale', NaN)), ...
                                helper_numberText(min(values)), helper_numberText(max(values)), ...
                                helper_numberText(mean(values)), helper_numberText(std(values, 0)));
end

end

function text = helper_propertyText(object, propertyName, defaultText)

if isprop(object, propertyName)
    text = helper_toText(object.(propertyName));
else
    text = defaultText;
end

end

function value = helper_componentField(component, fieldName, defaultValue)

if isfield(component, fieldName)
    value = component.(fieldName);
else
    value = defaultValue;
end

end

function text = helper_joinText(values)

if isempty(values)
    text = 'none';
elseif ischar(values)
    text = values;
else
    text = helper_toText(values{1});
    for i = 2:numel(values)
        text = [text, ', ', helper_toText(values{i})]; %#ok<AGROW>
    end
end

end

function text = helper_joinNumbers(values)

if isempty(values)
    text = 'none';
else
    text = num2str(values, '%g, ');
    text = regexprep(strtrim(text), ',\s*$', '');
end

end

function text = helper_numberText(value)

if isempty(value) || any(~isfinite(value(:)))
    text = 'N.A.';
else
    text = num2str(value, '%.6g');
end

end

function text = helper_toText(value)

if isstring(value) && isscalar(value)
    text = char(value);
elseif ischar(value)
    text = value;
elseif isnumeric(value) || islogical(value)
    text = num2str(value);
else
    text = class(value);
end

end

function text = helper_latexText(value)

text = helper_toText(value);
text = strrep(text, '\', '\textbackslash{}');
text = strrep(text, '&', '\&');
text = strrep(text, '%', '\%');
text = strrep(text, '$', '\$');
text = strrep(text, '#', '\#');
text = strrep(text, '_', '\_');

end
