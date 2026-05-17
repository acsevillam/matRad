function test_suite = test_scenarioConsumers

test_functions = localfunctions();

initTestSuite;

function test_planAnalysisAcceptsCompactScenarioModel
figureCleaner = onCleanup(@() close('all'));
ct.numOfCtScen = 2;
cst = helper_createCst();
pln = helper_createPlan(ct);
stf = struct();
resultGUI.physicalDose = [1; 2];
resultGUI.physicalDose_scen1 = [1; 3];
resultGUI.physicalDose_scen2 = [2; 4];
resultGUI.physicalDose_scen3 = [3; 5];

resultGUI = matRad_planAnalysis(resultGUI, ct, cst, stf, pln, ...
                                'refGy', 1, 'refVol', 50, ...
                                'showDVH', false, 'showQI', false);

assertTrue(isfield(resultGUI, 'dvh'));
assertTrue(isfield(resultGUI, 'qi'));

function test_latexReportBuildsCompactScenarioParameters
ct.numOfCtScen = 2;
pln = helper_createPlan(ct);

[parameterLines, modelLines, componentLines] = ...
    matRad_buildScenarioReportParameters(pln.multScen);
reportText = helper_linesToText([parameterLines; modelLines; componentLines]);

assertTrue(~isempty(strfind(reportText, 'compact-realization')));
assertTrue(~isempty(strfind(reportText, 'ct, setup, range, gantry, couch')));
assertTrue(~isempty(strfind(reportText, 'gantry.beam1')));
assertTrue(~isempty(strfind(reportText, 'couch.beam2')));
assertTrue(isempty(strfind(reportText, 'Max setup shift')));

function test_latexReportBuildsLegacyScenarioParameters
ct.numOfCtScen = 2;
model = matRad_NominalScenario(ct);

[parameterLines, modelLines, componentLines] = ...
    matRad_buildScenarioReportParameters(model);
reportText = helper_linesToText([parameterLines; modelLines; componentLines]);

assertTrue(~isempty(strfind(reportText, 'legacy-grid')));
assertTrue(~isempty(strfind(reportText, 'ct, setup, range')));
assertTrue(~isempty(strfind(reportText, 'setup.x')));
assertTrue(~isempty(strfind(reportText, 'range.absolute')));

function test_latexReportTemplateUsesGenericScenarioTables
matRad_cfg = MatRad_Config.instance();
templatePath = fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', ...
                        'samplingAnalysis', 'main_template.tex');
templateText = fileread(templatePath);

assertTrue(~isempty(strfind(templateText, 'scenarioModelSummary.tex')));
assertTrue(~isempty(strfind(templateText, 'scenarioComponentSummary.tex')));
assertTrue(isempty(strfind(templateText, '\numOfShiftScen')));
assertTrue(isempty(strfind(templateText, '\shiftCombType')));
assertTrue(isempty(strfind(templateText, '\rangeCombType')));
assertTrue(isempty(strfind(templateText, '\shiftSD')));
assertTrue(isempty(strfind(templateText, '\rangeAbsSD')));
assertTrue(isempty(strfind(templateText, '\rangeRelSD')));

function test_samplingSourceUsesCanonicalScenarioMetadata
matRad_cfg = MatRad_Config.instance();
samplingPath = fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', ...
                        'matRad_sampling.m');
samplingHelperPath = fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', ...
                              'sampling', 'matRad_calculateSamplingScenario.m');
samplingText = [fileread(samplingPath) newline fileread(samplingHelperPath)];

assertTrue(~isempty(strfind(samplingText, 'scenarioId')));
assertTrue(~isempty(strfind(samplingText, 'ctScenId')));
assertTrue(~isempty(strfind(samplingText, 'scenario')));
assertTrue(~isempty(strfind(samplingText, 'doseMapping')));
assertTrue(isempty(strfind(samplingText, 'caSampRes(i).isoShift')));
assertTrue(isempty(strfind(samplingText, 'caSampRes(i).relRangeShift')));
assertTrue(isempty(strfind(samplingText, 'caSampRes(i).absRangeShift')));
assertTrue(isempty(strfind(samplingText, 'sampleResult.bioModel')));
assertTrue(isempty(strfind(samplingText, 'bioParam.quantityOpt')));
assertTrue(isempty(strfind(samplingText, 'bioParam.quantityVis')));

function test_samplingScenarioWeightsUseCanonicalProbabilities
ct.numOfCtScen = 2;
pln = helper_createPlan(ct);

scenWeights = matRad_getSamplingScenarioWeights(pln, pln.multScen.numScenarios());

assertElementsAlmostEqual(scenWeights, [0.2; 0.3; 0.5], 'absolute', 1e-12);

function test_samplingScenarioWeightsAcceptExplicitWeights
ct.numOfCtScen = 2;
pln = helper_createPlan(ct);

scenWeights = matRad_getSamplingScenarioWeights(pln, pln.multScen.numScenarios(), [], [2 1 1]);

assertElementsAlmostEqual(scenWeights, [0.5; 0.25; 0.25], 'absolute', 1e-12);

function test_showDVHFromSamplingBuildsScaledTrustbandWithExplicitWeights
figureCleaner = onCleanup(@() close('all'));
fig = figure('Visible', 'off');
ct.numOfCtScen = 2;
pln = helper_createPlan(ct);
cst = helper_createCst();
caSamp = helper_createSamplingDvh();
scale = 2;
doseWindow = [0 20];
scenWeights = [2 1 1];

matRad_showDVHFromSampling(caSamp, scale, cst, pln, 1:numel(caSamp), ...
                           doseWindow, 'trustband', 1, 1, ...
                           'scenWeights', scenWeights);

assertTrue(ishandle(fig));
helper_assertTrustbandPatch(gca, caSamp, 1:numel(caSamp), scenWeights, scale);
assertElementsAlmostEqual(xlim(gca), doseWindow, 'absolute', 1e-12);
helper_assertTrustbandLinesOnTop(gca);

function test_showDVHFromSamplingUsesSelectedScenarioWeights
figureCleaner = onCleanup(@() close('all'));
fig = figure('Visible', 'off');
ct.numOfCtScen = 2;
pln = helper_createPlan(ct);
cst = helper_createCst();
caSamp = helper_createSamplingDvh();
scenarios = [1 3];
scenWeights = [0.25 0.25 0.5];

matRad_showDVHFromSampling(caSamp, 1, cst, pln, scenarios, ...
                           [0 10], 'trustband', 1, 1, ...
                           'scenWeights', scenWeights);

assertTrue(ishandle(fig));
helper_assertTrustbandPatch(gca, caSamp, scenarios, scenWeights, 1);
helper_assertTrustbandLinesOnTop(gca);

function cst = helper_createCst()
cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'testVoi';
cst{1, 3} = 'OAR';
cst{1, 4} = {[1; 2]};
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {};

function pln = helper_createPlan(ct)
pln.numOfFractions = 1;
pln.bioModel = 'none';
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.quantityVis = 'physicalDose';
pln.multScen = helper_createCompactScenarioModel(ct);

function caSamp = helper_createSamplingDvh()
doseGrid = [0 5 10];
volumePoints = [100 50 0; ...
                100 60 0; ...
                100 40 0];

for i = 1:size(volumePoints, 1)
    caSamp(i).dvh(1).name = 'testVoi';
    caSamp(i).dvh(1).doseGrid = doseGrid;
    caSamp(i).dvh(1).volumePoints = volumePoints(i, :);
end

function model = helper_createCompactScenarioModel(ct)
model = matRad_NominalScenario(ct);
activeDimensions = {'ct', 'setup', 'range', 'gantry', 'couch'};
model.numOfBeams = 2;
model.scenarioDimensionActive = activeDimensions;
components = matRad_createScenarioComponents([1 1 1], 1, 1, ...
                                             activeDimensions, 2);
scenarioValues = zeros(3, numel(components));
scenarioValues(:, strcmp({components.name}, 'gantry.beam1')) = [0; 1; -1];
scenarioValues(:, strcmp({components.name}, 'couch.beam2')) = [0; -2; 2];
ctScenIds = [2; 1; 2];
scenProb = [0.2; 0.3; 0.5];
scenWeight = scenProb;
scenForProb = [ctScenIds scenarioValues];
storageSubscripts = [(1:3)' ones(3, 1)];
storageSize = [3 1];

model.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                              scenProb, scenWeight, scenForProb, ...
                              storageSubscripts, storageSize, ...
                              'compact-realization');

function helper_assertTrustbandPatch(ax, caSamp, scenarios, scenWeights, scale)
patchHandles = findobj(ax, 'Type', 'patch');
assertEqual(numel(patchHandles), 1);

[expectedDoseGrid, expectedLower, expectedUpper] = ...
    helper_expectedTrustband(caSamp, scenarios, scenWeights, scale);
expectedXData = [expectedDoseGrid fliplr(expectedDoseGrid)];
expectedYData = [expectedLower fliplr(expectedUpper)];

assertElementsAlmostEqual(get(patchHandles, 'XData'), expectedXData(:), 'absolute', 1e-12);
assertElementsAlmostEqual(get(patchHandles, 'YData'), expectedYData(:), 'absolute', 1e-12);

function helper_assertTrustbandLinesOnTop(ax)
patchHandle = findobj(ax, 'Type', 'patch');
lineHandles = findobj(ax, 'Type', 'line');
assertEqual(numel(patchHandle), 1);
assertTrue(numel(lineHandles) >= 2);

children = get(ax, 'Children');
patchChildIx = find(children == patchHandle, 1);
lineChildIx = arrayfun(@(h) find(children == h, 1), lineHandles);
assertTrue(all(lineChildIx < patchChildIx));

function [doseGrid, lowerBound, upperBound] = ...
    helper_expectedTrustband(caSamp, scenarios, scenWeights, scale)
doseGrid = caSamp(1).dvh(1).doseGrid * scale;
allDvh = zeros(numel(scenarios), numel(doseGrid));

for s = 1:numel(scenarios)
    allDvh(s, :) = caSamp(scenarios(s)).dvh(1).volumePoints;
end

if numel(scenWeights) == numel(caSamp)
    currScenWeights = scenWeights(scenarios);
else
    currScenWeights = scenWeights;
end

weightedMean = helper_weightedMean(allDvh, currScenWeights);
weightedStd = sqrt(helper_weightedMean(bsxfun(@minus, allDvh, weightedMean).^2, ...
                                       currScenWeights));
lowerBound = max(weightedMean - weightedStd, 0);
upperBound = min(weightedMean + weightedStd, 100);

function meanValue = helper_weightedMean(values, weights)
weights = weights(:);
meanValue = weights' * values ./ sum(weights);

function text = helper_linesToText(lines)
text = sprintf('%s\n', lines{:});
