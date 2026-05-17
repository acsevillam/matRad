function test_suite = test_calcQualityIndicators

test_functions = localfunctions();

initTestSuite;

function test_covUsesPerFractionReferenceDose
[cst, pln] = helper_qiFixture(78, 39, 10);
doseCube = 2 * ones(10, 1);

qi = matRad_calcQualityIndicators(cst, pln, doseCube, [], []);

assertElementsAlmostEqual(qi.COV1, 1);
assertElementsAlmostEqual(qi.referenceDose, 2);

function test_covThresholdsAreFractional
[cst, pln] = helper_qiFixture(78, 39, 10);
doseCube = 0.99 * 2 * ones(10, 1);

qi = matRad_calcQualityIndicators(cst, pln, doseCube, [], []);

assertElementsAlmostEqual(qi.COV_99, 1);
assertElementsAlmostEqual(qi.COV1, 0);

function test_totalDisplayDoesNotChangeCanonicalQi
[cst, pln] = helper_qiFixture(78, 39, 10);
ct.numOfCtScen = 1;
stf = struct();
resultGUI.physicalDose = 2 * ones(10, 1);

resultGUI = matRad_planAnalysis(resultGUI, ct, cst, stf, pln, ...
                                'showDVH', false, 'showQI', false, 'evaluationMode', 'total');

assertElementsAlmostEqual(resultGUI.qi.COV1, 1);
assertElementsAlmostEqual(resultGUI.qi.referenceDose, 2);
assertEqual(resultGUI.evaluationModeBase, 'perFraction');
assertEqual(resultGUI.evaluationMode, 'total');
assertElementsAlmostEqual(resultGUI.displayQi.referenceDose, 78);
assertTrue(isfield(resultGUI.qi, 'V_2Gy'));
assertFalse(isfield(resultGUI.qi, 'V_40Gy'));
assertElementsAlmostEqual(resultGUI.displayDvh(1).doseGrid, ...
                          resultGUI.dvh(1).doseGrid * pln.numOfFractions);

function test_defaultReferenceDosesIgnoreNanOutsideVoi
[cst, pln] = helper_qiFixture(78, 39, 2);
doseCube = [1; 2; NaN];

qi = matRad_calcQualityIndicators(cst, pln, doseCube, [], []);

assertTrue(isfield(qi, 'V_2Gy'));
assertFalse(isfield(qi, 'V_NaNGy'));

function test_analysisCodeDoesNotUseDeprecatedQuantityProperties
matRad_cfg = MatRad_Config.instance();
files = { ...
         fullfile(matRad_cfg.matRadRoot, 'matRad', 'matRad_planAnalysis.m'), ...
         fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', 'matRad_sampling.m'), ...
         fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', 'samplingAnalysis', 'matRad_samplingAnalysis.m'), ...
         fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', 'samplingAnalysis', 'matRad_latexReport.m'), ...
         fullfile(matRad_cfg.matRadRoot, 'matRad', 'doseCalc', '+ScenarioBatch', ...
                  '+Quantity', 'matRad_getDefaultScenarioDoseQuantity.m')};

for i = 1:numel(files)
    fileText = fileread(files{i});
    deprecatedNeedles = {['bioModel' '.quantityOpt'], ...
                         ['bioModel' '.quantityVis'], ...
                         ['bioParam' '.quantityOpt'], ...
                         ['bioParam' '.quantityVis']};
    for j = 1:numel(deprecatedNeedles)
        assertTrue(isempty(strfind(fileText, deprecatedNeedles{j})));
    end
end

function [cst, pln] = helper_qiFixture(prescriptionDose, numOfFractions, numOfVoxels)
pln = struct();
pln.numOfFractions = numOfFractions;
pln.bioModel = 'none';
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.quantityVis = 'physicalDose';

objective = DoseObjectives.matRad_SquaredDeviation(1, prescriptionDose);

cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'TARGET';
cst{1, 3} = 'TARGET';
cst{1, 4}{1} = 1:numOfVoxels;
cst{1, 5} = struct('Visible', true, 'Priority', 1);
cst{1, 6}{1} = objective;
