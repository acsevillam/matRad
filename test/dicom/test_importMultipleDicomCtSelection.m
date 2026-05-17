function test_suite = test_importMultipleDicomCtSelection

test_functions = localfunctions();

initTestSuite;

function test_singlePatientSelectsAllCtSeries
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P1', 'CTB'; ...
                                 'RTSTRUCT', 'P1', 'RT1'});

[selectedPatientIDs, ctSeriesUIDs, rtssRowsByScenario] = ...
    matRad_resolveMultipleDicomCtImportSelection(fileList, {'P1'}, struct());

assertEqual(selectedPatientIDs, {'P1'});
assertEqual(ctSeriesUIDs, {'CTA'; 'CTB'});
assertEqual(rtssRowsByScenario{1}, fileList(3, :));
assertEqual(rtssRowsByScenario{2}, fileList(3, :));

function test_multiplePatientsRequireSelection
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P2', 'CTB'});

assertExceptionThrown(@() matRad_resolveMultipleDicomCtImportSelection( ...
                                                                       fileList, {'P1', 'P2'}, struct()), 'matRad:Error');

function test_explicitCtSeriesCanSpanPatients
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P2', 'CTB'; ...
                                 'CT', 'P2', 'CTC'});
metadata = struct('ctSeriesUIDs', {{'CTA', 'CTB'}});

[selectedPatientIDs, ctSeriesUIDs, rtssRowsByScenario] = ...
    matRad_resolveMultipleDicomCtImportSelection(fileList, {'P1', 'P2'}, metadata);

assertEqual(selectedPatientIDs, {});
assertEqual(ctSeriesUIDs, {'CTA'; 'CTB'});
assertTrue(isempty(rtssRowsByScenario{1}));
assertTrue(isempty(rtssRowsByScenario{2}));

function test_singleRtssUidIsReusedForAllCtScenarios
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P1', 'CTB'; ...
                                 'RTSTRUCT', 'P1', 'RT1'; ...
                                 'RTSTRUCT', 'P1', 'RT2'});
metadata = struct('rtssUID', 'RT1');

[~, ~, rtssRowsByScenario] = ...
    matRad_resolveMultipleDicomCtImportSelection(fileList, {'P1'}, metadata);

assertEqual(rtssRowsByScenario{1}, fileList(3, :));
assertEqual(rtssRowsByScenario{2}, fileList(3, :));

function test_rtssUidsMatchCtScenariosOneToOne
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P1', 'CTB'; ...
                                 'RTSTRUCT', 'P1', 'RT1'; ...
                                 'RTSTRUCT', 'P1', 'RT2'});
metadata = struct('rtssUIDs', {{'RT1', 'RT2'}});

[~, ~, rtssRowsByScenario] = ...
    matRad_resolveMultipleDicomCtImportSelection(fileList, {'P1'}, metadata);

assertEqual(rtssRowsByScenario{1}, fileList(3, :));
assertEqual(rtssRowsByScenario{2}, fileList(4, :));

function test_rtssUidsMustMatchNumberOfCtScenarios
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'CT', 'P1', 'CTB'; ...
                                 'RTSTRUCT', 'P1', 'RT1'});
metadata = struct('rtssUIDs', {{'RT1'}});

assertExceptionThrown(@() matRad_resolveMultipleDicomCtImportSelection( ...
                                                                       fileList, {'P1'}, metadata), 'matRad:Error');

function test_ambiguousRtstructRequiresExplicitSelection
fileList = helper_dicomFileList({'CT', 'P1', 'CTA'; ...
                                 'RTSTRUCT', 'P1', 'RT1'; ...
                                 'RTSTRUCT', 'P1', 'RT2'});

assertExceptionThrown(@() matRad_resolveMultipleDicomCtImportSelection( ...
                                                                       fileList, {'P1'}, struct()), 'matRad:Error');

function test_cstMetadataCompatibilityAcceptsMatchingScenarios
cstScenarios = {helper_cstFixture('PTV', 'TARGET'), helper_cstFixture('PTV', 'TARGET')};

assertTrue(matRad_hasCompatibleCstMetadata(cstScenarios));

function test_cstMetadataCompatibilityRejectsMismatchedScenarios
cstScenarios = {helper_cstFixture('PTV', 'TARGET'), helper_cstFixture('CTV', 'TARGET')};

assertFalse(matRad_hasCompatibleCstMetadata(cstScenarios));

function test_finalizeImportCropsCommonOverlapAndMergesScenarios
tmpCtOriginal = {helper_ctFixture(1:3), helper_ctFixture(2:4)};
tmpCstOriginal = {helper_cstWithIndices([sub2ind([3 3 2], 2, 1, 1); sub2ind([3 3 2], 2, 2, 1)]), ...
                  helper_cstWithIndices([sub2ind([3 3 2], 2, 1, 1); sub2ind([3 3 2], 2, 3, 1)])};
ctIn.numOfCtScen = 2;

[ct, cst] = matRad_finalizeMultipleDicomCtImport(ctIn, tmpCtOriginal, tmpCstOriginal);

expectedIndex = sub2ind([3 2 2], 2, 1, 1);
assertEqual(ct.numOfCtScen, 2);
assertEqual(ct.cubeDim, [3 2 2]);
assertEqual(ct.x, [2 3]);
assertEqual(numel(ct.cubeHU), 2);
assertEqual(size(ct.cubeHU{1}), [3 2 2]);
assertEqual(size(ct.cubeHU{2}), [3 2 2]);
assertEqual(cst{1, 4}{1}, expectedIndex);
assertEqual(cst{1, 4}{2}, expectedIndex);

function fileList = helper_dicomFileList(rows)
numRows = size(rows, 1);
fileList = cell(numRows, 5);
for rowIx = 1:numRows
    fileList{rowIx, 1} = ['file', num2str(rowIx), '.dcm'];
    fileList{rowIx, 2} = rows{rowIx, 1};
    fileList{rowIx, 3} = rows{rowIx, 2};
    fileList{rowIx, 4} = rows{rowIx, 3};
    fileList{rowIx, 5} = num2str(rowIx);
end

function cst = helper_cstFixture(name, type)
cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = name;
cst{1, 3} = type;
cst{1, 4}{1} = [1 2];
cst{1, 5} = struct('Visible', true, 'Priority', 1);
cst{1, 6} = {};

function ct = helper_ctFixture(xAxis)
ct.x = xAxis;
ct.y = 1:3;
ct.z = 1:2;
ct.cubeHU = {reshape(1:18, [3 3 2])};
ct.cubeDim = [3 3 2];
ct.dicomInfo = struct('SeriesInstanceUID', ['CT', num2str(xAxis(1))]);
ct.dicomMeta = struct('SeriesInstanceUID', ['CT', num2str(xAxis(1))]);

function cst = helper_cstWithIndices(indices)
cst = helper_cstFixture('PTV', 'TARGET');
cst{1, 4}{1} = indices(:);
