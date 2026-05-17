function test_suite = test_resolveStructureSelection

test_functions = localfunctions();

initTestSuite;

function test_allTargets
cst = helper_selectionFixture();

ix = matRad_resolveStructureSelection(cst, 'all', [], 'TARGET');

assertEqual(ix, [1 2]);

function test_allOars
cst = helper_selectionFixture();

ix = matRad_resolveStructureSelection(cst, 'all', [], 'OAR');

assertEqual(ix, 3);

function test_allTargetsAndOars
cst = helper_selectionFixture();

ix = matRad_resolveStructureSelection(cst, 'all', [], {'TARGET', 'OAR'});

assertEqual(ix, [1 2 3]);

function test_includeByNameAndIndex
cst = helper_selectionFixture();

ix = matRad_resolveStructureSelection(cst, 'include', {'PTV', 2}, 'TARGET');

assertEqual(ix, [1 2]);

function test_excludeByName
cst = helper_selectionFixture();

ix = matRad_resolveStructureSelection(cst, 'exclude', {'PTV'}, 'TARGET');

assertEqual(ix, 2);

function test_rejectsTypeNotAllowed
cst = helper_selectionFixture();

assertExceptionThrown(@() matRad_resolveStructureSelection(cst, 'include', {'RECTUM'}, 'TARGET'), ...
                      'matRad:Error');

function test_rejectsAmbiguousName
cst = helper_selectionFixture();
cst{4, 1} = 4;
cst{4, 2} = 'PTV';
cst{4, 3} = 'TARGET';
cst{4, 4}{1} = 4;
cst{4, 5} = struct();
cst{4, 6} = {};

assertExceptionThrown(@() matRad_resolveStructureSelection(cst, 'include', {'PTV'}, 'TARGET'), ...
                      'matRad:Error');

function test_rejectsEmptyIncludeList
cst = helper_selectionFixture();

assertExceptionThrown(@() matRad_resolveStructureSelection(cst, 'include', [], 'TARGET'), ...
                      'matRad:Error');

function test_rejectsInvalidMode
cst = helper_selectionFixture();

assertExceptionThrown(@() matRad_resolveStructureSelection(cst, 'selected', [], 'TARGET'), ...
                      'matRad:Error');

function cst = helper_selectionFixture()
cst = cell(3, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', 1);
cst = helper_addStructure(cst, 2, 'CTV', 'TARGET', 2);
cst = helper_addStructure(cst, 3, 'RECTUM', 'OAR', 3);

function cst = helper_addStructure(cst, ix, name, type, voxels)
cst{ix, 1} = ix;
cst{ix, 2} = name;
cst{ix, 3} = type;
cst{ix, 4}{1} = voxels;
cst{ix, 5} = struct();
cst{ix, 6} = {};
