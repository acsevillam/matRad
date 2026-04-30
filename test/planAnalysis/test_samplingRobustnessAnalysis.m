function test_suite = test_samplingRobustnessAnalysis

test_functions = localfunctions();

initTestSuite;

function test_default_uses_all_targets
    [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture();

    robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,[]);

    assertEqual({robustnessAnalysis.targets.name},{'PTV','CTV'});
    assertElementsAlmostEqual(robustnessAnalysis.index1.robPassRate,50);
    assertEqual(robustnessAnalysis.structureSelection.selectedCstIndices,[1 2]);

function test_include_uses_only_requested_target
    [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture();

    robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,[], ...
        'robustnessTargetMode','include','robustnessTargets',{'PTV'});

    assertEqual({robustnessAnalysis.targets.name},{'PTV'});
    assertElementsAlmostEqual(robustnessAnalysis.index1.robPassRate,100);
    assertEqual(robustnessAnalysis.structureSelection.selectedNames,{'PTV'});

function test_exclude_omits_requested_target
    [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture();

    robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,[], ...
        'robustnessTargetMode','exclude','robustnessTargets',{'CTV'});

    assertEqual({robustnessAnalysis.targets.name},{'PTV'});
    assertElementsAlmostEqual(robustnessAnalysis.index1.robPassRate,100);
    assertEqual(robustnessAnalysis.structureSelection.selectedCstIndices,1);

function test_rejects_oar_in_robustness_target_selection
    [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture();

    assertExceptionThrown(@() matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,[], ...
        'robustnessTargetMode','include','robustnessTargets',{'RECTUM'}), ...
        'matRad:Error');

function test_partial_sample_mask_skips_uncovered_targets
    [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture();
    sampleMask = false(ct.cubeDim);
    sampleMask([1 2]) = true;

    robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,[], ...
        'sampleMask',sampleMask);

    assertTrue(robustnessAnalysis.targets(1).isEvaluable);
    assertFalse(robustnessAnalysis.targets(2).isEvaluable);
    assertEqual(robustnessAnalysis.targets(2).samplingStatus,'partialSamplingCoverage');
    assertEqual(robustnessAnalysis.targets(2).numUnsampledVoxels,2);
    assertElementsAlmostEqual(robustnessAnalysis.index1.robPassRate,100);

function [meanCube,stdCube,criteria,ct,cst,pln] = robustnessFixture()
    ct = struct();
    ct.cubeDim = [5 1 1];
    ct.refScen = 1;

    meanCube = zeros(ct.cubeDim);
    meanCube([1 2]) = 2;
    meanCube([3 4]) = 10;
    stdCube = zeros(ct.cubeDim);
    criteria = [5 5];

    pln = struct();
    pln.numOfFractions = 1;
    pln.bioParam = struct('model','none');

    cst = cell(3,6);
    cst = addStructure(cst,1,'PTV','TARGET',[1 2],2);
    cst = addStructure(cst,2,'CTV','TARGET',[3 4],2);
    cst = addStructure(cst,3,'RECTUM','OAR',5,[]);

function cst = addStructure(cst,ix,name,type,voxels,prescriptionDose)
    cst{ix,1} = ix;
    cst{ix,2} = name;
    cst{ix,3} = type;
    cst{ix,4}{1} = voxels;
    cst{ix,5} = struct();
    if isempty(prescriptionDose)
        cst{ix,6} = {};
    else
        cst{ix,6}{1} = DoseObjectives.matRad_SquaredDeviation(1,prescriptionDose);
    end
