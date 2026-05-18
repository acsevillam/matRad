%% Example: 4D photon robust model comparison
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will
% (i)   create a small artificial 4D photon phantom
% (ii)  optimize nominal, COWC, c-COWC, PROB2, INTERVAL2, and INTERVAL3 plans
% (iii) sample all plans on the same 4D scenario model
% (iv)  compare DVH trustbands and sampled dose standard deviation

%% Set matRad runtime configuration
matRad_rc;

%% Create an artificial CT image series
xDim = 150;
yDim = 150;
zDim = 50;

ct.cubeDim = [xDim yDim zDim];
ct.resolution.x = 2; % mm
ct.resolution.y = 2; % mm
ct.resolution.z = 3; % mm
ct.numOfCtScen = 1;
ct.cubeHU{1} = ones(ct.cubeDim) * -1024;

%% Create the VOI data for the phantom
ixOAR = 1;
ixPTV = 2;

cst{ixOAR, 1} = 0;
cst{ixOAR, 2} = 'contour';
cst{ixOAR, 3} = 'OAR';
cst{ixPTV, 1} = 1;
cst{ixPTV, 2} = 'target';
cst{ixPTV, 3} = 'TARGET';

cst{ixOAR, 5}.TissueClass = 1;
cst{ixOAR, 5}.alphaX = 0.1000;
cst{ixOAR, 5}.betaX = 0.0500;
cst{ixOAR, 5}.Priority = 2;
cst{ixOAR, 5}.Visible = 1;
cst{ixOAR, 6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10, 30));

cst{ixPTV, 5}.TissueClass = 1;
cst{ixPTV, 5}.alphaX = 0.1000;
cst{ixPTV, 5}.betaX = 0.0500;
cst{ixPTV, 5}.Priority = 1;
cst{ixPTV, 5}.Visible = 1;
cst{ixPTV, 6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(50, 60));

%% Create a cylindrical OAR and spherical PTV
cubeHelper = zeros(ct.cubeDim);
centerXOAR = round(xDim / 2);
centerYOAR = round(yDim / 2);
radiusOAR = xDim / 6;
zLowOAR = round(zDim / 2 - zDim / 4);
zHighOAR = round(zDim / 2 + zDim / 4);

for x = 1:xDim
    for y = 1:yDim
        if (x - centerXOAR)^2 + (y - centerYOAR)^2 < radiusOAR^2
            for z = zLowOAR:zHighOAR
                cubeHelper(x, y, z) = 1;
            end
        end
    end
end
cst{ixOAR, 4}{1} = find(cubeHelper);

cubeHelper = zeros(ct.cubeDim);
radiusPTV = xDim / 14;
for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPos = [x y z] - round(ct.cubeDim ./ 2);
            if sqrt(sum(currPos.^2)) < radiusPTV
                cubeHelper(x, y, z) = 1;
            end
        end
    end
end
cst{ixPTV, 4}{1} = find(cubeHelper);

vIxOAR = cst{ixOAR, 4}{1};
vIxPTV = cst{ixPTV, 4}{1};
ct.cubeHU{1}(vIxOAR) = 300;
ct.cubeHU{1}(vIxPTV) = 0;

%% Add motion and create the 4D CT
amplitude = [5 0 0]; % [voxels]
numOfCtScen = 5;
motionPeriod = 2.5; % [s]

[ct, cst] = matRad_addMovement(ct, cst, motionPeriod, numOfCtScen, amplitude, 'visBool', false);
ct.refScen = 1;
ct.dvfMetadata.dvfUnits = 'mm';
ct.dvfMetadata.refScen = ct.refScen;
ct.dvfMetadata.referenceCtScen = ct.refScen;

clear x y z xDim yDim zDim centerXOAR centerYOAR radiusOAR zHighOAR zLowOAR vIxOAR vIxPTV cubeHelper currPos radiusPTV;

%% Treatment plan
pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.bioModel = 'none';

quantityOpt = 'physicalDose';
prescriptionDose = 60;

pln.numOfFractions = 20;
pln.propStf.gantryAngles = 0:60:300;
pln.propStf.couchAngles = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth = 5;
pln.propStf.numOfBeams = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter = ones(pln.propStf.numOfBeams, 1) * matRad_getIsoCenter(cst, ct, 0);
pln.propOpt.quantityOpt = quantityOpt;
pln.propOpt.quantityVis = quantityOpt;
pln.propOpt.runDAO = 0;
pln.propOpt.runSequencing = 0;

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

pln.multScen = matRad_createScenarioModel(ct, 'nomScen');

%% Generate beam geometry and dose influence
stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% Nominal optimization
plnNominal = pln;
cstNominal = cst;
resultGUINominal = matRad_fluenceOptimization(dij, cstNominal, plnNominal);

%% COWC optimization
plnCOWC = pln;
plnCOWC.propOpt.scen4D = 'all';
cstCOWC = cst;
cstCOWC{ixPTV, 6}{1}.robustness = 'COWC';
cstCOWC{ixOAR, 6}{1}.robustness = 'COWC';
resultGUICOWC = matRad_fluenceOptimization(dij, cstCOWC, plnCOWC);

%% c-COWC optimization
plnCCOWC = pln;
plnCCOWC.propOpt.scen4D = 'all';
plnCCOWC.propOpt.p1 = 1;
plnCCOWC.propOpt.p2 = min(3, pln.multScen.numScenarios());
ccowcLabel = sprintf('c-COWC (p1=%d, p2=%d)', ...
                     plnCCOWC.propOpt.p1, plnCCOWC.propOpt.p2);
cstCCOWC = cst;
cstCCOWC{ixPTV, 6}{1}.robustness = 'c-COWC';
cstCCOWC{ixOAR, 6}{1}.robustness = 'c-COWC';
resultGUICCOWC = matRad_fluenceOptimization(dij, cstCCOWC, plnCCOWC);

%% Common scenario-batch precomputation setup
scenarioBatchCfg = struct();
scenarioBatchCfg.Quantity = quantityOpt;
scenarioBatchCfg.refScen = ct.refScen;
scenarioBatchCfg.ProgressLevel = 'summary';
scenarioBatchCfg.targetStructSel = cst(ixPTV, 2);
scenarioBatchCfg.OARStructSel = cst(ixOAR, 2);

plnScenarioBatchInput = pln;
plnScenarioBatchInput.propOpt.scen4D = 'all';

%% PROB2 optimization
cstProb2 = cst;
cstProb2{ixPTV, 6}{1}.robustness = 'PROB2';
cstProb2{ixOAR, 6}{1}.robustness = 'PROB2';
meanVarianceWeight = 30;
meanVariancePenalty = cstProb2{ixPTV, 6}{1}.penalty * meanVarianceWeight;
cstProb2{ixPTV, 6}{end + 1} = struct(DoseObjectives.matRad_MeanVariance(meanVariancePenalty));
cstProb2{ixPTV, 6}{end}.robustness = 'PROB2';
prob2Label = sprintf('PROB2 (MeanVariance weight=%g)', meanVarianceWeight);
[plnProb2, dijProbContext] = matRad_calculateProbabilisticQuantities( ...
                                                                     ct, cstProb2, stf, plnScenarioBatchInput, dij, scenarioBatchCfg);
resultGUIProb2 = matRad_fluenceOptimization(dijProbContext, cstProb2, plnProb2);

%% INTERVAL2 optimization
cstInterval2 = cst;
cstInterval2{ixPTV, 6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation(50, prescriptionDose));
interval2Cfg = scenarioBatchCfg;
interval2Cfg.IntervalMode = 'INTERVAL2';
[plnInterval2, dijInterval2Context] = matRad_calcDoseInterval( ...
                                                              ct, cstInterval2, stf, plnScenarioBatchInput, dij, interval2Cfg);

plnInterval2.propOpt.theta1 = 30;
interval2Label = sprintf('INTERVAL2 (theta1=%g)', plnInterval2.propOpt.theta1);
cstInterval2{ixPTV, 6}{1}.robustness = 'INTERVAL2';
cstInterval2{ixOAR, 6}{1}.robustness = 'INTERVAL2';
resultGUIInterval2 = matRad_fluenceOptimization(dijInterval2Context, cstInterval2, plnInterval2);

%% INTERVAL3 optimization
cstInterval3 = cst;
cstInterval3{ixPTV, 6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation(50, prescriptionDose));
interval3Cfg = scenarioBatchCfg;
interval3Cfg.IntervalMode = 'INTERVAL3';
[plnInterval3, dijInterval3Context] = matRad_calcDoseInterval( ...
                                                              ct, cstInterval3, stf, plnScenarioBatchInput, dij, interval3Cfg);

plnInterval3.propOpt.theta1 = 30;
plnInterval3.propOpt.theta2 = 1.0;
interval3Label = sprintf('INTERVAL3 (theta1=%g, theta2=%g)', ...
                         plnInterval3.propOpt.theta1, plnInterval3.propOpt.theta2);
cstInterval3{ixPTV, 6}{1}.robustness = 'INTERVAL3';
cstInterval3{ixOAR, 6}{1}.robustness = 'INTERVAL3';
resultGUIInterval3 = matRad_fluenceOptimization(dijInterval3Context, cstInterval3, plnInterval3);

%% Sampling setup
plane = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1, :), ct);
slice = slice(3);
structSel = {};
evaluationMode = 'total';
[dvhScale, ~] = matRad_convertToEvaluationMode(1, pln, evaluationMode);
dvhDoseWindow = [0 1.6 * prescriptionDose];
dvhSamplingDoseWindow = matRad_convertFromEvaluationMode(dvhDoseWindow, pln, evaluationMode);
dvhSamplingArgs = {'dvhDoseWindow', dvhSamplingDoseWindow};
samplingMultScen = pln.multScen;
doseUnitLabel = 'Gy';
if strncmp(quantityOpt, 'RBExD', 5) || strncmp(quantityOpt, 'RBExDose', 8)
    doseUnitLabel = 'Gy(RBE)';
end
stdColorBarLabel = ['Standard deviation [' doseUnitLabel ']'];

planSamples(1).label = 'Nominal';
planSamples(1).cst = cstNominal;
planSamples(1).w = resultGUINominal.w;
planSamples(2).label = 'COWC';
planSamples(2).cst = cstCOWC;
planSamples(2).w = resultGUICOWC.w;
planSamples(3).label = ccowcLabel;
planSamples(3).cst = cstCCOWC;
planSamples(3).w = resultGUICCOWC.w;
planSamples(4).label = prob2Label;
planSamples(4).cst = cstProb2;
planSamples(4).w = resultGUIProb2.w;
planSamples(5).label = interval2Label;
planSamples(5).cst = cstInterval2;
planSamples(5).w = resultGUIInterval2.w;
planSamples(6).label = interval3Label;
planSamples(6).cst = cstInterval3;
planSamples(6).w = resultGUIInterval3.w;

%% Perform sampling
for i = 1:numel(planSamples)
    [planSamples(i).caSamp, planSamples(i).mSampDose, ...
        planSamples(i).plnSamp, planSamples(i).resultGUINomScen] = ...
        matRad_sampling(ct, stf, planSamples(i).cst, pln, planSamples(i).w, ...
                        structSel, samplingMultScen, dvhSamplingArgs{:});

    [planSamples(i).cstStat, planSamples(i).doseStat, planSamples(i).meta, ...
        planSamples(i).figures] = ...
        matRad_samplingAnalysis(ct, planSamples(i).cst, planSamples(i).plnSamp, ...
                                planSamples(i).caSamp, planSamples(i).mSampDose, ...
                                planSamples(i).resultGUINomScen, 'plane', plane);

    if ~isempty(planSamples(i).doseStat.robustnessAnalysis.index1.robustnessIndex)
        fprintf('%s robustness index1/index2 pass fractions: %.4f / %.4f\n', ...
                planSamples(i).label, ...
                planSamples(i).doseStat.robustnessAnalysis.index1.robustnessIndex, ...
                planSamples(i).doseStat.robustnessAnalysis.index2.robustnessIndex);
    end
    fprintf('%s gamma analysis status: %s\n', planSamples(i).label, ...
            planSamples(i).doseStat.gammaAnalysis.status);
    fprintf('%s expected dose difference status: %s, max expected over/under dose difference: %.4f / %.4f\n', ...
            planSamples(i).label, ...
            planSamples(i).doseStat.expectedDoseDifferenceAnalysis.status, ...
            planSamples(i).doseStat.expectedDoseDifferenceAnalysis.summary.maxOverReferenceExpectedDoseDifference, ...
            planSamples(i).doseStat.expectedDoseDifferenceAnalysis.summary.maxUnderReferenceExpectedDoseDifference);
end

%% Compare DVH trustbands
numPlotRows = 2;
numPlotCols = ceil(numel(planSamples) / numPlotRows);

figure('Name', 'DVH trustbands from 4D sampling');
for i = 1:numel(planSamples)
    subplot(numPlotRows, numPlotCols, i);
    matRad_showDVHFromSampling(planSamples(i).caSamp, dvhScale, planSamples(i).cst, ...
                               planSamples(i).plnSamp, 1:numel(planSamples(i).caSamp), ...
                               dvhDoseWindow, 'trustband', 1, 1, ...
                               'scenWeights', planSamples(i).meta.scenWeights);
    title(planSamples(i).label);
    ylim([0 105]);
end

%% Compare sampling robustness index
figure('Name', 'Robustness index 1 from 4D sampling');
for i = 1:numel(planSamples)
    subplot(numPlotRows, numPlotCols, i);
    matRad_plotSamplingRobustnessAnalysis(planSamples(i).doseStat.robustnessAnalysis, ...
                                          ct, planSamples(i).cst, slice, ...
                                          'axesHandle', gca, 'method', 'index1', ...
                                          'plane', plane, 'contourColorMap', colorcube);
    title([planSamples(i).label ' robustness index 1']);
end

%% Compare expected dose difference over/under nominal
figure('Name', 'Expected dose difference over/under nominal');
for i = 1:numel(planSamples)
    subplot(numPlotRows, numPlotCols, i);
    matRad_plotExpectedDoseDifferenceAnalysis( ...
                                              planSamples(i).doseStat.expectedDoseDifferenceAnalysis, ...
                                              ct, planSamples(i).cst, slice, ...
                                              'axesHandle', gca, 'plane', plane, ...
                                              'contourColorMap', colorcube);
    title([planSamples(i).label ' E[D - D_{nom}]']);
end

%% Compare sampled dose standard deviation
stdDoseValues = [];
for i = 1:numel(planSamples)
    stdDoseValues = [stdDoseValues; planSamples(i).doseStat.stdCube(:)];
end
stdDoseValues = stdDoseValues(isfinite(stdDoseValues));
if isempty(stdDoseValues) || max(stdDoseValues) <= 0
    stdDoseWindow = [0 1];
else
    stdDoseWindow = [0 max(stdDoseValues)];
end

figure('Name', 'Sampled dose standard deviation');
for i = 1:numel(planSamples)
    subplot(numPlotRows, numPlotCols, i);
    matRad_plotSlice(ct, 'axesHandle', gca, 'cst', planSamples(i).cst, ...
                     'cubeIdx', 1, 'dose', planSamples(i).doseStat.stdCube, ...
                     'plane', plane, 'slice', slice, 'contourColorMap', colorcube, ...
                     'doseWindow', stdDoseWindow, 'colorBarLabel', stdColorBarLabel);
    title(planSamples(i).label);
end
