%% Example: 4D INTERVAL3 robust treatment planning with photons
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
% (ii)  calculate photon dose influence on all CT phases
% (iii) optimize a nominal plan and an INTERVAL3 robust plan
% (iv)  sample both plans and compare DVH trustbands

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

[ct, cst] = matRad_addMovement(ct, cst, motionPeriod, numOfCtScen, amplitude, 'visBool', true);
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
resultGUI = matRad_fluenceOptimization(dij, cst, pln);

%% INTERVAL3 dose interval precomputation
cstInterval3 = cst;
cstInterval3{ixPTV, 6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation(50, prescriptionDose));

intervalCfg = struct();
intervalCfg.IntervalMode = 'INTERVAL3';
intervalCfg.Quantity = quantityOpt;
intervalCfg.refScen = ct.refScen;
intervalCfg.ProgressLevel = 'detailed';
intervalCfg.targetStructSel = cstInterval3(ixPTV, 2);
intervalCfg.OARStructSel = cstInterval3(ixOAR, 2);

plnIntervalInput = pln;
plnIntervalInput.propOpt.scen4D = 'all';
[plnInterval3, dijIntervalContext] = matRad_calcDoseInterval( ...
                                                             ct, cstInterval3, stf, plnIntervalInput, dij, intervalCfg);

%% INTERVAL3 robust optimization
cstInterval3{ixPTV, 6}{1}.robustness = 'INTERVAL3';
cstInterval3{ixOAR, 6}{1}.robustness = 'INTERVAL3';

plnInterval3.propOpt.theta1 = 30;
plnInterval3.propOpt.theta2 = 1.0;
interval3Label = sprintf('INTERVAL3 (theta1=%g, theta2=%g)', ...
                         plnInterval3.propOpt.theta1, plnInterval3.propOpt.theta2);

resultGUIInterval3 = matRad_fluenceOptimization(dijIntervalContext, cstInterval3, plnInterval3);

scenarioIds = pln.multScen.scenarioIds();
for i = 1:numel(scenarioIds)
    dijScenarioIx = pln.multScen.getDijScenarioIndex(scenarioIds(i));
    resultGUIScenario = matRad_calcCubes(resultGUIInterval3.w, dij, dijScenarioIx);
    resultGUIInterval3 = matRad_appendResultGUI(resultGUIInterval3, ...
                                                resultGUIScenario, false, sprintf('scen%d', i));
end

resultGUI = matRad_appendResultGUI(resultGUI, resultGUIInterval3, false, 'interval3');

%% Calculate accumulated 4D dose
totalPhaseMatrix = ones(dij.totalNumOfBixels, ct.numOfCtScen) / ct.numOfCtScen;
totalPhaseMatrix = bsxfun(@times, totalPhaseMatrix, resultGUI.w);
[resultGUI4D, ~] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, totalPhaseMatrix);

totalPhaseMatrix = ones(dij.totalNumOfBixels, ct.numOfCtScen) / ct.numOfCtScen;
totalPhaseMatrix = bsxfun(@times, totalPhaseMatrix, resultGUIInterval3.w);
[resultGUIInterval4D, ~] = matRad_calc4dDose(ct, pln, dij, stf, cstInterval3, ...
                                             resultGUIInterval3, totalPhaseMatrix);

%% Visualize nominal and INTERVAL3 doses
plane = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1, :), ct);
slice = slice(3);
maxDose = max([max(max(resultGUI.(quantityOpt)(:, :, slice))) ...
               max(max(resultGUIInterval3.(quantityOpt)(:, :, slice))) ...
               max(max(resultGUI4D.(quantityOpt)(:, :, slice))) ...
               max(max(resultGUI4D.accPhysicalDose(:, :, slice))) ...
               max(max(resultGUIInterval4D.(quantityOpt)(:, :, slice))) ...
               max(max(resultGUIInterval4D.accPhysicalDose(:, :, slice)))]) + 1e-4;
doseIsoLevels = linspace(0.1 * maxDose, maxDose, 10);

figure;
subplot(1, 2, 1);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, ...
                 'dose', resultGUI.(quantityOpt), 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', [0 maxDose], ...
                 'doseIsoLevels', doseIsoLevels);
title('nominal plan');
subplot(1, 2, 2);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cstInterval3, 'cubeIdx', 1, ...
                 'dose', resultGUIInterval3.(quantityOpt), 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', [0 maxDose], ...
                 'doseIsoLevels', doseIsoLevels);
title(interval3Label);

figure;
subplot(1, 2, 1);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, ...
                 'dose', resultGUI4D.accPhysicalDose, 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', [0 maxDose], ...
                 'doseIsoLevels', doseIsoLevels);
title('nominal accumulated 4D dose');
subplot(1, 2, 2);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cstInterval3, 'cubeIdx', 1, ...
                 'dose', resultGUIInterval4D.accPhysicalDose, 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', [0 maxDose], ...
                 'doseIsoLevels', doseIsoLevels);
title([interval3Label ' accumulated 4D dose']);

%% Interactive scenario dose view
scenarioIx = 1;
numScenarios = pln.multScen.numScenarios();
maxScenarioDose = max(max(resultGUIInterval3.([quantityOpt '_scen' num2str(scenarioIx)])(:, :, slice))) + 0.2;
doseWindow = [0 maxScenarioDose];
doseIsoLevels = linspace(0.1 * maxScenarioDose, maxScenarioDose, 10);
scenarioFigure = figure;
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cstInterval3, 'cubeIdx', 1, ...
                 'dose', resultGUIInterval3.([quantityOpt '_scen' num2str(scenarioIx)]), ...
                 'plane', plane, 'slice', slice, 'contourColorMap', colorcube, ...
                 'doseWindow', doseWindow, 'doseIsoLevels', doseIsoLevels);
title([interval3Label ' scenario ' num2str(scenarioIx)]);

[env, envver] = matRad_getEnvironment();
if numScenarios > 1 && (strcmp(env, 'MATLAB') || str2double(envver(1)) >= 5)
    sliderStep = 1 / (numScenarios - 1);
    scenarioSlider = uicontrol('Parent', scenarioFigure, 'Style', 'slider', ...
                               'Position', [50 5 419 23], 'value', scenarioIx, 'min', 1, ...
                               'max', numScenarios, 'SliderStep', [sliderStep sliderStep]);
    set(scenarioSlider, 'Callback', @(es, ed) ...
        matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cstInterval3, ...
                         'cubeIdx', round(get(es, 'Value')), ...
                         'dose', resultGUIInterval3.([quantityOpt '_scen' ...
                                                      num2str(round(get(es, 'Value')))]), ...
                         'plane', plane, 'slice', slice, 'contourColorMap', colorcube, ...
                         'doseWindow', doseWindow, 'doseIsoLevels', doseIsoLevels));
end

%% Sampling and DVH trustbands
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

[caSamp, mSampDose, plnSamp, resultGUINomScen] = ...
    matRad_sampling(ct, stf, cst, pln, resultGUI.w, structSel, samplingMultScen, ...
                    dvhSamplingArgs{:});
[~, resultGUISamp, meta] = ...
    matRad_samplingAnalysis(ct, cst, plnSamp, caSamp, mSampDose, resultGUINomScen, ...
                            'plane', plane);

[caSampInterval, mSampDoseInterval, plnSampInterval, resultGUINomScen] = ...
    matRad_sampling(ct, stf, cstInterval3, pln, resultGUIInterval3.w, ...
                    structSel, samplingMultScen, dvhSamplingArgs{:});
[~, resultGUISampInterval, metaInterval] = ...
    matRad_samplingAnalysis(ct, cstInterval3, plnSampInterval, caSampInterval, ...
                            mSampDoseInterval, resultGUINomScen, 'plane', plane);

if ~isempty(resultGUISamp.robustnessAnalysis.index1.robustnessIndex)
    fprintf('Nominal robustness index1 pass fraction: %.4f\n', ...
            resultGUISamp.robustnessAnalysis.index1.robustnessIndex);
end
if ~isempty(resultGUISampInterval.robustnessAnalysis.index1.robustnessIndex)
    fprintf('%s robustness index1 pass fraction: %.4f\n', interval3Label, ...
            resultGUISampInterval.robustnessAnalysis.index1.robustnessIndex);
end
fprintf('Nominal expected dose difference status: %s, max expected over/under dose difference: %.4f / %.4f\n', ...
        resultGUISamp.expectedDoseDifferenceAnalysis.status, ...
        resultGUISamp.expectedDoseDifferenceAnalysis.summary.maxOverReferenceExpectedDoseDifference, ...
        resultGUISamp.expectedDoseDifferenceAnalysis.summary.maxUnderReferenceExpectedDoseDifference);
fprintf('%s expected dose difference status: %s, max expected over/under dose difference: %.4f / %.4f\n', ...
        interval3Label, resultGUISampInterval.expectedDoseDifferenceAnalysis.status, ...
        resultGUISampInterval.expectedDoseDifferenceAnalysis.summary.maxOverReferenceExpectedDoseDifference, ...
        resultGUISampInterval.expectedDoseDifferenceAnalysis.summary.maxUnderReferenceExpectedDoseDifference);

figure('Name', 'Robustness index 1 comparison');
subplot(1, 2, 1);
matRad_plotSamplingRobustnessAnalysis(resultGUISamp.robustnessAnalysis, ct, cst, slice, ...
                                      'axesHandle', gca, 'method', 'index1', ...
                                      'plane', plane, 'contourColorMap', colorcube);
title('nominal robustness index 1');
subplot(1, 2, 2);
matRad_plotSamplingRobustnessAnalysis(resultGUISampInterval.robustnessAnalysis, ...
                                      ct, cstInterval3, slice, ...
                                      'axesHandle', gca, 'method', 'index1', ...
                                      'plane', plane, 'contourColorMap', colorcube);
title([interval3Label ' robustness index 1']);

figure('Name', 'Expected dose difference comparison');
subplot(1, 2, 1);
matRad_plotExpectedDoseDifferenceAnalysis(resultGUISamp.expectedDoseDifferenceAnalysis, ...
                                          ct, cst, slice, 'axesHandle', gca, ...
                                          'plane', plane, 'contourColorMap', colorcube);
title('nominal E[D - D_{nom}]');
subplot(1, 2, 2);
matRad_plotExpectedDoseDifferenceAnalysis(resultGUISampInterval.expectedDoseDifferenceAnalysis, ...
                                          ct, cstInterval3, slice, 'axesHandle', gca, ...
                                          'plane', plane, 'contourColorMap', colorcube);
title([interval3Label ' E[D - D_{nom}]']);

figure;
subplot(1, 2, 1);
matRad_showDVHFromSampling(caSamp, dvhScale, cst, plnSamp, 1:numel(caSamp), ...
                           dvhDoseWindow, 'trustband', 1, 1, 'scenWeights', meta.scenWeights);
title('nominal DVH trustband');
ylim([0 105]);
subplot(1, 2, 2);
matRad_showDVHFromSampling(caSampInterval, dvhScale, cstInterval3, plnSampInterval, ...
                           1:numel(caSampInterval), dvhDoseWindow, 'trustband', 1, 1, ...
                           'scenWeights', metaInterval.scenWeights);
title([interval3Label ' DVH trustband']);
ylim([0 105]);

stdDoseValues = [resultGUISamp.stdCube(:); resultGUISampInterval.stdCube(:)];
stdDoseValues = stdDoseValues(isfinite(stdDoseValues));
if isempty(stdDoseValues) || max(stdDoseValues) <= 0
    stdDoseWindow = [0 1];
else
    stdDoseWindow = [0 max(stdDoseValues)];
end

figure;
subplot(1, 2, 1);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, ...
                 'dose', resultGUISamp.stdCube, 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', stdDoseWindow, ...
                 'colorBarLabel', stdColorBarLabel);
title('nominal sampling std');
subplot(1, 2, 2);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cstInterval3, 'cubeIdx', 1, ...
                 'dose', resultGUISampInterval.stdCube, 'plane', plane, 'slice', slice, ...
                 'contourColorMap', colorcube, 'doseWindow', stdDoseWindow, ...
                 'colorBarLabel', stdColorBarLabel);
title([interval3Label ' sampling std']);
