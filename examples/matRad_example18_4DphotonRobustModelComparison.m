%% Example: 4D robust model comparison with photons
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
% (ii)  optimize nominal, COWC, c-COWC, and INTERVAL3 plans
% (iii) sample all plans on the same 4D CT scenario model
% (iv)  compare DVH trustbands and robustness index 1 in 2x2 figures

%% Set matRad runtime configuration
matRad_rc

%% Create an artificial CT image series
xDim = 150;
yDim = 150;
zDim = 50;

ct.cubeDim      = [xDim yDim zDim];
ct.resolution.x = 2; % mm
ct.resolution.y = 2; % mm
ct.resolution.z = 3; % mm
ct.numOfCtScen  = 1;
ct.cubeHU{1}    = ones(ct.cubeDim) * -1024;

%% Create the VOI data for the phantom
ixOAR = 1;
ixPTV = 2;

cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';
cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';

cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 2;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));

cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(50,60));

%% Create a cylindrical OAR and spherical PTV
cubeHelper = zeros(ct.cubeDim);
centerXOAR = round(xDim/2);
centerYOAR = round(yDim/2);
radiusOAR  = xDim/6;
zLowOAR    = round(zDim/2 - zDim/4);
zHighOAR   = round(zDim/2 + zDim/4);

for x = 1:xDim
   for y = 1:yDim
      if (x - centerXOAR)^2 + (y - centerYOAR)^2 < radiusOAR^2
         for z = zLowOAR:1:zHighOAR
            cubeHelper(x,y,z) = 1;
         end
      end
   end
end
cst{ixOAR,4}{1} = find(cubeHelper);

cubeHelper = zeros(ct.cubeDim);
radiusPTV = xDim/14;
for x = 1:xDim
   for y = 1:yDim
      for z = 1:zDim
         currPost = [x y z] - round(ct.cubeDim./2);
         if sqrt(sum(currPost.^2)) < radiusPTV
            cubeHelper(x,y,z) = 1;
         end
      end
   end
end
cst{ixPTV,4}{1} = find(cubeHelper);

vIxOAR = cst{ixOAR,4}{1};
vIxPTV = cst{ixPTV,4}{1};
ct.cubeHU{1}(vIxOAR) = 300;
ct.cubeHU{1}(vIxPTV) = 0;

%% Add motion and create the 4D CT
amplitude    = [5 0 0]; % [voxels]
numOfCtScen  = 5;
motionPeriod = 2.5; % [s]

[ct,cst] = matRad_addMovement(ct,cst,motionPeriod,numOfCtScen,amplitude,'visBool',true);
ct.refScen = 1;
ct.dvfMetadata.dvfUnits = 'mm';
ct.dvfMetadata.refScen = ct.refScen;
ct.dvfMetadata.referenceCtScen = ct.refScen;

clear x y z xDim yDim zDim centerXOAR centerYOAR radiusOAR zHighOAR zLowOAR vIxOAR vIxPTV cubeHelper currPost radiusPTV

%% Treatment plan
pln.radiationMode = 'photons';
pln.machine       = 'Generic';

modelName    = 'none';
quantityOpt  = 'physicalDose';

pln.numOfFractions        = 20;
pln.propStf.gantryAngles  = 0:60:300;
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
pln.multScen = matRad_createScenarioModel(ct,'nomScen');

%% Generate beam geometry and dose influence
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Nominal optimization
plnNominal = pln;
cstNominal = cst;
resultGUINominal = matRad_fluenceOptimization(dij,cstNominal,plnNominal);

%% COWC optimization
plnCOWC = pln;
plnCOWC.propOpt.scen4D = 'all';
cstCOWC = cst;
cstCOWC{ixPTV,6}{1}.robustness = 'COWC';
cstCOWC{ixOAR,6}{1}.robustness = 'COWC';
resultGUICOWC = matRad_fluenceOptimization(dij,cstCOWC,plnCOWC);

%% c-COWC optimization
plnCCOWC = pln;
plnCCOWC.propOpt.scen4D = 'all';
plnCCOWC.propOpt.p1 = 1;
plnCCOWC.propOpt.p2 = min(3,pln.multScen.totNumScen);
ccowcLabel = sprintf('c-COWC (p1=%d, p2=%d)', ...
    plnCCOWC.propOpt.p1,plnCCOWC.propOpt.p2);
cstCCOWC = cst;
cstCCOWC{ixPTV,6}{1}.robustness = 'c-COWC';
cstCCOWC{ixOAR,6}{1}.robustness = 'c-COWC';
resultGUICCOWC = matRad_fluenceOptimization(dij,cstCCOWC,plnCCOWC);

%% INTERVAL3 optimization
cstInterval3 = cst;
cstInterval3{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation(50,60));

intervalCfg = struct();
intervalCfg.Quantity = quantityOpt;
intervalCfg.refScen = ct.refScen;
intervalCfg.ProgressLevel = 'summary';
intervalCfg.targetStructSel = cstInterval3(ixPTV,2);
intervalCfg.OARStructSel = cstInterval3(ixOAR,2);
[plnInterval3,dijIntervalContext] = matRad_calcDoseInterval3(ct,cstInterval3,stf,pln,dij,intervalCfg);

plnInterval3.propOpt.scen4D = 'all';
plnInterval3.propOpt.theta1 = 30;
plnInterval3.propOpt.theta2 = 1.0;
interval3Label = sprintf('INTERVAL3 (theta1=%g, theta2=%g)', ...
    plnInterval3.propOpt.theta1,plnInterval3.propOpt.theta2);
cstInterval3{ixPTV,6}{1}.robustness = 'INTERVAL3';
cstInterval3{ixOAR,6}{1}.robustness = 'INTERVAL3';
resultGUIInterval3 = matRad_fluenceOptimization(dijIntervalContext,cstInterval3,plnInterval3);

%% Sampling setup
plane = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);

structSel = {};
targetDoseInfo = matRad_getTargetReferenceDoses(cstNominal,pln);
targetRefDose = [targetDoseInfo.refDose];
targetRefDose = targetRefDose(isfinite(targetRefDose));
dvhDoseWindow = [0 1.6*max(targetRefDose)];
dvhSamplingArgs = {'dvhDoseWindow',dvhDoseWindow};

robustnessCriteria = [5 5];
samplingAnalysisArgs = {'robustnessCriteria',robustnessCriteria};
samplingMultScen = pln.multScen;

planSamples(1).label = 'Nominal';
planSamples(1).cst = cstNominal;
planSamples(1).w = resultGUINominal.w;
planSamples(2).label = 'COWC';
planSamples(2).cst = cstCOWC;
planSamples(2).w = resultGUICOWC.w;
planSamples(3).label = ccowcLabel;
planSamples(3).cst = cstCCOWC;
planSamples(3).w = resultGUICCOWC.w;
planSamples(4).label = interval3Label;
planSamples(4).cst = cstInterval3;
planSamples(4).w = resultGUIInterval3.w;

%% Perform sampling
for i = 1:numel(planSamples)
    [planSamples(i).caSamp,planSamples(i).mSampDose, ...
        planSamples(i).plnSamp,planSamples(i).resultGUInomScen] = ...
        matRad_sampling(ct,stf,planSamples(i).cst,pln,planSamples(i).w, ...
        structSel,samplingMultScen,dvhSamplingArgs{:});

    [planSamples(i).cstStat,planSamples(i).doseStat,planSamples(i).meta] = ...
        matRad_samplingAnalysis(ct,planSamples(i).cst,planSamples(i).plnSamp, ...
        planSamples(i).caSamp,planSamples(i).mSampDose, ...
        planSamples(i).resultGUInomScen,samplingAnalysisArgs{:});
end

%% Compare DVH trustbands
figure('Name','DVH trustbands from 4D sampling')
for i = 1:numel(planSamples)
    subplot(2,2,i);
    matRad_showDVHFromSampling(planSamples(i).caSamp,[],planSamples(i).cst, ...
        planSamples(i).plnSamp,1:planSamples(i).plnSamp.multScen.totNumScen, ...
        dvhDoseWindow,'trustband',1,1,'scenWeights',planSamples(i).meta.scenWeights);
    title(planSamples(i).label);
    ylim([0 105]);
end

%% Compare robustness index 1
maxRob = 5.01;
robustnessWindow = [0 maxRob];
mMap1 = round(1/maxRob*256);
mMap2 = 256 - mMap1;
colormap1 = [linspace(0.40,1,mMap1)',linspace(1,1,mMap1)',linspace(0.40,1,mMap1)'];
colormap2 = matRad_getColormap('gammaIndex',2*mMap2);
robustnessColorMap = [colormap1; colormap2(mMap2+1:end-1,:)];

figure('Name','Robustness index 1 from 4D sampling')
for i = 1:numel(planSamples)
    subplot(2,2,i);
    robustnessAnalysis = planSamples(i).doseStat.robustnessAnalysis;

    if ~isempty(robustnessAnalysis.index1.robustnessCube)
        targetMaskCube = NaN(ct.cubeDim);
        targetMaskCube(robustnessAnalysis.evaluableTargetMask) = 1;
        robustnessCube = robustnessAnalysis.index1.robustnessCube .* targetMaskCube;
        [robustnessCube,plotColorMap] = matRad_prepareLimitedIndexPlot( ...
            robustnessCube,robustnessWindow,robustnessColorMap);
        matRad_plotSliceWrapper(gca,ct,planSamples(i).cst,1,robustnessCube, ...
            plane,slice,[],[],colorcube,plotColorMap,robustnessWindow,[],[], ...
            'Delta Index',false,'LineWidth',1.5);

        robustnessIndex = robustnessAnalysis.index1.robustnessIndex;
        if isempty(robustnessIndex) || ~isfinite(robustnessIndex)
            title([planSamples(i).label ' - RI n/a']);
        else
            title([planSamples(i).label ' - RI ' num2str(robustnessIndex,'%.3f')]);
        end
    else
        imagesc(NaN(ct.cubeDim(1),ct.cubeDim(2)));
        axis image;
        title([planSamples(i).label ' - RI n/a']);
    end
end
