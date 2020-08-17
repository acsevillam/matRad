%% Example: Photon Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize beamlet intensities
% (iv) how to visually and quantitatively evaluate the result

%% set matRad runtime configuration
matRad_rc

%% Import patient data multi-scenario 
load('images/patient_2/patient_2_scen_1_resized.mat');

%%
% The file TG119.mat contains two Matlab variables. Let's check what we 
% have just imported. First, the 'ct' variable comprises the ct cube along
%with some meta information describing properties of the ct cube (cube 
% dimensions, resolution, number of CT scenarios). Please note that 
%multiple ct cubes (e.g. 4D CT) can be stored in the cell array ct.cube{}
if param.logLevel == 1
    display(ct);
end

%%
% The 'cst' cell array defines volumes of interests along with information 
% required for optimization. Each row belongs to one certain volume of 
% interest (VOI), whereas each column defines different properties. 
% Specifically, the second and third column  show the name and the type of 
% the structure. The type can be set to OAR, TARGET or IGNORED. The fourth 
% column contains a linear index vector that lists all voxels belonging to 
% a certain VOI.
if param.logLevel == 1
    display(cst);
end

%% plot CT slice
if param.logLevel == 1
   
figure('Renderer', 'painters', 'Position', [10 10 300*ct.numOfCtScen 400]);

    for scen_iterator = 1:ct.numOfCtScen
        plane      = 1;
        slice      = 128;
        subplot(2,ct.numOfCtScen,scen_iterator); camroll(90);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
        plane      = 3;
        slice      = 30;
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
    end

end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = 128;
    numScen    = 1;

    matRad_geoSliceWrapper(gca,ct,cst,numScen,plane,slice,[],[],colorcube,[],[]);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
          'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_geoSliceWrapper(gca,ct,cst,round(es.Value),plane,slice,[],[],colorcube,[],[]);
end
clear  numScen plane slice ans f b;

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% The sixth column contains optimization information such as objectives and
% constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for
% individual structures. Here, we have defined a squared deviation 
% objective making it 'expensive/costly' for the optimizer to over- and 
% underdose the target structure (both are equally important). 

% Body
cst{1,6}(1).Priority    = 3;
cst{1,6}(1).type        = 'square overdosing';
cst{1,6}(1).dose        = 20;
cst{1,6}(1).penalty     = 200;
cst{1,6}(1).EUD         = NaN;
cst{1,6}(1).volume      = NaN;
cst{1,6}(1).robustness  = 'none';
 
% % Contralateral Lung
cst{2,5}.priority       = 2;
cst{2,6}(1).type        = 'max DVH objective';
cst{2,6}(1).dose        = 50;
cst{2,6}(1).penalty     = 80;
cst{2,6}(1).EUD         = NaN;
cst{2,6}(1).volume      = 5;
cst{2,6}(1).robustness  = 'none';
 
% Ipsilateral Lung 
cst{3,5}.priority       = 2;
cst{3,6}(1).type        = 'max DVH objective';
cst{3,6}(1).dose        = 20;
cst{3,6}(1).penalty     = 400;
cst{3,6}(1).EUD         = NaN;
cst{3,6}(1).volume      = 0;
cst{3,6}(1).robustness  = 'none';

% % Heart
cst{4,5}.priority       = 2;
cst{4,6}(1).type        = 'mean';
cst{4,6}(1).dose        = 4;
cst{4,6}(1).penalty     = 250;
cst{4,6}(1).EUD         = NaN;
cst{4,6}(1).volume      = NaN;
cst{4,6}(1).robustness  = 'none';

% Contralateral Breast
cst{5,5}.priority       = 2;
cst{5,6}(1).type        = 'mean';
cst{5,6}(1).dose        = 8;
cst{5,6}(1).penalty     = 80;
cst{5,6}(1).EUD         = NaN;
cst{5,6}(1).volume      = 1;
cst{5,6}(1).robustness  = 'none';

ixTarget = 6;

% PTV
cst{ixTarget,5}.Priority    = 1;
cst{ixTarget,6}.type        = 'square deviation';
cst{ixTarget,6}.dose        = 42.56;
cst{ixTarget,6}.penalty     = 800;
cst{ixTarget,6}.robustness  = 'none';

display(cst{ixTarget,6});

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this case
% we want to use photons. Then, we need to define a treatment machine to 
% correctly load the corresponding base data. matRad includes base data for
% generic photon linear accelerator called 'Generic'. By this means matRad 
% will look for 'photons_Generic.mat' in our root directory and will use 
% the data provided in there for dose calculation
pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

%%
% Define the flavor of optimization along with the quantity that should be
% used for optimization. Possible quantities used for optimization are: 
% physicalDose: physical dose based optimization; 
% effect: biological effect based optimization;
% RBExD: RBE weighted dose based optimzation;
% Possible biological models are:
% none:        use no specific biological model
% constRBE:    use a constant RBE
% MCN:         use the variable RBE McNamara model for protons
% WED:         use the variable RBE Wedenberg model for protons
% LEM:         use the biophysical variable RBE Local Effect model for carbons
% As we are  using photons, we simply set the parameter to 'physicalDose' and
% and 'none'
quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

%%
% Now we have to set some beam parameters. We can define multiple beam 
% angles for the treatment and pass these to the plan as a vector. matRad 
% will then interpret the vector as multiple beams. In this case, we define
% linear spaced beams from 0 degree to 359 degree in 40 degree steps. This 
% results in 9 beams. All corresponding couch angles are set to 0 at this 
% point. Moreover, we set the bixelWidth to 5, which results in a beamlet 
% size of 5 x 5 mm in the isocenter plane. The number of fractions is set 
% to 30. Internally, matRad considers the fraction dose for optimization, 
% however, objetives and constraints are defined for the entire treatment.
pln.numOfFractions         = 16;
pln.propStf.gantryAngles   = [15 45 75 105 135 315 345];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% define nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

%%
% Obtain the number of beams and voxels from the existing variables and 
% calculate the iso-center which is per default the center of gravity of 
% all target voxels.
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 0;
pln.propOpt.runDAO        = 0;

%%
% and et voila our treatment plan structure is ready. Lets have a look:
if param.logLevel == 1
    display(pln);
end

%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln,param);

%%
% Let's display the beam geometry information of the 6th beam
if param.logLevel == 1
    display(stf(1));
end

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst,param);

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil 
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to 
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln,param);
%matRadGUI;

%% Plot the results
% Let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet);
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

%% Import patient data multi-scenario 
load('images/patient_2/patient_2_multiScen_contourned.mat');

%% plot CT slice
if param.logLevel == 1
   
figure('Renderer', 'painters', 'Position', [10 10 300*ct.numOfCtScen 400]);

    for scen_iterator = 1:ct.numOfCtScen
        plane      = 1;
        slice      = 128;
        subplot(2,ct.numOfCtScen,scen_iterator); camroll(90);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
        plane      = 3;
        slice      = 30;
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
    end

end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = 128;
    numScen    = 1;

    matRad_geoSliceWrapper(gca,ct,cst,numScen,plane,slice,[],[],colorcube,[],[]);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
          'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_geoSliceWrapper(gca,ct,cst,round(es.Value),plane,slice,[],[],colorcube,[],[]);
end

%%
% The sixth column contains optimization information such as objectives and
% constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for
% individual structures. Here, we have defined a squared deviation 
% objective making it 'expensive/costly' for the optimizer to over- and 
% underdose the target structure (both are equally important). 

% Body
cst{1,6}(1).Priority    = 3;
cst{1,6}(1).type        = 'square overdosing';
cst{1,6}(1).dose        = 20;
cst{1,6}(1).penalty     = 200;
cst{1,6}(1).EUD         = NaN;
cst{1,6}(1).volume      = NaN;
cst{1,6}(1).robustness  = 'none';
 
% % Contralateral Lung
cst{2,5}.priority       = 2;
cst{2,6}(1).type        = 'max DVH objective';
cst{2,6}(1).dose        = 50;
cst{2,6}(1).penalty     = 80;
cst{2,6}(1).EUD         = NaN;
cst{2,6}(1).volume      = 5;
cst{2,6}(1).robustness  = 'none';
 
% Ipsilateral Lung 
cst{3,5}.priority       = 2;
cst{3,6}(1).type        = 'max DVH objective';
cst{3,6}(1).dose        = 20;
cst{3,6}(1).penalty     = 400;
cst{3,6}(1).EUD         = NaN;
cst{3,6}(1).volume      = 0;
cst{3,6}(1).robustness  = 'none';

% % Heart
cst{4,5}.priority       = 2;
cst{4,6}(1).type        = 'mean';
cst{4,6}(1).dose        = 4;
cst{4,6}(1).penalty     = 250;
cst{4,6}(1).EUD         = NaN;
cst{4,6}(1).volume      = NaN;
cst{4,6}(1).robustness  = 'none';

% Contralateral Breast
cst{5,5}.priority       = 2;
cst{5,6}(1).type        = 'mean';
cst{5,6}(1).dose        = 8;
cst{5,6}(1).penalty     = 80;
cst{5,6}(1).EUD         = NaN;
cst{5,6}(1).volume      = 1;
cst{5,6}(1).robustness  = 'none';

ixTarget = 6;

% PTV
cst{ixTarget,5}.Priority    = 1;
cst{ixTarget,6}.type        = 'square deviation';
cst{ixTarget,6}.dose        = 42.56;
cst{ixTarget,6}.penalty     = 800;
cst{ixTarget,6}.robustness  = 'none';

display(cst{ixTarget,6});

%%
% retrieve 9 worst case scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'wcScen');

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst,param);

%% Trigger robust optimization
% Make the objective to a composite worst case objective
cst{1,6}(1).robustness  = 'COWC';
cst{2,6}(1).robustness  = 'COWC';
cst{3,6}(1).robustness  = 'COWC';
cst{4,6}(1).robustness  = 'COWC';
cst{5,6}(1).robustness  = 'COWC';
cst{ixTarget,6}.robustness  = 'COWC';

%% Inverse Optimization for IMRT
resultGUIrobust = matRad_fluenceOptimization(dij,cst,pln,param);

% add resultGUIrobust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUIrobust,0,'robust');
%matRadGUI;

%% Plot results
  
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose       = max([max(resultGUI.([quantityOpt])(:,:,slice)) max(resultGUIrobust.([quantityOpt])(:,:,slice))])+1e-4;
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
figure,
subplot(121),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.([quantityOpt])      ,plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('conventional plan')
subplot(122),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);title('robust plan')

[~,~] = matRad_indicatorWrapper(cst,pln,resultGUI);
[~,~] = matRad_indicatorWrapper(cst,pln,resultGUIrobust);

%% Multiple scenario treatment evaluation

numScen=size(dij.([quantityOpt]));

resultGUI2 = cell(numScen(1),numScen(2),numScen(3));
resultGUIrobust2 = cell(numScen(1),numScen(2),numScen(3));

for ctScen=1:numScen(1)
    for shiftScen=1:numScen(2)
        for shiftRangeScen=1:numScen(3)
            if(isempty(dij.physicalDose{ctScen,shiftScen,shiftRangeScen})==false)
                metadata.ctScen=ctScen;
                metadata.shiftScen=shiftScen;
                metadata.shiftRangeScen=1;
                resultGUI2{ctScen,shiftScen,shiftRangeScen} = matRad_calcCubes(resultGUI.wUnsequenced,dij,metadata);
                resultGUIrobust2{ctScen,shiftScen,shiftRangeScen} = matRad_calcCubes(resultGUIrobust.wUnsequenced,dij,metadata);
            end
        end
    end
end

%% create an interactive plot to slide through individual ct scenarios
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose = 0;
for ctScen = 1:ct.numOfCtScen
    if (max(max(resultGUI2{ctScen,1,1}.([quantityOpt])(:,:,slice)))+0.2>maxDose)
        maxDose = max(max(resultGUI2{ctScen,1,1}.([quantityOpt])(:,:,slice)))+0.2;
    end
end

f      = figure; title('individual ct scenarios');
ctScen = 1;

doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI2{ctScen,1,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);
b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
      'value',ctScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(es.Value),resultGUI2{round(es.Value),1,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);

%% create an interactive plot to slide through individual shift scenarios
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose = 0;
for shiftScen = 1:pln.multScen.totNumShiftScen
    if (max(max(resultGUI2{1,shiftScen,1}.([quantityOpt])(:,:,slice)))+0.2>maxDose)
        maxDose = max(max(resultGUI2{1,shiftScen,1}.([quantityOpt])(:,:,slice)))+0.2;
    end
end

f      = figure; title('individual shift scenarios');
shiftScen = 1;

doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI2{1,shiftScen,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);
b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
      'value',1, 'min',1, 'max',pln.multScen.totNumShiftScen,'SliderStep', [1/(pln.multScen.totNumShiftScen-1) , 1/(pln.multScen.totNumShiftScen-1)]);
b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI2{1,round(es.Value),1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);

%% create an interactive plot to slide through individual ct scenarios
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose = 0;
for ctScen = 1:ct.numOfCtScen
    if (max(max(resultGUIrobust2{ctScen,1,1}.([quantityOpt])(:,:,slice)))+0.2>maxDose)
        maxDose = max(max(resultGUIrobust2{ctScen,1,1}.([quantityOpt])(:,:,slice)))+0.2;
    end
end

f      = figure; title('individual ct scenarios');
ctScen = 1;

doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust2{ctScen,1,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);
b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
      'value',ctScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(es.Value),resultGUIrobust2{round(es.Value),1,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);

%% create an interactive plot to slide through individual shift scenarios
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
maxDose = 0;
for shiftScen = 1:pln.multScen.totNumShiftScen
    if (max(max(resultGUIrobust2{1,shiftScen,1}.([quantityOpt])(:,:,slice)))+0.2>maxDose)
        maxDose = max(max(resultGUIrobust2{1,shiftScen,1}.([quantityOpt])(:,:,slice)))+0.2;
    end
end

f      = figure; title('individual shift scenarios');
shiftScen = 1;

doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust2{1,shiftScen,1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);
b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
      'value',1, 'min',1, 'max',pln.multScen.totNumShiftScen,'SliderStep', [1/(pln.multScen.totNumShiftScen-1) , 1/(pln.multScen.totNumShiftScen-1)]);
b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrobust2{1,round(es.Value),1}.([quantityOpt]),plane,slice,[],[],colorcube,[],[0 maxDose],doseIsoLevels);

%% Plot results

matRad_multScenDVHWrapper(cst,pln,resultGUI2);
matRad_multScenDVHWrapper(cst,pln,resultGUIrobust2);

%% Perform sampling
% select structures to include in sampling; leave empty to sample dose for all structures
   
% sampling does not know on which scenario sampling should be performed
structSel = {}; % structSel = {'PTV','OAR1'};
param.logLevel = 2;
[caSamp, mSampDose, plnSamp, resultGUInomScen]          = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,[],[]);
[cstStat, resultGUISamp, param]                         = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen,[]);
   
[caSampRob, mSampDoseRob, plnSampRob, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln,resultGUIrobust.w,structSel,[],[]);
[cstStatRob, resultGUISampRob, paramRob]                = matRad_samplingAnalysis(ct,cst,plnSampRob,caSampRob, mSampDoseRob, resultGUInomScen,[]);
 
   
%% Plot conventional and robust planning stdDev
   
figure,title('std dose cube based on sampling - conventional')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.stdCube(:))],[]);
   
figure,title('std dose cube based on sampling - robust')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.stdCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISampRob.stdCube(:))],[]);
