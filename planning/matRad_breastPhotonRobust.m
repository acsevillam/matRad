%% Example: Photon Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% % terms in the LICENSE file found in the top-level directory of this
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

%%
clear;
clc;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the
% matRad root directory with all its subdirectories is added to the Matlab
% search path.

load('patient3_5mm.mat');

%% plot CT slice
if param.logLevel == 1
    
    figure('Renderer', 'painters', 'Position', [10 10 300*ct.numOfCtScen 400]);
    
    isocenter = matRad_getIsoCenter(cst,ct,0);
    
    for scen_iterator = 1:ct.numOfCtScen
        plane      = 1;
        slice      = round(isocenter(2)./ct.resolution.y);
        subplot(2,ct.numOfCtScen,scen_iterator); camroll(90);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
        plane      = 3;
        slice      = round(isocenter(3)./ct.resolution.z);
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_geoSliceWrapper(gca,ct,cst,scen_iterator,plane,slice,[],[],colorcube,[],[]);
        
    end
    
end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = round(isocenter(3)./ct.resolution.z);
    numScen    = 1;
    
    matRad_geoSliceWrapper(gca,ct,cst,numScen,plane,slice,[],[],colorcube,[],[]);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_geoSliceWrapper(gca,ct,cst,round(es.Value),plane,slice,[],[],colorcube,[],[]);
end
clear  numScen plane slice ans f b;

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target

% define optimization parameter for both VOIs

% Body
cst{1,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));
cst{1,6}{1}.robustness  = 'none';
 
% Contralateral Lung
cst{2,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(80,50,5));
cst{2,6}{1}.robustness  = 'none';
 
% Ipsilateral Lung
cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(400,20,0));
cst{3,6}{1}.robustness  = 'none';

% Heart
cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(250,4));
cst{4,6}{1}.robustness  = 'none';

% Contralateral Breast
cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{5,6}{1} = struct(DoseObjectives.matRad_MeanDose(80,8));
cst{5,6}{1}.robustness  = 'none';

% CTV
ixCTV = 6;
cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,42.56));
cst{ixCTV,6}{1}.robustness  = 'none';

display(cst{ixCTV,6});
%%
% The file TG119.mat contains two Matlab variables. Let's check what we
% have just imported. First, the 'ct' variable comprises the ct cube along
%with some meta information describing properties of the ct cube (cube
% dimensions, resolution, number of CT scenarios). Please note that
%multiple ct cubes (e.g. 4D CT) can be stored in the cell array ct.cube{}
display(ct);
%%
% The 'cst' cell array defines volumes of interests along with information
% required for optimization. Each row belongs to one certain volume of
% interest (VOI), whereas each column defines different properties.
% Specifically, the second and third column  show the name and the type of
% the structure. The type can be set to OAR, TARGET or IGNORED. The fourth
% column contains a linear index vector that lists all voxels belonging to
% a certain VOI.
display(cst);

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This
% matlab structure requires input from the treatment planner and defines
% the most important cornerstones of your treatment plan.

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

%%
% Obtain the number of beams and voxels from the existing variables and
% calculate the iso-center which is per default the center of gravity of
% all target voxels.
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%% dose calculation settings
% set resolution of dose calculation and optimization
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 0;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

%%
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

%%
% and et voila our treatment plan structure is ready. Lets have a look:
display(pln);

%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geometry information of the 6th beam
display(stf(6));

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Export dij matrix
matRad_exportDij('dij0.bin',dij,stf);

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%matRadGUI;

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice

plane      = 3;
slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI.physicalDose(:)])];
figure
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);

%% Obtain dose statistics
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

% Print evaluation indexes
D95 = sprintf('D95: %.4f Gy',qi(4).mean*16); disp(D95);

%% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'rndScen');
pln.multScen.numOfShiftScen = [10 10 10];
pln.multScen.shiftSD = [4 6 8];
pln.multScen.shiftGenType = 'sampled';
pln.multScen.shiftCombType = 'combined';
pln.multScen.numOfRangeShiftScen = 0;
pln.multScen.includeNomScen = true;

%% Export structures voxels indexes
matRad_exportStructures(cst);

%% Export beam positions
matRad_exportBeamPositions(pln,stf)

%% Dose calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Export dij matrix
for shiftScen = 1:pln.multScen.totNumShiftScen
    metadata.numScen=shiftScen;
    dij_filename=append('dij',num2str(shiftScen),'.bin');
    matRad_exportDij(dij_filename,dij,stf,metadata);
end

%% Export scen probabilities 
matRad_exportScenProb('scenProb.txt',pln);

%% Plot dose per voxel histos 
wOnes = ones(dij.totalNumOfBixels,1);
matRad_plotDoseHistos(dij,pln,cst,ixCTV,wOnes);

%% Dose interval calculation 
[dij_interval, pln_interval] = matRad_calcDoseInterval(dij,pln,cst{ixCTV,4}{1,1},20,80);

%% Export interval dij matrix
metadata.numScen=1;
matRad_exportDij('dij_min20.bin',dij_interval,stf,metadata);

metadata.numScen=2;
matRad_exportDij('dij_max80.bin',dij_interval,stf,metadata);

%% Calculate interval dose vector
wOnes = ones(dij.totalNumOfBixels,1);
% Min interval dose
tmpResultGUI1 = matRad_calcCubes(wOnes,dij_interval,1);
% Max interval dose
tmpResultGUI2 = matRad_calcCubecst{ixCTV,4}(wOnes,dij_interval,2);
% Center interval dose
tmpResultGUIcenter=[];
tmpResultGUIcenter.physicalDose=(tmpResultGUI2.physicalDose+tmpResultGUI1.physicalDose)/2;

%% Plot the resulting dose slice
% Let's plot the transversal iso-center dose slice
plane      = 3;
slice      = round(pln_interval.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([tmpResultGUI2.physicalDose(:)])];

% Min interval dose 
figure
matRad_plotSliceWrapper(gca,ct,cst,1,tmpResultGUI1.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
% Center interval dose
figure
matRad_plotSliceWrapper(gca,ct,cst,1,tmpResultGUIcenter.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
% Max interval dose
figure
matRad_plotSliceWrapper(gca,ct,cst,1,tmpResultGUI2.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);

%% Trigger robust optimization
% Make the objective to a composite worst case objective

% Body
cst{1,6}{1}.robustness  = 'INTERVAL';

% CTV
cst{ixCTV,6}{1}.robustness  = 'INTERVAL';

%Activate 4D Optimization
pln.propOpt.scen4D = 'all';

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI_robust = matRad_fluenceOptimization(dij_interval,cst,pln_interval);

%% add resultGUIrobust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUI_robust,0,'robust');
%matRadGUI;

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice

plane      = 3;
slice      = round(pln_interval.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI_robust.physicalDose(:)])];
figure
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_robust.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);

%% Obtain dose statistics
[dvh,qi] = matRad_indicatorWrapper(cst,pln_interval,resultGUI_robust);

% Print evaluation indexes
D95 = sprintf('D95: %.4f Gy',qi(4).mean*16); disp(D95);
