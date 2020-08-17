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

%% In this plan we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize beamlet intensities
% (iv) how to visually and quantitatively evaluate the result

%% set matRad runtime configuration
matRad_rc

%% Patient Data Import
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
cst{1,6}(1).dose        = 30;
cst{1,6}(1).penalty     = 400;
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
pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve 9 worst case scenarios for dose calculation and optimziation
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

%% Export dij matrix
matRad_exportDij('dij.bin',dij,stf);

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
