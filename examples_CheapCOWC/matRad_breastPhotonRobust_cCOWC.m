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

%% Set output folder
description_folder = 'breast';
run_config.robustness = 'c-COWC';
run_config.mode = 'wcScen';
run_config.sampling_mode = 'impScen';
run_config.p = 1;
run_config.beta = run_config.p/13;

output_folder = ['output' filesep description_folder filesep run_config.robustness filesep num2str(run_config.beta) filesep run_config.mode filesep datestr(datetime)];

%Set up parent export folder and full file path
if ~(isfolder(output_folder))
    mkdir(matRad_cfg.matRadRoot, output_folder);
end

folderPath = [matRad_cfg.matRadRoot filesep output_folder];

%%
diary([folderPath filesep 'diary.log']) 
diary on

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

savefig([folderPath filesep 'ct.fig']);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs

% Body
cst{1,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
cst{1,6}{1}.robustness  = 'none';
 
% Contralateral Lung
cst{2,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(80,50,5));
cst{2,6}{1}.robustness  = 'none';
 
% Ipsilateral Lung
cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(400,20,20));
cst{3,6}{1}.robustness  = 'none';
%cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));


% Heart
cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(250,4));
cst{4,6}{1}.robustness  = 'none';

% Contralateral Breast
cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,5,0));
cst{5,6}{1}.robustness  = 'none';

% CTV
p=42.56;
ixCTV = 6;
cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
cst{ixCTV,6}{1}.robustness  = 'none';
cst{ixCTV,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));

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

%% retrieve scenarios for dose calculation and optimziation
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

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%matRadGUI;

%% Obtain dose statistics
[dvh,dqi] = matRad_indicatorWrapper(cst,pln,resultGUI);
savefig([folderPath filesep 'dvh_nominal.fig']);

%%
% retrieve 13 scenarios for dose calculation and optimziation
 pln_robust=pln;
if(run_config.mode=="wcScen")
    multScen = matRad_multScen(ct,'wcScen'); 
    multScen.wcFactor=1.5;
    multScen.shiftSD = [4 6 8];
end
%%
% retrieve 13 scenarios for dose calculation and optimziation
if(run_config.mode=="impScen")
    multScen = matRad_multScen(ct,'impScen'); 
    multScen.wcFactor=1.5;
    multScen.shiftSD = [4 6 8];
    multScen.numOfShiftScen = [4 4 4];
    multScen.numOfRangeShiftScen=4;
    multScen.includeNomScen=true;
end
%%
pln_robust.multScen=multScen;

%% Dose calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
now1 = tic();
dij_robust = matRad_calcPhotonDose(ct,stf,pln_robust,cst);
DCTime_robust = toc(now1);
time1=sprintf('DCTime_robust: %.2f\n',DCTime_robust); disp(time1);

%% Trigger robust optimization
% Make the objective to a composite worst case objective

% CTV
cst{ixCTV,6}{1}.robustness  = run_config.robustness;
cst{ixCTV,8}{1}.beta = run_config.beta;
cst{ixCTV,8}{1}.p = run_config.p;

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
now2 = tic();
resultGUI_robust = matRad_fluenceOptimization(dij_robust,cst,pln_robust);
OPTTime_robust = toc(now2);
time2=sprintf('OPTTime_robust: %.2f\n',OPTTime_robust); disp(time2);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice

plane      = 3;
slice      = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI_robust.physicalDose(:)*pln_robust.numOfFractions])];
figure
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_robust.physicalDose*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,[],[],'Dose [Gy]');
savefig([folderPath filesep 'dose_robust.fig']);

%% Obtain dose statistics
[dvh_robust,dqi_robust] = matRad_indicatorWrapper(cst,pln_robust,resultGUI_robust);
savefig([folderPath filesep 'dvh_robust.fig']);

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {}; % structSel = {'PTV','OAR1'};

if(run_config.sampling_mode=="rndScen")
    multScen = matRad_multScen(ct,'rndScen'); % 'impSamp' or 'wcSamp'
    multScen.numOfShiftScen = matRad_cfg.defaults.samplingScenarios * ones(3,1);
    multScen.shiftSD = [4 6 8];
    multScen.numOfRangeShiftScen = matRad_cfg.defaults.samplingScenarios;
end

if(run_config.sampling_mode=="impScen")
    multScen = matRad_multScen(ct,'impScen'); 
    multScen.wcFactor=1.5;
    multScen.shiftSD = [4 6 8];
    multScen.numOfShiftScen = [8 8 8];
    multScen.numOfRangeShiftScen=8;
    multScen.includeNomScen=true;
end

%% Perform sampling
[caSamp, mSampDose, plnSamp, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln_robust,resultGUI_robust.w,structSel,multScen);

%% Perform sampling analysis
varargin.GammaCriterion = [3 3]; % [%  mm]
[cstStat, resultGUISamp, meta] = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen);

%% Multi-scenario dose volume histogram (DVH)
figure,set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVH_sampledScen(caSamp,dvh_robust,cst,plnSamp,[1:25]);
savefig([folderPath filesep 'dvh_robust_multiscen.fig']);

%% Dose volume histogram (DVH)
resultGUISamp_ul=[];
resultGUISamp_ul.physicalDose=resultGUI_robust.physicalDose;
resultGUISamp_ul.physicalDose_lower=resultGUISamp.meanCube-resultGUISamp.stdCube;
resultGUISamp_ul.physicalDose_upper=resultGUISamp.meanCube+resultGUISamp.stdCube;
[dvh_sampled,dqi_sampled] = matRad_indicatorWrapper_sampled(cst,pln,resultGUISamp_ul,[20,50]/pln.numOfFractions,[2,5,10,20,30,40,50,60,70,80,90,95,98]);
savefig([folderPath filesep 'dvh_robust_trustband.fig']);

%% STD dose based on sampling
figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.stdCube*pln.numOfFractions,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.stdCube(:)*pln.numOfFractions)],[],[],'Dose uncertainty [Gy]');
savefig([folderPath filesep 'uncertainty.fig']);

%% Uncertainty volume histogram (UVH)
resultGUISamp.physicalDose=resultGUISamp.stdCube;
[uvh,uqi] = matRad_indicatorWrapper_UVH(cst,pln,resultGUI,resultGUISamp);
savefig([folderPath filesep 'uvh.fig']);

%% Gamma index based on sampling
figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.gammaAnalysis.gammaCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.gammaAnalysis.gammaCube(:))],[],[],'Gamma index');
savefig([folderPath filesep 'gamma.fig']);

%% Gamma index based on sampling
figure;
resultGUISamp.gammaAnalysis.robustnessViolationCube = (resultGUISamp.gammaAnalysis.gammaCube>power(1,1/2));
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.gammaAnalysis.robustnessViolationCube,plane,slice,[],[],colorcube,[],[0 max(resultGUISamp.gammaAnalysis.gammaCube(:))],[],[],'Gamma index');
savefig([folderPath filesep 'robustness_violation.fig']);

%% Print evaluation indexes
% CTV
disp('CTV evaluation');
DMean=sprintf('DMean: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).mean*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).mean*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).mean*pln.numOfFractions); disp(DMean);
D98=sprintf('D98: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).D_98*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_98*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_98*pln.numOfFractions); disp(D98);
D50=sprintf('D50: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).D_50*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_50*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_50*pln.numOfFractions); disp(D50);
D2=sprintf('D2: %.2f [%.2f - %.2f] Gy \n',dqi_sampled{1,1}(ixCTV).D_2*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_2*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_2*pln.numOfFractions); disp(D2);

%%
% robustness indexes

disp('Robustness evaluation (Method 1)');

U2=sprintf('U2: %.2f Gy',uqi(ixCTV).D_2*pln.numOfFractions); disp(U2);
U50=sprintf('U50: %.2f Gy',uqi(ixCTV).D_50*pln.numOfFractions); disp(U50);
U95=sprintf('U95: %.2f Gy',uqi(ixCTV).D_95*pln.numOfFractions); disp(U95);
U98=sprintf('U98: %.2f Gy',uqi(ixCTV).D_95*pln.numOfFractions); disp(U98);
UI=sprintf('UI (>=0): %.2f \n',uqi(ixCTV).D_2.*pln.numOfFractions/p); disp(UI);

disp('Robustness evaluation (Method 2)');

UncertaintyIndex=nnz(resultGUISamp.gammaAnalysis.robustnessViolationCube(cst{ixCTV,4}{1,1}))/numel(cst{ixCTV,4}{1,1});

UI2=sprintf('UI (0-1): %.2f',UncertaintyIndex); disp(UI2);

RobustnessIndex=1-UncertaintyIndex;
RI2=sprintf('RI (0-1): %.2f \n',RobustnessIndex); disp(RI2);

h1=histogram(resultGUISamp.gammaAnalysis.gammaCube(cst{6,4}{1,1}));
savefig([folderPath filesep 'gamma_histo.fig']);

%%
% Evaluating nominal and robust solutions

w_nominal=resultGUI.w;
f_nominal = matRad_calcObjectiveFunction(w_nominal,dij,cst,pln);

w_robust=resultGUI_robust.w;
f_robust = matRad_calcObjectiveFunction(w_robust,dij,cst,pln);

%%
% robustness price indexes

disp('Robustness price evaluation (Method 1)');
BodyTotal=dqi(1).mean*numel(cst{1,4}{1,1});
CTVTotal=dqi(ixCTV).mean*numel(cst{ixCTV,4}{1,1});
OARMean_nominal=(BodyTotal-CTVTotal)/(numel(cst{1,4}{1,1})-numel(cst{ixCTV,4}{1,1}))*pln.numOfFractions;

BodyTotal_robust=dqi_sampled{1,1}(1).mean*numel(cst{1,4}{1,1});
CTVTotal_robust=dqi_sampled{1,1}(ixCTV).mean*numel(cst{ixCTV,4}{1,1});
OARMean_robust=(BodyTotal_robust-CTVTotal_robust)/(numel(cst{1,4}{1,1})-numel(cst{ixCTV,4}{1,1}))*pln.numOfFractions;

RPI=sprintf('RPI: %.2f\n',(OARMean_robust-OARMean_nominal)/OARMean_nominal); disp(RPI);

disp('Robustness price evaluation (Method 2)');

F_nominal=sprintf('F(nominal): %.2f',f_nominal); disp(F_nominal); % 1.6149103490694448e+00
F_robust=sprintf('F(robust): %.2f \n',f_robust); disp(F_robust);

RPI2=sprintf('RPI: %.2f \n',(f_robust-f_nominal)/f_nominal); disp(RPI2);

%% Print evaluation indexes
% Contralateral Lung
disp('Contralateral lung evaluation');
V50ConLung=sprintf('V50: %.2f [%.2f - %.2f] %%',dqi_sampled{1,1}(2).V_3_12Gy*100,dqi_sampled{1,2}(2).V_3_12Gy*100,dqi_sampled{1,3}(2).V_3_12Gy*100); disp(V50ConLung);
D5ConLung=sprintf('D5: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(2).D_5*pln.numOfFractions,dqi_sampled{1,2}(2).D_5*pln.numOfFractions,dqi_sampled{1,3}(2).D_5*pln.numOfFractions); disp(D5ConLung);

%% Print evaluation indexes
% Ipsilateral Lung
disp('Ipsilateral lung evaluation');
V20IpsLung=sprintf('V20: %.2f [%.2f - %.2f] %%',dqi_sampled{1,1}(3).V_1_25Gy*100,dqi_sampled{1,2}(3).V_1_25Gy*100,dqi_sampled{1,3}(3).V_1_25Gy*100); disp(V20IpsLung);
D20IpsLung=sprintf('D20: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(3).D_20*pln.numOfFractions,dqi_sampled{1,2}(3).D_20*pln.numOfFractions,dqi_sampled{1,3}(3).D_20*pln.numOfFractions); disp(D20IpsLung);

%% Print evaluation indexes
% Heart
disp('Heart evaluation');
DMeanHeart=sprintf('DMean: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(4).mean*pln.numOfFractions,dqi_sampled{1,2}(4).mean*pln.numOfFractions,dqi_sampled{1,3}(4).mean*pln.numOfFractions); disp(DMeanHeart);

%% Print evaluation indexes
% Contralateral Breast
disp('Contralateral Breast evaluation');
DMaxConBreast=sprintf('DMax: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(5).max*pln.numOfFractions,dqi_sampled{1,2}(5).max*pln.numOfFractions,dqi_sampled{1,3}(5).max*pln.numOfFractions); disp(DMaxConBreast);

%%
diary off