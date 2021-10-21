function matRad_photonRobust_cCOWC2(description,plan_objectives,plan_target,plan_beams,shiftSD,robustness,beam_shapping_mode,scen_mode,wcFactor,rootPath,sampling,sampling_mode,sampling_wcFactor,p1,p2)
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

%% Clear variables
clearvars -except description plan_objectives plan_target plan_beams shiftSD robustness beam_shapping_mode scen_mode wcFactor rootPath sampling sampling_mode sampling_wcFactor p1 p2 ;
clc;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Set examples parameters
if ~exist('description','var') || isempty(description)
    run_config.description = 'prostate';
else
    run_config.description = description;
end
run_config.resolution = '5x5x5';

if ~exist('plan_objectives','var') || isempty(plan_objectives)
    run_config.plan_objectives = '1';
else
    run_config.plan_objectives = plan_objectives;
end

if ~exist('plan_target','var') || isempty(plan_target)
    run_config.plan_target = 'CTV';
else
    run_config.plan_target = plan_target;
end

if ~exist('plan_beams','var') || isempty(plan_beams)
    run_config.plan_beams = '9F';
else
    run_config.plan_beams = plan_beams;
end

if ~exist('shiftSD','var') || isempty(shiftSD)
    run_config.shiftSD = [5 5 5]; % mm
else
    run_config.shiftSD = shiftSD;
end

if ~exist('robustness','var') || isempty(robustness)
    run_config.robustness = 'COWC';
else
    run_config.robustness = robustness;
end

if ~exist('beam_shapping_mode','var') || isempty(beam_shapping_mode)
    run_config.beam_shapping_mode = 'wcScen';
else
    run_config.beam_shapping_mode = beam_shapping_mode;
end

if ~exist('scen_mode','var') || isempty(scen_mode)
    run_config.scen_mode = 'wcScen';
else
    run_config.scen_mode = scen_mode;
end

if ~exist('wcFactor','var') || isempty(wcFactor)
    run_config.wcFactor = 1.5;
else
    run_config.wcFactor = wcFactor;
end

if run_config.robustness == "c-COWC"
    
    switch run_config.scen_mode
        case "wcScen"
            run_config.numScens = 7;
        case "impScen"
            run_config.numScens = 13;
        case "impScen_permuted_truncated"
            run_config.numScens = 33;
    end
    
    if ~exist('p1','var') || isempty(p1)
        run_config.p1 = 1;
    else
        run_config.p1 = p1;
    end

    run_config.beta1 = run_config.p1/run_config.numScens;

    if ~exist('p2','var') || isempty(p2)
        run_config.p2 = 1;
    else
        run_config.p2 = p2;
    end

    run_config.beta2 = run_config.p2/run_config.numScens;
end

if ~exist('sampling','var') || isempty(sampling)
    run_config.sampling = false;
else
    run_config.sampling = sampling;
end

if ~exist('sampling_mode','var') || isempty(sampling_mode)
    run_config.sampling_mode = 'impScen_permuted_truncated';
else
    run_config.sampling_mode = sampling_mode;
end

if ~exist('sampling_wcFactor','var') || isempty(sampling_wcFactor)
    run_config.sampling_wcFactor = 2.0;
else
    run_config.sampling_wcFactor = sampling_wcFactor;
end

run_config.UncertaintyCriterion = 0.03;
run_config.GammaCriterion = [3 3];
run_config.sampling_size = 50;


if ~exist('rootPath','var') || isempty(rootPath)
    run_config.rootPath = matRad_cfg.matRadRoot;
else
    run_config.rootPath = rootPath;
end

if run_config.robustness == "c-COWC"
    output_folder = ['output' filesep run_config.description filesep run_config.robustness filesep num2str(run_config.beta1) '_to_' num2str(run_config.beta2) filesep run_config.scen_mode filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
else
    output_folder = ['output' filesep run_config.description filesep run_config.robustness filesep run_config.scen_mode filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
end

%Set up parent export folder and full file path
if ~(isfolder(output_folder))
    mkdir(run_config.rootPath, output_folder);
end

folderPath = [run_config.rootPath filesep output_folder];

%% Initiallize diary log
diary([folderPath filesep 'diary.log']) 
diary on

%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the
% matRad root directory with all its subdirectories is added to the Matlab
% search path.
[ct,cst] = matRad_loadGeometry(run_config);

%% Print run config
display(run_config);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
[cst,ixTarget,p,ixBody,ixCTV] = matRad_loadObjectives(run_config,'CTV',cst);
cst_nominal=cst;

%% Print target objectives
display(cst{ixTarget,2});
for i=1:length(cst{ixTarget,6})
    display(cst{ixTarget,6}{i});
end

%% Plot CT slice
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

%% Set plot and histograms window
run_config.doseWindow = [0 p*1.25];
run_config.doseWindow_dvh = [0 p*1.6];
run_config.doseWindow_uncertainty = [0 p*0.5];
run_config.doseWindow_relative_uncertainty1 = [0 100];
run_config.doseWindow_relative_uncertainty2 = [0 0.5];
run_config.doseWindow_uvh = [0 p*0.5];
run_config.gammaWindow = [0 1];

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

%% Load beams and couch setup
[pln] = matRad_loadBeams(run_config,pln,ct,cst);

%% Dose calculation settings
% set resolution of dose calculation and optimization
switch run_config.resolution
    case '3x3x3'
        pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
        pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
        pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
    case '5x5x5'
        pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
        pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
        pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
end

%% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 0;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

%% Generate dummy scenarios for beam shaping
[multScen] = matRad_multiScenGenerator(run_config.beam_shapping_mode,run_config,'optimization',ct);

%% Save multi-scenarios to plan
pln.multScen=multScen;

%%
% and et voila our treatment plan structure is ready. Lets have a look:
display(pln);

%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geometry information of the 6th beam
display(stf(1));

%% Retrieve nominal scenario for dose calculation and optimziation reference
pln.multScen = matRad_multScen(ct,'nomScen');

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

%% Evaluating nominal solution
w_nominal=resultGUI.w;
f_nominal = matRad_calcObjectiveFunction(w_nominal,dij,cst_nominal,pln);

%% Plot nominal fluence
matRad_visSpotWeights(stf,resultGUI.w);
savefig([folderPath filesep 'fluence_nominal.fig']);

%% Plot dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'dose3d_nominal.fig']);

%% Create target mask
target_mask = zeros(size(resultGUI.physicalDose));
target_mask(cst{ixCTV,4}{1,1}) = 1;

%% Create OAR mask
OAR_mask = zeros(size(resultGUI.physicalDose));
OAR_mask(cst{ixBody,4}{1,1}) = 1;
OAR_mask(cst{ixCTV,4}{1,1}) = 0;

%% Plot target dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose.*target_mask*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'target_dose3d_nominal.fig']);

%% Plot OAR dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose.*OAR_mask*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'OAR_dose3d_nominal.fig']);

%% Obtain dose statistics
[dvh,dqi] = matRad_indicatorWrapper(cst,pln,resultGUI.physicalDose*pln.numOfFractions,[10,20,30,40,50,60,70,80]/pln.numOfFractions,[2,5,10,20,30,40,50,60,70,80,90,95,98],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_nominal.fig']);

%% Copy reference plan
pln_robust=pln;

%% retrieve scenarios for dose calculation and optimziation
[multScen] = matRad_multiScenGenerator(run_config.scen_mode,run_config,'optimization',ct);

%% save multi scenarios to plan
pln_robust.multScen=multScen;

%% Dose calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
now1 = tic();
dij_robust = matRad_calcPhotonDose(ct,stf,pln_robust,cst);
DCTime_robust = toc(now1);
time1=sprintf('DCTime_robust: %.2f\n',DCTime_robust); disp(time1);
results.performance.DCTime_robust=DCTime_robust;

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
[cst,ixTarget,p,ixBody,ixCTV] = matRad_loadObjectives(run_config,run_config.plan_target,cst);

%% Print target objectives
display(cst{ixTarget,2});
for i=1:length(cst{ixTarget,6})
    display(cst{ixTarget,6}{i});
end

%% Trigger robust optimization
% Make the objective to a composite worst case objective

% Target
switch run_config.robustness
    case 'c-COWC'
        cst{ixTarget,8}{1}.beta1 = run_config.beta1;
        cst{ixTarget,8}{1}.p1 = run_config.p1;
        cst{ixTarget,8}{1}.beta2 = run_config.beta2;
        cst{ixTarget,8}{1}.p2 = run_config.p2;
        for i=1:length(cst{ixTarget,6})
            cst{ixTarget,6}{i}.robustness  = run_config.robustness;
        end
    otherwise
        for i=1:length(cst{ixTarget,6})
            cst{ixTarget,6}{i}.robustness  = run_config.robustness;
        end
end

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
results.performance.OPTTime_robust=OPTTime_robust;

%% Evaluating robust solution
w_robust=resultGUI_robust.w;
f_robust = matRad_calcObjectiveFunction(w_robust,dij,cst_nominal,pln);

%% Plot robust fluence
matRad_visSpotWeights(stf,resultGUI_robust.w);
savefig([folderPath filesep 'fluence_robust.fig']);

%% Create an interactive plot to slide through axial slices
for numScen=1:pln_robust.multScen.totNumScen
    
    if(pln_robust.multScen.totNumScen>1)
        quantityMap=[quantityOpt '_' num2str(round(numScen))];
    else
        quantityMap=quantityOpt;
    end
    
    plane      = 3;
    doseWindow = [0 max([max([resultGUI_robust.physicalDose(:)*pln_robust.numOfFractions]) run_config.doseWindow(2)])];
    doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
    f = figure;
    title([quantityMap]);
    set(gcf,'position',[10,10,550,400]);
    numSlices = ct.cubeDim(3);
    slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]');
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
          'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]');
    
    savefig([folderPath filesep 'dose_robust_' quantityMap '.fig']);   
end

%% Plot dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI_robust.physicalDose*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'dose3d_robust.fig']);

%% Create target mask
target_mask = zeros(size(resultGUI.physicalDose));
target_mask(cst{ixCTV,4}{1,1}) = 1;

%% Create OAR mask
OAR_mask = zeros(size(resultGUI.physicalDose));
OAR_mask(cst{ixBody,4}{1,1}) = 1;
OAR_mask(cst{ixCTV,4}{1,1}) = 0;

%% Plot target dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI_robust.physicalDose.*target_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'target_dose3d_robust.fig']);

%% Plot OAR dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI_robust.physicalDose.*OAR_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
savefig([folderPath filesep 'OAR_dose3d_robust.fig']);

%% Obtain dose statistics
[dvh_robust,dqi_robust] = matRad_indicatorWrapper(cst,pln_robust,resultGUI_robust.physicalDose*pln_robust.numOfFractions,[],[],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_robust.fig']);

%% check sampling option is activated
if ~run_config.sampling
    return;
end

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {}; % structSel = {'PTV','OAR1'};
[multScen] = matRad_multiScenGenerator(run_config.sampling_mode,run_config,'sampling',ct);

%% Perform sampling
[caSamp, mSampDose, plnSamp, resultGUInomScen] = matRad_sampling(ct,stf,cst,pln_robust,resultGUI_robust.w,structSel,multScen);

%% Perform sampling analysis
varargin.GammaCriterion = run_config.GammaCriterion; % [%  mm]
[cstStat, resultGUISamp, meta] = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen, varargin);

%% Multi-scenario dose volume histogram (DVH)
figure,set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVH_sampledScen(caSamp,dvh_robust,cst,plnSamp,[1:multScen.totNumShiftScen],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_robust_multiscen.fig']);

%% Dose volume histogram (DVH)
resultGUISamp_ul=[];
resultGUISamp_ul.physicalDose=resultGUI_robust.physicalDose;
resultGUISamp_ul.physicalDose_lower=resultGUISamp.meanCube-resultGUISamp.stdCube;
resultGUISamp_ul.physicalDose_upper=resultGUISamp.meanCube+resultGUISamp.stdCube;
[dvh_sampled,dqi_sampled] = matRad_indicatorWrapper_sampled(cst,pln,resultGUISamp_ul,[10,20,30,40,50,60,70,80]/pln.numOfFractions,[2,5,10,20,30,40,50,60,70,80,90,95,98],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_robust_trustband.fig']);

%% Create an interactive plot to slide through axial slices
quantityMap='stdCube';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.stdCube(:)*pln_robust.numOfFractions) run_config.doseWindow_uncertainty(2)])];
doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]');

savefig([folderPath filesep 'uncertainty.fig']);   

%% Plot uncertainty distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.stdCube*pln.numOfFractions,run_config.doseWindow_uncertainty,[0.005 0.00005],colorcube,[],'Dose uncertainty [Gy]');
savefig([folderPath filesep 'uncertainty3d.fig']);

%% Plot target uncertainty distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.stdCube.*target_mask*pln.numOfFractions,run_config.doseWindow_uncertainty,[0.1 0.00005],colorcube,[],'Dose uncertainty [Gy]');
savefig([folderPath filesep 'target_uncertainty3d_robust.fig']);

%% Calculate relative uncertainty
resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube = (resultGUISamp.stdCube./resultGUISamp.meanCube);

%% Create an interactive plot to slide through axial slices
analysis='uncertaintyAnalysis';
quantityMap='relativeUncertaintyCube';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.(analysis).(quantityMap)(:)*pln_robust.numOfFractions) run_config.doseWindow_relative_uncertainty1(2)])];
doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Relative uncertainty');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Relative uncertainty');

savefig([folderPath filesep 'relative_uncertainty.fig']);   

%% Plot uncertainty distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube*pln.numOfFractions,run_config.doseWindow_relative_uncertainty1,[0.005 0.00005],colorcube,[],'Relative uncertainty');
savefig([folderPath filesep 'relative_uncertainty3d.fig']);

%% Calculate robustness violation cube
resultGUISamp.uncertaintyAnalysis.robustnessViolationCube = resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube>run_config.UncertaintyCriterion;
resultGUISamp.uncertaintyAnalysis.robustnessViolationCube(isnan(resultGUISamp.uncertaintyAnalysis.robustnessViolationCube))=max(resultGUISamp.uncertaintyAnalysis.robustnessViolationCube(:));

%% Calculate robustness compliance cube
resultGUISamp.uncertaintyAnalysis.robustnessComplianceCube = resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube<=run_config.UncertaintyCriterion;
resultGUISamp.uncertaintyAnalysis.robustnessComplianceCube(isnan(resultGUISamp.uncertaintyAnalysis.robustnessComplianceCube))=max(resultGUISamp.uncertaintyAnalysis.robustnessComplianceCube(:));

%% Create an interactive plot to slide through axial slices
analysis='uncertaintyAnalysis';
quantityMap='robustnessViolationCube';
plane      = 3;
doseWindow = [0 1.01];
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');

savefig([folderPath filesep 'robustness_violation.fig']);  

%% uncertainty index distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.uncertaintyAnalysis.robustnessViolationCube,[0 1],[0.05 0.00005],colorcube,[],'Relative uncertainty violation');
savefig([folderPath filesep 'robustness_violation3d.fig']);

%% Create an interactive plot to slide through axial slices
analysis='uncertaintyAnalysis';
quantityMap='robustnessViolationCube';
plane      = 3;
doseWindow = [0 1.01];
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap).*target_mask,plane,slice,[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap).*target_mask,plane,round(es.Value),[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');

savefig([folderPath filesep 'target_robustness_violation.fig']);  

%% Plot target uncertainty index distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.uncertaintyAnalysis.robustnessViolationCube.*target_mask,[0 1],[0.5 0.00005],colorcube,[],'Relative uncertainty violation');
savefig([folderPath filesep 'target_robustness_violation3d.fig']);

%% Create an interactive plot to slide through axial slices
analysis='uncertaintyAnalysis';
quantityMap='robustnessComplianceCube';
plane      = 3;
doseWindow = [0 2.01];
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap),plane,slice,[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap),plane,round(es.Value),[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');

savefig([folderPath filesep 'robustness_compliance.fig']);

%% uncertainty index distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUISamp.uncertaintyAnalysis.robustnessComplianceCube,[0 2],[0.05 0.00005],colorcube,[],'Relative uncertainty violation');
savefig([folderPath filesep 'robustness_compliance3d.fig']);

%% Create an interactive plot to slide through axial slices
analysis='uncertaintyAnalysis';
quantityMap='robustnessComplianceCube';
plane      = 3;
doseWindow = [0 2.01];
f = figure;
title([quantityMap]);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);
slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap).*target_mask,plane,slice,[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
      'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(analysis).(quantityMap).*target_mask,plane,round(es.Value),[],[],colorcube,[],doseWindow,[],[],'Relative uncertainty');

savefig([folderPath filesep 'target_robustness_compliance.fig']);  

%% Uncertainty volume histogram (UVH)
resultGUISamp.physicalDose=resultGUISamp.stdCube;
[uvh,uqi] = matRad_indicatorWrapper_UVH(cst,pln,resultGUI,resultGUISamp,[],[],run_config.doseWindow_uvh);
savefig([folderPath filesep 'uvh.fig']);

%% Print evaluation indexes
% CTV
disp('CTV evaluation');

results.structures{ixCTV,1}=cst{ixCTV,2};
DMean=sprintf('DMean: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).mean*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).mean*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).mean*pln.numOfFractions); disp(DMean);
results.structures{ixCTV,2}.D_mean.nom=dqi_sampled{1,1}(ixCTV).mean*pln.numOfFractions;
results.structures{ixCTV,2}.D_mean.min=dqi_sampled{1,2}(ixCTV).mean*pln.numOfFractions;
results.structures{ixCTV,2}.D_mean.max=dqi_sampled{1,3}(ixCTV).mean*pln.numOfFractions;

D98=sprintf('D98: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).D_98*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_98*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_98*pln.numOfFractions); disp(D98);
results.structures{ixCTV,2}.D_98.nom=dqi_sampled{1,1}(ixCTV).D_98*pln.numOfFractions;
results.structures{ixCTV,2}.D_98.min=dqi_sampled{1,2}(ixCTV).D_98*pln.numOfFractions;
results.structures{ixCTV,2}.D_98.max=dqi_sampled{1,3}(ixCTV).D_98*pln.numOfFractions;

D95=sprintf('D95: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).D_95*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_95*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_95*pln.numOfFractions); disp(D95);
results.structures{ixCTV,2}.D_95.nom=dqi_sampled{1,1}(ixCTV).D_95*pln.numOfFractions;
results.structures{ixCTV,2}.D_95.min=dqi_sampled{1,2}(ixCTV).D_95*pln.numOfFractions;
results.structures{ixCTV,2}.D_95.max=dqi_sampled{1,3}(ixCTV).D_95*pln.numOfFractions;

D50=sprintf('D50: %.2f [%.2f - %.2f] Gy',dqi_sampled{1,1}(ixCTV).D_50*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_50*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_50*pln.numOfFractions); disp(D50);
results.structures{ixCTV,2}.D_50.nom=dqi_sampled{1,1}(ixCTV).D_50*pln.numOfFractions;
results.structures{ixCTV,2}.D_50.min=dqi_sampled{1,2}(ixCTV).D_50*pln.numOfFractions;
results.structures{ixCTV,2}.D_50.max=dqi_sampled{1,3}(ixCTV).D_50*pln.numOfFractions;

D2=sprintf('D2: %.2f [%.2f - %.2f] Gy \n',dqi_sampled{1,1}(ixCTV).D_2*pln.numOfFractions,dqi_sampled{1,2}(ixCTV).D_2*pln.numOfFractions,dqi_sampled{1,3}(ixCTV).D_2*pln.numOfFractions); disp(D2);
results.structures{ixCTV,2}.D_2.nom=dqi_sampled{1,1}(ixCTV).D_50*pln.numOfFractions;
results.structures{ixCTV,2}.D_2.min=dqi_sampled{1,2}(ixCTV).D_50*pln.numOfFractions;
results.structures{ixCTV,2}.D_2.max=dqi_sampled{1,3}(ixCTV).D_50*pln.numOfFractions;

% robustness indexes

disp('Robustness evaluation (Method 1)');

UMean=sprintf('UMean: %.2f Gy',uqi(ixCTV).mean*pln.numOfFractions); disp(UMean);
results.structures{ixCTV,2}.U_mean=uqi(ixCTV).mean*pln.numOfFractions;

U98=sprintf('U98: %.2f Gy',uqi(ixCTV).D_98*pln.numOfFractions); disp(U98);
results.structures{ixCTV,2}.U_98=uqi(ixCTV).D_98*pln.numOfFractions;

U95=sprintf('U95: %.2f Gy',uqi(ixCTV).D_95*pln.numOfFractions); disp(U95);
results.structures{ixCTV,2}.U_95=uqi(ixCTV).D_95*pln.numOfFractions;

U50=sprintf('U50: %.2f Gy',uqi(ixCTV).D_50*pln.numOfFractions); disp(U50);
results.structures{ixCTV,2}.U_50=uqi(ixCTV).D_50*pln.numOfFractions;

U2=sprintf('U2: %.2f Gy',uqi(ixCTV).D_2*pln.numOfFractions); disp(U2);
results.structures{ixCTV,2}.U_2=uqi(ixCTV).D_2*pln.numOfFractions;

UI=sprintf('UI (>=0): %.2f \n',uqi(ixCTV).D_2.*pln.numOfFractions/p); disp(UI);
results.robustness.UI1=uqi(ixCTV).D_2*pln.numOfFractions;

disp('Robustness evaluation (Method 2)');

UncertaintyIndex2=nnz(resultGUISamp.uncertaintyAnalysis.robustnessViolationCube(cst{ixCTV,4}{1,1}))/numel(cst{ixCTV,4}{1,1});

UI2=sprintf('UI (0-1): %.2f',UncertaintyIndex2); disp(UI2);
results.robustness.UI2=UncertaintyIndex2;

RobustnessIndex2=1-UncertaintyIndex2;
RI2=sprintf('RI (0-1): %.2f \n',RobustnessIndex2); disp(RI2);
results.robustness.RI2=RobustnessIndex2;

%% plot relative uncertainty distribution
figure;
edges = [0:0.025:run_config.doseWindow_relative_uncertainty1(2)];
h1=histogram(resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube,edges,'Normalization','probability', 'DisplayStyle','stairs','LineWidth', 2);
xlabel('Relative uncertainty');
ylabel('[%]');
ylim([0 0.05]);
savefig([folderPath filesep 'relative_uncertainty_histo.fig']);

%% plot target relative uncertainty distribution
figure;
edges = [0:0.0025:run_config.doseWindow_relative_uncertainty2(2)];
h2=histogram(resultGUISamp.uncertaintyAnalysis.relativeUncertaintyCube(cst{ixCTV,4}{1,1}),edges,'Normalization','probability', 'DisplayStyle','stairs','LineWidth', 2);
xlabel('Relative uncertainty');
ylabel('[%]');
ylim([0 0.5]);
savefig([folderPath filesep 'target_relative_uncertainty_histo.fig']);

%% robustness price indexes
disp('Robustness price evaluation (Method 1)');
OARMean_nominal=mean(resultGUI.physicalDose.*OAR_mask,'all')*pln.numOfFractions;
OARMean_robust=mean(resultGUI_robust.physicalDose.*OAR_mask,'all')*pln.numOfFractions;

RPI1=sprintf('RPI: %.2f\n',(OARMean_robust-OARMean_nominal)/OARMean_nominal); disp(RPI1);
results.robustness.RPI1=(OARMean_robust-OARMean_nominal)/OARMean_nominal;

disp('Robustness price evaluation (Method 2)');

F_nominal=sprintf('F(nominal): %.2f',f_nominal); disp(F_nominal);
results.robustness.f_nominal=f_nominal;

F_robust=sprintf('F(robust): %.2f \n',f_robust); disp(F_robust);
results.robustness.f_robust=f_robust;

RPI2=sprintf('RPI: %.2f \n',(f_robust-f_nominal)/f_nominal); disp(RPI2);
results.robustness.RPI2=(f_robust-f_nominal)/f_nominal;

%% print results
[results] = matRad_printResults(run_config,results,cst,pln,dqi_sampled);

%% Save outputs
save([folderPath filesep 'resultGUI.mat'],'resultGUI','resultGUI_robust');
save([folderPath filesep 'plan.mat'],'pln','pln_robust','ct','cst','stf','run_config');
save([folderPath filesep 'sampling.mat'],'caSamp', 'mSampDose', 'plnSamp', 'resultGUInomScen','cstStat','resultGUISamp','meta','dvh','dvh_robust');
save([folderPath filesep 'results.mat'],'results');

%%
diary off