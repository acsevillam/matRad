function matRad_photonRobust4D(radiationMode,description,varargin)

%% Example: 4D robust Treatment Planning with photons
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will
% (i)   import a 4D CT into a multiscenario ct and cst struct
% (ii)  create a photon treatment plan with seven beams
% (iii) perform dose calculation on each 4D CT
% (iv)  perform first a fluence optimization on the first CT scenario and then secondly
%       another fluence optimization using the composite worst case approach
%       considering all 4D CTs
% (v)   visualise all individual dose scenarios
% (vi)  sample discrete scenarios from Gaussian uncertainty assumptions

%% Clear variables
clearvars -except radiationMode description varargin ;
clc;
close 'all';

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;
s = settings;
s.matlab.general.matfile.SaveFormat.TemporaryValue = 'v7.3';

%% Set function parameters

validRadiationModes = {'photons','protons'};
validDescriptions = {'prostate','breast'};
validPatientIDs = {'3482','3648','3782','3790','3840','3477','3749','3832','3833','3929'};
validPlanObjectives = {'1','2','3','4','5','6'};
validPlanTargets = {'CTV','PTV'};
validPlanBeams = {'5F','7F','9F'};
validRobustness = {'none','STOCH','COWC','c-COWC','INTERVAL1','INTERVAL2','INTERVAL3'};
validScenModes = {'nomScen','wcScen','impScen','impScen_permuted','impScen_permuted_truncated','random','random_truncated'};

defaultPatientID = '3482';
defaultPlanObjective = '4';
defaultPlanTarget = 'CTV';
defaultPlanBeams = '9F';
defaultShiftSD = [5 5 5]; % mm
defaultRobustness = 'COWC';
defaultScenMode = 'wcScen';
defaultWCFactor = 1.0;
defaultP1 = 1;
defaultP2 = 1;
defaultTheta1 = 1.0;
defaultK = 10;
defaultTheta2 = 0.0;
defaultLoadDij = true;
defaultSampling = true;
defaultSamplingMode = 'impScen_permuted_truncated';
defaultSamplingWCFactor = 1.5;
defaultRootPath = matRad_cfg.matRadRoot;
defaultNCores = feature('numcores');

parser = inputParser;

addRequired(parser,'radiationMode',@(x) any(validatestring(x,validRadiationModes)));
addRequired(parser,'description',@(x) any(validatestring(x,validDescriptions)));
addParameter(parser,'caseID',defaultPatientID,@(x) any(validatestring(x,validPatientIDs)));
addParameter(parser,'plan_objectives',defaultPlanObjective,@(x) any(validatestring(x,validPlanObjectives)));
addParameter(parser,'plan_target',defaultPlanTarget,@(x) any(validatestring(x,validPlanTargets)));
addParameter(parser,'plan_beams',defaultPlanBeams,@(x) any(validatestring(x,validPlanBeams)));
addParameter(parser,'shiftSD',defaultShiftSD,@(x) numel(x) == 3 && isnumeric(x) && all(x > 0));
addParameter(parser,'robustness',defaultRobustness,@(x) any(validatestring(x,validRobustness)));
addParameter(parser,'scen_mode',defaultScenMode,@(x) any(validatestring(x,validScenModes)));
addParameter(parser,'wcFactor',defaultWCFactor,@(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(parser,'p1',defaultP1,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));
addParameter(parser,'p2',defaultP2,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));
addParameter(parser,'theta1',defaultTheta1,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','positive'}));
addParameter(parser,'k',defaultK,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));
addParameter(parser,'theta2',defaultTheta2,@(x) validateattributes(x,{'numeric'},...
            {'nonempty'}));
addParameter(parser,'loadDij',defaultLoadDij,@islogical);
addParameter(parser,'sampling',defaultSampling,@islogical);
addOptional(parser,'sampling_mode',defaultSamplingMode,@(x) any(validatestring(x,validScenModes)));
addOptional(parser,'sampling_wcFactor',defaultSamplingWCFactor,@(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(parser,'rootPath',defaultRootPath,@isfolder);
addParameter(parser,'n_cores',defaultNCores,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));

parse(parser,radiationMode,description,varargin{:});

run_config.radiationMode = parser.Results.radiationMode;
run_config.description = parser.Results.description;
run_config.caseID = parser.Results.caseID;
run_config.plan_objectives = parser.Results.plan_objectives;
run_config.plan_target = parser.Results.plan_target;
run_config.plan_beams = parser.Results.plan_beams;
run_config.shiftSD = parser.Results.shiftSD;
run_config.robustness = parser.Results.robustness;
run_config.scen_mode = parser.Results.scen_mode;
run_config.wcFactor = parser.Results.wcFactor;
run_config.sampling = parser.Results.sampling;
run_config.loadDij = parser.Results.loadDij;
run_config.sampling_mode = parser.Results.sampling_mode;
run_config.sampling_wcFactor = parser.Results.sampling_wcFactor;
run_config.rootPath = parser.Results.rootPath;
run_config.n_cores = parser.Results.n_cores;

switch run_config.robustness
    case "c-COWC"
        switch run_config.scen_mode
            case "wcScen"
                run_config.numScens = 7;
            case "impScen"
                run_config.numScens = 13;
            case "impScen_permuted_truncated"
                run_config.numScens = 33;
        end
        run_config.p1 = parser.Results.p1;
        run_config.p2 = parser.Results.p2;
        run_config.beta1 = run_config.p1/run_config.numScens;
        run_config.beta2 = run_config.p2/run_config.numScens;
        output_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.beta1) '_to_' num2str(run_config.beta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];

    case "INTERVAL2"
        run_config.theta1 = parser.Results.theta1;
        output_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.theta1) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
        dij_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval2.mat'];
        if run_config.loadDij && isfile(dij_file)
            load(dij_file,'pln_robust','dij_robust','dij_interval');
        end
    case "INTERVAL3"
        run_config.theta1 = parser.Results.theta1;
        run_config.k = parser.Results.k;
        run_config.theta2 = parser.Results.theta2;
        output_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.theta1) filesep num2str(run_config.theta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
        dij_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3.mat'];
        if run_config.loadDij && isfile(dij_file)
            load(dij_file,'pln_robust','dij_robust','dij_interval');
        end
    otherwise
        output_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
end

run_config.resolution = '5x5x5';
run_config.GammaCriterion = [3 3];
run_config.robustnessCriterion = [10 10];
run_config.sampling_size = 50;

%Set up parent export folder and full file path
if ~(isfolder(output_folder))
    mkdir(run_config.rootPath, output_folder);
end

folderPath = [run_config.rootPath filesep output_folder];

%% Initiallize diary log
diary([folderPath filesep 'diary.log'])
diary on

%% Set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Import 3D CT
[ct,cst] = matRad_loadGeometry(run_config);

%% Import 4D CT
metadata.resolution = [3 3 3]; % [mm]
[ct,cst] = matRad_importMultipleDicomCt('/jobs/images/breast/patient2_4D/dicom',metadata);
clear 'metadata';

%%
for  it = 1:size(cst,1)
    switch cst{it,2}
        case 'Skin'
            cst{it,2}='BODY';
        case 'PTV'
            cst{it,2}='PTV';
        case 'CTV'
            cst{it,2}='CTV';
        case 'CORAZON'
            cst{it,2}='HEART';
        case 'PULMON IZQUIERDO'
            cst{it,2}='LEFT LUNG';
        case 'PULMON DERECHO'
            cst{it,2}='RIGTH LUNG';
    end
    if isempty(cst{it,4}{1,1}) == true
        fprintf(' %s is an empty structure. \n',cst{it,2});
        fprintf('Deleting %s structure. \n',cst{it,2});
        cst(it,:) = [];
    end
end


%% Instantiate elastic registration
metadata.nItera = 100;
metadata.dvfType = 'pull';
register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
clear 'metadata';

%% Calculate deformation vector field
[ct,cst] = register.calcDVF();

%% Print run config
disp(run_config);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
[cst,ixTarget,p,ixBody,ixCTV,OARStructSel] = matRad_loadObjectives(run_config,'CTV',cst);

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
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice);
        
        plane      = 3;
        slice      = round(isocenter(3)./ct.resolution.z);
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice);
        
    end
    
end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = round(isocenter(2)./ct.resolution.y);
    numScen    = 1;
    matRad_plotSliceWrapper(gca,ct,cst,numScen,[],plane,slice);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(es.Value),[],plane,slice);
end
clear  numScen plane slice ans f b;

%% Plot CT slice
%{
if param.logLevel == 1
    
    figure('Renderer', 'painters', 'Position', [10 10 300*ct.numOfCtScen 400]);
    
    isocenter = matRad_getIsoCenter(cst,ct,0);
    
    for scen_iterator = 1:ct.numOfCtScen
        plane      = 1;
        slice      = round(isocenter(2)./ct.resolution.y);
        subplot(2,ct.numOfCtScen,scen_iterator); camroll(90);
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice,[],[],[],[],[],[],[],[],[],'LineWidth',1.2);
        
        plane      = 3;
        slice      = round(isocenter(3)./ct.resolution.z);
        subplot(2,ct.numOfCtScen,scen_iterator+ct.numOfCtScen);
        matRad_plotSliceWrapper(gca,ct,cst,scen_iterator,[],plane,slice,[],[],[],[],[],[],[],[],[],'LineWidth',1.2);
    end
    
end
clear  scen_iterator plane slice ans;

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios'); camroll(90);
    plane      = 1;
    slice      = round(isocenter(2)./ct.resolution.y);
    numScen    = 1;
    matRad_plotSliceWrapper(gca,ct,cst,numScen,[],plane,slice,[],[],[],[],[],[],[],[],[],'LineWidth',1.2);
    b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
        'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
    b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,round(es.Value),[],plane,slice,[],[],[],[],[],[],[],[],[],'LineWidth',1.2);
end
clear  numScen plane slice ans f b;
%}

savefig([folderPath filesep 'ct.fig']);

%% Set plot and histograms window
run_config.doseWindow = [0 p*1.25];
run_config.doseWindow_dvh = [0 p*1.6];
run_config.doseWindow_uncertainty = [0 p*0.5];
run_config.doseWindow_relative_uncertainty1 = [0 1];
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

if run_config.radiationMode == "photons"
    pln.radiationMode = 'photons';
    pln.machine       = 'Generic';
    quantityOpt    = 'physicalDose';
    modelName      = 'none';
    
end

if run_config.radiationMode == "protons"
    pln.radiationMode = 'protons';
    pln.machine       = 'Generic';
    quantityOpt   = 'RBExD';            % either  physicalDose / effect / RBExD
    modelName     = 'constRBE';         % none: for photons, protons, carbon                                    
    % constRBE: constant RBE model
    % MCN: McNamara-variable RBE model for protons
    % WED: Wedenberg-variable RBE model for protons
    % LEM: Local Effect Model for carbon ions
    % calculate LET distribution
    pln.propDoseCalc.calcLET = 0;
end

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

%% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

%% Retrieve nominal scenario for dose calculation and optimziation reference
multScen = matRad_multScen(ct,'nomScen');

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

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
if run_config.radiationMode == "photons"
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
end

if run_config.radiationMode == "protons"
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

%% Inverse Optimization  for IMPT based on RBE-weighted dose
% The goal of the fluence optimization is to find a set of bixel/spot
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
%
resultGUI_nominal = matRad_fluenceOptimization(dij,cst,pln);

% add resultGUI_nominal dose cubes to resultGUI structure to allow the visualization in the GUI
resultGUI = resultGUI_nominal;

%% Plot nominal fluence
matRad_visSpotWeights(stf,resultGUI.w);
savefig([folderPath filesep 'fluence_nominal.fig']);

%% Plot dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([folderPath filesep 'dose3d_nominal.fig']);

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
%savefig([folderPath filesep 'target_dose3d_nominal.fig']);

%% Plot OAR dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose.*OAR_mask*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([folderPath filesep 'OAR_dose3d_nominal.fig']);

%% Indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI_nominal,[],pln.numOfFractions,[],[],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_nominal.fig']);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
cst_robust=cst;
[cst_robust,ixTarget,p,ixBody,ixCTV,OARStructSel] = matRad_loadObjectives(run_config,run_config.plan_target,cst_robust);

%% Print target objectives
display(cst_robust{ixTarget,2});
for i=1:length(cst_robust{ixTarget,6})
    display(cst_robust{ixTarget,6}{i});
end

%% Copy reference plan
pln_robust=pln;

%% retrieve scenarios for dose calculation and optimziation
[multScen] = matRad_multiScenGenerator(run_config.scen_mode,run_config,'optimization',ct);

%% save multi scenarios to plan
pln_robust.multScen=multScen;

%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf_robust = matRad_generateStf(ct,cst_robust,pln_robust);

%%
% Let's display the beam geometry information of the 6th beam
display(stf_robust(1));

%% Dose calculation
% Let's generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows subsequent inverse optimization.
now1 = tic();
if ~exist('dij_robust','var') || isempty(dij_robust)
    if run_config.radiationMode == "photons"
        dij_robust = matRad_calcPhotonDose(ct,stf_robust,pln_robust,cst_robust);
    end

    if run_config.radiationMode == "protons"
        dij_robust = matRad_calcParticleDose(ct,stf_robust,pln_robust,cst_robust);
    end
end
DCTime_robust = toc(now1);
time1=sprintf('DCTime_robust: %.2f\n',DCTime_robust); disp(time1);
results.performance.DCTime_robust=DCTime_robust;

%% Dose interval pre-computing

switch run_config.robustness
    case 'INTERVAL2'
        targetStructSel = {'CTV'};
        now2 = tic();
        [dij_robust,pln_robust,dij_interval] = matRad_calcDoseInterval2(ct,cst,stf_robust,pln_robust,dij_robust,targetStructSel);
        dij_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval2.mat'];
        save(dij_file,'dij_robust','pln_robust','dij_interval', '-v7.3');
        %load('dij_interval2.mat');
        IDCTime_robust = toc(now2);
        time2=sprintf('IDCTime_robust: %.2f\n',IDCTime_robust); disp(time2);
        results.performance.IDCTime_robust=IDCTime_robust;
    case 'INTERVAL3'
        targetStructSel = {'CTV'};
        now2 = tic();
        if ~exist('dij_interval','var') || isempty(dij_interval)
            [dij_robust,pln_robust,dij_interval] = matRad_calcDoseInterval3(ct,cst,stf_robust,pln_robust,dij_robust,targetStructSel,OARStructSel,run_config.k);
            dij_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3.mat'];
            save(dij_file,'dij_robust','pln_robust','dij_interval', '-v7.3');
        end
        IDCTime_robust = toc(now2);
        time2=sprintf('IDCTime_robust: %.2f\n',IDCTime_robust); disp(time2);
        results.performance.IDCTime_robust=IDCTime_robust;
end


%% Trigger robust optimization
% Make the objective to a composite worst case objective

% Target
switch run_config.robustness
    case 'c-COWC'
        
        pln_robust.propOpt.p1=run_config.p1;
        pln_robust.propOpt.p2=run_config.p2;
        pln_robust.propOpt.useMaxApprox='cheapCOWC';
        
        for i=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{i}.robustness  = 'COWC';
        end
        
    case 'INTERVAL2'
        
        pln_robust.propOpt.dij_interval=dij_interval;
        pln_robust.propOpt.theta1=run_config.theta1;
        cst_robust{ixCTV,6}=[];
        cst_robust{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation2(800,p));
        cst_robust{ixCTV,6}{1}.robustness  = 'INTERVAL2';
    
    case 'INTERVAL3'
        
        pln_robust.propOpt.dij_interval=dij_interval;
        pln_robust.propOpt.theta1=run_config.theta1;
        pln_robust.propOpt.theta2=run_config.theta2;
        
        cst_robust{ixCTV,6}=[];
        cst_robust{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation2(800,p));
        cst_robust{ixCTV,6}{1}.robustness  = 'INTERVAL3';
        
        for i=1:size(cst,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst{i,2})
                    for k = 1:numel(cst{i,6})
                        cst_robust{i,6}{k}.robustness  = 'INTERVAL3';
                    end
                end
            end
        end

end

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
profile_master = parallel.importProfile('profile1.mlsettings');
pool=parpool(profile_master,run_config.n_cores);

profile on;
now3 = tic();
resultGUI_robust = matRad_fluenceOptimization(dij_robust,cst_robust,pln_robust);
% add resultGUI_robust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
resultGUI = matRad_appendResultGUI(resultGUI,resultGUI_robust,0,'robust');
OPTTime_robust = toc(now3);
time3=sprintf('OPTTime_robust: %.2f\n',OPTTime_robust); disp(time3);
results.performance.OPTTime_robust=OPTTime_robust;
profiler = profile('info');
save([folderPath filesep 'profiler'],'profiler');

delete(gcp('nocreate'))
parallel.internal.ui.MatlabProfileManager.removeProfile('profile1');

%% Plot robust fluence
matRad_visSpotWeights(stf_robust,resultGUI_robust.w);
savefig([folderPath filesep 'fluence_robust.fig']);

%% Create an interactive plot to slide through axial slices

if run_config.radiationMode == "photons" && false
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
        matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
        b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
            'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
        b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
        
        savefig([folderPath filesep 'dose_robust_' quantityMap '.fig']);
    end
end

if run_config.radiationMode == "protons" || true
    
    quantityMap=quantityOpt;
    
    plane      = 3;
    doseWindow = [0 max([max([resultGUI_robust.physicalDose(:)*pln_robust.numOfFractions]) run_config.doseWindow(2)])];
    doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
    f = figure;
    title([quantityMap]);
    set(gcf,'position',[10,10,550,400]);
    numSlices = ct.cubeDim(3);
    slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
    matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
    
    savefig(f,[folderPath filesep 'dose_robust_' quantityMap '.fig']);
    
end

%% Plot dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust.physicalDose*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([folderPath filesep 'dose3d_robust.fig']);

%% Create target mask
target_mask = zeros(size(resultGUI.physicalDose));
target_mask(cst_robust{ixCTV,4}{1,1}) = 1;

%% Create OAR mask
OAR_mask = zeros(size(resultGUI.physicalDose));
OAR_mask(cst_robust{ixBody,4}{1,1}) = 1;
OAR_mask(cst_robust{ixCTV,4}{1,1}) = 0;

%% Plot target dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust.physicalDose.*target_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([folderPath filesep 'target_dose3d_robust.fig']);

%% Plot OAR dose distribution
figure;
matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust.physicalDose.*OAR_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([folderPath filesep 'OAR_dose3d_robust.fig']);

%% Indicator calculation and show DVH and QI
[dvh_robust,dqi_robust] = matRad_indicatorWrapper(cst,pln,resultGUI_robust,[],pln_robust.numOfFractions,[],[],run_config.doseWindow_dvh);
savefig([folderPath filesep 'dvh_robust.fig']);

%% check sampling option is activated
if ~run_config.sampling
    return;
end

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {};
[multScen] = matRad_multiScenGenerator(run_config.sampling_mode,run_config,'sampling',ct);

%% Perform sampling for nominal optimization results
[caSamp, mSampDose, plnSamp, resultGUInomScen,resultGUIsampledScen] = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen);

%% Perform sampling analysis
phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
robustnessCriterion = run_config.robustnessCriterion;
GammaCriterion = run_config.GammaCriterion;
slice = round(isocenter(3)./ct.resolution.z);
[cstStat, resultGUISamp, meta] = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen,phaseProb,'robustnessCriterion',robustnessCriterion,'GammaCriterion',GammaCriterion,'slice',slice);
results.robustnessAnalysis_nominal=resultGUISamp.robustnessAnalysis;

savefig([folderPath filesep 'sampling_analysis_nominal.fig']);

%% Create an mean dose interactive plot to slide through axial slices
quantityMap='meanCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.(quantityMap)(:)) p*1.25])];
maxDose       = max(resultGUISamp.(quantityMap)(:));
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for nominal optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

savefig([folderPath filesep 'mean_dose_nominal.fig']);

%% Create an std dose interactive plot to slide through axial slices
quantityMap='stdCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISamp.(quantityMap)(:)) p*0.5])];
maxDose       = max(resultGUISamp.(quantityMap)(:));
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for nominal optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

savefig([folderPath filesep 'std_dose_nominal.fig']);

%% Multi-scenario dose volume histogram (DVH)
f = figure;
set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSamp,pln.numOfFractions,cst,plnSamp,[1:plnSamp.multScen.totNumScen],run_config.doseWindow_dvh,'trustband',1);
title('Multi-scenario DVH for nominal optimization results');

savefig([folderPath filesep 'dvh_trustband_nominal.fig']);

%% Perform sampling for robust optimization results
[caSampRob, mSampDoseRob, plnSampRob, resultGUIRobNomScen,resultGUIsampledScenRob] = matRad_sampling(ct,stf_robust,cst,pln_robust,resultGUI_robust.w,structSel,multScen);

%% Perform sampling analysis
phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
robustnessCriterion = run_config.robustnessCriterion;
GammaCriterion = run_config.GammaCriterion; 
slice = round(isocenter(3)./ct.resolution.z);
[cstStatRob, resultGUISampRob, metaRob] = matRad_samplingAnalysis(ct,cst,plnSampRob,caSampRob, mSampDoseRob, resultGUIRobNomScen,phaseProb,'robustnessCriterion',robustnessCriterion,'GammaCriterion',GammaCriterion,'slice',slice);
results.robustnessAnalysis_robust=resultGUISampRob.robustnessAnalysis;

savefig([folderPath filesep 'sampling_analysis_robust.fig']);

%% Create an mean dose interactive plot to slide through axial slices
quantityMap='meanCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p*1.25])];
maxDose       = max(resultGUISampRob.(quantityMap)(:));
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for robust optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

savefig([folderPath filesep 'mean_dose_robust.fig']);

%% Create an std dose interactive plot to slide through axial slices
quantityMap='stdCubeW';
plane      = 3;
doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p*0.5])];
maxDose       = max(resultGUISampRob.(quantityMap)(:));
doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
f = figure;
title([quantityMap 'for robust optimization results']);
set(gcf,'position',[10,10,550,400]);
numSlices = ct.cubeDim(3);

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

savefig([folderPath filesep 'std_dose_robust.fig']);

%% Multi-scenario dose volume histogram (DVH)
f = figure;
set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSampRob,pln_robust.numOfFractions,cst,plnSampRob,[1:plnSampRob.multScen.totNumScen],run_config.doseWindow_dvh,'trustband',1);
title('Multi-scenario DVH for robust optimization results');

savefig([folderPath filesep 'dvh_trustband_robust.fig']);

%% Perform price of robustness analysis in nominal scenario
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
[resultGUISampRob] = matRad_priceOfRobustnessIndex(resultGUISampRob,resultGUI_nominal.(quantityOpt),resultGUIRobNomScen.(quantityOpt),ct,cst,pln_robust,[],[],[],[-5 5],'relative',slice);
savefig([folderPath filesep 'price_in_nominal.fig']);

%% Perform price of robustness analysis using mean dose
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
[resultGUISampRob] = matRad_priceOfRobustnessIndex(resultGUISampRob,resultGUISamp.meanCubeW,resultGUISampRob.meanCubeW,ct,cst,pln_robust,[],[],[],[-5 5],'relative',slice);
results.priceOfRobustnesAnalysis=resultGUISampRob.priceOfRobustnesAnalysis;

savefig([folderPath filesep 'price_in_mean.fig']);

%% print results
disp('Performance:');
disp(results.performance);
disp('Robustness analysis (nominal plan):');
disp(results.robustnessAnalysis_nominal);
disp('Robustness analysis (robust plan):');
disp(results.robustnessAnalysis_robust);
disp('Price of Robustness analysis (robust plan):');
disp(results.priceOfRobustnesAnalysis.priceOfRobustness);

%% Save outputs
save([folderPath filesep 'resultGUI.mat'],'resultGUI','resultGUI_robust');
save([folderPath filesep 'plan.mat'],'ct','cst','cst_robust','pln','pln_robust','stf','stf_robust','run_config');
save([folderPath filesep 'sampling.mat'],'caSamp', 'mSampDose', 'plnSamp', 'resultGUInomScen', 'resultGUISampRob','cstStat','resultGUISamp','meta','dvh','dvh_robust');
save([folderPath filesep 'results.mat'],'results');

%%
diary off