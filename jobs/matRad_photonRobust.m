function matRad_photonRobust(radiationMode,description,varargin)

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
validPatientIDs = {'3482','3648','3782','3790','3840','3477','3749','3832','3833','3929','4136','4155','4203','4357','4390','4428','4494','4531','4585','4681','1758'};
validAcquisitionTypes = {'mat','dicom'};
validPlanObjectives = {'1','2','3','4','5','6'};
validDosePulling1Targets = {'CTV','PTV'};
validDosePulling1Criteria = {'COV_95','COV_98','COV_99','COV1'};
validDosePulling2Criteria = {'minQiTarget','meanQiTarget'};
validPlanTargets = {'CTV','PTV'};
validPlanBeams = {'5F','7F','9F'};
validRobustness = {'none','STOCH','COWC','COWC2','c-COWC','c-COWC2','INTERVAL1','INTERVAL2','INTERVAL3'};
validKdin = {'dinamic','static'};
validScenModes = {'nomScen','wcScen','impScen5','impScen7','impScen_permuted5','impScen_permuted7','impScen_permuted_truncated5','impScen_permuted_truncated7','random','random_truncated'};

defaultPatientID = '3482';
defaultDoseResolution = [5 5 5]; % mm
defaultAcquisitionType = 'dicom';
defaultPlanObjective = '4';
defaultDosePulling1 = false;
defaultDosePulling1Target = 'CTV';
defaultDosePulling1Criteria = 'COV1';
defaultDosePulling1Limit = 0.98;
defaultDosePulling1Start = 0;
defaultScaleFactor = 1.0;
defaultDosePulling2 = false;
defaultDosePulling2Criteria = 'meanQiTarget';
defaultDosePulling2Limit = 0.80;
defaultDosePulling2Start = 0;
defaultPlanTarget = 'CTV';
defaultPlanBeams = '9F';
defaultShiftSD = [5 10 5]; % mm
defaultRobustness = 'COWC';
defaultScenMode = 'wcScen';
defaultWCFactor = 1.0;
defaultP1 = 1;
defaultP2 = 1;
defaultTheta1 = 1.0;
defaultKdin = 'dinamic';
defaultKmax = 10;
defaultRetentionThreshold = 0.95;
defaultTheta2 = 1.0;
defaultLoadDij = true;
defaultSampling = true;
defaultSamplingMode = 'impScen_permuted_truncated5';
defaultSamplingWCFactor = 1.5;
defaultRootPath = matRad_cfg.matRadRoot;
defaultNCores = feature('numcores');

parser = inputParser;

addRequired(parser,'radiationMode',@(x) any(validatestring(x,validRadiationModes)));
addRequired(parser,'description',@(x) any(validatestring(x,validDescriptions)));
addParameter(parser,'caseID',defaultPatientID,@(x) any(validatestring(x,validPatientIDs)));
addParameter(parser,'doseResolution',defaultDoseResolution,@(x) numel(x) == 3 && isnumeric(x) && all(x > 0));
addParameter(parser,'AcquisitionType',defaultAcquisitionType,@(x) any(validatestring(x,validAcquisitionTypes)));
addParameter(parser,'plan_objectives',defaultPlanObjective,@(x) any(validatestring(x,validPlanObjectives)));
addParameter(parser,'dose_pulling1',defaultDosePulling1,@islogical);
addParameter(parser,'dose_pulling1_target',defaultDosePulling1Target,@(x) numel(x) >= 1 && all(ismember(x,validDosePulling1Targets)));
addParameter(parser,'dose_pulling1_criteria',defaultDosePulling1Criteria,@(x) numel(x) >= 1 && all(ismember(x,validDosePulling1Criteria)));
addOptional(parser,'dose_pulling1_limit',defaultDosePulling1Limit,@(x) numel(x) >= 1 && isnumeric(x) && all(x > 0));
addOptional(parser,'dose_pulling1_start',defaultDosePulling1Start,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','nonnegative'}));
addParameter(parser,'scale_factor',defaultScaleFactor,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','nonnegative'}));
addParameter(parser,'dose_pulling2',defaultDosePulling2,@islogical);
addParameter(parser,'dose_pulling2_criteria',defaultDosePulling2Criteria,@(x) numel(x) >= 1 && all(ismember(x,validDosePulling2Criteria)));
addOptional(parser,'dose_pulling2_limit',defaultDosePulling2Limit,@(x) numel(x) >= 1 && isnumeric(x) && all(x > 0));
addOptional(parser,'dose_pulling2_start',defaultDosePulling2Start,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','nonnegative'}));
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
addParameter(parser,'theta1',defaultTheta1,@(x) numel(x) >= 1 && isnumeric(x) && all(x >= 0));
addParameter(parser,'kdin',defaultKdin,@(x) any(validatestring(x,validKdin)));
addParameter(parser,'kmax',defaultKmax,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));
addParameter(parser,'retentionThreshold',defaultRetentionThreshold,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','nonnegative'}));
addParameter(parser,'theta2',defaultTheta2,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','nonnegative'}));
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
run_config.doseResolution = parser.Results.doseResolution;
run_config.AcquisitionType = parser.Results.AcquisitionType;
run_config.plan_objectives = parser.Results.plan_objectives;
run_config.dose_pulling1 = parser.Results.dose_pulling1;
run_config.dose_pulling1_target = parser.Results.dose_pulling1_target;
run_config.dose_pulling1_criteria = parser.Results.dose_pulling1_criteria;
run_config.dose_pulling1_limit = parser.Results.dose_pulling1_limit;
run_config.dose_pulling1_start = parser.Results.dose_pulling1_start;
run_config.scale_factor = parser.Results.scale_factor;
run_config.dose_pulling2 = parser.Results.dose_pulling2;
run_config.dose_pulling2_criteria = parser.Results.dose_pulling2_criteria;
run_config.dose_pulling2_limit = parser.Results.dose_pulling2_limit;
run_config.dose_pulling2_start = parser.Results.dose_pulling2_start;
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
    case {"c-COWC","c-COWC2"}
        switch run_config.scen_mode
            case "wcScen"
                run_config.numScens = 7;
            case "impScen5"
                run_config.numScens = 13;
            case "impScen_permuted_truncated5"
                run_config.numScens = 33;
        end
        run_config.p1 = parser.Results.p1;
        run_config.p2 = parser.Results.p2;
        run_config.beta1 = run_config.p1/run_config.numScens;
        run_config.beta2 = run_config.p2/run_config.numScens;
        num_plans=1;
        output_folder=cell(1,1);
        output_folder{1} = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep num2str(run_config.shiftSD(1)) 'x' num2str(run_config.shiftSD(2)) 'x' num2str(run_config.shiftSD(3)) filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.beta1) '_to_' num2str(run_config.beta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];

    case "INTERVAL2"
        run_config.theta1 = parser.Results.theta1;
        num_plans=1;
        output_folder=cell(1,1);
        output_folder{1} = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep num2str(run_config.shiftSD(1)) 'x' num2str(run_config.shiftSD(2)) 'x' num2str(run_config.shiftSD(3)) filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.theta1) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
        dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval2.mat'];
        if run_config.loadDij && isfile(dij_interval_file)
            load(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval');
        end
        dij_robust_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_robust_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '.mat'];
        if run_config.loadDij && isfile(dij_robust_file)
            load(dij_robust_file,'dij_robust');
        end
    case "INTERVAL3"
        run_config.theta1 = parser.Results.theta1;
        run_config.kdin = parser.Results.kdin;
        run_config.kmax = parser.Results.kmax;
        run_config.retentionThreshold = parser.Results.retentionThreshold;
        num_plans=length(run_config.theta1);
        run_config.theta2 = parser.Results.theta2;
        root_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep num2str(run_config.shiftSD(1)) 'x' num2str(run_config.shiftSD(2)) 'x' num2str(run_config.shiftSD(3)) filesep run_config.scen_mode filesep num2str(run_config.wcFactor)];
        output_folder=cell(1, length(run_config.theta1));
        if length(run_config.theta1)>1
            root_folder = [root_folder  filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
            for planIx = 1:num_plans
                output_folder{planIx} = [root_folder filesep num2str(run_config.theta1(planIx)) filesep num2str(run_config.retentionThreshold) filesep num2str(run_config.theta2)];
            end
        else
            output_folder{1} = [root_folder filesep num2str(run_config.theta1(1)) filesep num2str(run_config.retentionThreshold) filesep num2str(run_config.theta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
            root_folder = output_folder{1};
        end
        if(isequal(run_config.kdin,'dinamic'))
            dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_' num2str(run_config.retentionThreshold) '.mat'];
        elseif(isequal(run_config.kdin,'static'))
            dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_k_' num2str(run_config.kmax) '.mat'];
        end 
        if run_config.loadDij && isfile(dij_interval_file)
            load(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval');
        end
        dij_robust_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_robust_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '.mat'];
        if run_config.loadDij && isfile(dij_robust_file)
            load(dij_robust_file,'dij_robust');
        end
    otherwise
        num_plans=1;
        output_folder=cell(1, 1);
        output_folder{1} = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep num2str(run_config.shiftSD(1)) 'x' num2str(run_config.shiftSD(2)) 'x' num2str(run_config.shiftSD(3)) filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
end

if ~exist('root_folder','var') || isempty(root_folder)
    root_folder=output_folder;
end

run_config.resolution = [3 3 3];
run_config.GammaCriteria = [3 3];
run_config.robustnessCriteria = [5 5];
run_config.sampling_size = 50;

%Set up parent export folder and full file path
if ~(isfolder(root_folder))
    mkdir(run_config.rootPath, root_folder);
end
for planIx = 1:num_plans
    if ~(isfolder(output_folder{planIx}))
        mkdir(run_config.rootPath, output_folder{planIx});
    end
end

rootPath = [run_config.rootPath filesep root_folder];
folderPath=cell(1,num_plans);
for planIx = 1:num_plans
    folderPath{planIx} = [run_config.rootPath filesep output_folder{planIx}];
end

%% Initiallize diary log
diary([rootPath filesep 'diary.log'])
diary on

%% Set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Import CT and rename structures
[ct,cst] = matRad_loadGeometry(run_config);
cst = matRad_renameStructures(cst,run_config);

%% Calculate deformation vector field
if (ct.numOfCtScen>1)
    metadata.nItera = 100;
    metadata.dvfType = 'pull';
    register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
    [ct] = register.calcDVF();
    clear metadata;
end

%% Print run config
disp(run_config);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
run_config_tmp=run_config;
run_config_tmp.plan_objectives='4';
dp_start=[run_config_tmp.dose_pulling1_start 0];
[cst,ixTarget,p,ixBody,ixCTV,~] = matRad_loadObjectives(run_config_tmp,run_config_tmp.dose_pulling1_target{end},dp_start,cst);

%% Check target visibility 
cst{ixTarget,5}.Visible=true;

%% Add 0 - 20 mm ring 
vInnerMargin.x=0;
vInnerMargin.y=0;
vInnerMargin.z=0;

vOuterMargin.x=20;
vOuterMargin.y=20;
vOuterMargin.z=20;

metadata.name='RING 0 - 20 mm';
metadata.type='OAR';
metadata.visibleColor=[0,1,0.501960784313726];
[cst,ixRing1] = matRad_createRing(ixTarget,ixBody,cst,ct,vOuterMargin,vInnerMargin,metadata);
clear metadata;

% Add 20 - 50 mm ring 
vInnerMargin.x=20;
vInnerMargin.y=20;
vInnerMargin.z=20;

vOuterMargin.x=50;
vOuterMargin.y=50;
vOuterMargin.z=50;

metadata.name='RING 20 - 50 mm';
metadata.type='OAR';
metadata.visibleColor=[0,1,0.501960784313726];
[cst,ixRing2] = matRad_createRing(ixTarget,ixBody,cst,ct,vOuterMargin,vInnerMargin,metadata);
clear metadata;

%% Define ring objectives
cst{ixRing1,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixRing1,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.10,0));
cst{ixRing1,6}{1}.robustness  = 'none';
cst{ixRing1,6}{1}.dosePulling  = false;

cst{ixRing2,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
cst{ixRing2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
cst{ixRing2,6}{1}.robustness  = 'none';
cst{ixRing2,6}{1}.dosePulling  = false;

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

savefig([rootPath filesep 'ct.fig']);

if (ct.numOfCtScen>1)
    f          = figure; title('individual scenarios');
    plane      = 3;
    slice      = round(isocenter(2)./ct.resolution.y);
    cubeIdx1=1;
    cubeIdx2=2;
    tmp_cube1=ct.cubeHU{cubeIdx1};
    tmp_cube2=ct.cubeHU{cubeIdx2};
    tmp_cube2_moved=imwarp(tmp_cube2, permute(ct.dvf{cubeIdx2},[2 3 4 1]));
    subplot(2,2,1);
    matRad_plotSliceWrapper(gca,ct,cst,cubeIdx1,[],plane,slice);
    title(['Ct scenario No. ' int2str(cubeIdx1)]);%camroll(90);
    subplot(2,2,2);
    matRad_plotSliceWrapper(gca,ct,cst,cubeIdx2,[],plane,slice);
    title(['Ct scenario No. ' int2str(cubeIdx2)]);%camroll(90);
    subplot(2,2,3);
    matRad_compareCtSlice(gca,tmp_cube1,tmp_cube2,plane,slice);
    title('comparison without correction');%camroll(90);
    subplot(2,2,4);
    matRad_compareCtSlice(gca,tmp_cube1,tmp_cube2_moved,plane,slice);
    title('comparison with correction');%camroll(90);
    savefig([rootPath filesep 'image_registration_1to2.fig']);
   
end

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
pln.propDoseCalc.doseGrid.resolution.x = run_config.doseResolution(1); % [mm]
pln.propDoseCalc.doseGrid.resolution.y = run_config.doseResolution(2); % [mm]
pln.propDoseCalc.doseGrid.resolution.z = run_config.doseResolution(3); % [mm]

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

%% Indicator calculation and show DVH and QI
[~,qi] = matRad_indicatorWrapper(cst,pln,resultGUI_nominal,[],pln.numOfFractions,[],[],run_config.doseWindow_dvh);
savefig([rootPath filesep 'dvh_nominal.fig']);

%% Dose pulling Step 1

ixTargetQi=zeros(size(run_config.dose_pulling1_target));
for iX=1:size(run_config.dose_pulling1_target,2)
    for j=1:size(qi,2)
        if strcmp(qi(j).name,run_config.dose_pulling1_target(iX))
            ixTargetQi(iX)=j;
        end
    end
end

numIteration=run_config.dose_pulling1_start+1;

while(run_config.dose_pulling1 && numIteration<=100 && any(arrayfun(@(ixTarget,criteria,limit) qi(ixTarget).(criteria)<limit,ixTargetQi,run_config.dose_pulling1_criteria,run_config.dose_pulling1_limit)))

    [cst,optimization_flag] = matRad_pullDose(cst,1);

    if(optimization_flag)
        % optimize using the new objectives and penalties
        resultGUI_nominal = matRad_fluenceOptimization(dij,cst,pln);

        % Indicator calculation and show DVH and QI
        [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI_nominal,[],pln.numOfFractions,[],[],run_config.doseWindow_dvh);
        savefig([rootPath filesep 'dvh_nominal.fig']);

        % add resultGUI_nominal dose cubes to resultGUI structure to allow the visualization in the GUI
        resultGUI = resultGUI_nominal;

        % Create an interactive plot to slide through axial slices
        quantityMap=quantityOpt;
        plane      = 3;
        doseWindow = [0 max([max([resultGUI.physicalDose(:)*pln.numOfFractions]) run_config.doseWindow(2)])];
        doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
        f = figure;
        title([quantityMap]);
        set(gcf,'position',[10,10,550,400]);
        numSlices = ct.cubeDim(3);
        slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
        matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.(quantityMap)*pln.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
        b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
            'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
        b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.(quantityMap)*pln.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
        savefig(f,[rootPath filesep 'dose_' quantityMap '.fig']);

        disp(['!!--####################### ITERATION No. ' num2str(numIteration) ' (STEP 1) #######################--!!']);

        for iX=1:size(run_config.dose_pulling1_target,2)
            sprintf('%s (%s) = %5.3f %% ', run_config.dose_pulling1_criteria(iX),qi(ixTargetQi(iX)).name,qi(ixTargetQi(iX)).(run_config.dose_pulling1_criteria(iX)) )
        end

        % Print target objectives
        for  structure = 1:size(cst,1)
            display(cst{structure,2});
            for iX=1:length(cst{structure,6})
                display(cst{structure,6}{iX});
            end
        end

        numIteration=numIteration+1;

    else
        break;
    end

end

%% Plot nominal fluence
matRad_visSpotWeights(stf,resultGUI.w);
savefig([rootPath filesep 'fluence_nominal.fig']);

%% Plot dose distribution
%figure;
%matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([rootPath filesep 'dose3d_nominal.fig']);

%% Create target mask
%target_mask = zeros(size(resultGUI.physicalDose));
%target_mask(cst{ixCTV,4}{1,1}) = 1;

%% Create OAR mask
%OAR_mask = zeros(size(resultGUI.physicalDose));
%OAR_mask(cst{ixBody,4}{1,1}) = 1;
%OAR_mask(cst{ixCTV,4}{1,1}) = 0;

%% Plot target dose distribution
%figure;
%matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose.*target_mask*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([rootPath filesep 'target_dose3d_nominal.fig']);

%% Plot OAR dose distribution
%figure;
%matRad_geo3DWrapper(gca,ct,cst,resultGUI.physicalDose.*OAR_mask*pln.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%savefig([rootPath filesep 'OAR_dose3d_nominal.fig']);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
cst_robust=cst;

%% Load objectives
dp_start=[run_config_tmp.dose_pulling1_start run_config_tmp.dose_pulling2_start];
[cst_robust,ixTarget,p,ixBody,ixCTV,OARStructSel] = matRad_loadObjectives(run_config,run_config.plan_target,dp_start,cst_robust);

%% Copy OAR objectives from nominal dose pulling and relax them by scale factor%

for itOARStructure = 1:size(OARStructSel,2)
    for  itStructure = 1:size(cst_robust,1)
        if(strcmp(cst{itStructure,2},OARStructSel{itOARStructure}))
            cst_robust{itStructure,6}=cst{itStructure,6};
        end
    end
end

cst_robust = matRad_scaleDoseObjectives(cst_robust,OARStructSel,run_config.scale_factor);

%% Define ring objectives
cst_robust{ixRing1,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
cst_robust{ixRing1,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.10,0));
cst_robust{ixRing1,6}{1}.robustness  = 'none';
cst_robust{ixRing1,6}{1}.dosePulling  = false;
%OARStructSel{length(OARStructSel)+1}='RING 0 - 20 mm';

cst_robust{ixRing2,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
cst_robust{ixRing2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
cst_robust{ixRing2,6}{1}.robustness  = 'none';
cst_robust{ixRing2,6}{1}.dosePulling  = false;
%OARStructSel{length(OARStructSel)+1}='RING 20 - 50 mm';

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
if (~exist('dij_robust','var') || isempty(dij_robust)) && (~exist('dij_interval','var') || isempty(dij_interval))
    if run_config.radiationMode == "photons"
        dij_robust = matRad_calcPhotonDose(ct,stf_robust,pln_robust,cst_robust);
    end

    if run_config.radiationMode == "protons"
        dij_robust = matRad_calcParticleDose(ct,stf_robust,pln_robust,cst_robust);
    end
    dij_robust_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_robust_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '.mat'];
    save(dij_robust_file,'dij_robust', '-v7.3');
end
DCTime_robust = toc(now1);
time1=sprintf('DCTime_robust: %.2f\n',DCTime_robust); disp(time1);
results.performance.DCTime_robust=DCTime_robust;

%% Dose interval pre-computing

switch run_config.robustness
    case 'INTERVAL2'
        targetStructSel = {'CTV'};
        now2 = tic();
        if ~exist('dij_interval','var') || isempty(dij_interval)
            [dij_dummy, pln_dummy,dij_robust,pln_robust,dij_interval] = matRad_calcDoseInterval2b(ct,cst,stf_robust,pln_robust,dij_robust,targetStructSel,OARStructSel);
            dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval2.mat'];
            save(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval', '-v7.3');
        end
        dij_robust=dij_dummy;
        pln_robust=pln_dummy;
        IDCTime_robust = toc(now2);
        time2=sprintf('IDCTime_robust: %.2f\n',IDCTime_robust); disp(time2);
        results.performance.IDCTime_robust=IDCTime_robust;
    case 'INTERVAL3'
        targetStructSel = {'CTV'};
        now2 = tic();
        if ~exist('dij_interval','var') || isempty(dij_interval)
            [dij_dummy, pln_dummy,dij_robust,pln_robust,dij_interval] = matRad_calcDoseInterval3e(ct,cst,stf_robust,pln_robust,dij_robust,targetStructSel,OARStructSel,run_config.kdin,run_config.kmax,run_config.retentionThreshold);
            if(isequal(run_config.kdin,'dinamic'))
                dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_' num2str(run_config.retentionThreshold) '.mat'];
            elseif(isequal(run_config.kdin,'static'))
                dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_k_' num2str(run_config.kmax) '.mat'];
            end
            save(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval', '-v7.3');
        end
        dij_robust=dij_dummy;
        pln_robust=pln_dummy;
        IDCTime_robust = toc(now2);
        time2=sprintf('IDCTime_robust: %.2f\n',IDCTime_robust); disp(time2);
        results.performance.IDCTime_robust=IDCTime_robust;
end


%% Trigger robust optimization
% Make the objective to a composite worst case objective

dp_target_factor=1.0;
switch run_config.plan_objectives
    case '1'
        dp_target_factor=8.0;
    case '2'
        dp_target_factor=4.0;
    case '3'
        dp_target_factor=2.0;
    case '4'
        dp_target_factor=1.0;
    case '5'
        dp_target_factor=0.5;
end

% Target
switch run_config.robustness
    case 'STOCH'
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'STOCH';
        end

    case 'STOCH2'        
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'STOCH';
        end
        for iX=1:size(cst_robust,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst_robust{iX,2})
                    for k = 1:numel(cst_robust{iX,6})
                        cst_robust{iX,6}{k}.robustness  = 'STOCH';
                    end
                end
            end
        end

    case 'COWC'
        pln_robust.propOpt.useMaxApprox='logsumexp';
        
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'COWC';
        end

    case 'COWC2'
        pln_robust.propOpt.useMaxApprox='logsumexp';
        
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'COWC';
        end
        for iX=1:size(cst_robust,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst_robust{iX,2})
                    for k = 1:numel(cst_robust{iX,6})
                        cst_robust{iX,6}{k}.robustness  = 'COWC';
                    end
                end
            end
        end

    case 'c-COWC'
        
        pln_robust.propOpt.p1=run_config.p1;
        pln_robust.propOpt.p2=run_config.p2;
        pln_robust.propOpt.useMaxApprox='cheapCOWC';
        
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'COWC';
        end

    case 'c-COWC2'
        
        pln_robust.propOpt.p1=run_config.p1;
        pln_robust.propOpt.p2=run_config.p2;
        pln_robust.propOpt.useMaxApprox='cheapCOWC';
        
        for iX=1:length(cst_robust{ixTarget,6})
            cst_robust{ixTarget,6}{iX}.robustness  = 'COWC';
        end
        for iX=1:size(cst,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst{iX,2})
                    for k = 1:numel(cst{iX,6})
                        cst_robust{iX,6}{k}.robustness  = 'COWC';
                    end
                end
            end
        end

    case 'INTERVAL2'
        
        pln_robust.propOpt.dij_interval=dij_interval;
        clear('dij_interval');
        cst_robust{ixCTV,6}=[];
        cst_robust{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation2(30*dp_target_factor,p));
        cst_robust{ixCTV,6}{1}.robustness  = 'INTERVAL2';

        for iX=1:size(cst,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst{iX,2})
                    for k = 1:numel(cst{iX,6})
                        cst_robust{iX,6}{k}.robustness  = 'INTERVAL2';
                    end
                end
            end
        end
    
    case 'INTERVAL3'
        
        pln_robust.propOpt.dij_interval=dij_interval;
        clear('dij_interval');
        pln_robust.propOpt.theta2=run_config.theta2;        
        cst_robust{ixCTV,6}=[];
        cst_robust{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredBertoluzzaDeviation2(30*dp_target_factor,p));
        cst_robust{ixCTV,6}{1}.robustness  = 'INTERVAL3';
        
        for iX=1:size(cst,1)
            for j = 1:numel(OARStructSel)
                if strcmp(OARStructSel{j}, cst{iX,2})
                    for k = 1:numel(cst{iX,6})
                        cst_robust{iX,6}{k}.robustness  = 'INTERVAL3';
                    end
                end
            end
        end

end

%% Print target objectives
for  structure = 1:size(cst_robust,1)
    display(cst_robust{structure,2});
    for iX=1:length(cst_robust{structure,6})
        display(cst_robust{structure,6}{iX});
    end
end

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation
% treatment. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
profile on;
now3 = tic();
% Target

resultGUI_robust=cell(1, length(run_config.theta1));
for planIx = 1:num_plans
    tic
    switch run_config.robustness
        case {'INTERVAL2','INTERVAL3'}
            pln_robust.propOpt.theta1 = run_config.theta1(planIx);
    end
    resultGUI_robust{planIx} = matRad_fluenceOptimization(dij_robust,cst_robust,pln_robust);
    toc
end

%% Delete
switch run_config.robustness
    case {'INTERVAL2','INTERVAL3'}
        pln_robust.propOpt.dij_interval=[];
end

%% Indicator calculation and show DVH and QI

for planIx = 1:num_plans
    for scenIt = 1:pln_robust.multScen.totNumScen
        qi_robust{scenIt} = matRad_calcQualityIndicators(cst_robust,pln_robust,resultGUI_robust{planIx}.([pln_robust.bioParam.quantityVis '_' num2str(scenIt)])*pln_robust.numOfFractions,[],[]);
    end
    
    for iX=1:size(run_config.dose_pulling1_target,2)
        for scenIt = 1:pln_robust.multScen.totNumScen
            qiTarget{scenIt,iX}=qi_robust{scenIt}(ixTarget(iX)).([run_config.dose_pulling1_criteria(iX)]);
        end
        minQiTarget{iX} = min(cell2mat(qiTarget(:,iX)),[],1);
        meanQiTarget{iX} = sum(cell2mat(qiTarget(:,iX)).*multScen.scenProb);
        results.dosePulling.meanQiTarget{1,iX}=meanQiTarget{iX};
        results.dosePulling.minQiTarget{1,iX}=minQiTarget{iX};
    end
    
    
    % print results
    for iX=1:size(run_config.dose_pulling1_target,2)
        fprintf('meanQiTarget = %5.3f %% \n', meanQiTarget{end,iX}*100);
        fprintf('minQiTarget = %5.3f %% \n', minQiTarget{end,iX}*100);
    end
end

%% Dose pulling Step 2

numIteration=run_config.dose_pulling2_start+1;
criteria_array = eval(run_config.dose_pulling2_criteria);

for planIx = 1:num_plans
    switch run_config.robustness
        case {'INTERVAL2','INTERVAL3'}
            pln_robust.propOpt.theta1 = run_config.theta1(planIx);
    end
    while(run_config.dose_pulling2 && numIteration<=100 && any(arrayfun(@(criteria,limit) criteria<limit,cell2mat(criteria_array),run_config.dose_pulling2_limit)))
    
        [cst_robust,optimization_flag] = matRad_pullDose(cst_robust,2);
    
        if(optimization_flag)

            if(planIx>1)
                resultGUI_robust{planIx} = matRad_fluenceOptimization(dij_robust,cst_robust,pln_robust,resultGUI_robust{planIx-1}.w);
            else
                resultGUI_robust{planIx} = matRad_fluenceOptimization(dij_robust,cst_robust,pln_robust);
            end
        
            for scenIt = 1:pln_robust.multScen.totNumScen
                qi_robust{scenIt} = matRad_calcQualityIndicators(cst_robust,pln_robust,resultGUI_robust{planIx}.([pln_robust.bioParam.quantityVis '_' num2str(scenIt)])*pln_robust.numOfFractions,[],[]);
            end
            
            for iX=1:size(run_config.dose_pulling1_target,2)
                for scenIt = 1:pln_robust.multScen.totNumScen
                    qiTarget{scenIt,iX}=qi_robust{scenIt}(ixTarget(iX)).([run_config.dose_pulling1_criteria(iX)]);
                end
                meanQiTarget{iX}=sum(cell2mat( qiTarget(:,iX) ).*multScen.scenProb);
                minQiTarget{iX} = min(cell2mat( qiTarget(:,iX) ),[],1);
                results.dosePulling.meanQiTarget{end+1,iX}=meanQiTarget{iX};
                results.dosePulling.minQiTarget{end+1,iX}=minQiTarget{iX};
            end
            
            disp(['!!--####################### ITERATION No. ' num2str(numIteration) ' (STEP 2) #######################--!!']);
    
            for iX=1:size(run_config.dose_pulling1_target,2)
                fprintf('meanQiTarget = %5.3f %% \n', meanQiTarget{end,iX}*100);
                fprintf('minQiTarget = %5.3f %% \n', minQiTarget{end,iX}*100);
            end
    
            % Print target objectives
            for  structure = 1:size(cst_robust,1)
                display(cst_robust{structure,2});
                for iX=1:length(cst_robust{structure,6})
                    display(cst_robust{structure,6}{iX});
                end
            end
    
            numIteration=numIteration+1;
            criteria_array = eval(run_config.dose_pulling2_criteria);
    
        else
            break;
        end
    
    end
    
    % Indicator calculation and show DVH and QI
    [dvh_robust,qi_robust] = matRad_indicatorWrapper(cst,pln,resultGUI_robust{planIx},[],pln_robust.numOfFractions,[],[],run_config.doseWindow_dvh);
    savefig([folderPath{planIx} filesep 'dvh_robust.fig']);
    
    % add resultGUI_robust dose cubes to the existing resultGUI structure to allow the visualization in the GUI
    resultGUI = matRad_appendResultGUI(resultGUI,resultGUI_robust{planIx},0,'robust');
    OPTTime_robust = toc(now3);
    time3=sprintf('OPTTime_robust: %.2f\n',OPTTime_robust); disp(time3);
    results.performance.OPTTime_robust=OPTTime_robust;
    profiler = profile('info');
    save([folderPath{planIx} filesep 'profiler'],'profiler');

end

delete(gcp('nocreate'));

%% Plot robust fluence
for planIx = 1:num_plans
    matRad_visSpotWeights(stf_robust,resultGUI_robust{planIx}.w);
    savefig([folderPath{planIx} filesep 'fluence_robust.fig']);
end

%% Save outputs
save([rootPath filesep 'resultGUI.mat'],'resultGUI','resultGUI_robust');

%% Create an interactive plot to slide through axial slices
for planIx = 1:num_plans
%   for numScen=1:pln_robust.multScen.totNumScen
%
%       if(pln_robust.multScen.totNumScen>1)
%           quantityMap=[quantityOpt '_' num2str(round(numScen))];
%       else
%           quantityMap=quantityOpt;
%       end
%
%       plane      = 3;
%       doseWindow = [0 max([max([resultGUI_robust{planIx}.physicalDose(:)*pln_robust.numOfFractions]) run_config.doseWindow(2)])];
%       doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
%       f = figure;
%       title([quantityMap]);
%       set(gcf,'position',[10,10,550,400]);
%       numSlices = ct.cubeDim(3);
%       slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
%       matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust{planIx}.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);
%       b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
%           'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
%       b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust{planIx}.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]',[],'LineWidth',1.2);

%       savefig([folderPath{planIx} filesep 'dose_robust_' quantityMap '.fig']);
%   end

%   % Plot dose distribution
%   figure;
%   matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust{planIx}.physicalDose*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%   savefig([folderPath{planIx} filesep 'dose3d_robust.fig']);

%   % Create target mask
%   target_mask = zeros(size(resultGUI.physicalDose));
%   target_mask(cst_robust{ixCTV,4}{1,1}) = 1;

%   % Create OAR mask
%   OAR_mask = zeros(size(resultGUI.physicalDose));
%   OAR_mask(cst_robust{ixBody,4}{1,1}) = 1;
%   OAR_mask(cst_robust{ixCTV,4}{1,1}) = 0;

%   % Plot target dose distribution
%   figure;
%   matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust{planIx}.physicalDose.*target_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%   savefig([folderPath{planIx} filesep 'target_dose3d_robust.fig']);

%   % Plot OAR dose distribution
%   figure;
%   matRad_geo3DWrapper(gca,ct,cst_robust,resultGUI_robust{planIx}.physicalDose.*OAR_mask*pln_robust.numOfFractions,run_config.doseWindow,[0.002 0.00005],colorcube,[],'Dose [Gy]');
%   savefig([folderPath{planIx} filesep 'OAR_dose3d_robust.fig']);

end

%% check sampling option is activated
if ~run_config.sampling
    return;
end

if(run_config.plan_target=="PTV")
    cst{ixTarget,3}='OAR';
    cst_robust{ixTarget,3}='OAR';
end

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {};
[multScen] = matRad_multiScenGenerator(run_config.sampling_mode,run_config,'sampling',ct);

%% Perform sampling for robust optimization results
for planIx = 1:num_plans
    [caSampRob, mSampDoseRob, plnSampRob, resultGUIRobNomScen,resultGUIsampledScenRob] = matRad_sampling(ct,stf_robust,cst_robust,pln_robust,resultGUI_robust{planIx}.w,structSel,multScen);

    %% Perform sampling analysis
    phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
    robustnessCriteria = run_config.robustnessCriteria;
    GammaCriteria = run_config.GammaCriteria; 
    slice = round(isocenter(3)./ct.resolution.z);
    [cstStatRob, resultGUISampRob, metaRob, gammaFig, robustnessFig1, robustnessFig2] = matRad_samplingAnalysis(ct,cst_robust,plnSampRob,caSampRob, mSampDoseRob, resultGUIRobNomScen,phaseProb,'robustnessCriteria',robustnessCriteria,'GammaCriteria',GammaCriteria,'slice',slice);
    results.robustnessAnalysis_robust=resultGUISampRob.robustnessAnalysis;
    
    savefig(gammaFig, [folderPath{planIx} filesep 'gamma_analysis_robust.fig']);
    savefig(robustnessFig1, [folderPath{planIx} filesep 'robustness_analysis1_robust.fig']);
    savefig(robustnessFig2, [folderPath{planIx} filesep 'robustness_analysis2_robust.fig']);
    
    %%
    if (ct.numOfCtScen>1)
        f          = figure; title('individual scenarios'); camroll(90);
        plane      = 1;
        slice      = round(isocenter(2)./ct.resolution.y);
        numScen    = 1;
        matRad_plotSliceWrapper(gca,ct,cst_robust,numScen,[],plane,slice);
        b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
            'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
        b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,round(es.Value),[],plane,slice);
    end
    clear  numScen plane slice ans f b;
    
    %% 
    if (ct.numOfCtScen>1)
        f          = figure; title('individual scenarios'); %camroll(90);
        plane      = 3;
        slice      = round(isocenter(2)./ct.resolution.y);
        numScen    = 1;
        matRad_plotSliceWrapper(gca,ct,cst_robust,numScen, resultGUIRobNomScen.([quantityMap '_' int2str(numScen)])*pln_robust.numOfFractions,plane,slice);
        b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
            'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
        b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,round(es.Value),resultGUIRobNomScen.([quantityMap '_' int2str(round(es.Value))])*pln_robust.numOfFractions,plane,slice);
    end
    clear  numScen plane slice ans f b;
    
    %%
    if (ct.numOfCtScen>1)
        f          = figure; title('individual scenarios'); %camroll(90);
        plane      = 3;
        slice      = round(isocenter(2)./ct.resolution.y);
        numScen    = 1;
        matRad_plotSliceWrapper(gca,ct,cst_robust,numScen, resultGUIsampledScenRob.(quantityMap){numScen}*pln_robust.numOfFractions,plane,slice);
        b             = uicontrol('Parent',f,'Style','slider','Position',[50,5,419,23],...
            'value',numScen, 'min',1, 'max',ct.numOfCtScen,'SliderStep', [1/(ct.numOfCtScen-1) , 1/(ct.numOfCtScen-1)]);
        b.Callback    = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUIsampledScenRob.(quantityMap){round(es.Value)}*pln_robust.numOfFractions,plane,slice);
    end
    clear  numScen plane slice ans f b;
    
    %% Create an mean dose interactive plot to slide through axial slices
    quantityMap='meanCubeW';
    plane      = 3;
    doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p*1.25])];
    maxDose       = max(resultGUISampRob.(quantityMap)(:))*pln_robust.numOfFractions;
    doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
    meanDoseFig1 = figure;
    title([quantityMap 'for robust optimization results']);
    set(gcf,'position',[10,10,550,400]);
    set(gcf,'Color',[1 1 1]);
    numSlices = ct.cubeDim(3);
    
    slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
    matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Expected Dose [Gy]',[],'LineWidth',1.5);
    b = uicontrol('Parent',meanDoseFig1,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback = @(es,ed) matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Expected Dose [Gy]',[],'LineWidth',1.5);
    
    savefig(meanDoseFig1,[folderPath{planIx} filesep 'mean_dose_robust.fig']);
    
    %% Create an std dose interactive plot to slide through axial slices
    quantityMap='stdCubeW';
    plane      = 3;
    doseWindow = [0 max([max(resultGUISampRob.(quantityMap)(:)) p*0.5])];
    maxDose       = max(resultGUISampRob.(quantityMap)(:))*pln_robust.numOfFractions;
    doseIsoLevels = linspace(0.1 * maxDose,maxDose,10);
    stdDoseFig1 = figure;
    title([quantityMap 'for robust optimization results']);
    set(gcf,'position',[10,10,550,400]);
    set(gcf,'Color',[1 1 1]);
    numSlices = ct.cubeDim(3);
    
    slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
    matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
    b = uicontrol('Parent',stdDoseFig1,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback = @(es,ed) matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUISampRob.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
    
    savefig(stdDoseFig1,[folderPath{planIx} filesep 'std_dose_robust.fig']);
    
    %% Calculate nominal and expected dose difference
    quantityMap1='stdCube';
    quantityMap2='physicalDose';
    resultGUISampRob.diffCube=resultGUISampRob.(quantityMap1)-resultGUIRobNomScen.(quantityMap2);
    quantityMap1='stdCubeW';
    quantityMap2='physicalDose';
    resultGUISampRob.diffCubeW=resultGUISampRob.(quantityMap1)-resultGUIRobNomScen.(quantityMap2);
    
    %% Create an high and low doses expectation interactive plot to slide through axial slices
    quantityMap='diffCubeW';
    quantityMapRef='physicalDose';
    maxDose = 20.01; % [%]
    doseWindow = [-maxDose maxDose];
    plane      = 3;
    doseIsoLevels = linspace(-maxDose,maxDose,10);
    f = figure;
    title([quantityMap 'for nominal optimization results']);
    set(gcf,'position',[10,10,550,400]);
    numSlices = ct.cubeDim(3);
    
    mMap1=10;
    colormap1 = [linspace(0.20,1,mMap1)',linspace(0.20,1,mMap1)', linspace(1,1,mMap1)'];
    colormap2 = [linspace(1,1,mMap1)',linspace(1,0.20,mMap1)', linspace(1,0.20,mMap1)'];
    myColormap = [colormap1; colormap2];
    
    slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)./resultGUIRobNomScen.(quantityMapRef)*pln.numOfFractions,plane,slice,[],[],colorcube,myColormap,doseWindow,doseIsoLevels,[],'Relative Dose Difference [%]',[],'LineWidth',1.2);
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISampRob.(quantityMap)./resultGUIRobNomScen.(quantityMapRef)*pln.numOfFractions,plane,round(es.Value),[],[],colorcube,myColormap,doseWindow,doseIsoLevels,[],'Relative Dose Difference [%]',[],'LineWidth',1.2);
    
    savefig([folderPath{planIx} filesep 'dose_difference_robust.fig']);
    
    %% Multi-scenario dose volume histogram (DVH)
    f = figure;
    set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
    matRad_showDVHFromSampling(caSampRob,pln_robust.numOfFractions,cst_robust,plnSampRob,[1:plnSampRob.multScen.totNumScen],run_config.doseWindow_dvh,'multiscenario',1);
    title('Multi-scenario DVH for robust optimization results');
    
    savefig([folderPath{planIx} filesep 'dvh_multiscenario_robust.fig']);
    
    %% Trust band dose volume histogram (DVH)
    f = figure;
    set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
    matRad_showDVHFromSampling(caSampRob,pln_robust.numOfFractions,cst_robust,plnSampRob,[1:plnSampRob.multScen.totNumScen],run_config.doseWindow_dvh,'trustband',1);
    title('Trust band DVH for robust optimization results');
    
    savefig([folderPath{planIx} filesep 'dvh_trustband_robust.fig']);
    
    % print results
    disp('Robustness analysis (robust plan):');
    disp(results.robustnessAnalysis_robust);
    % Save outputs
    save([folderPath{planIx}  filesep 'results.mat'],'results');
    save([folderPath{planIx} filesep 'sampling.mat'],'caSampRob', 'mSampDoseRob', 'plnSampRob', 'resultGUIRobNomScen', 'resultGUISampRob','cstStatRob','metaRob','qi','qi_robust');

end

% print results
disp('Performance:');
disp(results.performance);
% save outputs 
save([rootPath filesep 'plan.mat'],'ct','cst','cst_robust','pln','pln_robust','stf','stf_robust','run_config');

diary off
