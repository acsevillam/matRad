function matRad_photonDosePulling(radiationMode,description,varargin)

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
validPatientIDs = {'3482','3648','3782','3790','3840','3477','3749','3832','3833','3929','1758'};
validAcquisitionTypes = {'mat','dicom'};
validPlanObjectives = {'1','2','3','4','5','6'};
validDosePullingTargets = {'CTV','PTV'};
validDosePullingCriteria = {'COV_95','COV_98','COV_99','COV1'};
validPlanTargets = {'CTV','PTV'};
validPlanBeams = {'5F','7F','9F'};
validRobustness = {'none'};
validScenModes = {'nomScen','wcScen','impScen5','impScen7','impScen_permuted5','impScen_permuted7','impScen_permuted_truncated5','impScen_permuted_truncated7','random','random_truncated'};

defaultPatientID = '3482';
defaultAcquisitionType = 'dicom';
defaultPlanObjective = '4';
defaultDosePulling = false;
defaultDosePullingTarget = 'CTV';
defaultDosePullingCriteria = 'COV1';
defaultDosePullingLimit = 0.98;
defaultDosePullingStart = 0;
defaultPlanTarget = 'CTV';
defaultPlanBeams = '9F';
defaultShiftSD = [5 10 5]; % mm
defaultRobustness = 'none';
defaultScenMode = 'nomScen';
defaultSampling = true;
defaultSamplingMode = 'impScen_permuted_truncated5';
defaultSamplingWCFactor = 1.5;
defaultRootPath = matRad_cfg.matRadRoot;
defaultNCores = feature('numcores');

parser = inputParser;

addRequired(parser,'radiationMode',@(x) any(validatestring(x,validRadiationModes)));
addRequired(parser,'description',@(x) any(validatestring(x,validDescriptions)));
addParameter(parser,'caseID',defaultPatientID,@(x) any(validatestring(x,validPatientIDs)));
addParameter(parser,'AcquisitionType',defaultAcquisitionType,@(x) any(validatestring(x,validAcquisitionTypes)));
addParameter(parser,'plan_objectives',defaultPlanObjective,@(x) any(validatestring(x,validPlanObjectives)));
addParameter(parser,'dose_pulling',defaultDosePulling,@islogical);
addOptional(parser,'dose_pulling_target',defaultDosePullingTarget,@(x) numel(x) >= 1 && all(ismember(x,validDosePullingTargets)));
addOptional(parser,'dose_pulling_criteria',defaultDosePullingCriteria,@(x) numel(x) >= 1 && all(ismember(x,validDosePullingCriteria)));
addOptional(parser,'dose_pulling_limit',defaultDosePullingLimit,@(x) numel(x) >= 1 && isnumeric(x) && all(x > 0));
addOptional(parser,'dose_pulling_start',defaultDosePullingStart,@(x) validateattributes(x,{'numeric'},...
            {'nonempty','integer','nonnegative'}));
addParameter(parser,'plan_target',defaultPlanTarget,@(x) any(validatestring(x,validPlanTargets)));
addParameter(parser,'plan_beams',defaultPlanBeams,@(x) any(validatestring(x,validPlanBeams)));
addParameter(parser,'shiftSD',defaultShiftSD,@(x) numel(x) == 3 && isnumeric(x) && all(x > 0));
addParameter(parser,'robustness',defaultRobustness,@(x) any(validatestring(x,validRobustness)));
addParameter(parser,'scen_mode',defaultScenMode,@(x) any(validatestring(x,validScenModes)));
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
run_config.AcquisitionType = parser.Results.AcquisitionType;
run_config.plan_objectives = parser.Results.plan_objectives;
run_config.dose_pulling = parser.Results.dose_pulling;
run_config.dose_pulling_target = parser.Results.dose_pulling_target;
run_config.dose_pulling_criteria = parser.Results.dose_pulling_criteria;
run_config.dose_pulling_limit = parser.Results.dose_pulling_limit;
run_config.dose_pulling_start = parser.Results.dose_pulling_start;
run_config.plan_target = parser.Results.plan_target;
run_config.plan_beams = parser.Results.plan_beams;
run_config.shiftSD = parser.Results.shiftSD;
run_config.robustness = parser.Results.robustness;
run_config.scen_mode = parser.Results.scen_mode;
run_config.sampling = parser.Results.sampling;
run_config.sampling_mode = parser.Results.sampling_mode;
run_config.sampling_wcFactor = parser.Results.sampling_wcFactor;
run_config.rootPath = parser.Results.rootPath;

output_folder = ['output' filesep run_config.radiationMode filesep run_config.description filesep run_config.caseID filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep num2str(run_config.shiftSD(1)) 'x' num2str(run_config.shiftSD(2)) 'x' num2str(run_config.shiftSD(3)) filesep run_config.scen_mode filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];

run_config.resolution = [3 3 3];
run_config.doseResolution = [3 3 3];
run_config.GammaCriteria = [3 3];
run_config.robustnessCriteria = [5 5];
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

%% Import CT and rename structures
[ct,cst] = matRad_loadGeometry(run_config);
cst = matRad_renameStructures(cst,run_config);

%% Print run config
disp(run_config);

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
% define optimization parameter for both VOIs
dp_start=[run_config.dose_pulling_start 0];
[cst,ixTarget,p,ixBody,~,~] = matRad_loadObjectives(run_config,run_config.plan_target,dp_start,cst);

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
savefig([folderPath filesep 'dvh_nominal.fig']);

%% Dose pulling

ixTargetQi=zeros(size(run_config.dose_pulling_target));
for i=1:size(run_config.dose_pulling_target,2)
    for j=1:size(qi,2)
        if strcmp(qi(j).name,run_config.dose_pulling_target(i))
            ixTargetQi(i)=j;
        end
    end
end

numIteration=run_config.dose_pulling_start+1;

while(run_config.dose_pulling && numIteration<=100 && any(arrayfun(@(ixTarget,criteria,limit) qi(ixTarget).([criteria])<limit,ixTargetQi,run_config.dose_pulling_criteria,run_config.dose_pulling_limit)))

    [cst,optimization_flag] = matRad_pullDose(cst,1);

    if(optimization_flag)
        % optimize using the new objectives and penalties
        resultGUI_nominal = matRad_fluenceOptimization(dij,cst,pln);

        % Indicator calculation and show DVH and QI
        [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI_nominal,[],pln.numOfFractions,[],[],run_config.doseWindow_dvh);
        savefig([folderPath filesep 'dvh_nominal.fig']);

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
        savefig(f,[folderPath filesep 'dose_' quantityMap '.fig']);

        disp(['!!--####################### ITERATION No. ' num2str(numIteration) ' #######################--!!']);

        for i=1:size(run_config.dose_pulling_target,2)
            sprintf('%s (%s) = %5.3f %% ', run_config.dose_pulling_criteria(i),qi(ixTargetQi(i)).name,qi(ixTargetQi(i)).(run_config.dose_pulling_criteria(i)) )
        end

        % Print target objectives
        for  structure = 1:size(cst,1)
            display(cst{structure,2});
            for i=1:length(cst{structure,6})
                display(cst{structure,6}{i});
            end
        end

        numIteration=numIteration+1;

    else
        break;
    end

end

%% Plot nominal fluence
matRad_visSpotWeights(stf,resultGUI.w);
savefig([folderPath filesep 'fluence_nominal.fig']);

%% check sampling option is activated
if ~run_config.sampling
    return;
end

if(run_config.plan_target=="PTV")
    cst{ixTarget,3}='OAR';
end

%% Define sampling parameters
% select structures to include in sampling; leave empty to sample dose for all structures
% sampling does not know on which scenario sampling should be performed
structSel = {};
[multScen] = matRad_multiScenGenerator(run_config.sampling_mode,run_config,'sampling',ct);

%% Perform sampling for nominal optimization results
delete(gcp('nocreate'));
[caSamp, mSampDose, plnSamp, resultGUInomScen,~] = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen);

%% Perform sampling analysis
phaseProb = ones(1,ct.numOfCtScen)/ct.numOfCtScen;
robustnessCriteria = run_config.robustnessCriteria;
GammaCriteria = run_config.GammaCriteria;
slice = round(isocenter(3)./ct.resolution.z);
[cstStat, resultGUISamp, meta, gammaFig, robustnessFig1, robustnessFig2] = matRad_samplingAnalysis(ct,cst,plnSamp,caSamp, mSampDose, resultGUInomScen,phaseProb,'robustnessCriteria',robustnessCriteria,'GammaCriteria',GammaCriteria,'slice',slice);
results.robustnessAnalysis_nominal=resultGUISamp.robustnessAnalysis;

savefig(gammaFig, [folderPath filesep 'gamma_analysis_nominal.fig']);
savefig(robustnessFig1, [folderPath filesep 'robustness_analysis1_nominal.fig']);
savefig(robustnessFig2, [folderPath filesep 'robustness_analysis2_nominal.fig']);

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
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

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
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);
b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
    'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst,1,resultGUISamp.(quantityMap)*pln.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose uncertainty [Gy]',[],'LineWidth',1.2);

savefig([folderPath filesep 'std_dose_nominal.fig']);

%% Multi-scenario dose volume histogram (DVH)
f = figure;
set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSamp,pln.numOfFractions,cst,plnSamp,[1:plnSamp.multScen.totNumScen],run_config.doseWindow_dvh,'multiscenario',1);
title('Multi-scenario DVH for nominal optimization results');

savefig([folderPath filesep 'dvh_multiscenario_nominal.fig']);

%% Trust band dose volume histogram (DVH)
f = figure;
set(gcf,'Color',[1 1 1],'position',[10,10,600,400]);
matRad_showDVHFromSampling(caSamp,pln.numOfFractions,cst,plnSamp,[1:plnSamp.multScen.totNumScen],run_config.doseWindow_dvh,'trustband',1);
title('Trust band DVH for nominal optimization results');

savefig([folderPath filesep 'dvh_trustband_nominal.fig']);

%% print results
disp('Robustness analysis (nominal plan):');
disp(results.robustnessAnalysis_nominal);

%% Save outputs
save([folderPath filesep 'resultGUI.mat'],'resultGUI');
save([folderPath filesep 'plan.mat'],'ct','cst','pln','stf','run_config');
save([folderPath filesep 'sampling.mat'],'caSamp', 'mSampDose', 'plnSamp', 'resultGUInomScen','cstStat','resultGUISamp','meta','dvh');
save([folderPath filesep 'results.mat'],'results');

%%
diary off

%%
delete(gcp('nocreate'));