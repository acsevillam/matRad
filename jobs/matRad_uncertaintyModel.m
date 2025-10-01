function matRad_photonRobust2(radiationMode,description,varargin)

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

%% Parallel pool test
%matRad_testSLURMParpool();
%matRad_testSLURMParpool();
%ok = matRad_testSLURMParpool();
%if ~ok
%    error('Parallel pool test failed.');
%end

%% Set function parameters

validRadiationModes = {'photons','protons'};
validDescriptions = {'prostate','breast'};
validPatientIDs = {'1_mct','2_mct','3482','3648','3782','3790','3840','1758_mct','3477','3749','3832','3833','3929','4136','4136_mct','4155','4203','4357','4390','4428','4494','4531','4585','4585_mct','4681'};
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
defaultFilesPrefix = '';
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
addParameter(parser,'filesPrefix',defaultFilesPrefix,@(x) ischar(x));
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
run_config.filesPrefix = parser.Results.filesPrefix;
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
        dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID run_config.filesPrefix '_dij_interval2_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3))  '.mat'];
        if run_config.loadDij && isfile(dij_interval_file)
            load(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval');
        end
        dij_robust_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID run_config.filesPrefix '_dij_robust_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '.mat'];
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
        if(isequal(run_config.kdin,'dinamic'))
            if length(run_config.theta1)>1
                root_folder = [root_folder  filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
                for planIx = 1:num_plans
                    output_folder{planIx} = [root_folder filesep num2str(run_config.theta1(planIx)) filesep num2str(run_config.retentionThreshold) filesep num2str(run_config.theta2)];
                end
            else
                output_folder{1} = [root_folder filesep num2str(run_config.theta1(1)) filesep num2str(run_config.retentionThreshold) filesep num2str(run_config.theta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
                root_folder = output_folder{1};
            end
        elseif(isequal(run_config.kdin,'static'))
            if length(run_config.theta1)>1
                root_folder = [root_folder  filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
                for planIx = 1:num_plans
                    output_folder{planIx} = [root_folder filesep num2str(run_config.theta1(planIx)) filesep 'k_' num2str(run_config.kmax) filesep num2str(run_config.theta2)];
                end
            else
                output_folder{1} = [root_folder filesep num2str(run_config.theta1(1)) filesep 'k_' num2str(run_config.kmax) filesep num2str(run_config.theta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
                root_folder = output_folder{1};
            end
        end 
        if(isequal(run_config.kdin,'dinamic'))
            dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_' num2str(run_config.retentionThreshold) run_config.filesPrefix '.mat'];
        elseif(isequal(run_config.kdin,'static'))
            dij_interval_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID '_dij_interval3_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '_k_' num2str(run_config.kmax) run_config.filesPrefix '.mat'];
        end 
        if run_config.loadDij && isfile(dij_interval_file)
            load(dij_interval_file,'dij_dummy','pln_dummy','pln_robust','dij_interval');
        end
        dij_robust_file = [run_config.rootPath  filesep 'jobs' filesep 'images' filesep run_config.description filesep run_config.caseID run_config.filesPrefix '_dij_robust_' num2str(run_config.doseResolution(1)) '_' num2str(run_config.doseResolution(2)) '_' num2str(run_config.doseResolution(3)) '.mat'];
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
    metadata.dvfType = 'push';
    register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
    if(ct.numOfCtScen>numel(cst{1,4}))
        [ct,cst] = register.propContours();
    end
    metadata.dvfType = 'pull';
    register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
    [ct] = register.calcDVF();
    clear metadata;
end

%% Hiperparámetros para TPS-RMP
opts = struct( ...
    'fps', 600, ...
    'max_iters_outer', 70, ...
    'nsinkhorn', 6, ...
    'anneal_rate', 0.85, ...
    'sigma2_final', 1e-4, ...
    'lambda_init', 2e-2, ...
    'lambda_final', 1e-4, ...
    'lambda_anneal', 0.95, ...
    'outlier_scale', 1.0, ...
    'weight_anchors', 2.0 ...
);

%%
structIx  = 1;
Scell = cell(ct.numOfCtScen,1);
Xcell = cell(ct.numOfCtScen,1);
Vcell = cell(ct.numOfCtScen,1);
Fcell = cell(ct.numOfCtScen,1);
for k = 1:ct.numOfCtScen
    Scell{k} = cst_to_surface3d(cst, ct, k, structIx, opts.fps);
    % Malla completa (para SDF, visualización, etc.):
    figure;
    trisurf(Scell{k}.F, Scell{k}.V(:,1), Scell{k}.V(:,2), Scell{k}.V(:,3), ...
            'FaceColor',[0.2 0.6 1.0],'EdgeColor','none');
    axis equal; camlight; lighting gouraud;
    Xcell{k} = Scell{k}.PC;
    Vcell{k} = Scell{k}.V;
    Fcell{k} = Scell{k}.F;
end

Xavg = average_shape_per_patient1(Xcell, opts);

%%
[Xavg2, avg_raw] = postprocess_soft_average(Xavg, Xcell{1}, struct( ...
    'lambda_s', 1e-2, ...          % suavizado TPS
    'keep_percent', 10, ...         % descarta 10% de vértices menos confiables
    'alpha_shape', [] ...           % o un valor ~30-60 si quieres proyección
));

%% Calcular límites globales de todos los puntos
allPts = Xavg;  % empezamos con la media
for k = 1:ct.numOfCtScen
    allPts = [allPts; Xcell{k}]; %#ok<AGROW>
end

xmin = min(allPts(:,1)); xmax = max(allPts(:,1));
ymin = min(allPts(:,2)); ymax = max(allPts(:,2));
zmin = min(allPts(:,3)); zmax = max(allPts(:,3));

% Graficar cada escenario con la misma escala
for k = 1:ct.numOfCtScen
    figure; hold on; axis equal vis3d; grid on; view(3);
    [~, model_k] = tps_rpm_3d(Xcell{1}, Xcell{k}, opts);
    Xw = tps3d_warp_with_model(Xcell{1}, model_k);
    Xw2 = bsxfun(@minus, Xw, model_k.pre_shift);
    scatter3(Xw2(:,1), Xw2(:,2), Xw2(:,3), 6, 'blue', 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('x'); ylabel('y'); zlabel('z');
    legend({sprintf('Scan %d',k)});
    xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);
end

% Gráfica del promedio con la misma escala
figure; hold on; axis equal vis3d; grid on; view(3);
scatter3(Xavg(:,1), Xavg(:,2), Xavg(:,3), 6, 'red', 'filled', 'MarkerFaceAlpha',0.25);
xlabel('x'); ylabel('y'); zlabel('z');
legend({'Promedio'});
xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);

% Gráfica del promedio con la misma escala
figure; hold on; axis equal vis3d; grid on; view(3);
scatter3(Xavg2(:,1), Xavg2(:,2), Xavg2(:,3), 6, 'red', 'filled', 'MarkerFaceAlpha',0.25);
xlabel('x'); ylabel('y'); zlabel('z');
legend({'Promedio'});
xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);

diary off
