function matRad_processCDVH(metadata)
    % Initialize environment and log level
    clc; close all;
    matRad_rc;
    param.logLevel = 1;

    % Build the output folder dynamically depending on robustness strategy
    output_folder = buildOutputFolder(metadata);
    foldername = fullfile(metadata.defaultRootPath, metadata.jobPath, ...
        metadata.job_folder, output_folder);

    % List folder content and find latest subfolder (typically timestamped)
    listing = dir(foldername);
    if isempty(listing)
        error('No subdirectories found in: %s', foldername);
    end
    lastJobFolder = fullfile(foldername, listing(end).name);

    % Define expected file paths
    filename1 = fullfile(lastJobFolder, ['dvh_trustband_' metadata.robustness_approach '.fig']);
    filename2 = fullfile(lastJobFolder, ['dvh_' metadata.robustness_approach '.fig']);
    filename3 = fullfile(lastJobFolder, 'results.mat');
    filename4 = fullfile(lastJobFolder, 'plan.mat');

    % Validate all files exist
    requiredFiles = {filename1, filename2, filename3, filename4};
    for i = 1:length(requiredFiles)
        if ~isfile(requiredFiles{i})
            error('File not found: %s', requiredFiles{i});
        end
    end

    % Load figures and result data
    fig1 = openfig(filename1);
    fig2 = openfig(filename2);
    load(filename3,'results');
    load(filename4,'run_config');

    % Start logging to a .log file
    diaryname = fullfile(lastJobFolder, 'results.log');
    if exist(diaryname, 'file') == 2
        delete(diaryname);
    end
    diary(diaryname);
    diary on;

    % Print metadata and configuration for traceability
    fprintf('Description: \t %s \n', run_config.description);
    fprintf('CaseID: \t %s \n', run_config.caseID);
    fprintf('Resolution: \t [%.2f %.2f %.2f] \n', run_config.resolution);
    fprintf('Objectives: \t %s \n', run_config.plan_objectives);
    fprintf('Target: \t %s \n', run_config.plan_target);
    fprintf('Beam setup: \t %s \n', run_config.plan_beams);
    fprintf('shiftSD: \t [%.2f %.2f %.2f] \n', run_config.shiftSD);
    fprintf('Robustness: \t %s \n', run_config.robustness);
    fprintf('Scen mode: \t %s \n', run_config.scen_mode);
    if isfield(run_config, 'wcFactor'), fprintf('wcFactor: \t %.2f \n', run_config.wcFactor); end
    if isfield(run_config, 'beta1'), fprintf('beta1: \t %.2f \n', run_config.beta1); end
    if isfield(run_config, 'beta2'), fprintf('beta2: \t %.2f \n', run_config.beta2); end
    if isfield(run_config, 'theta1'), fprintf('theta1: \t %.2f \n', run_config.theta1); end
    if isfield(run_config, 'kdin'), fprintf('kdin: \t %s \n', run_config.kdin); end
    if isfield(run_config, 'kmax'), fprintf('kmax: \t %.2f \n', run_config.kmax); end
    if isfield(run_config, 'theta2'), fprintf('theta2: \t %.2f \n', run_config.theta2); end
    fprintf('Samp. mode: \t %s \n', run_config.sampling_mode);
    fprintf('Samp. wcFactor:  %.2f \n\n\n', run_config.sampling_wcFactor);

    % Analyze selected DVH structures
    dvhSpecs = getDVHStructureSpecs(metadata.description);
    for i = 1:size(dvhSpecs, 1)
        structName = dvhSpecs{i, 1};
        metric = dvhSpecs{i, 2};
        analyzeStructureDVH(fig1, fig2, structName, metric);
    end

    analyzeCTV(fig1, results, metadata.robustness_approach);

    % Stop logging
    diary off;
end

function output_folder = buildOutputFolder(metadata)
    % Returns the appropriate output folder path depending on the robustness model
    switch metadata.robustness
        case {'COWC', 'COWC2', 'c-COWC', 'c-COWC2', 'INTERVAL2'}
            % Advanced robust optimization: includes shiftSD and scenario details
            output_folder = fullfile('output', metadata.radiationMode, metadata.description, ...
                metadata.caseID, metadata.robustness, metadata.plan_target, ...
                metadata.plan_beams, metadata.plan_objectives, metadata.shiftSD, ...
                metadata.scen_mode, num2str(metadata.wcFactor), ...
                [num2str(metadata.beta1) '_to_' num2str(metadata.beta2)]);
        case {'INTERVAL3'}
            % Advanced robust optimization: includes shiftSD and scenario details
            output_folder = fullfile('output', metadata.radiationMode, metadata.description, ...
                metadata.caseID, metadata.robustness, metadata.plan_target, ...
                metadata.plan_beams, metadata.plan_objectives, metadata.shiftSD, ...
                metadata.scen_mode, num2str(metadata.wcFactor), ...
                num2str(metadata.theta1),num2str(metadata.retentionThreshold),num2str(metadata.theta2));
        case 'none'
            % Nominal plan or PTV-based plan with simpler folder structure
            output_folder = fullfile('output', metadata.radiationMode, metadata.description, ...
                metadata.caseID, metadata.robustness, metadata.plan_target, ...
                metadata.plan_beams, metadata.plan_objectives);
        otherwise
            error('Unsupported robustness type: %s', metadata.robustness);
    end
end

function dvhSpecs = getDVHStructureSpecs(description)
% Returns a list of structures to analyze with metric type and reference value

switch lower(description)
    case 'breast'
        dvhSpecs = {
            'LEFT LUNG', struct('type', 'Vx', 'value', 20);
            'HEART',     struct('type', 'Dmean', 'value', []);
            'SKIN',      struct('type', 'Dx', 'value', 5)
        };
    case 'prostate'
        dvhSpecs = {
            'BLADDER', struct('type', 'Vx', 'value', 60);
            'RECTUM',  struct('type', 'Vx', 'value', 40)
        };
    otherwise
        warning('No DVH specs defined for "%s".', description);
        dvhSpecs = {};
end
end


function analyzeStructureDVH(fig1, fig2, structName, metric)
% Analyze a DVH structure using either Vx, Dx or Dmean metric
% Inputs:
% - structName: structure name
% - metric: struct with fields 'type' and 'value'

fprintf('>>> %s <<<\n', upper(structName));

switch metric.type
    case 'Vx'
        % Analyze volume receiving at least x Gy
        h = findobj(fig2, 'Type', 'line', 'DisplayName', structName);
        if isempty(h), warning('Structure "%s" not found.', structName); return; end
        x = h(1).XData;
        y = h(1).YData;
        Vx = interp1(x, y, metric.value, 'linear', 'extrap');
        fprintf('%s: V%.1f = %.2f %%\n', structName, metric.value, Vx);
        str = sprintf('V%.1f (%s) = %.2f %%', metric.value, structName, Vx);

    case 'Dx'
        % Analyze dose to x% of volume
        h = findobj(fig2, 'Type', 'line', 'DisplayName', structName);
        if isempty(h), warning('Structure "%s" not found.', structName); return; end
        x = h(1).XData;
        y = h(1).YData;
        [y_unique, iy] = unique(y);
        x_unique = x(iy);
        Dx = interp1(y_unique, x_unique, metric.value, 'linear', 'extrap');
        fprintf('%s: D%.1f%% = %.2f Gy\n', structName, metric.value, Dx);
        str = sprintf('D%.0f%% (%s) = %.2f Gy', metric.value, structName, Dx);

    case 'Dmean'
        % Compute mean dose
        h = findobj(fig2, 'Type', 'line', 'DisplayName', structName);
        if isempty(h), warning('Structure "%s" not found.', structName); return; end
        x = h(1).XData;
        y = h(1).YData;
        dx = x(2) - x(1);
        d = -diff(y) / dx;
        x_d = x(1:end-1);
        f = @(z) interp1(x_d, d, z, 'linear', 'extrap');
        mean_dose = x_d * f(x_d)' / sum(f(x_d));
        fprintf('%s: Dmean = %.2f Gy\n', structName, mean_dose);
        str = sprintf('Dmean (%s) = %.2f Gy', structName, mean_dose);

    otherwise
        warning('Unknown metric type: %s', metric.type);
        return;
end

% Add annotation to figure
dim = [0.3 0.2 + 0.1*rand() 0.3 0.1]; % Randomized vertical position
annotation(fig2, 'textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 8, ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');
end

function analyzeCTV(fig1, results, robustness_approach)
% Analyze CTV robustness metrics and annotate the c-DVH plot.
%
% Inputs:
% - fig1: Handle to c-DVH figure
% - results: Struct loaded from results.mat
% - robustness_approach: string key for accessing robustness data

% Robustness index
ri = results.(['robustnessAnalysis_' robustness_approach]).robustnessIndex1;
fprintf('CTV: Robustness Index (RI) = %.4f\n', ri);

% Find the CTV confidence band minimum line
h_patch = findobj(fig1, 'Type', 'patch', 'DisplayName', 'CTV');
if isempty(h_patch)
    warning('CTV confidence band not found.');
    return;
end

% Extract lower envelope of the confidence band
x_min = h_patch(1).XData(1:1000);
y_min = h_patch(1).YData(1:1000);
f_min = @(z) interp1(x_min, y_min, z, 'linear', 'extrap');

% Evaluate coverage at prescription dose (e.g., 42.56 Gy)
x0 = 42.56;
coverage = f_min(x0);
fprintf('CTV: V%.2f min = %.2f %%\n', x0, coverage);
fprintf('CTV: RCvI = %.4f\n', coverage / 100);

% Annotate RI on figure
dim = [0.35 0.5 0.3 0.1];
str = sprintf('RI = %.4f', ri);
annotation(fig1, 'textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 8, ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');
end
