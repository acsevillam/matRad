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
    [lastJobFolder,lastDataFolder] = getFinalJobFolder(foldername, metadata);

    % Validate all required files exist
    requiredFiles = {
        fullfile(lastJobFolder, 'plan.mat');
        fullfile(lastDataFolder, ['dvh_trustband_' metadata.robustness_approach '.fig']);
        fullfile(lastDataFolder, ['dvh_' metadata.robustness_approach '.fig']);
        fullfile(lastDataFolder, 'results.mat');
    };
    
    % Check existence of each file
    missingFiles = requiredFiles(~cellfun(@isfile, requiredFiles));
    
    if ~isempty(missingFiles)
        fprintf('\n[ERROR] The following required files were not found:\n');
        for i = 1:length(missingFiles)
            fprintf('  • %s\n', missingFiles{i});
        end
        error('Missing required DVH-related files. Execution stopped.');
    end

    % Load figures and result data
    load(requiredFiles{1},'run_config');
    fig1 = openfig(requiredFiles{2});
    fig2 = openfig(requiredFiles{3});
    load(requiredFiles{4},'results');

    % Start logging to a .log file
    diaryname = fullfile(lastDataFolder, 'results.log');
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
    if isfield(run_config, 'theta1'), fprintf('theta1: \t %.2f \n', metadata.theta1); end
    if isfield(run_config, 'kdin'), fprintf('kdin: \t %s \n', run_config.kdin); end
    if isfield(run_config, 'kmax'), fprintf('kmax: \t %.2f \n', run_config.kmax); end
    if isfield(run_config, 'theta2'), fprintf('theta2: \t %.2f \n', run_config.theta2); end
    fprintf('Samp. mode: \t %s \n', run_config.sampling_mode);
    fprintf('Samp. wcFactor:  %.2f \n\n\n', run_config.sampling_wcFactor);

    % Custimize DVH
    customizeDVHGraphics(fig1,fig2);
   
    dvhSpecs = getDVHStructureSpecs(metadata.description);
    
    % --- Nominal DVH Section ---
    fprintf('\n\n');
    fprintf('!!!------------- Nominal DVH -------------!!!\n');
    fprintf('\n');

    % Analyze selected DVH structures
    for i = 1:size(dvhSpecs, 1)
        structName = dvhSpecs{i, 1};
        metric = dvhSpecs{i, 3};
        analyzeNominalDVH(fig2, structName, metric);
    end

    % --- Robust Confidence DVH Section ---
    fprintf('\n\n');
    fprintf('!!!------------- c-DVH -------------!!!\n');
    fprintf('\n');

    % Analyze selected DVH structures
    for i = 1:size(dvhSpecs, 1)
        structName = dvhSpecs{i, 1};
        structLabel = dvhSpecs{i, 2};
        metric = dvhSpecs{i, 3};
        dim = dvhSpecs{i, 4};
        analyzeCDVH(fig1, structName, structLabel, metric, dim);
    end    

    analyzeCTVCDVH(fig1, results, metadata.robustness_approach);

    % Stop logging
    diary off;

    metadata.title=getTitleFromRobustness(metadata);

    cleanAndExportFigures(fig1, fig2, lastDataFolder, metadata);

    close all;

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
                metadata.scen_mode, num2str(metadata.wcFactor));
        case 'none'
            % Nominal plan or PTV-based plan with simpler folder structure
            output_folder = fullfile('output', metadata.radiationMode, metadata.description, ...
                metadata.caseID, metadata.robustness, metadata.plan_target, ...
                metadata.plan_beams, metadata.plan_objectives);
        otherwise
            error('Unsupported robustness type: %s', metadata.robustness);
    end
end

function [lastJobFolder,lastDataFolder] = getFinalJobFolder(foldername_base, metadata)
% Returns the final job folder depending on robustness strategy
% Applies optional suffixes for INTERVAL3 models

% List subdirectories
listing = dir(foldername_base);
listing = listing([listing.isdir] & ~startsWith({listing.name}, '.'));

if isempty(listing)
    error('No subdirectories found in: %s', foldername_base);
end

if isfield(metadata, 'selected_folder')
    selectedFolder = metadata.selected_folder;
else
    selectedFolder = listing(end).name;
end

switch metadata.robustness
    case {'INTERVAL3'}
        % Extract suffix components safely
        if isfield(metadata, 'theta1')
            theta1 = metadata.theta1;
        else
            theta1 = 0;
        end

        if isfield(metadata, 'retentionThreshold')
            rT = metadata.retentionThreshold;
        else
            rT = 0;
        end

        if isfield(metadata, 'theta2')
            theta2 = metadata.theta2;
        else
            theta2 = 0;
        end
        suffix = fullfile(num2str(theta1), num2str(rT), num2str(theta2));

        lastJobFolder = fullfile(foldername_base, selectedFolder);
        lastDataFolder = fullfile(foldername_base, selectedFolder, suffix);

    case {'none','COWC', 'COWC2', 'c-COWC', 'c-COWC2', 'INTERVAL2'}
        lastJobFolder = fullfile(foldername_base, selectedFolder);
        lastDataFolder = fullfile(foldername_base, selectedFolder);

    otherwise
        error('Unsupported robustness strategy: %s', metadata.robustness);
end
end

function dvhSpecs = getDVHStructureSpecs(description)
% Returns a list of structures to analyze with metric type and reference value

switch lower(description)
    case 'breast'
        dvhSpecs = {
            'LEFT LUNG', 'Left Lung', struct('type', 'Vx', 'value', 20), [0.15 0.325 .3 .3];
            'HEART', 'Heart', struct('type', 'Dmean', 'value', []), [0.2 0.2 .3 .3] ;
            'SKIN', 'Skin', struct('type', 'Dx', 'value', 5) [0.66 0.13 .3 .3], 
        };
    case 'prostate'
        dvhSpecs = {
            'BLADDER', 'Bladder', struct('type', 'Vx', 'value', 60, [0.5 0.1 .3 .3]);
            'RECTUM', 'Rectum',  struct('type', 'Vx', 'value', 40, [0.4 0.2 .3 .3])
        };
    otherwise
        warning('No DVH specs defined for "%s".', description);
        dvhSpecs = {};
end
end

function analyzeNominalDVH(fig2, structName, metric)
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
        fprintf('%s: DMEAN = %.2f Gy\n', structName, mean_dose);

    otherwise
        warning('Unknown metric type: %s', metric.type);
        return;
end

end

function analyzeCDVH(fig1, structName, structLabel, metric, dim)
% Analyze robust c-DVH for a given structure using central, lower, and upper bands
% Inputs:
% - fig1: Figure handle to the c-DVH plot
% - structName: structure name (string)
% - metric: struct with fields 'type' and 'value' (e.g. D5 → value = 5)

fprintf('>>> %s <<<\n', upper(structName));

% Get central curve
h_central = findobj(fig1, 'Type', 'line', 'LineStyle', ':', 'DisplayName', structName);
if isempty(h_central)
    warning('Structure "%s" central curve not found.', structName);
    return;
end
x_central = h_central(1).XData;
y_central = h_central(1).YData;

% Get confidence band
h_patch = findobj(fig1, 'Type', 'patch', 'DisplayName', structName);
if isempty(h_patch)
    warning('Structure "%s" patch not found.', structName);
    return;
end

% Validate expected number of points
if length(h_patch(1).XData) < 3200
    warning('Patch data for "%s" does not contain enough points.', structName);
    return;
end

% Extract lower and upper envelope
x_min = h_patch(1).XData(1:1000);
y_min = h_patch(1).YData(1:1000);
x_max = h_patch(1).XData(1601:3200);
y_max = h_patch(1).YData(1601:3200);

% Compute based on metric type
switch metric.type
    case 'Vx'
        x0 = metric.value;

        f_central = @(z) interp1(x_central, y_central, z, 'linear', 'extrap');
        f_min = @(z) interp1(x_min, y_min, z, 'linear', 'extrap');
        f_max = @(z) interp1(x_max, y_max, z, 'linear', 'extrap');

        val_c = f_central(x0);
        val_l = f_min(x0);
        val_u = f_max(x0);

        fprintf('%s: V%.1f = %.2f [ %.2f - %.2f ] %%\n', structName, x0, val_c, val_l, val_u);
        str = sprintf('V%.0f (%s) = %.2f [ %.2f - %.2f ] %%', x0, structLabel, val_c, val_l, val_u);

    case 'Dx'
        y0 = metric.value;

        % Interpolation from DVH to Dose for each band
        [y_u, iy] = unique(y_central); x_u = x_central(iy);
        f_central = @(z) interp1(y_u, x_u, z, 'linear', 'extrap');
        [y_u, iy] = unique(y_min);     x_u = x_min(iy);
        f_min = @(z) interp1(y_u, x_u, z, 'linear', 'extrap');
        [y_u, iy] = unique(y_max);     x_u = x_max(iy);
        f_max = @(z) interp1(y_u, x_u, z, 'linear', 'extrap');

        val_c = f_central(y0);
        val_l = f_min(y0);
        val_u = f_max(y0);

        fprintf('%s: D%.1f%% = %.2f [ %.2f - %.2f ] Gy\n', structName, y0, val_c, val_l, val_u);

        % Build formatted string
        label = sprintf('D%.0f%% (%s)', y0, structLabel);
        value = sprintf('%.2f [ %.2f - %.2f ] Gy', val_c, val_l, val_u);
        
        % Combine them with a line break if label is too long
        if strlength(value) > 15  % umbral ajustable según tu layout
            str = sprintf('%s =\n%s', label, value);  % dos líneas
        else
            str = sprintf('%s = %s', label, value);  % una línea
        end

    case 'Dmean'
        % Central
        dx = x_central(2) - x_central(1);
        d = -diff(y_central) / dx;
        x_d = x_central(1:end-1);
        Dmean_c = sum(x_d .* d) / sum(d);

        % Lower
        dx = x_min(2) - x_min(1);
        d = -diff(y_min) / dx;
        x_d = x_min(1:end-1);
        Dmean_l = sum(x_d .* d) / sum(d);

        % Upper
        dx = x_max(2) - x_max(1);
        d = -diff(y_max) / dx;
        x_d = x_max(1:end-1);
        Dmean_u = sum(x_d .* d) / sum(d);

        fprintf('%s: DMEAN = %.2f [ %.2f - %.2f ] Gy\n', structName, Dmean_c, Dmean_l, Dmean_u);
        str = sprintf('DMEAN (%s) = %.2f [ %.2f - %.2f ] Gy', structLabel, Dmean_c, Dmean_l, Dmean_u);

    otherwise
        warning('Unknown metric type: %s', metric.type);
        return;
end

% Annotate c-DVH figure
annotation(fig1, 'textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 8, ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');
end


function analyzeCTVCDVH(fig1, results, robustness_approach)
% Analyze CTV robustness metrics and annotate the c-DVH plot.
%
% Inputs:
% - fig1: Handle to c-DVH figure
% - results: Struct loaded from results.mat
% - robustness_approach: string key for accessing robustness data

fprintf('>>> CTV <<<\n');

% Robustness index
ri = results.(['robustnessAnalysis_' robustness_approach]).robustnessIndex1;
fprintf('CTV: Robustness Index (RI) = %.4f\n', ri);

% Find the CTV confidence band minimum line
findobj(fig1,'LineStyle',':','Type','line','DisplayName','CTV');
findobj(fig1,'LineStyle','-','Type','patch','DisplayName','CTV');

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
dim = [0.3 0.5 .3 .3];
str = sprintf('RI = %.4f', ri);
annotation(fig1, 'textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 8, ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');
end

function titleStr = getTitleFromRobustness(metadata)
% Returns a formatted title string based on the robustness strategy
%
% Inputs:
%   metadata - struct with field 'robustness' and optionally 'p2'
%
% Output:
%   titleStr - formatted string to use as plot title

switch metadata.robustness
    case 'none'
        switch metadata.plan_target
            case 'CTV'
                titleStr = 'Nominal';
            case 'PTV'
                titleStr = 'PTV (Margin-based)';
        end
    case {'COWC', 'COWC2'}
        titleStr = sprintf('Minimax');
    case {'c-COWC', 'c-COWC2'}
        if isfield(metadata, 'p2')
            titleStr = sprintf('Cheap-Minimax (%d of 13)', metadata.p2);
        else
            titleStr = sprintf('Cheap-Minimax');
        end
    case 'INTERVAL2'
        if isfield(metadata, 'theta1')
            titleStr = sprintf('Interval2 (theta = %d)', metadata.theta1);
        else
            titleStr = sprintf('Interval2');
        end
    case 'INTERVAL3'
        if isfield(metadata, 'theta1')
            titleStr = sprintf('Interval3 (theta = %d)', metadata.theta1);
        else
            titleStr = sprintf('Interval3');
        end
    otherwise
        titleStr = sprintf('%s', metadata.robustness);
end
end


function cleanAndExportFigures(fig1, fig2, datafolder, metadata)
% Cleans up specified DVH figures and exports them to PDF
%
% Inputs:
%   fig1                - Figure handle to trustband c-DVH
%   fig2                - Figure handle to nominal DVH
%   foldername          - Path to job folder

    % List of structures to delete from both figures
    delStructs = {'PTV', 'RING 0 - 20 mm', 'RING 20 - 50 mm', ...
                  'BULB', 'CTV Min', 'CTV Max', ...
                  'CONTRALATERAL BREAST', 'SPINAL CORD'};

    % Remove overlays from fig1 and fig2
    for i = 1:numel(delStructs)
        delete(findobj(fig1, 'DisplayName', delStructs{i}));
        delete(findobj(fig2, 'DisplayName', delStructs{i}));
    end

    % Format fig1
    fig1.Children(1).FontSize = 9;
    fig1.Position = [120 100 540 230];
    fig1.Children(2).Title.String = metadata.title;

    % Export fig1
    set(fig1, 'PaperOrientation', 'landscape');
    set(fig1, 'PaperPositionMode', 'auto');
    set(fig1, 'PaperSize', [5.8 4.0]);
    filename1 = fullfile(datafolder, ...
        sprintf('dvh_trustband_%s_%s.pdf', metadata.robustness_approach, metadata.plan_target));
    print(fig1, filename1, '-dpdf', '-r0', '-fillpage');

end

function customizeDVHGraphics(fig1, fig2)
% Applies custom styling and label adjustments to DVH figures
%
% - Assigns display names to confidence bands
% - Applies specific color styling (e.g. for SKIN)
% - Designed to be optional and safe if elements are not found

    % === Patch label copying for central lines ===
    h_lines = findobj(fig1, 'LineStyle', ':');
    h_patches = findobj(fig1, 'FaceAlpha', 0.2);

    for i = 1:min(numel(h_lines), numel(h_patches))
        h_patches(i).DisplayName = h_lines(i).DisplayName;
    end

    % === Custom color styling for SKIN ===
    skinColor = [0 0 1];  % Blue
    
    % Set color for SKIN in various line types
    applyColor(fig1, 'SKIN', '-', 'line', skinColor);
    applyColor(fig1, 'SKIN', ':', 'line', skinColor);
    applyColor(fig2, 'SKIN', '-', 'line', skinColor);
    
    % Set color for SKIN patch
    applyColor(fig1, 'SKIN', '-', 'patch', skinColor);
end

function applyColor(fig, structName, lineStyle, objType, color)
% Utility to apply color styling to matching graphic objects
    h = findobj(fig, 'LineStyle', lineStyle, 'Type', objType, 'DisplayName', structName);
    for i = 1:numel(h)
        try
            if strcmp(objType, 'patch')
                h(i).EdgeColor = color;
                h(i).FaceColor = color;
            else
                h(i).Color = color;
            end
        end
    end
end
