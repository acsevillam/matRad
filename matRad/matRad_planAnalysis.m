function resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln,varargin)
% matRad plan analysis function
% This function performs analysis on radiation therapy plans, including DVH (Dose-Volume Histogram) and quality indicators.
% It optionally displays these analyses based on input parameters.
%
% input
%   resultGUI:              matRad resultGUI struct containing the analysis results
%   ct:                     matRad ct struct with computed tomography data
%   cst:                    matRad cst cell array with structure definitions
%   stf:                    matRad stf struct with beam information
%   pln:                    matRad pln struct with plan information
%   name / value pairs:     Optional parameters for analysis customization
%   quantity: (optional)    resultGUI dose quantity to analyse
%   displayDoseMode:(optional) 'perFraction' or 'total' for figures/tables
%   doseWindow: (optional)  dose axis window for DVH display
%   refGy: (optional)       Per-fraction dose values for V_XGy calculation
%   refVol:(optional)       Volume percentages for D_X calculation (default: [2 5 95 98])
%
% output
%   resultGUI:              Updated resultGUI with analysis data

% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize input parser for optional parameters
p = inputParser();

% Define required inputs again for clarity
p.addRequired('ct',@isstruct);
p.addRequired('cst',@iscell);
p.addRequired('stf',@isstruct);
p.addRequired('pln',@isstruct);

% Define optional parameters with default values
p.addParameter('refGy',[],@isnumeric); % Per-fraction reference dose values for V_XGy calculation
p.addParameter('refVol',[2 5 95 98],@isnumeric); % Reference volume percentages for D_X calculation
p.addParameter('showDVH',true,@islogical); % Flag to show or hide the DVH plot
p.addParameter('showQI',true,@islogical); % Flag to show or hide the Quality Indicators plot
p.addParameter('quantity','',@(x) ischar(x) || isstring(x)); % Dose quantity to analyze
p.addParameter('displayDoseMode','perFraction',@(x) ischar(x) || isstring(x)); % Display dose mode
p.addParameter('doseWindow',[],@(x) isempty(x) || (isnumeric(x) && numel(x) == 2)); % DVH dose axis window

% Parse input arguments to extract values
p.parse(ct,cst,stf,pln,varargin{:});

% Assign parsed values to variables
ct = p.Results.ct;
cst = p.Results.cst;
stf = p.Results.stf;
pln = p.Results.pln;
refGy = p.Results.refGy;
refVol = p.Results.refVol;
showDVH = p.Results.showDVH;
showQI = p.Results.showQI;
quantity = p.Results.quantity;
doseWindow = p.Results.doseWindow;
displayDoseMode = p.Results.displayDoseMode;

if isstring(quantity)
    quantity = char(quantity);
end

[displayDoseScale,displayDoseMode] = matRad_getDisplayDoseScale(pln,displayDoseMode);

if ~isempty(doseWindow)
    doseWindow = doseWindow(:)';
end

% Determine which dose cube to use based on resultGUI structure
if ~isempty(quantity)
    visQ = quantity;
elseif isfield(pln,'bioParam') && isobject(pln.bioParam) && isprop(pln.bioParam,'quantityVis')
    visQ = pln.bioParam.quantityVis;
elseif isfield(pln,'bioParam') && isstruct(pln.bioParam) && isfield(pln.bioParam,'quantityVis')
    visQ = pln.bioParam.quantityVis;
elseif isfield(resultGUI,'RBExD')
    visQ = 'RBExD';
else
    visQ = 'physicalDose';
end

if ~isfield(resultGUI,visQ)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Unknown quantity ''%s'' to analyse!',visQ);
end

doseCube = resultGUI.(visQ);

% Calculate DVH and quality indicators
resultGUI.dvh = matRad_calcDVH(cst,doseCube,'cum'); % Calculate cumulative DVH
resultGUI.qi  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol); % Calculate quality indicators
resultGUI.analysisQuantity = visQ;
resultGUI.analysisDoseMode = 'perFraction';
resultGUI.displayDoseMode = displayDoseMode;
resultGUI.displayDoseScale = displayDoseScale;
resultGUI.displayDvh = scaleDvhForDisplay(resultGUI.dvh,displayDoseScale);
resultGUI.displayQi = scaleQiForDisplay(resultGUI.qi,displayDoseScale);

dvhScen = {};
if isfield(pln,'multScen') && pln.multScen.totNumScen > 1
    for i = 1:pln.multScen.totNumScen
        scenFieldName = sprintf('%s_scen%d',visQ,i);
        if isfield(resultGUI,scenFieldName)
            dvhScen{i} = scaleDvhForDisplay( ...
                matRad_calcDVH(cst,resultGUI.(scenFieldName),'cum'), ...
                displayDoseScale); % Calculate cumulative scenario DVH
        end
    end
end


if showDVH || showQI
    % Configuration for GUI appearance
    matRad_cfg = MatRad_Config.instance();

    % Create figure for plots with background color from configuration
    hF = figure('Color',matRad_cfg.gui.backgroundColor);

    colorSpec = {'Color',matRad_cfg.gui.elementColor,...
            'XColor',matRad_cfg.gui.textColor,...
            'YColor',matRad_cfg.gui.textColor,...
            'GridColor',matRad_cfg.gui.textColor,...
            'MinorGridColor',matRad_cfg.gui.backgroundColor};

    % Determine subplot layout based on flags
    if showDVH && showQI
        hDVHax = subplot(2,1,1,colorSpec{:}); % DVH plot area
        hQIax = subplot(2,1,2,colorSpec{:}); % Quality Indicators plot area
    elseif showDVH
        hDVHax = subplot(1,1,1,colorSpec{:}); % Only DVH plot
    elseif showQI
        hQIax = subplot(1,1,1,colorSpec{:}); % Only Quality Indicators plot
    end
end

% Display DVH if enabled
if showDVH
    matRad_showDVH(resultGUI.displayDvh,cst,pln,'axesHandle',hDVHax,'LineWidth',3); % Show DVH plot

    for i = 1:numel(dvhScen)
        matRad_showDVH(dvhScen{i},cst,pln,'axesHandle',hDVHax,'LineWidth',0.5,'plotLegend',false,'LineStyle','--'); % Show DVH plot
    end

    if ~isempty(doseWindow)
        xlim(hDVHax,doseWindow);
    end
end

% Display Quality Indicators if enabled
if showQI
    matRad_showQualityIndicators(hQIax,resultGUI.displayQi); % Show Quality Indicators plot
end

end

function dvh = scaleDvhForDisplay(dvh,scale)
if scale == 1 || isempty(dvh)
    return;
end

for i = 1:numel(dvh)
    if isfield(dvh(i),'doseGrid') && ~isempty(dvh(i).doseGrid)
        dvh(i).doseGrid = dvh(i).doseGrid * scale;
    end
end
end

function qi = scaleQiForDisplay(qi,scale)
if scale == 1 || isempty(qi)
    return;
end

doseFields = {'mean','std','max','min','referenceDose'};
for i = 1:numel(qi)
    for j = 1:numel(doseFields)
        fieldName = doseFields{j};
        if isfield(qi(i),fieldName) && isnumeric(qi(i).(fieldName))
            qi(i).(fieldName) = qi(i).(fieldName) * scale;
        end
    end

    fields = fieldnames(qi(i));
    for j = 1:numel(fields)
        fieldName = fields{j};
        if strncmp(fieldName,'D_',2) && isnumeric(qi(i).(fieldName))
            qi(i).(fieldName) = qi(i).(fieldName) * scale;
        end
    end
end
end
