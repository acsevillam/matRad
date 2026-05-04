% Configuration script for uncertainty sampling.
%
% This template assumes that ct, cst, stf, pln, and resultGUI are available
% in the MATLAB workspace. Use the PatientMatFile name-value argument in
% matRad_calcStudy for batch runs from a saved patient file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% examine only specific structures? improves calculation time
examineStructures = {}; % e.g. examinedStructures = {'CTV', 'OAR'};

if exist('ct','var')
    multScen = matRad_createScenarioModel(ct,'impScen');
else
    multScen = matRad_createScenarioModel([],'impScen');
end

% Define active public uncertainty dimensions.
multScen.scenarioDimensionActive = {'ct','setup','range'};

% Define gridded setup and range scenario counts.
multScen.numOfSetupGridPoints = 3;
multScen.numOfRangeGridPoints = 9;

% Define scenario combination behavior.
multScen.combinations = 'none'; % 'none', 'shift', or 'all'
multScen.combineRange = true;   % combine absolute and relative range shifts

% Define uncertainty standard deviations.
multScen.shiftSD    = [1.8 1.8 1.8]; % [mm]
multScen.rangeAbsSD = 1.8;           % [mm]
multScen.rangeRelSD = 0;             % [%]

%% path for output pdf and mat
param.outputPath = pwd;

%% report parameters
param.operator = 'matrad';

% default set of percentiles for scenario analysis; uncomment to override
% param.percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
% default criteria for gamma analysis between nominal and mean dose
% param.criteria = [3 3]; %%% [X % Y mm]

%% start calculation
matRad_calcStudy(multScen, ...
    'SelectStructures',examineStructures, ...
    'OutputPath',param.outputPath, ...
    'OperatorName',param.operator);

% exit;
