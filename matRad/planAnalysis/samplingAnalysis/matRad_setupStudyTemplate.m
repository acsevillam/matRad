% configuration script for uncertainty sampling
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

multScen = matRad_ScenarioModel.create('impScen', []);

% a) define setup uncertainty scenarios
multScen.numOfSetupGridPoints = 1;        % only the nominal setup scenario
multScen.combinations         = 'none';   % independent setup and range realizations

% b) define range error scenarios
multScen.numOfRangeGridPoints = 30;       % number of range grid points
multScen.combineRange         = true;     % combine absolute and relative range shifts

% define standard deviation of normal distribution - important for probabilistic treatment planning
multScen.rangeRelSD           = 4.5;             % given in [%]
multScen.rangeAbsSD           = 2;               % given in [mm]
multScen.shiftSD              = [4.5 4.5 4.5];   % given in [mm]

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
                 'SelectStructures', examineStructures, ...
                 'OutputPath', param.outputPath, ...
                 'OperatorName', param.operator);

% exit;
