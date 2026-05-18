function samplingDoseFig = matRad_plotSamplingDoseCubeAnalysis(doseAnalysis, doseCube, ct, cst, slice, varargin)
% matRad_plotSamplingDoseCubeAnalysis plots a sampling dose analysis cube
%
% call
%   samplingDoseFig = matRad_plotSamplingDoseCubeAnalysis(doseAnalysis,doseCube,ct,cst,slice)
%
% input
%   doseAnalysis:    analysis metadata struct from matRad_samplingAnalysis
%   doseCube:        per-fraction dose cube
%   ct:              matRad ct struct
%   cst:             matRad cst cell array
%   slice:           CT slice to plot
%   varargin:        optional Name/Value pairs:
%                    - 'axesHandle': axes used for plotting
%                    - 'plane': matRad plane index
%
% output
%   samplingDoseFig: figure handle when a new figure is created
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5 || isempty(slice)
    samplingDoseFig = [];
    return
end

p = inputParser;
p.CaseSensitive = false;
p.addParameter('axesHandle', [], @(axesHandle) isempty(axesHandle) || ...
               strcmp(get(axesHandle, 'type'), 'axes'));
p.addParameter('plane', 3, @(plane) isnumeric(plane) && isscalar(plane) && ...
               any(plane == [1 2 3]));
parse(p, varargin{:});

axesHandle = p.Results.axesHandle;
if isempty(axesHandle)
    samplingDoseFig = figure('Name', doseAnalysis.title);
    set(samplingDoseFig, 'Color', [1 1 1]);
    axesHandle = axes('Parent', samplingDoseFig);
else
    samplingDoseFig = ancestor(axesHandle, 'figure');
end

displayCube = doseCube .* doseAnalysis.displayScale;
refScen = matRad_getSamplingDoseReferenceScenario(ct);

matRad_plotSlice(ct, 'axesHandle', axesHandle, 'cst', cst, ...
                 'cubeIdx', refScen, 'dose', displayCube, ...
                 'plane', p.Results.plane, 'slice', slice, ...
                 'doseWindow', doseAnalysis.displayDoseWindow, ...
                 'colorBarLabel', doseAnalysis.colorBarLabel, ...
                 'LineWidth', 1.2);
title(axesHandle, doseAnalysis.title);

end

function refScen = matRad_getSamplingDoseReferenceScenario(ct)
if isstruct(ct) && isfield(ct, 'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end
end
