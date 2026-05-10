function dij = matRad_buildSyntheticDijForScenarioResolution(templateDij,pln,scenarioDijIx,matRad_cfg,calculationName)
% matRad_buildSyntheticDijForScenarioResolution builds metadata-only robust dij
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

if nargin < 5 || isempty(calculationName)
    calculationName = 'streaming scenario dose';
end

dij = templateDij;
if ~isfield(templateDij,'doseGrid') || ~isfield(templateDij.doseGrid,'numOfVoxels') || ...
        ~isfield(templateDij,'totalNumOfBixels')
    matRad_cfg.dispError('%s template dij is missing doseGrid or totalNumOfBixels metadata.', ...
        calculationName);
end

containerSize = pln.multScen.getDijContainerSize();
fields = fieldnames(templateDij);
for f = 1:numel(fields)
    fieldName = fields{f};
    if ~iscell(templateDij.(fieldName))
        continue;
    end

    sample = getFirstNonEmptyCell(templateDij.(fieldName));
    if isempty(sample) || ~isnumeric(sample) || ~ismatrix(sample)
        continue;
    end

    cells = cell(containerSize);
    for ix = scenarioDijIx(:)'
        cells{ix} = spalloc(templateDij.doseGrid.numOfVoxels, ...
            templateDij.totalNumOfBixels,0);
    end
    dij.(fieldName) = cells;
end
end

function value = getFirstNonEmptyCell(cells)
value = [];
firstIx = find(~cellfun(@isempty,cells(:)),1,'first');
if ~isempty(firstIx)
    value = cells{firstIx};
end
end
