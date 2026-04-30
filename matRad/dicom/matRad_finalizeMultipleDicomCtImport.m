function [ct,cst] = matRad_finalizeMultipleDicomCtImport(ct,tmpCtOriginal,tmpCstOriginal)
% matRad helper to crop and merge imported DICOM CT scenarios
%
% call
%   [ct,cst] = matRad_finalizeMultipleDicomCtImport(ct,tmpCtOriginal,tmpCstOriginal)
%
% input
%   ct:              partially initialized multi-scenario ct struct
%   tmpCtOriginal:   cell array with one imported single-scenario ct struct
%                    per CT scenario
%   tmpCstOriginal:  cell array with one imported single-scenario cst cell
%                    array per CT scenario
%
% output
%   ct:              matRad ct multi-scenario struct cropped to the common
%                    overlap
%   cst:             matRad cst multi-scenario struct cropped to the common
%                    overlap
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

matRad_cfg = MatRad_Config.instance();

try
    [tmpCt, ~, cropIdx] = matRad_cropCtToOverlap(tmpCtOriginal,1e-5);
catch ME
    matRad_cfg.dispError('Unable to crop CT scenarios to their common overlap: %s',ME.message);
end

try
    tmpCst = matRad_cropCstToOverlap(tmpCstOriginal,cropIdx,tmpCtOriginal,tmpCt);
catch ME
    matRad_cfg.dispError('Unable to crop structure masks to the common overlap: %s',ME.message);
end

if ~matRad_hasCompatibleCstMetadata(tmpCst)
    matRad_cfg.dispError(['The imported CT scenarios contain incompatible structure metadata. ' ...
        'Select a single compatible RTSTRUCT or harmonize structure names, types, properties, and objectives before importing.']);
end

[ct,cst] = matRad_mergeMultipleDicomCtScenarios(ct,tmpCt,tmpCst);

end
