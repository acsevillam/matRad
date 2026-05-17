function [ct, cst] = matRad_mergeMultipleDicomCtScenarios(ct, tmpCt, tmpCst)
% matRad merges cropped DICOM CT scenarios into multi-scenario ct and cst structs
%
% call
%   [ct,cst] = matRad_mergeMultipleDicomCtScenarios(ct,tmpCt,tmpCst)
%
% input
%   ct:        partially initialized multi-scenario ct struct
%   tmpCt:     cell array with cropped single-scenario ct structs
%   tmpCst:    cell array with cropped single-scenario cst structs
%
% output
%   ct:        multi-scenario ct struct
%   cst:       multi-scenario cst struct
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

matRadCfg = MatRad_Config.instance();

ct.x = matRad_getCommonCtField(tmpCt, 'x', matRadCfg);
ct.y = matRad_getCommonCtField(tmpCt, 'y', matRadCfg);
ct.z = matRad_getCommonCtField(tmpCt, 'z', matRadCfg);
ct.cubeDim = matRad_getCommonCtField(tmpCt, 'cubeDim', matRadCfg);

cst = tmpCst{1};
numOfStruct = size(cst, 1);

for structure = 1:numOfStruct
    cst{structure, 4} = cell(1, ct.numOfCtScen);
end

for ctPhase = 1:ct.numOfCtScen
    ct.cubeHU{ctPhase} = tmpCt{ctPhase}.cubeHU{1};
    ct.dicomInfo(ctPhase) = tmpCt{ctPhase}.dicomInfo;
    ct.dicomMeta(ctPhase) = tmpCt{ctPhase}.dicomMeta;

    for structure = 1:numOfStruct
        cst{structure, 4}{ctPhase} = tmpCst{ctPhase}{structure, 4}{1};
    end
end
end

function value = matRad_getCommonCtField(ctScenarios, fieldName, matRadCfg)
value = ctScenarios{1}.(fieldName);

for ctPhase = 2:numel(ctScenarios)
    if ~isequal(ctScenarios{ctPhase}.(fieldName), value)
        matRadCfg.dispError('There are scenarios with incompatible %s values.', fieldName);
    end
end
end
