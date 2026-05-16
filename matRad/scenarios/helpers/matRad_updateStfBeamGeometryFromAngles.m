function stf = matRad_updateStfBeamGeometryFromAngles(stf)
% matRad_updateStfBeamGeometryFromAngles refreshes world geometry from BEV.
%
% The helper recomputes source points, ray positions, target points, and
% photon beamlet corners after gantry/couch angles have changed. BEV fields
% are treated as the geometry source of truth.
%
% call
%   stf = matRad_updateStfBeamGeometryFromAngles(stf)
%
% input
%   stf:                matRad steering information struct with gantryAngle,
%                       couchAngle, and BEV ray geometry fields
%
% output
%   stf:                steering information with refreshed world-coordinate
%                       geometry fields
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

for i = 1:numel(stf)
    rotMatVectorsT = transpose(matRad_getRotationMatrix(stf(i).gantryAngle, stf(i).couchAngle));

    if isfield(stf(i), 'sourcePoint_bev')
        stf(i).sourcePoint = stf(i).sourcePoint_bev * rotMatVectorsT;
    end

    if ~isfield(stf(i), 'ray')
        continue
    end

    for j = 1:numel(stf(i).ray)
        if isfield(stf(i).ray(j), 'rayPos_bev')
            stf(i).ray(j).rayPos = stf(i).ray(j).rayPos_bev * rotMatVectorsT;
        end

        if isfield(stf(i).ray(j), 'targetPoint_bev')
            stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev * rotMatVectorsT;
        end

        if isfield(stf(i), 'bixelWidth') && isfield(stf(i).ray(j), 'rayPos_bev') && ...
                isfield(stf(i).ray(j), 'beamletCornersAtIso')
            rayPos = stf(i).ray(j).rayPos_bev;
            stf(i).ray(j).beamletCornersAtIso = [ ...
                                                 rayPos + [+stf(i).bixelWidth / 2, 0, +stf(i).bixelWidth / 2]; ...
                                                 rayPos + [-stf(i).bixelWidth / 2, 0, +stf(i).bixelWidth / 2]; ...
                                                 rayPos + [-stf(i).bixelWidth / 2, 0, -stf(i).bixelWidth / 2]; ...
                                                 rayPos + [+stf(i).bixelWidth / 2, 0, -stf(i).bixelWidth / 2]] * rotMatVectorsT;
        end
    end
end

end
