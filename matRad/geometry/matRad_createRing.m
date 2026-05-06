function [cst,ixRing] = matRad_createRing(ixBase,ixLimit,cst,ct,vOuterMargin,vInnerMargin,metadata)
% matRad function to create an outer-minus-inner isotropic ring VOI clipped to a limiting VOI
%
% call
%   [cst,ixRing] = matRad_createRing(ixBase,ixLimit,cst,ct,vOuterMargin,vInnerMargin,metadata)
%
% input
%   ixBase:        row index of the base VOI in the cst struct
%   ixLimit:       row index of the limiting VOI in the cst struct
%   cst:           matRad cst struct
%   ct:            matRad ct struct
%   vOuterMargin:  outer margin in mm, with fields x, y and z
%   vInnerMargin:  inner margin in mm, with fields x, y and z
%   metadata:      struct with fields name, type and visibleColor for the
%                  created ring VOI
%
% output
%   cst:           updated matRad cst struct
%   ixRing:        row index of the created ring VOI
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

VOIRing = cell(1, ct.numOfCtScen);
useDiagonalConnectivity = false;

for scen_iterator = 1:ct.numOfCtScen
    VOIBase = cst{ixBase,4}{scen_iterator};
    VOILimit = cst{ixLimit,4}{scen_iterator};

    geo_tmp = zeros(ct.cubeDim);
    geo_tmp(VOIBase) = 1;

    % Add margin to base volume
    VOIEnlargedOuter = find(matRad_addMargin(geo_tmp,cst,ct.resolution,vOuterMargin,useDiagonalConnectivity));
    VOIEnlargedInner = find(matRad_addMargin(geo_tmp,cst,ct.resolution,vInnerMargin,useDiagonalConnectivity));

    % Delete voxels from base volume
    VOIRing{scen_iterator} = setdiff(VOIEnlargedOuter,VOIEnlargedInner);

    % Delete voxels outside limit volume
    VOIRing{scen_iterator} = intersect(VOIRing{scen_iterator},VOILimit);
end

ixRing = size(cst,1) + 1;

cst{ixRing,1} = cst{end,1} + 1;
cst{ixRing,2} = metadata.name;
cst{ixRing,3} = metadata.type;
cst{ixRing,4} = VOIRing;
cst{ixRing,5} = cst{ixBase,5};
cst{ixRing,5}.visibleColor = metadata.visibleColor;

end
