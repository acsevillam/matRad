function [cst,ixSkin] = matRad_createSkin(ixPTV,ixCTV,cst,metadata)
% matRad add margin function
% 
% call
%   mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
%
% input
%   mVOI:           image stack in dimensions of X x Y x Z holding ones for
%                   object and zeros otherwise 
%   cst:            matRad cst struct
%   vResolution     ct resolution
%   vMargin:        margin in mm 
%   bDiaElem        if true 26-connectivity is used otherwise 6-connectivity
%
% output
%   mVOIEnlarged:   enlarged VOI
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete voxels from base volume
VOIPTV=cst{ixPTV,4}{1,1};
VOICTV=cst{ixCTV,4}{1,1};

VOISkin=setdiff(VOICTV,VOIPTV);

% Delete voxels outside limit volume
%VOIRing=intersect(VOIRing,VOILimit);

ixSkin=size(cst,1)+1;

cst{ixSkin,1}=cst{end,1}+1;
cst{ixSkin,2}=metadata.name;
cst{ixSkin,3}=metadata.type;
cst{ixSkin,4}{1,1}=VOISkin;
cst{ixSkin,5}=cst{ixPTV,5};
cst{ixSkin,5}.visibleColor=metadata.visibleColor;

end

