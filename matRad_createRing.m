function [cst,ixRing] = matRad_createRing(ixBase,ixLimit,cst,ct,vOuterMargin,vInnerMargin,metadata)
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

VOIBase=cst{ixBase,4}{1,1};
VOILimit=cst{ixLimit,4}{1,1};

geo_tmp=zeros(ct.cubeDim);
geo_tmp(VOIBase)=1;

% Add margin to base volume
VOIEnlargedOuter=find(matRad_addMargin(geo_tmp,cst,ct.resolution, vOuterMargin));
VOIEnlargedInner=find(matRad_addMargin(geo_tmp,cst,ct.resolution, vInnerMargin));

% Delete voxels from base volume
VOIRing=setdiff(VOIEnlargedOuter,VOIEnlargedInner);

% Delete voxels outside limit volume
VOIRing=intersect(VOIRing,VOILimit);

ixRing=size(cst,1)+1;

cst{ixRing,1}=cst{end,1}+1;
cst{ixRing,2}=metadata.name;
cst{ixRing,3}=metadata.type;
cst{ixRing,4}{1,1}=VOIRing;
cst{ixRing,5}=cst{ixBase,5};
cst{ixRing,5}.visibleColor=metadata.visibleColor;

end

