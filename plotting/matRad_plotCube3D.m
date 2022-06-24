function hvol3d = matRad_plotCube3D(axesHandle,doseCube,cMap,window,alpha)
% matRad function that plots isolines in 3d
%
% call
%   hpatch = matRad_plotIsoDose3D(axesHandle,xMesh,yMesh,zMesh,doseCube,isoLevels,cMap,window,alpha)
%
% input
%   axesHandle  handle to axes the slice should be displayed in
%   x/y/zMesh   meshs
%   doseCube    dose cube
%   isoLevels   levels for computation
%   cMap        colormap
%   window      window for dose display
%   alpha       transparency
%
% output
%   hpatch: handle to the patch object
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

matRad_cfg = MatRad_Config.instance();

if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.001;
end

if ~exist('window','var') || isempty(window)
    window = [min(doseCube(:)) max(doseCube(:))];
end

if nargin < 3
    cMap = jet(64);
end

colormap(axesHandle,cMap);
hvol3d = vol3d('Cdata',doseCube,'Parent',axesHandle,'alpha',alpha);

end

