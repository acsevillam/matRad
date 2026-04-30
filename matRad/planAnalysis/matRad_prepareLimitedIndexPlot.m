function [plotCube,doseColorMap] = matRad_prepareLimitedIndexPlot(plotCube,doseWindow,doseColorMap)
% matRad_prepareLimitedIndexPlot Prepare an index cube with an upper plot limit.
%
% call
%   [plotCube,doseColorMap] = matRad_prepareLimitedIndexPlot(plotCube,doseWindow,doseColorMap)
%
% input
%   plotCube:      cube to be shown as colorwash
%   doseWindow:    display range, values >= doseWindow(2) are clipped
%   doseColorMap:  colormap used for display
%
% output
%   plotCube:      display cube with values above the upper limit clipped
%   doseColorMap:  display colormap with the last color set to dark red if
%                  values exceed the upper limit
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

if nargin < 3
    doseColorMap = [];
end

if nargin < 2 || isempty(doseWindow) || numel(doseWindow) < 2 || ~isfinite(doseWindow(2))
    return;
end

exceededLimit = plotCube >= doseWindow(2);
if any(exceededLimit(:))
    plotCube(exceededLimit) = doseWindow(2) - eps(doseWindow(2));
    if isempty(doseColorMap)
        doseColorMap = jet(64);
    end
    doseColorMap(end,:) = [0.65 0 0];
end

end
