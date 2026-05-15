function tf = matRad_isInteractiveSession(varargin)
% matRad_isInteractiveSession checks whether interactive UI features are usable
%
% call
%   tf = matRad_isInteractiveSession()
%   tf = matRad_isInteractiveSession('requireDesktop',true)
%   tf = matRad_isInteractiveSession('requireFigureWindows',true)
%
% input
%   varargin:             optional Name-Value pair arguments
%   requireDesktop:       require the MATLAB desktop, for example for
%                         command-window callbacks
%   requireFigureWindows: require visible figure-window support
%
% output
%   tf:                   true when the current session can host the
%                         requested interactive UI feature
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

p = inputParser;
p.addParameter('requireDesktop',false, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.addParameter('requireFigureWindows',false, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.parse(varargin{:});

matRad_cfg = MatRad_Config.instance();
tf = false;

if matRad_cfg.disableGUI || isdeployed
    return;
end

if matRad_cfg.isOctave
    tf = true;
    if p.Results.requireFigureWindows
        try
            tf = ~strcmpi(graphics_toolkit(),'gnuplot');
        catch
            tf = false;
        end
    end
    if p.Results.requireDesktop
        tf = false;
    end
    return;
end

if p.Results.requireDesktop
    try
        if ~(usejava('desktop') && usejava('awt'))
            return;
        end
    catch
        return;
    end
end

if p.Results.requireFigureWindows
    try
        if exist('feature','builtin') == 5 && ...
                ~feature('ShowFigureWindows')
            return;
        end
    catch
        return;
    end
end

tf = true;
end
