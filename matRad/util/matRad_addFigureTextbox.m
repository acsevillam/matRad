function textboxHandle = matRad_addFigureTextbox(figHandle,textString,varargin)
% matRad_addFigureTextbox adds a styled annotation textbox to a figure
%
% call
%   textboxHandle = matRad_addFigureTextbox(figHandle,textString)
%   textboxHandle = matRad_addFigureTextbox(figHandle,textString,'position',position)
%
% input
%   figHandle:     handle to the target figure
%   textString:    text shown in the annotation textbox
%
% input (optional Name-Value pairs)
%   varargin:     optional Name-Value pairs
%   position:      normalized annotation position, defaults to lower left
%   edgeColor:     textbox edge color
%   backgroundColor: textbox background color
%   fontSize:      textbox font size
%   fontWeight:    textbox font weight
%   margin:        textbox margin
%   tag:           graphics object tag
%
% output
%   textboxHandle: handle to the annotation textbox, empty if creation fails
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

textboxHandle = [];

p = inputParser;
p.addParameter('position',[0.02 0.02 0.28 0.06],@isValidPosition);
p.addParameter('edgeColor',[0.25 0.25 0.25],@isValidColor);
p.addParameter('backgroundColor',[1 1 1],@isValidColor);
p.addParameter('fontSize',10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('fontWeight','bold',@(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('margin',3,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('tag','matRadFigureTextbox',@(x) ischar(x) || (isstring(x) && isscalar(x)));
p.parse(varargin{:});

if isempty(figHandle) || ~ishandle(figHandle)
    return;
end

try
    textboxHandle = annotation(figHandle,'textbox',p.Results.position, ...
        'String',formatTextboxString(textString), ...
        'FitBoxToText','on', ...
        'EdgeColor',p.Results.edgeColor, ...
        'BackgroundColor',p.Results.backgroundColor, ...
        'FontSize',p.Results.fontSize, ...
        'FontWeight',char(p.Results.fontWeight), ...
        'Margin',p.Results.margin, ...
        'Tag',char(p.Results.tag));
catch
    textboxHandle = [];
end
end

function tf = isValidPosition(position)
tf = isnumeric(position) && numel(position) == 4 && all(isfinite(position(:)));
end

function tf = isValidColor(colorValue)
tf = isnumeric(colorValue) && numel(colorValue) == 3 && ...
    all(isfinite(colorValue(:))) && all(colorValue(:) >= 0) && ...
    all(colorValue(:) <= 1);
end

function textString = formatTextboxString(textString)
if isstring(textString)
    textString = cellstr(textString);
    if isscalar(textString)
        textString = textString{1};
    end
end
end
