function [displayDoseScale,displayDoseMode] = matRad_getDisplayDoseScale(pln,displayDoseMode)
% matRad_getDisplayDoseScale Resolve dose display scaling.
%
% call
%   displayDoseScale = matRad_getDisplayDoseScale(pln,displayDoseMode)
%
% input
%   pln:              matRad pln struct
%   displayDoseMode:  'perFraction' or 'total'
%
% output
%   displayDoseScale: scalar multiplier for display copies
%   displayDoseMode:  normalized display dose mode string

if nargin < 2 || isempty(displayDoseMode)
    displayDoseMode = 'perFraction';
end

if isstring(displayDoseMode)
    displayDoseMode = char(displayDoseMode);
end

switch lower(displayDoseMode)
    case 'perfraction'
        displayDoseMode = 'perFraction';
        displayDoseScale = 1;
    case 'total'
        displayDoseMode = 'total';
        displayDoseScale = getNumOfFractions(pln);
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Unknown displayDoseMode "%s". Use "perFraction" or "total".',displayDoseMode);
end

end

function numOfFractions = getNumOfFractions(pln)

numOfFractions = 1;
if isfield(pln,'numOfFractions') && ~isempty(pln.numOfFractions)
    numOfFractions = pln.numOfFractions;
end

if ~(isnumeric(numOfFractions) && isscalar(numOfFractions) && ...
        isfinite(numOfFractions) && numOfFractions > 0)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('pln.numOfFractions must be a positive finite scalar.');
end

end
