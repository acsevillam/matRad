function specId = macroSpecId(site,particleType,caseID,robustness,samplingProfile)
% macroSpecId Build a local macro spec ID from wrapper selectors.

if nargin < 5 || isempty(samplingProfile)
    samplingProfile = 'default';
end

site = normalizeText(site,'site');
particleType = normalizeText(particleType,'particleType');
caseID = normalizeText(caseID,'caseID');
robustness = normalizeText(robustness,'robustness');
samplingProfile = normalizeText(samplingProfile,'samplingProfile');

switch samplingProfile
    case 'default'
        planAlias = robustness;
    case 'noAngles'
        planAlias = [robustness '_noAngles'];
    otherwise
        error('planWorkflow:macros:InvalidSamplingProfile', ...
            'Unsupported macro sampling profile "%s".',samplingProfile);
end

specId = strjoin({site,particleType,caseID,planAlias},'.');

end

function value = normalizeText(value,fieldName)

if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('planWorkflow:macros:InvalidMacroSpecSelector', ...
        '%s must be a non-empty text scalar.',fieldName);
end
value = char(string(value));
if isempty(strtrim(value))
    error('planWorkflow:macros:InvalidMacroSpecSelector', ...
        '%s must be a non-empty text scalar.',fieldName);
end

end
