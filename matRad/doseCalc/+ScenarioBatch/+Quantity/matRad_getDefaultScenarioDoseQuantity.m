function value = matRad_getDefaultScenarioDoseQuantity(pln)
% matRad_getDefaultScenarioDoseQuantity returns the plan default quantity
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

value = 'physicalDose';
if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'quantityOpt') && ...
        ~isempty(pln.propOpt.quantityOpt)
    value = pln.propOpt.quantityOpt;
elseif isfield(pln, 'bioModel') && isobject(pln.bioModel) && ...
        isprop(pln.bioModel, 'defaultReportQuantity') && ...
        ~isempty(pln.bioModel.defaultReportQuantity)
    value = pln.bioModel.defaultReportQuantity;
elseif isfield(pln, 'bioModel') && isstruct(pln.bioModel) && ...
        isfield(pln.bioModel, 'defaultReportQuantity') && ...
        ~isempty(pln.bioModel.defaultReportQuantity)
    value = pln.bioModel.defaultReportQuantity;
end

if isstring(value)
    value = char(value);
end
end
