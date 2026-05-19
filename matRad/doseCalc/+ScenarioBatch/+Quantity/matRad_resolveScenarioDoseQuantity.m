function quantity = matRad_resolveScenarioDoseQuantity(dij, pln, cfg, matRadCfg)
% matRad_resolveScenarioDoseQuantity resolves the linear scenario dose field
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

quantityName = matRad_normalizeText(cfg.Quantity, 'Quantity', matRadCfg);
quantityField = matRad_normalizeOptionalText(cfg.QuantityField, 'QuantityField', matRadCfg);
scale = 1;

if ~isempty(quantityField)
    field = quantityField;
elseif matRad_hasLinearDijField(dij, quantityName)
    field = quantityName;
elseif strcmpi(quantityName, 'physicalDose')
    field = 'physicalDose';
elseif strcmpi(quantityName, 'RBExDose')
    if isfield(dij, 'RBE') && isnumeric(dij.RBE) && isscalar(dij.RBE) && isfinite(dij.RBE)
        field = 'physicalDose';
        scale = dij.RBE;
    else
        matRadCfg.dispError(['RBExDose scenario dose calculation requires a scalar dij.RBE ' ...
                             'or an explicit linear QuantityField.']);
    end
elseif any(strcmpi(quantityName, {'effect', 'BED'}))
    matRadCfg.dispError('%s scenario dose calculation requires an explicit linear QuantityField.', ...
                        quantityName);
else
    matRadCfg.dispError('Could not resolve linear scenario dose quantity ''%s''.', quantityName);
end

if ~matRad_hasLinearDijField(dij, field)
    matRadCfg.dispError('dij.%s must be a cell array of dose influence matrices.', field);
end

quantity.name = quantityName;
quantity.field = field;
quantity.scale = scale;
quantity.matrixCell = dij.(field);
quantity.optimizationQuantity = ScenarioBatch.Quantity.matRad_getDefaultScenarioDoseQuantity(pln);
end

function textValue = matRad_normalizeText(value, name, matRadCfg)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || isempty(value)
    matRadCfg.dispError('%s must be a non-empty character vector.', name);
end
textValue = value;
end

function textValue = matRad_normalizeOptionalText(value, name, matRadCfg)
if isempty(value)
    textValue = '';
    return
end
textValue = matRad_normalizeText(value, name, matRadCfg);
end

function tf = matRad_hasLinearDijField(dij, fieldName)
tf = isfield(dij, fieldName) && iscell(dij.(fieldName));
end
