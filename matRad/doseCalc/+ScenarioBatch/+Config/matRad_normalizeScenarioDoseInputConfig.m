function cfg = matRad_normalizeScenarioDoseInputConfig(cfg, ct, pln, calculationMode, matRadCfg)
% matRad_normalizeScenarioDoseInputConfig normalizes scenario-dose inputs
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

if nargin < 1 || isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRadCfg.dispError('Scenario dose configuration must be a struct.');
end

matRad_normalizeText(calculationMode, 'calculationMode', matRadCfg);

cfg = matRad_setDefault(cfg, 'Quantity', ...
                        ScenarioBatch.Quantity.matRad_getDefaultScenarioDoseQuantity(pln));
cfg = matRad_setDefault(cfg, 'QuantityField', []);
cfg = matRad_setDefault(cfg, 'refScen', matRad_getDefaultReferenceScenario(ct));
cfg = matRad_setDefault(cfg, 'targetStructSel', []);
cfg = matRad_setDefault(cfg, 'OARStructSel', []);
cfg = matRad_setDefault(cfg, 'UseParallel', false);
cfg = matRad_setDefault(cfg, 'parallelOptions', []);
cfg = matRad_setDefault(cfg, 'CollectTiming', false);
cfg = matRad_setDefault(cfg, 'MemoryLimitMB', 256);
cfg = matRad_setDefault(cfg, 'BatchSize', []);
cfg = matRad_setDefault(cfg, 'ProgressLevel', 'summary');

cfg.refScen = matRad_validatePositiveInteger(cfg.refScen, 'refScen', matRadCfg);
if isfield(ct, 'numOfCtScen') && cfg.refScen > ct.numOfCtScen
    matRadCfg.dispError('refScen (%d) exceeds ct.numOfCtScen (%d).', ...
                        cfg.refScen, ct.numOfCtScen);
end

if isfield(ct, 'refScen') && ~isempty(ct.refScen) && ct.refScen ~= cfg.refScen
    matRadCfg.dispError('Requested refScen %d does not match ct.refScen (%d).', ...
                        cfg.refScen, ct.refScen);
end

cfg.UseParallel = matRad_logicalScalar(cfg.UseParallel, 'UseParallel', matRadCfg);
cfg.CollectTiming = matRad_logicalScalar(cfg.CollectTiming, 'CollectTiming', matRadCfg);
ScenarioBatch.Pool.matRad_doseParallelPoolOptions(cfg, matRadCfg, 'parallelOptions');

if ~isnumeric(cfg.MemoryLimitMB) || ~isscalar(cfg.MemoryLimitMB) || ...
        ~isfinite(cfg.MemoryLimitMB) || cfg.MemoryLimitMB <= 0
    matRadCfg.dispError('MemoryLimitMB must be a positive finite scalar.');
end

if ~isempty(cfg.BatchSize)
    cfg.BatchSize = matRad_validatePositiveInteger(cfg.BatchSize, 'BatchSize', matRadCfg);
end

cfg.ProgressLevel = matRad_validateProgressLevel(cfg.ProgressLevel, matRadCfg);
end

function progressLevel = matRad_validateProgressLevel(progressLevel, matRadCfg)
if isstring(progressLevel) && isscalar(progressLevel)
    progressLevel = char(progressLevel);
end

if ~ischar(progressLevel) || ...
        ~any(strcmpi(progressLevel, {'summary', 'detailed'}))
    matRadCfg.dispError('ProgressLevel must be ''summary'' or ''detailed''.');
end

progressLevel = lower(progressLevel);
end

function refScen = matRad_getDefaultReferenceScenario(ct)
if isfield(ct, 'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function cfg = matRad_setDefault(cfg, fieldName, value)
if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = value;
end
end

function value = matRad_validatePositiveInteger(value, name, matRadCfg)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || fix(value) ~= value
    matRadCfg.dispError('%s must be a positive integer scalar.', name);
end
value = double(value);
end

function value = matRad_logicalScalar(value, name, matRadCfg)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
    matRadCfg.dispError('%s must be a logical scalar.', name);
end
value = logical(value);
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
