function [cfg, intervalMode] = matRad_normalizeDoseIntervalConfig(cfg, intervalMode, matRadCfg)
% matRad_normalizeDoseIntervalConfig normalizes dose interval options
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

if isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRadCfg.dispError('Dose interval scenario-batch configuration must be a struct.');
end

intervalMode = matRad_normalizeText(intervalMode, 'intervalMode', matRadCfg);
intervalMode = upper(intervalMode);
if ~any(strcmp(intervalMode, {'INTERVAL2', 'INTERVAL3'}))
    matRadCfg.dispError('intervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end
cfg.IntervalMode = intervalMode;
cfg = matRad_setDefault(cfg, 'RadiusMode', 'std');
cfg = matRad_setDefault(cfg, 'KMode', 'dynamic');
cfg = matRad_setDefault(cfg, 'KMax', 10);
cfg = matRad_setDefault(cfg, 'RetentionThreshold', 1.0);

cfg.RadiusMode = matRad_validateRadiusMode(cfg.RadiusMode, matRadCfg);
if strcmp(intervalMode, 'INTERVAL3')
    cfg = matRad_normalizeInterval3Config(cfg, matRadCfg);
end

cfg = ScenarioBatch.Config.matRad_normalizeScenarioDosePrecomputeConfig(cfg, matRadCfg, ...
                                                                        ScenarioBatch.Cache.matRad_defaultScenarioDoseCacheRoot('doseInterval'));
end

function cfg = matRad_normalizeInterval3Config(cfg, matRadCfg)
cfg.KMax = matRad_validatePositiveInteger(cfg.KMax, 'KMax', matRadCfg);
if isstring(cfg.KMode) && isscalar(cfg.KMode)
    cfg.KMode = char(cfg.KMode);
end
if ~ischar(cfg.KMode) || ~any(strcmpi(cfg.KMode, {'dynamic', 'static'}))
    matRadCfg.dispError('KMode must be ''dynamic'' or ''static''.');
end
cfg.KMode = lower(cfg.KMode);
if ~isnumeric(cfg.RetentionThreshold) || ~isscalar(cfg.RetentionThreshold) || ...
        ~isfinite(cfg.RetentionThreshold) || cfg.RetentionThreshold <= 0 || ...
        cfg.RetentionThreshold > 1
    matRadCfg.dispError('RetentionThreshold must be in the interval (0,1].');
end
end

function radiusMode = matRad_validateRadiusMode(radiusMode, matRadCfg)
if isstring(radiusMode) && isscalar(radiusMode)
    radiusMode = char(radiusMode);
end

if ~ischar(radiusMode) || ...
        ~any(strcmpi(radiusMode, {'std', 'extreme'}))
    matRadCfg.dispError('RadiusMode must be ''std'' or ''extreme''.');
end

radiusMode = lower(radiusMode);
end

function value = matRad_validatePositiveInteger(value, name, matRadCfg)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || fix(value) ~= value
    matRadCfg.dispError('%s must be a positive integer scalar.', name);
end
value = double(value);
end

function cfg = matRad_setDefault(cfg, fieldName, value)
if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = value;
end
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
