function cfg = matRad_normalizeScenarioDosePrecomputeConfig(cfg, matRadCfg, defaultCacheRoot)
% matRad_normalizeScenarioDosePrecomputeConfig normalizes shared scenario-batch options
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
    matRadCfg.dispError('Scenario dose precomputation configuration must be a struct.');
end

cfg = matRad_setDefault(cfg, 'SecondPassStrategy', 'disk');
cfg.SecondPassStrategy = lower(matRad_normalizeText(cfg.SecondPassStrategy, ...
                                                    'SecondPassStrategy', matRadCfg));
if ~any(strcmp(cfg.SecondPassStrategy, {'disk', 'recompute'}))
    matRadCfg.dispError('cfg.SecondPassStrategy must be ''disk'' or ''recompute''.');
end

cfg = matRad_setDefault(cfg, 'KeepCache', false);
cfg.KeepCache = matRad_logicalScalar(cfg.KeepCache, 'KeepCache', matRadCfg);

if ~isfield(cfg, 'CacheRoot') || isempty(cfg.CacheRoot)
    cfg.CacheRoot = defaultCacheRoot;
else
    cfg.CacheRoot = matRad_normalizeText(cfg.CacheRoot, 'CacheRoot', matRadCfg);
end
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

function value = matRad_logicalScalar(value, name, matRadCfg)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
    matRadCfg.dispError('%s must be a logical scalar.', name);
end
value = logical(value);
end
