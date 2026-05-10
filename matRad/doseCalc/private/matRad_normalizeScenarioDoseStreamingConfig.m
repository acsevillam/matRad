function cfg = matRad_normalizeScenarioDoseStreamingConfig(cfg,matRad_cfg,defaultCacheRoot)
% matRad_normalizeScenarioDoseStreamingConfig normalizes shared streaming options
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
    matRad_cfg.dispError('Scenario dose streaming configuration must be a struct.');
end

cfg = setDefault(cfg,'SecondPassStrategy','disk');
cfg.SecondPassStrategy = lower(normalizeText(cfg.SecondPassStrategy, ...
    'SecondPassStrategy',matRad_cfg));
if ~any(strcmp(cfg.SecondPassStrategy,{'disk','recompute'}))
    matRad_cfg.dispError('cfg.SecondPassStrategy must be ''disk'' or ''recompute''.');
end

cfg = setDefault(cfg,'KeepCache',false);
cfg.KeepCache = logicalScalar(cfg.KeepCache,'KeepCache',matRad_cfg);

if ~isfield(cfg,'CacheRoot') || isempty(cfg.CacheRoot)
    cfg.CacheRoot = defaultCacheRoot;
else
    cfg.CacheRoot = normalizeText(cfg.CacheRoot,'CacheRoot',matRad_cfg);
end
end

function cfg = setDefault(cfg,fieldName,value)
if ~isfield(cfg,fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = value;
end
end

function textValue = normalizeText(value,name,matRad_cfg)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || isempty(value)
    matRad_cfg.dispError('%s must be a non-empty character vector.',name);
end
textValue = value;
end

function value = logicalScalar(value,name,matRad_cfg)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
    matRad_cfg.dispError('%s must be a logical scalar.',name);
end
value = logical(value);
end
