function [memoryLimitMB, memoryInfo] = matRad_resolveScenarioDoseMemoryLimitMB(cfg, matRadCfg, memoryInfo)
% matRad_resolveScenarioDoseMemoryLimitMB resolves scenario-dose memory guard
%
% call
%   memoryLimitMB = ScenarioBatch.Config.matRad_resolveScenarioDoseMemoryLimitMB(cfg,matRadCfg)
%
% input
%   cfg:        scenario-dose configuration with MemoryLimitMB,
%               MemoryLimitFraction, and MemoryLimitFallbackMB
%   matRadCfg: MatRad_Config instance for diagnostics
%   memoryInfo: optional memory information struct for tests
%
% output
%   memoryLimitMB: numeric effective memory limit in MB
%   memoryInfo:    memory information used for automatic resolution
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

if nargin < 2 || isempty(matRadCfg)
    matRadCfg = MatRad_Config.instance();
end
if nargin < 3
    memoryInfo = [];
end

if isnumeric(cfg.MemoryLimitMB)
    if ~isscalar(cfg.MemoryLimitMB) || ~isfinite(cfg.MemoryLimitMB) || ...
            cfg.MemoryLimitMB <= 0
        matRadCfg.dispError('MemoryLimitMB must be a positive finite scalar or ''auto''.');
    end
    memoryLimitMB = double(cfg.MemoryLimitMB);
    return
end

if ~matRad_isAutoMemoryLimit(cfg.MemoryLimitMB)
    matRadCfg.dispError('MemoryLimitMB must be a positive finite scalar or ''auto''.');
end

if isempty(memoryInfo)
    memoryInfo = matRad_getSystemMemoryInfo();
end

if ~isstruct(memoryInfo) || ~isfield(memoryInfo, 'usableBytes') || ...
        isempty(memoryInfo.usableBytes) || ~isfinite(memoryInfo.usableBytes) || ...
        memoryInfo.usableBytes <= 0
    memoryLimitMB = double(cfg.MemoryLimitFallbackMB);
    matRadCfg.dispWarning(['Could not detect a reliable scenario-batch memory ', ...
                           'budget from "%s". Falling back to MemoryLimitMB ', ...
                           '%.2f MB.\n'], matRad_memoryInfoSource(memoryInfo), ...
                          memoryLimitMB);
    return
end

memoryLimitMB = double(memoryInfo.usableBytes) * double(cfg.MemoryLimitFraction) / 1e6;
if ~isfinite(memoryLimitMB) || memoryLimitMB <= 0
    memoryLimitMB = double(cfg.MemoryLimitFallbackMB);
    matRadCfg.dispWarning(['Auto scenario-batch MemoryLimitMB resolved to an ', ...
                           'invalid value. Falling back to %.2f MB.\n'], ...
                          memoryLimitMB);
    return
end

matRadCfg.dispInfo(['matRad: Auto scenario-batch MemoryLimitMB resolved to ', ...
                    '%.2f MB from %s usable memory %.2f MB with fraction %.2f.\n'], ...
                   memoryLimitMB, matRad_memoryInfoSource(memoryInfo), ...
                   double(memoryInfo.usableBytes) / 1e6, ...
                   double(cfg.MemoryLimitFraction));
end

function tf = matRad_isAutoMemoryLimit(value)
tf = false;
if ischar(value)
    tf = strcmpi(strtrim(value), 'auto');
elseif isstring(value) && isscalar(value)
    tf = strcmpi(strtrim(char(value)), 'auto');
end
end

function source = matRad_memoryInfoSource(memoryInfo)
source = 'unknown';
if isstruct(memoryInfo) && isfield(memoryInfo, 'source') && ...
        ~isempty(memoryInfo.source)
    source = char(string(memoryInfo.source));
end
end
