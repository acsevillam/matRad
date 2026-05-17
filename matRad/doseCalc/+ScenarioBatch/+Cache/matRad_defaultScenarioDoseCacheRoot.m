function cacheRoot = matRad_defaultScenarioDoseCacheRoot(cacheName)
% matRad_defaultScenarioDoseCacheRoot returns the default scenario-batch cache root
%
% Scenario-dose scenario-batch writes temporary matrix blocks during disk-backed
% second passes. The default location lives outside the matRad checkout;
% callers can still override it explicitly through cfg.CacheRoot.
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

if nargin < 1 || isempty(cacheName)
    cacheName = 'scenarioDose';
end
if isstring(cacheName) && isscalar(cacheName)
    cacheName = char(cacheName);
end
if ~ischar(cacheName) || isempty(cacheName)
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('Scenario dose cache name must be a non-empty character vector.');
end

cacheRoot = fullfile(tempdir(), 'matRad', 'scenarioDose', cacheName);
end
