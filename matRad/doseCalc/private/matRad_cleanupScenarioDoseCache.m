function matRad_cleanupScenarioDoseCache(cacheContext,keepCache)
% matRad_cleanupScenarioDoseCache removes a streaming cache folder if needed
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

if isempty(cacheContext) || ~isfield(cacheContext,'created') || ...
        ~cacheContext.created || keepCache
    return;
end

if exist(cacheContext.runDir,'dir') == 7
    [success,~] = rmdir(cacheContext.runDir,'s');
    if ~success
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispWarning('Could not remove scenario dose cache folder "%s".\n', ...
            cacheContext.runDir);
    end
end
end
