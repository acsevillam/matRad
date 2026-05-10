function rows = matRad_readScenarioDoseCacheBlock(cacheContext,scenarioIx,kind,blockIx,expectedRows)
% matRad_readScenarioDoseCacheBlock reads and validates a row cache block
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

fileName = cacheBlockFile(cacheContext,scenarioIx,kind,blockIx);
matRad_cfg = MatRad_Config.instance();
if exist(fileName,'file') ~= 2
    matRad_cfg.dispError('Scenario dose cache block is missing: %s',fileName);
end

try
    data = load(fileName,'rowIx','rows');
catch ME
    matRad_cfg.dispError('Scenario dose cache block is unreadable: %s (%s)', ...
        fileName,ME.message);
end
if ~isfield(data,'rowIx') || ~isfield(data,'rows') || ...
        ~isequal(data.rowIx(:),expectedRows(:))
    matRad_cfg.dispError('Scenario dose cache block metadata is inconsistent: %s',fileName);
end
rows = data.rows;
end

function fileName = cacheBlockFile(cacheContext,scenarioIx,kind,blockIx)
fileName = fullfile(cacheContext.runDir, ...
    sprintf('scenario_%04d_%s_block_%04d.mat',scenarioIx,kind,blockIx));
end
