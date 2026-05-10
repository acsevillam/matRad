function matRad_writeScenarioDoseCacheBlock(cacheContext,scenarioIx,kind,blockIx,rowIx,rows)
% matRad_writeScenarioDoseCacheBlock writes a sparse row cache block
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
save(fileName,'scenarioIx','kind','blockIx','rowIx','rows','-v7');
end

function fileName = cacheBlockFile(cacheContext,scenarioIx,kind,blockIx)
fileName = fullfile(cacheContext.runDir, ...
    sprintf('scenario_%04d_%s_block_%04d.mat',scenarioIx,kind,blockIx));
end
