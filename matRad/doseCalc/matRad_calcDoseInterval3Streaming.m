function [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3Streaming(ct,cst,stf,pln,varargin)
% matRad_calcDoseInterval3Streaming calculates INTERVAL3 data scenario by scenario
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3Streaming(ct,cst,stf,pln)
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3Streaming(ct,cst,stf,pln,dij)
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3Streaming(ct,cst,stf,pln,cfg)
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3Streaming(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model
%   dij:    optional precomputed robust dose influence struct. If provided,
%           streaming uses dij scenario matrices instead of recalculating
%           scenario dose influence data
%   cfg:    optional configuration struct accepted by matRad_calcDoseInterval3 plus:
%           UseParallel: use safe available scenario parallelism for
%               scenario-reducible streaming passes when the Parallel
%               Computing Toolbox and enough workers/memory are available.
%               matRad may create or reduce the active parallel pool.
%               The OAR radius-factor pass keeps its batch/voxel parallel
%               path. Precomputed dij inputs fall back to serial streaming.
%           SecondPassStrategy: 'disk' (default) or 'recompute'
%           CacheRoot: root folder for disk blocks
%           KeepCache: keep the hash cache folder after the run (default false)
%
% output
%   pln_interval:        plan struct containing propOpt.dij_interval. The
%                        propOpt.dij_interval.streamingSize field summarizes
%                        compact result bytes and peak streaming auxiliary
%                        bytes used during precomputation.
%   dij_intervalContext: lightweight single-scenario dij context for
%                        interval fluence optimization
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

matRad_cfg = MatRad_Config.instance();
[~,cfg] = matRad_parseScenarioDoseStreamingArguments(matRad_cfg,varargin{:});

[pln_interval,dij_intervalContext] = ...
    matRad_calcDoseIntervalStreamingCore(ct,cst,stf,pln,cfg,'INTERVAL3');

end
