function [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalStreaming(ct,cst,stf,pln,varargin)
% matRad_calcDoseIntervalStreaming calculates interval data scenario by scenario
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalStreaming(ct,cst,stf,pln,cfg)
%   [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalStreaming(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model
%   dij:    optional precomputed robust dose influence struct. If provided,
%           streaming uses dij scenario matrices instead of recalculating
%           scenario dose influence data
%   cfg:    configuration struct. Required field:
%           IntervalMode: 'INTERVAL2' or 'INTERVAL3'
%           Optional fields accepted by matRad_calcDoseInterval2/3 plus:
%           SecondPassStrategy: 'disk' (default) or 'recompute'
%           CacheRoot: root folder for disk blocks
%           KeepCache: keep the hash cache folder after the run (default false)
%
% output
%   pln_interval:        plan struct containing propOpt.dij_interval
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
if ~isfield(cfg,'IntervalMode') || isempty(cfg.IntervalMode)
    matRad_cfg.dispError('cfg.IntervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end

[pln_interval,dij_intervalContext] = ...
    matRad_calcDoseIntervalStreamingCore(ct,cst,stf,pln,cfg,cfg.IntervalMode);

end
