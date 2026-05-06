function [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3(ct,cst,stf,pln,dij,cfg)
% matRad_calcDoseInterval3 calculates interval dose influence data with OAR radii
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3(ct,cst,stf,pln,dij)
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval3(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model; optional
%           pln.propOpt.scen4D selects CT scenarios used for the interval
%           calculation (default: 1; 'all' or positive integer CT
%           scenario ids)
%   dij:    matRad dose influence struct for the robust scenarios
%   cfg:    optional configuration struct with fields:
%           Quantity: optimized linear quantity (default: pln.bioParam.quantityOpt)
%           QuantityField: explicit linear dij field (default: [])
%           refScen: reference CT scenario (default: ct.refScen, otherwise 1)
%           targetStructSel: target structure selection (default: all TARGETs)
%           OARStructSel: OAR structure selection (default: all OARs)
%           UseParallel: parallelize INTERVAL3 OAR center/scenario SVD
%               batches when the Parallel Computing Toolbox is available;
%               target interval accumulation remains serial (default: false)
%           MemoryLimitMB: batch and OAR SVD memory limit, [] uses 256 MB internally
%           BatchSize: explicit voxel batch size (default: derived from MemoryLimitMB)
%           ProgressLevel: 'summary' (default) or 'detailed'
%           KMode: OAR SVD truncation mode, 'dynamic' (default) or 'static'
%           KMax: maximum retained OAR SVD rank (default: 10)
%           RetentionThreshold: dynamic retained singular-value energy in (0,1] (default: 1)
%           CollectTiming: collect internal timing diagnostics in
%               pln_interval.propOpt.dij_interval.timing (default: false)
%
% output
%   pln_interval:         plan struct; pln_interval.propOpt.dij_interval
%                         contains the full interval data
%   dij_intervalContext:  lightweight dij context passed to
%                         matRad_fluenceOptimization
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

if nargin < 6
    cfg = struct();
end

[pln_interval,dij_intervalContext] = ...
    matRad_calcDoseIntervalCore(ct,cst,stf,pln,dij,cfg,'INTERVAL3');

end
