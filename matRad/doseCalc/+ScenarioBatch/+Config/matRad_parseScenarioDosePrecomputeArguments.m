function [dij, cfg] = matRad_parseScenarioDosePrecomputeArguments(matRadCfg, varargin)
% matRad_parseScenarioDosePrecomputeArguments parses optional dij/cfg inputs
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

dij = [];
cfg = struct();

if isempty(varargin)
    return
end

if numel(varargin) > 2
    matRadCfg.dispError('Scenario dose precomputation accepts at most optional dij and cfg arguments.');
end

if numel(varargin) == 1
    value = varargin{1};
    if isempty(value)
        return
    elseif matRad_looksLikeDij(value)
        dij = value;
    else
        cfg = value;
    end
else
    dij = varargin{1};
    cfg = varargin{2};
end

if isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRadCfg.dispError('Scenario dose precomputation configuration must be a struct.');
end

if ~isempty(dij)
    if ~isstruct(dij)
        matRadCfg.dispError('Optional dij argument must be a dose influence struct.');
    end
    if isfield(cfg, 'PrecomputedDij') && ~isempty(cfg.PrecomputedDij)
        matRadCfg.dispError(['Pass precomputed dose influence either as the optional ', ...
                             'dij argument or as cfg.PrecomputedDij, not both.']);
    end
    cfg.PrecomputedDij = dij;
end
end

function tf = matRad_looksLikeDij(value)
tf = false;
if ~isstruct(value)
    return
end

knownMatrixFields = {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose'};
for i = 1:numel(knownMatrixFields)
    if matRad_hasDoseMatrixCell(value, knownMatrixFields{i})
        tf = true;
        return
    end
end

fieldNames = fieldnames(value);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if matRad_hasTextSuffix(fieldName, 'Dose') && matRad_hasDoseMatrixCell(value, fieldName)
        tf = true;
        return
    end
end
end

function tf = matRad_hasTextSuffix(textValue, suffix)
tf = numel(textValue) >= numel(suffix) && ...
    strcmp(textValue(end - numel(suffix) + 1:end), suffix);
end

function tf = matRad_hasDoseMatrixCell(value, fieldName)
tf = false;
if ~isfield(value, fieldName) || ~iscell(value.(fieldName))
    return
end

matrixCell = value.(fieldName);
for i = 1:numel(matrixCell)
    if isnumeric(matrixCell{i}) && ismatrix(matrixCell{i})
        tf = true;
        return
    end
end
end
