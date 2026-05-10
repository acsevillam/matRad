function [dij,cfg] = matRad_parseScenarioDoseStreamingArguments(matRad_cfg,varargin)
% matRad_parseScenarioDoseStreamingArguments parses optional dij/cfg inputs
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
    return;
end

if numel(varargin) > 2
    matRad_cfg.dispError('Streaming dose calculation accepts at most optional dij and cfg arguments.');
end

if numel(varargin) == 1
    value = varargin{1};
    if isempty(value)
        return;
    elseif looksLikeDij(value)
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
    matRad_cfg.dispError('Scenario dose streaming configuration must be a struct.');
end

if ~isempty(dij)
    if ~isstruct(dij)
        matRad_cfg.dispError('Optional dij argument must be a dose influence struct.');
    end
    if isfield(cfg,'PrecomputedDij') && ~isempty(cfg.PrecomputedDij)
        matRad_cfg.dispError(['Pass precomputed dose influence either as the optional ', ...
            'dij argument or as cfg.PrecomputedDij, not both.']);
    end
    cfg.PrecomputedDij = dij;
end
end

function tf = looksLikeDij(value)
tf = false;
if ~isstruct(value)
    return;
end

knownMatrixFields = {'physicalDose','mAlphaDose','mSqrtBetaDose'};
for i = 1:numel(knownMatrixFields)
    if hasDoseMatrixCell(value,knownMatrixFields{i})
        tf = true;
        return;
    end
end

fieldNames = fieldnames(value);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if hasTextSuffix(fieldName,'Dose') && hasDoseMatrixCell(value,fieldName)
        tf = true;
        return;
    end
end
end

function tf = hasTextSuffix(textValue,suffix)
tf = numel(textValue) >= numel(suffix) && ...
    strcmp(textValue(end - numel(suffix) + 1:end),suffix);
end

function tf = hasDoseMatrixCell(value,fieldName)
tf = false;
if ~isfield(value,fieldName) || ~iscell(value.(fieldName))
    return;
end

matrixCell = value.(fieldName);
for i = 1:numel(matrixCell)
    if isnumeric(matrixCell{i}) && ismatrix(matrixCell{i})
        tf = true;
        return;
    end
end
end
