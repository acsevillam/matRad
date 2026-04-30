function [selectedIndices,selection] = matRad_resolveStructureSelection(cst,selectionMode,structures,structureTypes)
% matRad_resolveStructureSelection resolves cst structure selections
%
% call
%   selectedIndices = matRad_resolveStructureSelection(cst,selectionMode,structures,structureTypes)
%   [selectedIndices,selection] = matRad_resolveStructureSelection(...)
%
% input
%   cst:                matRad cst cell array
%   selectionMode:      structure selection mode, 'all', 'include', or
%                       'exclude'
%   structures:         selected structures as cst row indices, names, or a
%                       cell array combining scalar indices and names
%   structureTypes:     optional allowed cst structure types, e.g. 'TARGET',
%                       'OAR', or {'TARGET','OAR'}. Empty means no type
%                       filter.
%
% output
%   selectedIndices:    selected cst row indices
%   selection:          struct describing the normalized selection
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

if nargin < 2 || isempty(selectionMode)
    selectionMode = 'all';
end

if nargin < 3
    structures = [];
end

if nargin < 4
    structureTypes = [];
end

matRad_cfg = MatRad_Config.instance();
validateCst(cst,matRad_cfg);

selectionMode = normalizeSelectionMode(selectionMode,matRad_cfg);
allowedTypes = normalizeStructureTypes(structureTypes,matRad_cfg);
allowedIndices = findAllowedStructureIndices(cst,allowedTypes);

switch selectionMode
    case 'all'
        if hasStructureRequest(structures)
            matRad_cfg.dispError('Structure selection mode ''all'' does not accept a structure list.');
        end
        selectedIndices = allowedIndices;
    case 'include'
        requestedIndices = resolveRequestedStructures(cst,structures,matRad_cfg);
        validateRequestedStructureTypes(cst,requestedIndices,allowedTypes,matRad_cfg);
        selectedIndices = allowedIndices(ismember(allowedIndices,requestedIndices));
    case 'exclude'
        requestedIndices = resolveRequestedStructures(cst,structures,matRad_cfg);
        validateRequestedStructureTypes(cst,requestedIndices,allowedTypes,matRad_cfg);
        selectedIndices = setdiff(allowedIndices,requestedIndices,'stable');
end

selectedIndices = selectedIndices(:)';
selection.mode = selectionMode;
selection.requestedStructures = structures;
selection.allowedTypes = allowedTypes;
selection.selectedCstIndices = selectedIndices;
selection.selectedNames = cst(selectedIndices,2)';

end

function validateCst(cst,matRad_cfg)
if ~iscell(cst) || size(cst,2) < 3
    matRad_cfg.dispError('Structure selection requires a valid cst cell array.');
end
end

function selectionMode = normalizeSelectionMode(selectionMode,matRad_cfg)
if isstring(selectionMode) && isscalar(selectionMode)
    selectionMode = char(selectionMode);
end

if ~ischar(selectionMode)
    matRad_cfg.dispError('Structure selection mode must be ''all'', ''include'', or ''exclude''.');
end

selectionMode = lower(selectionMode);
if ~any(strcmp(selectionMode,{'all','include','exclude'}))
    matRad_cfg.dispError('Structure selection mode must be ''all'', ''include'', or ''exclude''.');
end
end

function structureTypes = normalizeStructureTypes(structureTypes,matRad_cfg)
if isempty(structureTypes)
    structureTypes = {};
    return;
end

if ischar(structureTypes) || (isstring(structureTypes) && isscalar(structureTypes))
    structureTypes = {char(structureTypes)};
elseif isstring(structureTypes)
    structureTypes = cellstr(structureTypes(:))';
elseif iscell(structureTypes)
    validTypes = cellfun(@(typeValue) ischar(typeValue) || ...
        (isstring(typeValue) && isscalar(typeValue)),structureTypes);
    if ~all(validTypes)
        matRad_cfg.dispError('Structure types must be specified as char, string, or cell array.');
    end
    structureTypes = cellfun(@char,structureTypes(:)','UniformOutput',false);
else
    matRad_cfg.dispError('Structure types must be specified as char, string, or cell array.');
end

emptyTypes = cellfun(@isempty,structureTypes);
if any(emptyTypes)
    matRad_cfg.dispError('Structure types must not contain empty values.');
end

structureTypes = upper(structureTypes);
structureTypes = unique(structureTypes,'stable');
end

function allowedIndices = findAllowedStructureIndices(cst,allowedTypes)
allowedIndices = 1:size(cst,1);
if isempty(allowedTypes)
    return;
end

cstTypes = cellfun(@char,cst(:,3),'UniformOutput',false);
isAllowed = false(size(cstTypes));
for typeIx = 1:numel(allowedTypes)
    isAllowed = isAllowed | strcmpi(cstTypes,allowedTypes{typeIx});
end
allowedIndices = find(isAllowed)';
end

function tf = hasStructureRequest(structures)
tf = ~isempty(structures);
if isstring(structures)
    tf = any(strlength(structures) > 0);
elseif iscell(structures)
    tf = any(~cellfun(@isempty,structures));
end
end

function requestedIndices = resolveRequestedStructures(cst,structures,matRad_cfg)
if ~hasStructureRequest(structures)
    matRad_cfg.dispError('Structure selection mode requires a non-empty structure list.');
end

requestedIndices = [];
if isnumeric(structures)
    requestedIndices = validateNumericIndices(cst,structures,matRad_cfg);
elseif ischar(structures) || isstring(structures)
    requestedIndices = resolveStructureNames(cst,cellstr(structures),matRad_cfg);
elseif iscell(structures)
    for i = 1:numel(structures)
        item = structures{i};
        if isnumeric(item)
            requestedIndices = [requestedIndices validateNumericIndices(cst,item,matRad_cfg)]; %#ok<AGROW>
        elseif ischar(item) || (isstring(item) && isscalar(item))
            requestedIndices = [requestedIndices resolveStructureNames(cst,{char(item)},matRad_cfg)]; %#ok<AGROW>
        else
            matRad_cfg.dispError('Structure list entries must be names or scalar cst row indices.');
        end
    end
else
    matRad_cfg.dispError('Structures must be specified as names or cst row indices.');
end

requestedIndices = unique(requestedIndices,'stable');
end

function indices = validateNumericIndices(cst,indices,matRad_cfg)
validIndices = isnumeric(indices) && isvector(indices) && all(isfinite(indices(:))) && ...
    all(indices(:) == round(indices(:))) && all(indices(:) >= 1) && ...
    all(indices(:) <= size(cst,1));
if ~validIndices
    matRad_cfg.dispError('Structure indices must be finite integer cst row indices.');
end
indices = indices(:)';
end

function indices = resolveStructureNames(cst,names,matRad_cfg)
indices = zeros(1,numel(names));
cstNames = cst(:,2);
for nameIx = 1:numel(names)
    name = names{nameIx};
    if isempty(name)
        matRad_cfg.dispError('Structure names must not be empty.');
    end

    matches = find(strcmp(cstNames,name));
    if isempty(matches)
        matRad_cfg.dispError('Structure "%s" was not found in cst.',name);
    elseif numel(matches) > 1
        matRad_cfg.dispError('Structure name "%s" is ambiguous. Use the cst row index instead.',name);
    end
    indices(nameIx) = matches;
end
end

function validateRequestedStructureTypes(cst,requestedIndices,allowedTypes,matRad_cfg)
if isempty(allowedTypes)
    return;
end

cstTypes = cellfun(@char,cst(requestedIndices,3),'UniformOutput',false);
isAllowed = false(size(cstTypes));
for typeIx = 1:numel(allowedTypes)
    isAllowed = isAllowed | strcmpi(cstTypes,allowedTypes{typeIx});
end

if any(~isAllowed)
    badIndex = requestedIndices(find(~isAllowed,1,'first'));
    matRad_cfg.dispError('Structure "%s" has type "%s", but allowed types are: %s.', ...
        cst{badIndex,2},cst{badIndex,3},strjoin(allowedTypes,', '));
end
end
