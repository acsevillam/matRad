function [selectedIndices, selection] = matRad_resolveStructureSelection(cst, selectionMode, structures, structureTypes)
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

matRadCfg = MatRad_Config.instance();
matRad_validateCst(cst, matRadCfg);

selectionMode = matRad_normalizeSelectionMode(selectionMode, matRadCfg);
allowedTypes = matRad_normalizeStructureTypes(structureTypes, matRadCfg);
allowedIndices = matRad_findAllowedStructureIndices(cst, allowedTypes);

switch selectionMode
    case 'all'
        if matRad_hasStructureRequest(structures)
            matRadCfg.dispError('Structure selection mode ''all'' does not accept a structure list.');
        end
        selectedIndices = allowedIndices;
    case 'include'
        requestedIndices = matRad_resolveRequestedStructures(cst, structures, matRadCfg);
        matRad_validateRequestedStructureTypes(cst, requestedIndices, allowedTypes, matRadCfg);
        selectedIndices = allowedIndices(ismember(allowedIndices, requestedIndices));
    case 'exclude'
        requestedIndices = matRad_resolveRequestedStructures(cst, structures, matRadCfg);
        matRad_validateRequestedStructureTypes(cst, requestedIndices, allowedTypes, matRadCfg);
        selectedIndices = setdiff(allowedIndices, requestedIndices, 'stable');
end

selectedIndices = selectedIndices(:)';
selection.mode = selectionMode;
selection.requestedStructures = structures;
selection.allowedTypes = allowedTypes;
selection.selectedCstIndices = selectedIndices;
selection.selectedNames = cst(selectedIndices, 2)';

end

function matRad_validateCst(cst, matRadCfg)
if ~iscell(cst) || size(cst, 2) < 3
    matRadCfg.dispError('Structure selection requires a valid cst cell array.');
end
end

function selectionMode = matRad_normalizeSelectionMode(selectionMode, matRadCfg)
if isstring(selectionMode) && isscalar(selectionMode)
    selectionMode = char(selectionMode);
end

if ~ischar(selectionMode)
    matRadCfg.dispError('Structure selection mode must be ''all'', ''include'', or ''exclude''.');
end

selectionMode = lower(selectionMode);
if ~any(strcmp(selectionMode, {'all', 'include', 'exclude'}))
    matRadCfg.dispError('Structure selection mode must be ''all'', ''include'', or ''exclude''.');
end
end

function structureTypes = matRad_normalizeStructureTypes(structureTypes, matRadCfg)
if isempty(structureTypes)
    structureTypes = {};
    return
end

if ischar(structureTypes) || (isstring(structureTypes) && isscalar(structureTypes))
    structureTypes = {char(structureTypes)};
elseif isstring(structureTypes)
    structureTypes = cellstr(structureTypes(:))';
elseif iscell(structureTypes)
    validTypes = cellfun(@(typeValue) ischar(typeValue) || ...
                         (isstring(typeValue) && isscalar(typeValue)), structureTypes);
    if ~all(validTypes)
        matRadCfg.dispError('Structure types must be specified as char, string, or cell array.');
    end
    structureTypes = cellfun(@char, structureTypes(:)', 'UniformOutput', false);
else
    matRadCfg.dispError('Structure types must be specified as char, string, or cell array.');
end

emptyTypes = cellfun(@isempty, structureTypes);
if any(emptyTypes)
    matRadCfg.dispError('Structure types must not contain empty values.');
end

structureTypes = upper(structureTypes);
structureTypes = unique(structureTypes, 'stable');
end

function allowedIndices = matRad_findAllowedStructureIndices(cst, allowedTypes)
allowedIndices = 1:size(cst, 1);
if isempty(allowedTypes)
    return
end

cstTypes = cellfun(@char, cst(:, 3), 'UniformOutput', false);
isAllowed = false(size(cstTypes));
for typeIx = 1:numel(allowedTypes)
    isAllowed = isAllowed | strcmpi(cstTypes, allowedTypes{typeIx});
end
allowedIndices = find(isAllowed)';
end

function tf = matRad_hasStructureRequest(structures)
tf = ~isempty(structures);
if isstring(structures)
    tf = any(strlength(structures) > 0);
elseif iscell(structures)
    tf = any(~cellfun(@isempty, structures));
end
end

function requestedIndices = matRad_resolveRequestedStructures(cst, structures, matRadCfg)
if ~matRad_hasStructureRequest(structures)
    matRadCfg.dispError('Structure selection mode requires a non-empty structure list.');
end

requestedIndices = [];
if isnumeric(structures)
    requestedIndices = matRad_validateNumericIndices(cst, structures, matRadCfg);
elseif ischar(structures) || isstring(structures)
    requestedIndices = matRad_resolveStructureNames(cst, cellstr(structures), matRadCfg);
elseif iscell(structures)
    for i = 1:numel(structures)
        item = structures{i};
        if isnumeric(item)
            requestedIndices = [requestedIndices matRad_validateNumericIndices(cst, item, matRadCfg)]; %#ok<AGROW>
        elseif ischar(item) || (isstring(item) && isscalar(item))
            requestedIndices = [requestedIndices matRad_resolveStructureNames(cst, {char(item)}, matRadCfg)]; %#ok<AGROW>
        else
            matRadCfg.dispError('Structure list entries must be names or scalar cst row indices.');
        end
    end
else
    matRadCfg.dispError('Structures must be specified as names or cst row indices.');
end

requestedIndices = unique(requestedIndices, 'stable');
end

function indices = matRad_validateNumericIndices(cst, indices, matRadCfg)
validIndices = isnumeric(indices) && isvector(indices) && all(isfinite(indices(:))) && ...
    all(indices(:) == round(indices(:))) && all(indices(:) >= 1) && ...
    all(indices(:) <= size(cst, 1));
if ~validIndices
    matRadCfg.dispError('Structure indices must be finite integer cst row indices.');
end
indices = indices(:)';
end

function indices = matRad_resolveStructureNames(cst, names, matRadCfg)
indices = zeros(1, numel(names));
cstNames = cst(:, 2);
for nameIx = 1:numel(names)
    name = names{nameIx};
    if isempty(name)
        matRadCfg.dispError('Structure names must not be empty.');
    end

    matches = find(strcmp(cstNames, name));
    if isempty(matches)
        matRadCfg.dispError('Structure "%s" was not found in cst.', name);
    elseif numel(matches) > 1
        matRadCfg.dispError('Structure name "%s" is ambiguous. Use the cst row index instead.', name);
    end
    indices(nameIx) = matches;
end
end

function matRad_validateRequestedStructureTypes(cst, requestedIndices, allowedTypes, matRadCfg)
if isempty(allowedTypes)
    return
end

cstTypes = cellfun(@char, cst(requestedIndices, 3), 'UniformOutput', false);
isAllowed = false(size(cstTypes));
for typeIx = 1:numel(allowedTypes)
    isAllowed = isAllowed | strcmpi(cstTypes, allowedTypes{typeIx});
end

if any(~isAllowed)
    badIndex = requestedIndices(find(~isAllowed, 1, 'first'));
    matRadCfg.dispError('Structure "%s" has type "%s", but allowed types are: %s.', ...
                        cst{badIndex, 2}, cst{badIndex, 3}, strjoin(allowedTypes, ', '));
end
end
