function rtssRowsByScenario = matRad_selectMultipleDicomRtssRows(fileList, selectedPatientIDs, ctSeriesUIDs, metadata)
% matRad selects RTSTRUCT file rows for multiple DICOM CT scenarios
%
% call
%   rtssRowsByScenario = matRad_selectMultipleDicomRtssRows(fileList, ...
%       selectedPatientIDs,ctSeriesUIDs,metadata)
%
% input
%   fileList:             DICOM file list returned by matRad_scanDicomImportFolder
%   selectedPatientIDs:   PatientIDs used to filter the import. Empty means
%                         explicit CT series selection across all patients
%   ctSeriesUIDs:         selected CT SeriesInstanceUIDs
%   metadata:             struct with optional rtssUID and rtssUIDs selections
%
% output
%   rtssRowsByScenario:   cell array with selected RTSTRUCT file rows per CT
%                         scenario, or empty rows when no RTSTRUCT is present
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

matRadCfg = MatRad_Config.instance();

rtssMask = strcmp(fileList(:, 2), 'RTSTRUCT') & matRad_matchesPatientFilter(fileList(:, 3), selectedPatientIDs);
availableRtssRows = fileList(rtssMask, :);
rtssRowsByScenario = cell(1, numel(ctSeriesUIDs));

if isfield(metadata, 'rtssUIDs') && ~isempty(metadata.rtssUIDs)
    rtssUIDs = matRad_metadataValuesToCell(metadata.rtssUIDs, 'metadata.rtssUIDs', matRadCfg);
    if numel(rtssUIDs) ~= numel(ctSeriesUIDs)
        matRadCfg.dispError('metadata.rtssUIDs must contain one RTSTRUCT SOPInstanceUID per CT scenario.');
    end
    for ctPhase = 1:numel(ctSeriesUIDs)
        rtssRowsByScenario{ctPhase} = matRad_selectRtssRow(fileList, rtssMask, ...
                                                           rtssUIDs{ctPhase}, matRadCfg);
    end
elseif isfield(metadata, 'rtssUID') && ~isempty(metadata.rtssUID)
    rtssUID = matRad_metadataValueToChar(metadata.rtssUID, 'metadata.rtssUID', matRadCfg);
    rtssRow = matRad_selectRtssRow(fileList, rtssMask, rtssUID, matRadCfg);
    rtssRowsByScenario(:) = {rtssRow};
else
    rtssRowsByScenario = matRad_selectDefaultRtssRows(availableRtssRows, ...
                                                      rtssRowsByScenario, matRadCfg);
end
end

function rtssRowsByScenario = matRad_selectDefaultRtssRows(availableRtssRows, ...
                                                           rtssRowsByScenario, matRadCfg)
if isempty(availableRtssRows)
    rtssRowsByScenario(:) = {[]};
    return
end
if size(availableRtssRows, 1) ~= 1
    matRadCfg.dispError('Multiple RTSTRUCT files found. Specify metadata.rtssUID or metadata.rtssUIDs.');
end
rtssRowsByScenario(:) = {availableRtssRows};
end

function rtssRow = matRad_selectRtssRow(fileList, rtssMask, rtssUID, matRadCfg)
rtssMatch = rtssMask & strcmp(fileList(:, 4), rtssUID);
if ~any(rtssMatch)
    matRadCfg.dispError('RTSTRUCT SOPInstanceUID "%s" was not found for the selected PatientID filter.', rtssUID);
end
rtssRow = fileList(rtssMatch, :);
if size(rtssRow, 1) ~= 1
    matRadCfg.dispError('RTSTRUCT SOPInstanceUID "%s" matched %d files.', rtssUID, size(rtssRow, 1));
end
end

function mask = matRad_matchesPatientFilter(patientIDs, selectedPatientIDs)
if isempty(selectedPatientIDs)
    mask = true(size(patientIDs));
else
    mask = ismember(patientIDs, selectedPatientIDs);
end
end

function values = matRad_metadataValuesToCell(values, fieldName, matRadCfg)
if ischar(values)
    values = {values};
elseif isstring(values)
    values = cellstr(values(:));
elseif iscell(values)
    values = values(:);
else
    matRadCfg.dispError('%s must be a character vector, string, or cell array of character vectors.', fieldName);
end

for valueIx = 1:numel(values)
    values{valueIx} = matRad_metadataValueToChar(values{valueIx}, fieldName, matRadCfg);
end
end

function value = matRad_metadataValueToChar(value, fieldName, matRadCfg)
if isstring(value) && isscalar(value)
    value = char(value);
elseif iscell(value) && isscalar(value)
    value = matRad_metadataValueToChar(value{1}, fieldName, matRadCfg);
elseif ~ischar(value)
    matRadCfg.dispError('%s must be a character vector or string scalar.', fieldName);
end
end
