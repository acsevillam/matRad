function [selectedPatientIDs, ctSeriesUIDs, rtssRowsByScenario] = matRad_resolveMultipleDicomCtImportSelection(fileList, patientList, metadata)
% matRad resolves patient, CT series, and RTSTRUCT selections for multi-CT import
%
% call
%   [selectedPatientIDs,ctSeriesUIDs,rtssRowsByScenario] = ...
%       matRad_resolveMultipleDicomCtImportSelection(fileList,patientList,metadata)
%
% input
%   fileList:             DICOM file list returned by matRad_scanDicomImportFolder
%   patientList:          PatientID list returned by matRad_scanDicomImportFolder
%   metadata:             struct with optional patientID, ctSeriesUIDs,
%                         rtssUID, and rtssUIDs import selections
%
% output
%   selectedPatientIDs:   PatientIDs used to filter the import. Empty means
%                         explicit CT series selection across all patients
%   ctSeriesUIDs:         selected CT SeriesInstanceUIDs
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

selectedPatientIDs = matRad_selectPatientIds(fileList, patientList, metadata, matRadCfg);
ctSeriesUIDs = matRad_selectCtSeriesUids(fileList, selectedPatientIDs, metadata, matRadCfg);
rtssRowsByScenario = matRad_selectRtssRows(fileList, selectedPatientIDs, ...
                                           ctSeriesUIDs, metadata);

end

function selectedPatientIDs = matRad_selectPatientIds(fileList, patientList, metadata, matRadCfg)
if isfield(metadata, 'patientID') && ~isempty(metadata.patientID)
    selectedPatientID = matRad_metadataValueToChar(metadata.patientID, ...
                                                   'metadata.patientID', matRadCfg);
    if ~any(strcmp(fileList(:, 3), selectedPatientID))
        matRadCfg.dispError('No DICOM files found for PatientID "%s".', selectedPatientID);
    end
    selectedPatientIDs = {selectedPatientID};
elseif isscalar(patientList)
    selectedPatientIDs = {matRad_metadataValueToChar(patientList{1}, 'PatientID', matRadCfg)};
elseif isfield(metadata, 'ctSeriesUIDs') && ~isempty(metadata.ctSeriesUIDs)
    selectedPatientIDs = {};
else
    matRadCfg.dispError('Multiple PatientIDs found. Specify metadata.patientID or metadata.ctSeriesUIDs.');
end
end

function ctSeriesUIDs = matRad_selectCtSeriesUids(fileList, selectedPatientIDs, metadata, matRadCfg)
ctMask = strcmp(fileList(:, 2), 'CT') & matRad_matchesPatientFilter(fileList(:, 3), selectedPatientIDs);
availableCtSeriesUIDs = matRad_uniqueStableCell(fileList(ctMask, 4));

if isempty(availableCtSeriesUIDs)
    matRadCfg.dispError('No CT SeriesInstanceUID found for the selected PatientID filter.');
end

if isfield(metadata, 'ctSeriesUIDs') && ~isempty(metadata.ctSeriesUIDs)
    ctSeriesUIDs = matRad_metadataValuesToCell(metadata.ctSeriesUIDs, ...
                                               'metadata.ctSeriesUIDs', matRadCfg);
    for seriesIx = 1:numel(ctSeriesUIDs)
        if ~any(strcmp(availableCtSeriesUIDs, ctSeriesUIDs{seriesIx}))
            matRadCfg.dispError(['CT SeriesInstanceUID "%s" was not found ' ...
                                 'for the selected PatientID filter.'], ctSeriesUIDs{seriesIx});
        end
    end
else
    ctSeriesUIDs = availableCtSeriesUIDs;
end
end

function rtssRowsByScenario = matRad_selectRtssRows(fileList, selectedPatientIDs, ctSeriesUIDs, metadata)
rtssRowsByScenario = matRad_selectMultipleDicomRtssRows(fileList, ...
                                                        selectedPatientIDs, ctSeriesUIDs, metadata);
end

function mask = matRad_matchesPatientFilter(patientIDs, selectedPatientIDs)
if isempty(selectedPatientIDs)
    mask = true(size(patientIDs));
else
    mask = ismember(patientIDs, selectedPatientIDs);
end
end

function values = matRad_uniqueStableCell(values)
values = values(:);
uniqueValues = {};

for valueIx = 1:numel(values)
    currentValue = values{valueIx};
    if isempty(uniqueValues) || ~any(strcmp(uniqueValues, currentValue))
        uniqueValues{end + 1, 1} = currentValue; %#ok<AGROW>
    end
end

values = uniqueValues;
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
