function [selectedPatientIDs,ctSeriesUIDs,rtssRowsByScenario] = matRad_resolveMultipleDicomCtImportSelection(fileList,patientList,metadata)
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

matRad_cfg = MatRad_Config.instance();

selectedPatientIDs = selectPatientIDs(fileList,patientList,metadata,matRad_cfg);
ctSeriesUIDs = selectCtSeriesUIDs(fileList,selectedPatientIDs,metadata,matRad_cfg);
rtssRowsByScenario = selectRtssRows(fileList,selectedPatientIDs, ...
    ctSeriesUIDs,metadata,matRad_cfg);

end

function selectedPatientIDs = selectPatientIDs(fileList,patientList,metadata,matRad_cfg)
if isfield(metadata,'patientID') && ~isempty(metadata.patientID)
    selectedPatientID = metadataValueToChar(metadata.patientID, ...
        'metadata.patientID',matRad_cfg);
    if ~any(strcmp(fileList(:,3),selectedPatientID))
        matRad_cfg.dispError('No DICOM files found for PatientID "%s".',selectedPatientID);
    end
    selectedPatientIDs = {selectedPatientID};
elseif isscalar(patientList)
    selectedPatientIDs = {metadataValueToChar(patientList{1},'PatientID',matRad_cfg)};
elseif isfield(metadata,'ctSeriesUIDs') && ~isempty(metadata.ctSeriesUIDs)
    selectedPatientIDs = {};
else
    matRad_cfg.dispError('Multiple PatientIDs found. Specify metadata.patientID or metadata.ctSeriesUIDs.');
end
end

function ctSeriesUIDs = selectCtSeriesUIDs(fileList,selectedPatientIDs,metadata,matRad_cfg)
ctMask = strcmp(fileList(:,2),'CT') & matchesPatientFilter(fileList(:,3),selectedPatientIDs);
availableCtSeriesUIDs = uniqueStableCell(fileList(ctMask,4));

if isempty(availableCtSeriesUIDs)
    matRad_cfg.dispError('No CT SeriesInstanceUID found for the selected PatientID filter.');
end

if isfield(metadata,'ctSeriesUIDs') && ~isempty(metadata.ctSeriesUIDs)
    ctSeriesUIDs = metadataValuesToCell(metadata.ctSeriesUIDs, ...
        'metadata.ctSeriesUIDs',matRad_cfg);
    for seriesIx = 1:numel(ctSeriesUIDs)
        if ~any(strcmp(availableCtSeriesUIDs,ctSeriesUIDs{seriesIx}))
            matRad_cfg.dispError(['CT SeriesInstanceUID "%s" was not found ' ...
                'for the selected PatientID filter.'],ctSeriesUIDs{seriesIx});
        end
    end
else
    ctSeriesUIDs = availableCtSeriesUIDs;
end
end

function rtssRowsByScenario = selectRtssRows(fileList,selectedPatientIDs,ctSeriesUIDs,metadata,~)
rtssRowsByScenario = matRad_selectMultipleDicomRtssRows(fileList, ...
    selectedPatientIDs,ctSeriesUIDs,metadata);
end

function mask = matchesPatientFilter(patientIDs,selectedPatientIDs)
if isempty(selectedPatientIDs)
    mask = true(size(patientIDs));
else
    mask = ismember(patientIDs,selectedPatientIDs);
end
end

function values = uniqueStableCell(values)
values = values(:);
uniqueValues = {};

for valueIx = 1:numel(values)
    currentValue = values{valueIx};
    if isempty(uniqueValues) || ~any(strcmp(uniqueValues,currentValue))
        uniqueValues{end + 1,1} = currentValue; %#ok<AGROW>
    end
end

values = uniqueValues;
end

function values = metadataValuesToCell(values,fieldName,matRad_cfg)
if ischar(values)
    values = {values};
elseif isstring(values)
    values = cellstr(values(:));
elseif iscell(values)
    values = values(:);
else
    matRad_cfg.dispError('%s must be a character vector, string, or cell array of character vectors.',fieldName);
end

for valueIx = 1:numel(values)
    values{valueIx} = metadataValueToChar(values{valueIx},fieldName,matRad_cfg);
end
end

function value = metadataValueToChar(value,fieldName,matRad_cfg)
if isstring(value) && isscalar(value)
    value = char(value);
elseif iscell(value) && isscalar(value)
    value = metadataValueToChar(value{1},fieldName,matRad_cfg);
elseif ~ischar(value)
    matRad_cfg.dispError('%s must be a character vector or string scalar.',fieldName);
end
end
