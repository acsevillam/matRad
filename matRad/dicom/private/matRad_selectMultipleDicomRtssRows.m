function rtssRowsByScenario = matRad_selectMultipleDicomRtssRows(fileList,selectedPatientIDs,ctSeriesUIDs,metadata)
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

matRad_cfg = MatRad_Config.instance();

rtssMask = strcmp(fileList(:,2),'RTSTRUCT') & matchesPatientFilter(fileList(:,3),selectedPatientIDs);
availableRtssRows = fileList(rtssMask,:);
rtssRowsByScenario = cell(1,numel(ctSeriesUIDs));

if isfield(metadata,'rtssUIDs') && ~isempty(metadata.rtssUIDs)
    rtssUIDs = metadataValuesToCell(metadata.rtssUIDs,'metadata.rtssUIDs',matRad_cfg);
    if numel(rtssUIDs) ~= numel(ctSeriesUIDs)
        matRad_cfg.dispError('metadata.rtssUIDs must contain one RTSTRUCT SOPInstanceUID per CT scenario.');
    end
    for ctPhase = 1:numel(ctSeriesUIDs)
        rtssRowsByScenario{ctPhase} = selectRtssRow(fileList,rtssMask, ...
            rtssUIDs{ctPhase},matRad_cfg);
    end
elseif isfield(metadata,'rtssUID') && ~isempty(metadata.rtssUID)
    rtssUID = metadataValueToChar(metadata.rtssUID,'metadata.rtssUID',matRad_cfg);
    rtssRow = selectRtssRow(fileList,rtssMask,rtssUID,matRad_cfg);
    rtssRowsByScenario(:) = {rtssRow};
else
    rtssRowsByScenario = selectDefaultRtssRows(availableRtssRows, ...
        rtssRowsByScenario,matRad_cfg);
end
end

function rtssRowsByScenario = selectDefaultRtssRows(availableRtssRows, ...
    rtssRowsByScenario,matRad_cfg)
if isempty(availableRtssRows)
    rtssRowsByScenario(:) = {[]};
    return;
end
if size(availableRtssRows,1) ~= 1
    matRad_cfg.dispError('Multiple RTSTRUCT files found. Specify metadata.rtssUID or metadata.rtssUIDs.');
end
rtssRowsByScenario(:) = {availableRtssRows};
end

function rtssRow = selectRtssRow(fileList,rtssMask,rtssUID,matRad_cfg)
rtssMatch = rtssMask & strcmp(fileList(:,4),rtssUID);
if ~any(rtssMatch)
    matRad_cfg.dispError('RTSTRUCT SOPInstanceUID "%s" was not found for the selected PatientID filter.',rtssUID);
end
rtssRow = fileList(rtssMatch,:);
if size(rtssRow,1) ~= 1
    matRad_cfg.dispError('RTSTRUCT SOPInstanceUID "%s" matched %d files.',rtssUID,size(rtssRow,1));
end
end

function mask = matchesPatientFilter(patientIDs,selectedPatientIDs)
if isempty(selectedPatientIDs)
    mask = true(size(patientIDs));
else
    mask = ismember(patientIDs,selectedPatientIDs);
end
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
