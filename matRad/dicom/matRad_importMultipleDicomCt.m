function [ct, cst] = matRad_importMultipleDicomCt(filesDir, metadata)
% matRad wrapper function to import multiple DICOM CT series into one multi-scenario ct
%
% call
%   [ct,cst] = matRad_importMultipleDicomCt(filesDir,metadata)
%
% input
%   filesDir:       folder with DICOM CT and RT structure set files. Can be
%                   an absolute path or a path relative to matRadRoot
%   metadata:       struct with import options
%                   metadata.resolution:   [x y z] CT import resolution in mm
%                   metadata.useDoseGrid:  optional scalar logical to import
%                                          CT on the RT dose grid
%                   metadata.patientID:     optional PatientID to import
%                   metadata.ctSeriesUIDs:  optional CT SeriesInstanceUIDs
%                                          to import as scenarios. If
%                                          metadata.patientID is omitted,
%                                          explicit CT series may span
%                                          PatientIDs
%                   metadata.rtssUID:       optional RTSTRUCT SOPInstanceUID
%                                          used for all scenarios
%                   metadata.rtssUIDs:      optional RTSTRUCT SOPInstanceUIDs
%                                          matched one-to-one with
%                                          metadata.ctSeriesUIDs
%
% output
%   ct:             matRad ct multi-scenario struct
%   cst:            matRad cst multi-scenario struct
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

if ~exist('metadata', 'var') || ~isstruct(metadata) || ~isfield(metadata, 'resolution') || ...
        ~isnumeric(metadata.resolution) || numel(metadata.resolution) ~= 3
    matRad_cfg.dispError('metadata.resolution must be a numeric [x y z] vector.');
end

if isstring(filesDir)
    filesDir = char(filesDir);
end

if isfolder(filesDir)
    dicomFolder = filesDir;
else
    if startsWith(filesDir, filesep)
        filesDir = filesDir(2:end);
    end
    dicomFolder = fullfile(matRad_cfg.matRadRoot, filesDir);
end

if ~isfolder(dicomFolder)
    matRad_cfg.dispError('DICOM import folder not found: %s', dicomFolder);
end

% Scan the container directory using the canonical DICOM importer in dev.
dicomImporter = matRad_DicomImporter(dicomFolder);
fileList = dicomImporter.allfiles;
patientList = dicomImporter.patients;

[selectedPatientIDs, ctSeriesUIDs, rtssRowsByScenario] = ...
    matRad_resolveMultipleDicomCtImportSelection(fileList, patientList, metadata);

% Get CT scenario information
ct.numOfCtScen = numel(ctSeriesUIDs);
ct.timeStamp = datestr(clock);
ct.resolution.x = metadata.resolution(1);
ct.resolution.y = metadata.resolution(2);
ct.resolution.z = metadata.resolution(3);

files.resx = ct.resolution.x;
files.resy = ct.resolution.y;
files.resz = ct.resolution.z;

if isfield(metadata, 'useDoseGrid') && ~isempty(metadata.useDoseGrid)
    validUseDoseGrid = isscalar(metadata.useDoseGrid) && ...
        (islogical(metadata.useDoseGrid) || ...
         (isnumeric(metadata.useDoseGrid) && any(metadata.useDoseGrid == [0 1])));
    if ~validUseDoseGrid
        matRad_cfg.dispError('metadata.useDoseGrid must be a scalar logical value.');
    end
    files.useDoseGrid = logical(metadata.useDoseGrid);
else
    files.useDoseGrid = false;
end
files.rtdose = dicomImporter.importFiles.rtdose;
if files.useDoseGrid && isempty(files.rtdose)
    matRad_cfg.dispError('metadata.useDoseGrid requires at least one RTDOSE file in the DICOM import folder.');
end

% Initialize ct and cst structs for each scenario
tmpCtOriginal = cell(1, ct.numOfCtScen);
tmpCstOriginal = cell(1, ct.numOfCtScen);

for ctPhase = 1:ct.numOfCtScen

    ctMask = strcmp(fileList(:, 2), 'CT') & ...
        matRad_matchesPatientFilter(fileList(:, 3), selectedPatientIDs) & ...
        strcmp(fileList(:, 4), ctSeriesUIDs{ctPhase});
    files.ct = fileList(ctMask, :);
    files.rtss = rtssRowsByScenario{ctPhase};

    % Import DICOM information into ct and cst structs.
    [tmpCtOriginal{ctPhase}, tmpCstOriginal{ctPhase}] = ...
        matRad_importSingleCtScenario(dicomImporter, files);

end

[ct, cst] = matRad_finalizeMultipleDicomCtImport(ct, tmpCtOriginal, tmpCstOriginal);

end

function [ct, cst] = matRad_importSingleCtScenario(dicomImporter, files)
dicomImporter.importFiles = struct();
dicomImporter.importFiles.ct = files.ct;
dicomImporter.importFiles.rtss = files.rtss;
dicomImporter.importFiles.rtplan = {};
dicomImporter.importFiles.rtdose = files.rtdose;
dicomImporter.importFiles.resx = files.resx;
dicomImporter.importFiles.resy = files.resy;
dicomImporter.importFiles.resz = files.resz;
dicomImporter.importFiles.useImportGrid = files.useDoseGrid;

dicomImporter.ct = [];
dicomImporter.cst = [];
dicomImporter.stf = [];
dicomImporter.pln = [];
dicomImporter.resultGUI = [];
dicomImporter.importCT = [];
dicomImporter.importRtss = [];
dicomImporter.importRTDose = [];
dicomImporter.ImportGrid = [];

dicomImporter.matRad_importDicom();

ct = dicomImporter.ct;
cst = dicomImporter.cst;
end

function mask = matRad_matchesPatientFilter(patientIDs, selectedPatientIDs)
if isempty(selectedPatientIDs)
    mask = true(size(patientIDs));
else
    mask = ismember(patientIDs, selectedPatientIDs);
end
end
