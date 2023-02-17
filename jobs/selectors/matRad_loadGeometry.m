function [ct,cst] = matRad_loadGeometry(run_config)

AcquisitionType = run_config.AcquisitionType;

switch AcquisitionType
    case 'mat'
        % Import 3D CT
        load([run_config.caseID '.mat'],'ct','cst');

    case 'dicom'
        % Import 4D CT
        metadata.resolution = run_config.resolution;
        [ct,cst] = matRad_importMultipleDicomCt(['/jobs/images/breast/' run_config.caseID '/dicom'],metadata);
        clear 'metadata';

end

end