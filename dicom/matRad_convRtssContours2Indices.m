function indices = matRad_convRtssContours2Indices(structure,ct)
% matRad function to convert a polygon segmentation from an rt structure
% set into a binary segmentation as required within matRad's cst struct
% 
% call
%   indices = matRad_convRtssContours2Indices(contPoints,ct)
%
% input
%   structure:      information about a single structure
%   ct:             matRad ct struct where the binary segmentations will
%                   be aligned to
%
% output
%   indicies:       indices of voxels of the ct cube that are inside the
%                   contour
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voiCube = zeros(ct.cubeDim);

% loop over all closed contour items
for i = 1:size(structure.item,2)

    if ~isempty(structure.item(i).points)

        dicomCtSlicePos = unique(structure.item(i).points(:,3));
        
        if numel(dicomCtSlicePos) > 1 || isempty(dicomCtSlicePos)
            error('Contour defined over multiple planes!');
        end

        % Vector de posiciones z del CT
        zCT = ct.dicomInfo.SlicePositions(:);
        [zCT, sortIdx] = sort(zCT,'ascend');

        % ¿Hay SliceThickness por corte?
        hasPerSliceST = isfield(ct.dicomInfo,'SliceThickness') && ...
                        numel(ct.dicomInfo.SliceThickness) == numel(zCT);

        if hasPerSliceST
            stVec = ct.dicomInfo.SliceThickness(:);
            stVec = stVec(sortIdx);
        end

        % Espaciado medio entre cortes (dzMed), si se puede calcular
        dzMed = NaN;
        if numel(zCT) >= 2
            dz = abs(diff(zCT));
            dz = dz(dz > eps);
            if ~isempty(dz)
                dzMed = median(dz);
            end
        end

        % Corte más cercano a la z del contorno
        [dzMin, idxNear] = min(abs(zCT - dicomCtSlicePos));

        % Tolerancia: mínimo 0.5 mm; si hay dzMed, usa 20% de dzMed
        tol = 0.5;
        if ~isnan(dzMed) && dzMed > 0
            tol = max(0.5, 0.2 * dzMed);
        end

        if isempty(idxNear) || dzMin > tol
            error('Slice at z=%0.3f mm not found within tolerance (min Δ=%.3f > tol=%.3f).', ...
                   dicomCtSlicePos, dzMin, tol);
        end

        % Grosor del corte a usar
        if hasPerSliceST
            dicomCtSliceThickness = stVec(idxNear);
        else
            % Estimar grosor con vecinos si es posible; si no, usar dzMed
            if numel(zCT) >= 3 && idxNear > 1 && idxNear < numel(zCT)
                dicomCtSliceThickness = 0.5 * (zCT(idxNear+1) - zCT(idxNear-1));
            elseif numel(zCT) >= 2
                if idxNear == 1
                    dicomCtSliceThickness = (zCT(2) - zCT(1));
                else
                    dicomCtSliceThickness = (zCT(end) - zCT(end-1));
                end
            elseif ~isnan(dzMed)
                dicomCtSliceThickness = dzMed;
            else
                dicomCtSliceThickness = []; % sin forma de estimar
            end
            dicomCtSliceThickness = abs(dicomCtSliceThickness);
        end

        %Sanity check
        msg = checkSliceThickness(dicomCtSliceThickness);
        if ~isempty(msg)
            error('Slice Thickness of slice at %f could not be identified: %s',dicomCtSlicePos,msg);
        end
        
        slicesInMatradCt = find(dicomCtSlicePos+dicomCtSliceThickness/2 > ct.z & dicomCtSlicePos-dicomCtSliceThickness/2 <= ct.z);
        
        coords1 = interp1(ct.x,1:ct.cubeDim(2),structure.item(i).points(:,1),'linear','extrap');
        coords2 = interp1(ct.y,1:ct.cubeDim(1),structure.item(i).points(:,2),'linear','extrap');
        
        binIn = poly2mask(coords1,coords2,ct.cubeDim(1),ct.cubeDim(2));
        
        % loop over all slices in matRad ct
        for j = 1:numel(slicesInMatradCt)
            voiCube(:,:,slicesInMatradCt(j)) = voiCube(:,:,slicesInMatradCt(j)) | binIn;
        end
        
    end
    
end

indices = find(voiCube(:));

end

function msg = checkSliceThickness(dicomCtSliceThickness)
    if isempty(dicomCtSliceThickness)
        msg = 'Slice could not be identified (empty)';
    elseif ~isscalar(dicomCtSliceThickness)
        msg = 'Slice thickness not unique';
    elseif ~isnumeric(dicomCtSliceThickness)
        msg = 'unexpected value';
    else
        msg = '';
    end
end
