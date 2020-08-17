function  [ct,cst] = matRad_contourPropagation(ct,cst,pyramLevels,initialItera,smoothLevels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function non rigid registration to estimate the contour
% propagation of a tomography sequence, from ct and cst structure set.
% 
% call
%   [ct,cst] = matRad_ContourPropagation(ct,cst,pyramLevels,initialItera,smoothLevels)
%
% input
%   ct:            matRad ct structure 
%   cst:           matRad cst struct 
%   pyramLevels:   number of multi-resolution image pyramid levels to use,
%                  specified as the positive integer scalar
%   initialItera:  number of iterations, specified as a positive integer
%                  scalar 
%   smoothLevels:  smoothing applied at each iteration. This parameter
%                  controls the amount of diffusion-like regularization, it
%                  applies the standard deviation of the Gaussian
%                  smoothing. Values typically are in the range [0.5 , 3.0]
%
% output
%   ct:            matRad struct 
%   cst:           matRad struct of CTs segmented
% References
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set type
ct.dvfType = 'pull'; % push or pull

%   Non rigid registration demons-based. Calculates the DVF(Displacement
%   Vector Field) that models the transformation.
ct.dvf = cell(1,ct.numOfCtScen);
for Scen = 1:ct.numOfCtScen
    fprintf('Registrando escenario %d.\n',Scen);
    [ct.dvf{Scen},~] = imregdemons(ct.cubeHU{1},ct.cubeHU{Scen},initialItera(Scen),'PyramidLevels',pyramLevels(Scen),'AccumulatedFieldSmoothing',smoothLevels(Scen));
end
    

[num_Struct, ~] = size(cst);
for struct_iterator = 1:num_Struct

    if isempty(cst{struct_iterator,4}{1,1}) == false

%           Obtaining the fixed cubic structure from the linear indices
        cube_fixed_Scena = zeros(ct.cubeDim);
        struct_Fixed_cst = cst{struct_iterator,4}{1,1}; 
        [x,y,z] = ind2sub(ct.cubeDim,struct_Fixed_cst);

%           The HU value of the corresponding tomography is assigned to each position
        for j=1:length(x)
            cube_fixed_Scena(x(j),y(j),z(j)) = 1;%ct.cubeHU{1}(x(j),y(j),z(j));
        end

%           The DVF transformation is applied and the linear values are found
        for scen_moving = 2:ct.numOfCtScen
            fprintf('Propagando contornos de la estructura %d escenario %d. \n',struct_iterator,scen_moving);
            scen_estimated = imwarp(cube_fixed_Scena,ct.dvf{scen_moving}); 
            cst{struct_iterator,4}{1,scen_moving} = find(scen_estimated);
        end
    else
        fprintf('Estructura %d vacia. \n',struct_iterator);
        fprintf('Eliminando estructura %d. \n',struct_iterator);
        cst(struct_iterator,:) = [];
    end
end     


