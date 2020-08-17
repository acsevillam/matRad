
%%
clear
clc

%%

for scen = [1 2 3 4 5]
    
    
    switch scen
            case 1
                load('images/patient_2/patient_2_scen_1.mat');
            case 2
                load('images/patient_2/patient_2_scen_2.mat');
            case 3
                load('images/patient_2/patient_2_scen_3.mat');
            case 4
                load('images/patient_2/patient_2_scen_4.mat');
            case 5
                load('images/patient_2/patient_2_scen_5.mat');
    end
    
    [num_Struct, ~] = size(cst);
	for struct = 1:num_Struct
        if isempty(cst{struct,4}{1,1}) == true
            fprintf('Estructura %d vacia. \n',struct);
            fprintf('Eliminando estructura %d. \n',struct);
            cst(struct,:) = [];
        end
	end
        
    if(ct.cubeDim(3)>65)
    
        ct.z=ct.z(11:75);
        ct.cubeHU{1}=ct.cubeHU{1}(:,:,11:75);
        ct.dicomInfo(1).SlicePositions=ct.dicomInfo(1).SlicePositions(11:75);
        ct.dicomInfo(1).SliceThickness=ct.dicomInfo(1).SliceThickness(11:75);
        ct.dicomInfo(1).ImagePositionPatient = [-1.816830000000000e+02;-4.740810000000000e+02;-846];
        ct.cubeDim = [256 256 65];
        
        [num_Struct, ~] = size(cst);
        for struct = 1:num_Struct
            %Obtaining the fixed cubic structure from the linear indices
            cube_fixed_Scena = zeros([256 256 75]);
            struct_Fixed_cst = cst{struct,4}{1,1}; 
            [x,y,z] = ind2sub([256 256 75],struct_Fixed_cst);

            %The HU value of the corresponding tomography is assigned to each position
            for i=1:length(x)
                cube_fixed_Scena(x(i),y(i),z(i)) = 1;
            end

            cube_fixed_Scena=cube_fixed_Scena(:,:,11:75);
            cst{struct,4}{1,1} = find(cube_fixed_Scena);
        end
    
    end
    
	switch scen
        case 1
            save('images/patient_2/patient_2_scen_1_resized.mat','ct','cst');
        case 2
            save('images/patient_2/patient_2_scen_2_resized.mat','ct','cst');
        case 3
            save('images/patient_2/patient_2_scen_3_resized.mat','ct','cst');
        case 4
            save('images/patient_2/patient_2_scen_4_resized.mat','ct','cst');
        case 5
            save('images/patient_2/patient_2_scen_5_resized.mat','ct','cst');
    end
    
    clear cst ct cube_fixed_Scena i num_Struct struct struct_Fixed_cst x y z;
    
end

clear scen