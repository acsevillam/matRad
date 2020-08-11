%%
clear
clc

%%
load('patient_2_scen_1_resized.mat');
ct1 = ct;
cst1 = cst;

%%
for i = [2 3 4 5]
    
    switch i
        case 2
            load('patient_2_scen_2_resized.mat');
        case 3
            load('patient_2_scen_3_resized.mat');
        case 4
            load('patient_2_scen_4_resized.mat');
        case 5
            load('patient_2_scen_5_resized.mat');
    end
    
    ct1.cubeHU{i} = ct.cubeHU{1};
    ct1.dicomInfo(i) = ct.dicomInfo;
    ct1.dicomMeta(i) = ct.dicomMeta;  
    ct1.numOfCtScen = i;
end

%%
clear ct cst

%%
ct = ct1;
cst = cst1;

%%
clear ct1 cst1 i

%%
save('images/patient_2/patient_2_multiScen.mat','ct','cst');