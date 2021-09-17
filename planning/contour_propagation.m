%%
clear 
clc

%%
load('images/patient_2/patient_2_multiScen.mat');
ct_mama =ct;
cst_mama=cst;

%%
pyramLevels = [1 1 1 1 1 1 1 4 4 4];
initialItera = [100 100 100 100 100 100 100 300 100 100];
smoothLevels = [1.5 1.5 2.6 1.8 2.9 1.3 3.5 7.5 2.6 1.3];

%%
[ct_mama, cst_mama] = matRad_contourPropagation(ct_mama,cst_mama,pyramLevels,initialItera,smoothLevels);

%%
% Calculo DICE para los contornos
[num_Struct, ~] = size(cst_mama);
dice_value = cell(1,num_Struct);
for struct_iterator = 1:num_Struct
    dice_value{struct_iterator}=zeros(1,ct_mama.numOfCtScen);
    for scen = 1:ct_mama.numOfCtScen
        
        switch scen
            case 1
                load('patient_2_scen_1_resized.mat');
            case 2
                load('patient_2_scen_2_resized.mat');
            case 3
                load('patient_2_scen_3_resized.mat');
            case 4
                load('patient_2_scen_4_resized.mat');
            case 5
                load('patient_2_scen_5_resized.mat');
        end
        
        cube_Original = zeros(ct.cubeDim);       
        struct_Original_cst  = cst{struct_iterator,4}{1,1};         
        [x,y,z] = ind2sub(ct.cubeDim,struct_Original_cst);
        for i=1:length(x)
            cube_Original(x(i),y(i),z(i)) = 1;
        end
        
        cube_Estimated = zeros(ct_mama.cubeDim);       
        struct_Estimated_cst = cst_mama{struct_iterator,4}{1,scen};         
        [xe,ye,ze] = ind2sub(ct_mama.cubeDim,struct_Estimated_cst);
        for j=1:length(xe)
            cube_Estimated(xe(j),ye(j),ze(j)) = 1;
        end
        fprintf('Calculo del coeficiente DICE de la estructura %d en el escenario %d.\n',struct_iterator,scen);
        dice_value{struct_iterator}(scen)=matRad_dice3D(cube_Estimated,cube_Original);
        fprintf('%d. \n',dice_value{struct_iterator}(scen));
    end    
end
%% show the deformation vector field
slice = 30; %ctPhase = 9;   % select a specific slice and and ct phase to plot the vector field
[yDim,~,zDim,a] = size(ct_mama.dvf{2});

[mY,mZ]      = meshgrid(1:zDim,1:yDim);

figure,
for ctPhase = 2:ct_mama.numOfCtScen 
   clf;
   yVectorField = squeeze(ct_mama.dvf{ctPhase}(slice,:,:,1));  % retrieve the deformation vector field in x-direction of slice 25
   zVectorField = squeeze(ct_mama.dvf{ctPhase}(slice,:,:,3));  % retrieve the deformation vector field in y-direction of slice 25
   quiver(mZ,mY,zVectorField,yVectorField); title(['deformation vector field of phase ' num2str(ctPhase)]),
   
   set(gca,'XLim',[1 256]);set(gca,'YLim',[1 65]);
   pause(0.9);   
end

%%
% graficas de los valores Dice para cada escenario
figure('Renderer', 'painters', 'Position', [10 10 1200 600]);

subplot(2,3,1);
bar(1:5,dice_value{1,1});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('Skin');

subplot(2,3,2);
bar(1:5,dice_value{1,2});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('PULMON DERECHO');

subplot(2,3,3);
bar(1:5,dice_value{1,3});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('PULMON IZQUIERDO');

subplot(2,3,4);
bar(1:5,dice_value{1,5});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('SENO CONTRALATERAL');

subplot(2,3,5);
bar(1:5,dice_value{1,6});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('CTV-TARGET');

subplot(2,3,6);
bar(1:5,dice_value{1,4});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('CORAZON');

%%
ct = ct_mama;
cst = cst_mama;

%%
save('images/patient_2/patient_2_multiScen_contourned.mat','ct','cst');
