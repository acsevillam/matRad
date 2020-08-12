
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
for struct = 1:num_Struct
    dice_value{struct}=zeros(1,ct_mama.numOfCtScen);
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
        struct_Original_cst  = cst{struct,4}{1,1};         
        [x,y,z] = ind2sub(ct.cubeDim,struct_Original_cst);
        for i=1:length(x)
            cube_Original(x(i),y(i),z(i)) = 1;
        end
        
        cube_Estimated = zeros(ct_mama.cubeDim);       
        struct_Estimated_cst = cst_mama{struct,4}{1,scen};         
        [xe,ye,ze] = ind2sub(ct_mama.cubeDim,struct_Estimated_cst);
        for j=1:length(xe)
            cube_Estimated(xe(j),ye(j),ze(j)) = 1;%ct_mama.cubeHU{scen}(xe(j),ye(j),ze(j));
        end
        fprintf('Calculo del coeficiente DICE de la estructura %d en el escenario %d.\n',struct,scen);
%         cube_Estimated; 
%         common = (cube_Estimated & cube_Original);
%          a = sum(common(:));
%          b = sum(cube_Estimated (:));
%          c = sum(cube_Original(:));
%          dice_value{struct}(scen) = 2*a/(b+c);
        dice_value{struct}(scen)=Dice3D(cube_Estimated,cube_Original);
%         fprintf('%d %d %d %d.\n',a,b,c,dice_value{struct}(scen));
        fprintf('%d. \n',dice_value{struct}(scen));
    end    
end
%% show the deformation vector field
slice = 25; %ctPhase = 9;   % select a specific slice and and ct phase to plot the vector field
[xDim,yDim,~,a] = size(ct_mama.dvf{2});

[mX,mY]      = meshgrid(1:xDim,1:yDim);

figure,
for ctPhase = 2:ct_mama.numOfCtScen 
   clf;
   xVectorField = squeeze(ct_mama.dvf{ctPhase}(:,:,slice,1));  % retrieve the deformation vector field in x-direction of slice 25
   yVectorField = squeeze(ct_mama.dvf{ctPhase}(:,:,slice,2));  % retrieve the deformation vector field in y-direction of slice 25
   quiver(mX,mY,yVectorField,xVectorField); title(['deformation vector field of phase ' num2str(ctPhase)]),
   set(gca,'XLim',[60 130]);set(gca,'YLim',[60 130]);
   % flip y axis to be consistent with the previous plot
   ax = gca; ax.YDir = 'reverse';
   pause(0.9);  
end

%%
% graficas de los valores Dice para cada escenario
figure(1)
% subplot(2,2,1);
bar(1:5,dice_value{1,1});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('Skin');

figure(2)
% subplot(2,2,2);
bar(1:5,dice_value{1,2});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('PULMON DERECHO');

figure(3)
% subplot(2,2,3);
bar(1:5,dice_value{1,3});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('PULMON IZQUIERDO');

figure(4)
% subplot(2,2,1);
bar(1:5,dice_value{1,5});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('SENO CONTRALATERAL');

figure(5)
% subplot(2,2,2);
bar(1:5,dice_value{1,6});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('CTV-TARGET');

figure(6)
% subplot(2,2,3);
bar(1:5,dice_value{1,4});
xlabel('Escenario');ylabel('Coeficiente DICE');
title('CORAZON');

%%
ct = ct_mama;
cst = cst_mama;

%%
save('images/patient_2/patient_2_multiScen_contourned.mat','ct','cst');
