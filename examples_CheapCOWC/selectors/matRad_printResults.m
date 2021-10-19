function [results] = matRad_printResults(run_config,results,cst,pln,dqi_sampled)

patient=run_config.description;

switch patient
    case 'prostate'
        % Bladder
        disp('Bladder evaluation');
        results.structures{2,1}=cst{2,2};
        V60Bladder=sprintf('V60: %.2f [%.2f - %.2f] %%\n',dqi_sampled{1,1}(2).V_1_54Gy*100,dqi_sampled{1,2}(2).V_1_54Gy*100,dqi_sampled{1,3}(2).V_1_54Gy*100); disp(V60Bladder);
        results.structures{2,2}.V_60.nom=dqi_sampled{1,1}(2).V_1_54Gy*100;
        results.structures{2,2}.V_60.min=dqi_sampled{1,2}(2).V_1_54Gy*100;
        results.structures{2,2}.V_60.max=dqi_sampled{1,3}(2).V_1_54Gy*100;

        % Rectum
        disp('Rectum evaluation');
        results.structures{3,1}=cst{3,2};
        V40Rectum=sprintf('V40: %.2f [%.2f - %.2f] %%',dqi_sampled{1,1}(3).V_1_03Gy*100,dqi_sampled{1,2}(3).V_1_03Gy*100,dqi_sampled{1,3}(3).V_1_03Gy*100); disp(V40Rectum);
        results.structures{3,2}.V_40.nom=dqi_sampled{1,1}(3).V_1_03Gy*100;
        results.structures{3,2}.V_40.min=dqi_sampled{1,2}(3).V_1_03Gy*100;
        results.structures{3,2}.V_40.max=dqi_sampled{1,3}(3).V_1_03Gy*100;

        V50Rectum=sprintf('V50: %.2f [%.2f - %.2f] %%\n',dqi_sampled{1,1}(3).V_1_28Gy*100,dqi_sampled{1,2}(3).V_1_28Gy*100,dqi_sampled{1,3}(3).V_1_28Gy*100); disp(V50Rectum);
        results.structures{3,2}.V_50.nom=dqi_sampled{1,1}(3).V_1_28Gy*100;
        results.structures{3,2}.V_50.min=dqi_sampled{1,2}(3).V_1_28Gy*100;
        results.structures{3,2}.V_50.max=dqi_sampled{1,3}(3).V_1_28Gy*100;
        
	case 'breast'
        % Contralateral Lung
        disp('Contralateral lung evaluation');
        results.structures{2,1}=cst{2,2};
        V50ConLung=sprintf('V50: %.2f [%.2f - %.2f] %%',dqi_sampled{1,1}(2).V_3_12Gy*100,dqi_sampled{1,2}(2).V_3_12Gy*100,dqi_sampled{1,3}(2).V_3_12Gy*100); disp(V50ConLung);
        results.structures{2,2}.V_50.nom=dqi_sampled{1,1}(2).V_3_12Gy*100;
        results.structures{2,2}.V_50.min=dqi_sampled{1,2}(2).V_3_12Gy*100;
        results.structures{2,2}.V_50.max=dqi_sampled{1,3}(2).V_3_12Gy*100;

        D5ConLung=sprintf('D5: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(2).D_5*pln.numOfFractions,dqi_sampled{1,2}(2).D_5*pln.numOfFractions,dqi_sampled{1,3}(2).D_5*pln.numOfFractions); disp(D5ConLung);
        results.structures{2,2}.D_5.nom=dqi_sampled{1,1}(2).D_5*pln.numOfFractions;
        results.structures{2,2}.D_5.min=dqi_sampled{1,2}(2).D_5*pln.numOfFractions;
        results.structures{2,2}.D_5.max=dqi_sampled{1,3}(2).D_5*pln.numOfFractions;

        % Ipsilateral Lung
        disp('Ipsilateral lung evaluation');
        results.structures{3,1}=cst{3,2};
        V20IpsLung=sprintf('V20: %.2f [%.2f - %.2f] %%',dqi_sampled{1,1}(3).V_1_25Gy*100,dqi_sampled{1,2}(3).V_1_25Gy*100,dqi_sampled{1,3}(3).V_1_25Gy*100); disp(V20IpsLung);
        results.structures{3,2}.V_20.nom=dqi_sampled{1,1}(3).V_1_25Gy*100;
        results.structures{3,2}.V_20.min=dqi_sampled{1,2}(3).V_1_25Gy*100;
        results.structures{3,2}.V_20.max=dqi_sampled{1,3}(3).V_1_25Gy*100;

        D20IpsLung=sprintf('D20: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(3).D_20*pln.numOfFractions,dqi_sampled{1,2}(3).D_20*pln.numOfFractions,dqi_sampled{1,3}(3).D_20*pln.numOfFractions); disp(D20IpsLung);
        results.structures{3,2}.D_20.nom=dqi_sampled{1,1}(3).D_20*pln.numOfFractions;
        results.structures{3,2}.D_20.min=dqi_sampled{1,2}(3).D_20*pln.numOfFractions;
        results.structures{3,2}.D_20.max=dqi_sampled{1,3}(3).D_20*pln.numOfFractions;

        % Heart
        disp('Heart evaluation');
        results.structures{4,1}=cst{4,2};
        DMeanHeart=sprintf('DMean: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(4).mean*pln.numOfFractions,dqi_sampled{1,2}(4).mean*pln.numOfFractions,dqi_sampled{1,3}(4).mean*pln.numOfFractions); disp(DMeanHeart);
        results.structures{4,2}.D_mean.nom=dqi_sampled{1,1}(4).mean*pln.numOfFractions;
        results.structures{4,2}.D_mean.min=dqi_sampled{1,2}(4).mean*pln.numOfFractions;
        results.structures{4,2}.D_mean.max=dqi_sampled{1,3}(4).mean*pln.numOfFractions;

        % Contralateral Breast
        disp('Contralateral Breast evaluation');
        results.structures{5,1}=cst{5,2};
        DMaxConBreast=sprintf('DMax: %.2f [%.2f - %.2f] Gy\n',dqi_sampled{1,1}(5).max*pln.numOfFractions,dqi_sampled{1,2}(5).max*pln.numOfFractions,dqi_sampled{1,3}(5).max*pln.numOfFractions); disp(DMaxConBreast);
        results.structures{5,2}.D_mean.nom=dqi_sampled{1,1}(5).max*pln.numOfFractions;
        results.structures{5,2}.D_mean.min=dqi_sampled{1,2}(5).max*pln.numOfFractions;
        results.structures{5,2}.D_mean.max=dqi_sampled{1,3}(5).max*pln.numOfFractions;    
end

end

