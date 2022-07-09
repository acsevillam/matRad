function [cst,ixTarget,p,ixBody,ixCTV] = matRad_loadObjectives(run_config,target,cst)

patient=run_config.description;
plan_objectives = run_config.plan_objectives;
radiationMode = run_config.radiationMode;

switch radiationMode
    case 'protons'
        switch patient
            case 'prostate'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(600,76));
                cst{ixBody,6}{1}.robustness  = 'none';
                
                if(target=="PTV")
                    % CTV
                    p=78;
                    ixCTV = 7;
                    cst{ixCTV,3}  = 'TARGET';
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,95));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(200,81,5));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(200,83.5,0));
                    cst{ixTarget,6}{3}.robustness  = 'none';
                    %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                else
                    % CTV
                    p=78;
                    ixCTV = 7;
                    ixTarget = 7;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                    
                    %PTV
                    ixPTV = 8;
                    cst{ixPTV,3}  = 'OAR';
                end
                
                switch plan_objectives
                    case '1'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,20));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,25,25));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '2'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,25));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,30,30));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '3'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,55,30));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,35,35));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '4'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,35));
                        cst{5,6}{1}.robustness  = 'none';
                       
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,40,40));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '5'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,65,40));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,45));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '6'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,100));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,100));
                        cst{4,6}{1}.robustness  = 'none';
                        
                end
                               
        end
    case 'photons'
        switch patient
            case 'prostate'
                
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(600,76));
                cst{ixBody,6}{1}.robustness  = 'none';
                
                if(target=="PTV")
                    % CTV
                    p=78;
                    ixCTV = 7;
                    cst{ixCTV,3}  = 'TARGET';
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(600,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,95));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(200,81,5));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(200,83.5,0));
                    cst{ixTarget,6}{3}.robustness  = 'none';
                    %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                else
                    % CTV
                    p=78;
                    ixCTV = 7;
                    ixTarget = 7;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(600,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                    
                    %PTV
                    ixPTV = 8;
                    cst{ixPTV,3}  = 'OAR';
                end
                
                switch plan_objectives
                    case '1'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,20));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,25,25));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '2'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,25));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,30,30));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '3'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,55,30));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,35,35));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '4'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,35));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,40,40));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '5'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,65,40));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,45));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '6'
                        
                        % Bladder
                        cst{5,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,100));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Rectum
                        cst{4,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,100));
                        cst{4,6}{1}.robustness  = 'none';
                        
                end
                
            case 'breast'
                
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(600,40.56));
                cst{ixBody,6}{1}.robustness  = 'none';
                
                if(target=="PTV")
                    % CTV
                    p=42.56;
                    ixCTV = 7;
                    cst{ixCTV,3}  = 'TARGET';
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,95));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(200,p*1.04,5));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(200,p*1.07,0));
                    cst{ixTarget,6}{3}.robustness  = 'none';
                    %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                else
                    % CTV
                    p=42.56;
                    ixCTV = 7;
                    ixTarget = 7;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                    
                    %PTV
                    ixPTV = 8;
                    cst{ixPTV,3}  = 'OAR';
                end
                
                switch plan_objectives
                    
                    case '1'
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,14,14));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(10,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '2'
                        
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,16,16));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(10,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '3'
                        
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,18,18));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(10,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '4'
                        
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,20,20));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(2.5,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '5'
                        
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,22,22));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(1,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                    case '6'
                        
                        % Ipsilateral Lung
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,20,100));
                        cst{3,6}{1}.robustness  = 'none';
                        %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                        
                        % Heart
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(0.5,0));
                        cst{4,6}{1}.robustness  = 'none';
                        
                end
                
            case 'H&N'
                
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(600,48));
                cst{ixBody,6}{1}.robustness  = 'none';
                
                % Mandible
                cst{2,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,49,0));
                cst{2,6}{1}.robustness  = 'none';
                
                if(target=="PTV")
                    % CTV
                    p=50.00;
                    ixCTV = 8;
                    cst{ixCTV,3}  = 'TARGET';
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 9;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,95));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(200,p*1.04,5));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(200,p*1.07,0));
                    cst{ixTarget,6}{3}.robustness  = 'none';
                    %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                    
                else
                    % CTV
                    p=50;
                    ixCTV = 8;
                    ixTarget = 8;
                    cst{ixTarget,3}  = 'TARGET';
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,p,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                    
                    %PTV
                    ixPTV = 9;
                    cst{ixPTV,3}  = 'OAR';
                end
                
                switch plan_objectives
                    
                    case '1'
                        
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,8));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,8));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,20,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,33,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,20,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,12));
                        cst{7,6}{2}.robustness  = 'none';
                        
                    case '2'
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,9));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,9));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,24,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,35,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,24,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,14));
                        cst{7,6}{2}.robustness  = 'none';
                        
                    case '3'
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,10));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,10));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,28,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,37,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,28,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,16));
                        cst{7,6}{2}.robustness  = 'none';
                        
                    case '4'
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,11));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,11));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,32,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,39,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,32,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,18));
                        cst{7,6}{2}.robustness  = 'none';
                        
                    case '5'
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,13));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,13));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,36,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,41,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,36,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,20));
                        cst{7,6}{2}.robustness  = 'none';
                        
                    case '6'
                        % ParotidGland  Left
                        cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,100));
                        cst{3,6}{1}.robustness  = 'none';
                        
                        % ParotidGland  Right
                        cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(300,100));
                        cst{4,6}{1}.robustness  = 'none';
                        
                        % SpinalCord
                        cst{5,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{5,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,100,2));
                        cst{5,6}{1}.robustness  = 'none';
                        
                        % Brainstem
                        cst{6,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{6,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,100,2));
                        cst{6,6}{1}.robustness  = 'none';
                        
                        % Aux
                        cst{7,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                        cst{7,6}{1} = struct(DoseObjectives.matRad_MaxDVH(300,100,2));
                        cst{7,6}{1}.robustness  = 'none';
                        cst{7,6}{2} = struct(DoseObjectives.matRad_MeanDose(300,100));
                        cst{7,6}{2}.robustness  = 'none';
                        
                        
                end
        end
        
end

end

