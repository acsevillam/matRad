function [cst,ixTarget,p,ixBody,ixCTV] = matRad_loadObjectives(run_config,target,cst)

patient=run_config.description;
plan_objectives = run_config.plan_objectives;

switch patient
    case 'prostate'
        
        % Body
        ixBody = 1;
        cst{ixBody,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
        cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,76));
        cst{ixBody,6}{1}.robustness  = 'none';
        %cst{ixBody,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,70,0));
        %cst{ixBody,6}{2}.robustness  = 'none';
        
        % Bladder
        %cst{2,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,78,0));
        %cst{2,6}{2}.robustness  = 'none';
        
        % Rectum
        %cst{3,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,78,0));
        %cst{3,6}{2}.robustness  = 'none';

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
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,20));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,25,25));
                cst{3,6}{1}.robustness  = 'none';

            case '2'
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,25));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,30,30));
                cst{3,6}{1}.robustness  = 'none';

            case '3'
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,55,30));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,35,35));
                cst{3,6}{1}.robustness  = 'none';
                
            case '4'
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,35));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,40,40));
                cst{3,6}{1}.robustness  = 'none';
                
            case '5'
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,65,40));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,45,45));
                cst{3,6}{1}.robustness  = 'none';
                
            case '6'
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,100));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,100));
                cst{3,6}{1}.robustness  = 'none';

        end
        
    case 'breast'

        % Body
        ixBody = 1;
        cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
        cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
        cst{ixBody,6}{1}.robustness  = 'none';
        %cst{ixBody,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,70,0));
        %cst{ixBody,6}{2}.robustness  = 'none';
        
        % Ipsilateral Lung
        %cst{3,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,42.56,0));
        %cst{3,6}{2}.robustness  = 'none';
        
        % Heart
        %cst{4,6}{2} = struct(DoseObjectives.matRad_MaxDVH(600,42.56,0));
        %cst{4,6}{2}.robustness  = 'none';
        
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
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,2.5));
                cst{4,6}{1}.robustness  = 'none';

            case '2'

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,16,16));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,3.0));
                cst{4,6}{1}.robustness  = 'none';
                
            case '3'

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,18,18));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,3.5));
                cst{4,6}{1}.robustness  = 'none';
                
            case '4'

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,20,20));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,4));
                cst{4,6}{1}.robustness  = 'none';
                
            case '5'

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,22,22));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,5));
                cst{4,6}{1}.robustness  = 'none';
                
            case '6'

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,20,100));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(600,100));
                cst{4,6}{1}.robustness  = 'none';

        end  
end
    
end

