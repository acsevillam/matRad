function [cst,ixTarget,p,ixBody,ixCTV] = matRad_loadObjectives(run_config,target,cst)

patient=run_config.description;
plan_objectives = run_config.plan_objectives;

switch patient
    case 'prostate'
        switch plan_objectives
            case '1'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,76));
                cst{ixBody,6}{1}.robustness  = 'none';

                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,35));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,40,40));
                cst{3,6}{1}.robustness  = 'none';

                if(target=="PTV")
                    % CTV
                    p=78;
                    ixCTV = 7;
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
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
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                end

            case '2'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 4; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,76));
                cst{ixBody,6}{1}.robustness  = 'none';

                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,45));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,55));
                cst{3,6}{1}.robustness  = 'none';

                if(target=="PTV")
                    % CTV
                    p=78;
                    ixCTV = 7;
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
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
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                end

            case '3'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,76));
                cst{ixBody,6}{1}.robustness  = 'none';
                
                % Bladder
                cst{2,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{2,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,60,100));
                cst{2,6}{1}.robustness  = 'none';

                % Rectum
                cst{3,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(600,50,100));
                cst{3,6}{1}.robustness  = 'none';
                
                if(target=="PTV")
                    % CTV
                    p=78;
                    ixCTV = 7;
                    cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixCTV,6}{1}.robustness  = 'none';
                    cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixCTV,6}{2}.robustness  = 'none';
                    
                    % PTV
                    ixTarget = 8;
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
                    cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                    cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(200,78,99));
                    cst{ixTarget,6}{1}.robustness  = 'none';
                    cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                    cst{ixTarget,6}{2}.robustness  = 'none';
                    %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                end

        end
    case 'breast'
        switch plan_objectives
            case '1'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
                cst{ixBody,6}{1}.robustness  = 'none';

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(400,20,20));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(250,4));
                cst{4,6}{1}.robustness  = 'none';

                % CTV
                p=42.56;
                ixTarget = 6;
                ixCTV = 6;
                cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                cst{ixCTV,6}{1}.robustness  = 'none';
                cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MinDVH(400,p,95));
                cst{ixTarget,6}{2}.robustness  = 'none';
                cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(400,44,5));
                cst{ixTarget,6}{3}.robustness  = 'none';
                %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));

            case '2'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
                cst{ixBody,6}{1}.robustness  = 'none';

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(400,20,35));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(250,7));
                cst{4,6}{1}.robustness  = 'none';

                % CTV
                p=42.56;
                ixTarget = 6;
                ixCTV = 6;
                cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                cst{ixCTV,6}{1}.robustness  = 'none';
                cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MinDVH(400,p,95));
                cst{ixTarget,6}{2}.robustness  = 'none';
                cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(400,44,5));
                cst{ixTarget,6}{3}.robustness  = 'none';
                %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
                
            case '3'
                % Body
                ixBody = 1;
                cst{ixBody,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,40));
                cst{ixBody,6}{1}.robustness  = 'none';

                % Ipsilateral Lung
                cst{3,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(400,20,100));
                cst{3,6}{1}.robustness  = 'none';
                %cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));

                % Heart
                cst{4,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(250,100));
                cst{4,6}{1}.robustness  = 'none';
                
                % CTV
                p=42.56;
                ixTarget = 6;
                ixCTV = 6;
                cst{ixTarget,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,p));
                cst{ixCTV,6}{1}.robustness  = 'none';
                cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MinDVH(400,p,95));
                cst{ixTarget,6}{2}.robustness  = 'none';
                cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(400,44,5));
                cst{ixTarget,6}{3}.robustness  = 'none';
                %cst{ixTarget,6}{3} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
        end  
end
    
end

