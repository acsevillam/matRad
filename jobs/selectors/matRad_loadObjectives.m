function [cst,ixTarget,p,ixBody,ixCTV,OARStructSel] = matRad_loadObjectives(run_config,target,cst)

description=run_config.description;
plan_objectives = run_config.plan_objectives;
dp_target_factor=1.0;
switch plan_objectives
    case '1'
        dp_target_factor=8.0;
    case '2'
        dp_target_factor=4.0;
    case '3'
        dp_target_factor=2.0;
    case '4'
        dp_target_factor=1.0;
    case '5'
        dp_target_factor=0.5;
end

dp_start = run_config.dose_pulling_start;

for structure = 1:size(cst,1)
    cst{structure,6}=[];
end

switch description
    case 'prostate'

        for  it = 1:size(cst,1)
            switch cst{it,2}
                case 'BODY'
                    ixBody=it;
                case 'PTV'
                    ixPTV=it;
                case 'CTV'
                    ixCTV=it;
                case 'BLADDER'
                    ixBladder=it;
                case 'RECTUM'
                    ixRectum=it;
            end
        end

        p=78;
        OARStructSel = {'BLADDER','RECTUM'};

        % Body
        if exist('ixBody','var') && ixBody~=0
            cst{ixBody,5}.Priority = 5; % overlap priority for optimization - a lower number corresponds to a higher priority
            cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,0.5*p));
            cst{ixBody,6}{1}.robustness  = 'none';
            cst{ixBody,6}{1}.penaltyPullingRate  = 0.0;
            cst{ixBody,6}{1}.objectivePullingRate{1} = 0.0;
        end

        % CTV
        if exist('ixCTV','var') && ixCTV~=0
            cst{ixCTV,3}  = 'TARGET';
            cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
            cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(dp_target_factor*30,p,100));
            cst{ixCTV,6}{1}.robustness  = 'none';
            cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(0.01,p));
            cst{ixCTV,6}{2}.robustness  = 'none';
            cst{ixCTV,6}{3} = struct(DoseObjectives.matRad_MaxDVH(10,p*1.04,5));
            cst{ixCTV,6}{3}.robustness  = 'none';
            cst{ixCTV,6}{4} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.07,0));
            cst{ixCTV,6}{4}.robustness  = 'none';
            %PTV
            if exist('ixPTV','var') && ixPTV~=0
                % PTV = PTV + CTV
                for i=1:size(cst{ixPTV,4})
                    cst{ixPTV,4}{i} = union(cst{ixCTV,4}{i},cst{ixPTV,4}{i});
                end
            end
        end

        if(target=="PTV")
            % PTV
            if exist('ixPTV','var') && ixPTV~=0
                ixTarget = ixPTV;
                cst{ixTarget,3}  = 'TARGET';
                cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(dp_target_factor*30,p,100));
                cst{ixTarget,6}{1}.robustness  = 'none';
                cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(10,p*1.04,5));
                cst{ixTarget,6}{2}.robustness  = 'none';
                cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.07,0));
                cst{ixTarget,6}{3}.robustness  = 'none';
                %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
            end
        else
            %CTV
            if exist('ixCTV','var') && ixCTV~=0
                ixTarget = ixCTV;
            end
            %PTV
            if exist('ixPTV','var') && ixPTV~=0
                cst{ixPTV,3}  = 'OAR';
            end
        end

        switch plan_objectives
            case {'1','2','3','4','5'}

                % Bladder
                if exist('ixBladder','var') && ixBladder~=0
                    cst{ixBladder,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                    tmpPullingRate{2}  = +0.375;
                    cst{ixBladder,6}{1} = struct(DoseObjectives.matRad_MaxDVH(2,60,0+dp_start*tmpPullingRate{2}));
                    clear tmpPullingRate;
                    cst{ixBladder,6}{1}.penaltyPullingRate  = 0.0;
                    cst{ixBladder,6}{1}.objectivePullingRate{1}  = 0.0;
                    cst{ixBladder,6}{1}.objectivePullingRate{2}  = +0.375;
                    cst{ixBladder,6}{1}.robustness  = 'none';
                    cst{ixBladder,6}{2} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
                    cst{ixBladder,6}{2}.penaltyPullingRate  = 0.0;
                    cst{ixBladder,6}{2}.objectivePullingRate{1}  = 0.0;
                    cst{ixBladder,6}{2}.objectivePullingRate{2}  = 0.0;
                    cst{ixBladder,6}{2}.robustness  = 'none';
                end

                % Rectum
                if exist('ixRectum','var') && ixRectum~=0
                    cst{ixRectum,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                    tmpPullingRate{2}  = +0.5;
                    cst{ixRectum,6}{1} = struct(DoseObjectives.matRad_MaxDVH(2,40,dp_start*tmpPullingRate{2}));
                    clear tmpPullingRate;
                    cst{ixRectum,6}{1}.penaltyPullingRate  = 0.0;
                    cst{ixRectum,6}{1}.objectivePullingRate{1}  = 0.0;
                    cst{ixRectum,6}{1}.objectivePullingRate{2}  = +0.5;
                    cst{ixRectum,6}{1}.robustness  = 'none';
                    cst{ixRectum,6}{2} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
                    cst{ixRectum,6}{2}.penaltyPullingRate  = 0.0;
                    cst{ixRectum,6}{2}.objectivePullingRate{1}  = 0.0;
                    cst{ixRectum,6}{2}.objectivePullingRate{2}  = 0.0;
                    cst{ixRectum,6}{2}.robustness  = 'none';
                end

        end

    case 'breast'
        for  it = 1:size(cst,1)
            switch cst{it,2}
                case 'BODY'
                    ixBody=it;
                case 'PTV'
                    ixPTV=it;
                case 'CTV'
                    ixCTV=it;
                case 'HEART'
                    ixHeart=it;
                case 'LEFT LUNG'
                    ixLeftLung=it;
            end
        end

        p=42.56;
        OARStructSel = {'HEART','LEFT LUNG'};

        % Body
        if exist('ixBody','var') && ixBody~=0
            cst{ixBody,5}.Priority = 5; % overlap priority for optimization - a lower number corresponds to a higher priority
            cst{ixBody,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,0.5*p));
            cst{ixBody,6}{1}.robustness  = 'none';
            cst{ixBody,6}{1}.penaltyPullingRate  = 0.0;
            cst{ixBody,6}{1}.objectivePullingRate{1} = 0.0;
        end

        % CTV
        if exist('ixCTV','var') && ixCTV~=0
            cst{ixCTV,3}  = 'TARGET';
            cst{ixCTV,5}.Priority = 1; % overlap priority for optimization - a lower number corresponds to a higher priority
            cst{ixCTV,6}{1} = struct(DoseObjectives.matRad_MinDVH(dp_target_factor*30,p,100));
            cst{ixCTV,6}{1}.robustness  = 'none';
            cst{ixCTV,6}{2} = struct(DoseObjectives.matRad_SquaredDeviation(0.01,p));
            cst{ixCTV,6}{2}.robustness  = 'none';
            cst{ixCTV,6}{3} = struct(DoseObjectives.matRad_MaxDVH(10,p*1.04,5));
            cst{ixCTV,6}{3}.robustness  = 'none';
            cst{ixCTV,6}{4} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.07,0));
            cst{ixCTV,6}{4}.robustness  = 'none';
            %PTV
            if exist('ixPTV','var') && ixPTV~=0
                % PTV = PTV + CTV
                for i=1:size(cst{ixPTV,4})
                    cst{ixPTV,4}{i} = union(cst{ixCTV,4}{i},cst{ixPTV,4}{i});
                end
            end
        end

        if(target=="PTV")

            % PTV
            if exist('ixPTV','var') && ixPTV~=0
                ixTarget = ixPTV;
                cst{ixTarget,3}  = 'TARGET';
                cst{ixTarget,5}.Priority = 2; % overlap priority for optimization - a lower number corresponds to a higher priority
                cst{ixTarget,6}{1} = struct(DoseObjectives.matRad_MinDVH(dp_target_factor*30,p,100));
                cst{ixTarget,6}{1}.robustness  = 'none';
                cst{ixTarget,6}{2} = struct(DoseObjectives.matRad_MaxDVH(10,p*1.04,5));
                cst{ixTarget,6}{2}.robustness  = 'none';
                cst{ixTarget,6}{3} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.07,0));
                cst{ixTarget,6}{3}.robustness  = 'none';
                %cst{ixTarget,6}{4} = struct(DoseConstraints.matRad_MinMaxDVH(p,95,100));
            end

        else
            %CTV
            if exist('ixCTV','var') && ixCTV~=0
                ixTarget = ixCTV;
            end
            %PTV
            if exist('ixPTV','var') && ixPTV~=0
                cst{ixPTV,3}  = 'OAR';
            end
        end

        switch plan_objectives

            case {'1','2','3','4','5'}

                % Ipsilateral Lung
                if exist('ixLeftLung','var') && ixLeftLung~=0
                    cst{ixLeftLung,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                    tmpPullingRate{2}  = 0.5;
                    cst{ixLeftLung,6}{1} = struct(DoseObjectives.matRad_MaxDVH(2,20,dp_start*tmpPullingRate{2}));
                    clear tmpPullingRate;
                    cst{ixLeftLung,6}{1}.penaltyPullingRate  = 0.0;
                    cst{ixLeftLung,6}{1}.objectivePullingRate{1}  = 0.0;
                    cst{ixLeftLung,6}{1}.objectivePullingRate{2}  = 0.5;
                    cst{ixLeftLung,6}{1}.robustness  = 'none';
                    cst{ixLeftLung,6}{2} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
                    cst{ixLeftLung,6}{2}.penaltyPullingRate  = 0.0;
                    cst{ixLeftLung,6}{2}.objectivePullingRate{1}  = 0.0;
                    cst{ixLeftLung,6}{2}.objectivePullingRate{2}  = 0.0;
                    cst{ixLeftLung,6}{2}.robustness  = 'none';
                    %cst{ixLeftLung,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(20,0,20));
                end

                % Heart
                if exist('ixHeart','var') && ixHeart~=0
                    cst{ixHeart,5}.Priority = 3; % overlap priority for optimization - a lower number corresponds to a higher priority
                    tmpPullingRate{1}  = 0.1;
                    cst{ixHeart,6}{1} = struct(DoseObjectives.matRad_MeanDose(2,dp_start*tmpPullingRate{1},'Quadratic'));
                    clear tmpPullingRate;
                    cst{ixHeart,6}{1}.penaltyPullingRate  = 0.0;
                    cst{ixHeart,6}{1}.objectivePullingRate{1}  = 0.1;
                    cst{ixHeart,6}{1}.objectivePullingRate{2}  = 0.0;
                    cst{ixHeart,6}{1}.robustness  = 'none';
                    cst{ixHeart,6}{2} = struct(DoseObjectives.matRad_MaxDVH(100,p*1.00,0));
                    cst{ixHeart,6}{2}.penaltyPullingRate  = 0.0;
                    cst{ixHeart,6}{2}.objectivePullingRate{1}  = 0.0;
                    cst{ixHeart,6}{2}.objectivePullingRate{2}  = 0.0;
                    cst{ixHeart,6}{2}.robustness  = 'none';
                end

        end

end

end