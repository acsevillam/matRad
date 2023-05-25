function [cst] = matRad_scaleDoseObjectives(cst,structSel,scale_factor)

for itSelStructure = 1:size(structSel,2)
    for  itStructure = 1:size(cst,1)
        if(strcmp(cst{itStructure,2},structSel{itSelStructure}))
            for itObjective = 1:size(cst{itStructure,6},2)
                objective = cst{itStructure,6}{itObjective};
                if(objective.dosePulling)
                    for itObjParam = 1:size(objective.parameters,2)
                        if ~isempty(objective.objectivePullingRate{itObjParam}) && objective.objectivePullingRate{itObjParam}~=0
                            objective.parameters{itObjParam}=scale_factor*objective.parameters{itObjParam};
                        end
                    end
                end
                cst{itStructure,6}{itObjective} = objective;
            end
        end
    end
end

end

