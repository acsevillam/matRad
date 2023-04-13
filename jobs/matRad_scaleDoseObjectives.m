function [cst] = matRad_scaleDoseObjectives(cst,structSel,scale_factor)

for itSelStructure = 1:size(structSel,2)
    for  itStructure = 1:size(cst,1)
        if(strcmp(cst{itStructure,2},structSel{itSelStructure}))
            for itObjective = 1:size(cst{itStructure,6},2)
                for itObjParam = 1:size(cst{itStructure,6}{itObjective}.parameters,2)
                    if cst{itStructure,6}{itObjective}.objectivePullingRate{itObjParam}~=0
                        cst{itStructure,6}{itObjective}.parameters{itObjParam}=scale_factor*cst{itStructure,6}{itObjective}.parameters{itObjParam};
                    end
                end
            end
        end
    end
end

end

