function scenProb = matRad_getScenProb(pln,phaseProb)

scenProb = zeros(pln.multScen.totNumScen,1);

for l = 1:pln.multScen.totNumScen
 
    [ctScen,shiftScen,RangeScen] = deal(pln.multScen.linearMask(l,1),pln.multScen.linearMask(l,2),pln.multScen.linearMask(l,3));
    shiftScenMask = find(squeeze(pln.multScen.scenMask(1,:,:)));
    indProb = sub2ind([pln.multScen.totNumShiftScen pln.multScen.totNumRangeScen],shiftScen,RangeScen);
    
    numCtScen = nnz(pln.multScen.scenMask(:,shiftScen,RangeScen));
    if(numCtScen>1)
        scenProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb))*phaseProb(ctScen);
    else
        scenProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb));
    end
end

end
