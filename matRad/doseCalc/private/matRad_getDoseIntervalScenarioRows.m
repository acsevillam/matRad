function rows = matRad_getDoseIntervalScenarioRows(quantity,scenarioIx,scenarioMap,rowIx,matRad_cfg)
% Return influence rows, mapped to the reference CT scenario if required.

matrix = quantity.matrixCell{scenarioIx};
if ~scenarioMap.mapped
    rows = matrix(rowIx,:) .* quantity.scale;
else
    rows = mapInfluenceRowsToReference(matrix,rowIx,scenarioMap,matRad_cfg) .* ...
        quantity.scale;
end
end

function rows = mapInfluenceRowsToReference(matrix,rowIx,scenarioMap,matRad_cfg)
doseDim = scenarioMap.doseDim;
numBixels = size(matrix,2);
rows = sparse(numel(rowIx),numBixels);

[yRef,xRef,zRef] = ind2sub(doseDim,rowIx(:));
dvfX = scenarioMap.dvfX(rowIx);
dvfY = scenarioMap.dvfY(rowIx);
dvfZ = scenarioMap.dvfZ(rowIx);

xSrc = double(xRef) - dvfX(:);
ySrc = double(yRef) - dvfY(:);
zSrc = double(zRef) - dvfZ(:);

x0 = floor(xSrc);
y0 = floor(ySrc);
z0 = floor(zSrc);
x1 = x0 + 1;
y1 = y0 + 1;
z1 = z0 + 1;

wx1 = xSrc - x0;
wy1 = ySrc - y0;
wz1 = zSrc - z0;
wx0 = 1 - wx1;
wy0 = 1 - wy1;
wz0 = 1 - wz1;

rows = rows + weightedNeighborRows(matrix,doseDim,x0,y0,z0,wx0.*wy0.*wz0,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x1,y0,z0,wx1.*wy0.*wz0,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x0,y1,z0,wx0.*wy1.*wz0,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x1,y1,z0,wx1.*wy1.*wz0,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x0,y0,z1,wx0.*wy0.*wz1,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x1,y0,z1,wx1.*wy0.*wz1,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x0,y1,z1,wx0.*wy1.*wz1,numel(rowIx));
rows = rows + weightedNeighborRows(matrix,doseDim,x1,y1,z1,wx1.*wy1.*wz1,numel(rowIx));

if any(~isfinite(nonzeros(rows)))
    matRad_cfg.dispError('Mapped interval influence rows contain non-finite values.');
end
end

function rows = weightedNeighborRows(matrix,doseDim,x,y,z,weights,numRows)
rows = sparse(numRows,size(matrix,2));
valid = weights ~= 0 & x >= 1 & x <= doseDim(2) & ...
    y >= 1 & y <= doseDim(1) & z >= 1 & z <= doseDim(3);

if ~any(valid)
    return;
end

sourceIx = sub2ind(doseDim,y(valid),x(valid),z(valid));
rows(valid,:) = spdiags(weights(valid),0,nnz(valid),nnz(valid)) * matrix(sourceIx,:);
end
