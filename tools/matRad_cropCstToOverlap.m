function cst_out = matRad_cropCstToOverlap(cst_in, cropIdx, ct_original, ct_cropped)
% Recorta/remapea máscaras de estructuras al subvolumen común.
% Soporta dos layouts:
%  A) per-escenario: cst_in = 1×S cell; cst_in{1,s} = [N×6] cell; col 4 = {1×1 cell -> linIdx}
%  B) multisecenario en un solo cst: cst_in = [N×6] cell; col 4 = {1×S cell -> linIdx}

cst_out = cst_in;

% Helpers robustos para dims [Ny Nx Nz]
getDims = @(ct) ( ...
    isfield(ct,'cubeDim') && isnumeric(ct.cubeDim) && numel(ct.cubeDim)==3 ...
    ) .* double(ct.cubeDim(:)') + ...
    (~(isfield(ct,'cubeDim') && isnumeric(ct.cubeDim) && numel(ct.cubeDim)==3)) ...
    .* double(size(ct.cubeHU{1}));

nScen = numel(ct_original);
for s = 1:nScen
    if s > numel(cst_out), continue; end
    cst_scen = cst_out{1,s};
    if ~iscell(cst_scen) || isempty(cst_scen), continue; end

    % Índices de recorte
    Ix = cropIdx{s}.Ix(:)'; Iy = cropIdx{s}.Iy(:)'; Iz = cropIdx{s}.Iz(:)';

    % Dims originales y recortadas
    d0 = getDims(ct_original{s});  d0 = d0(:)';  % [Ny0 Nx0 Nz0]
    d1 = getDims(ct_cropped{s});   d1 = d1(:)';  % [Ny1 Nx1 Nz1]
    [Ny0,Nx0,Nz0] = deal(d0(1),d0(2),d0(3));
    [Ny1,Nx1,Nz1] = deal(d1(1),d1(2),d1(3));

    % Mapas -> índices relativos
    mapX = zeros(1, Nx0, 'uint32'); mapX(Ix) = uint32(1:numel(Ix));
    mapY = zeros(1, Ny0, 'uint32'); mapY(Iy) = uint32(1:numel(Iy));
    mapZ = zeros(1, Nz0, 'uint32'); mapZ(Iz) = uint32(1:numel(Iz));

    [numStruct, numCols] = size(cst_scen); %#ok<ASGLU>
    for st = 1:numStruct
        if size(cst_scen,2) < 4 || ~iscell(cst_scen{st,4}) || isempty(cst_scen{st,4})
            continue;
        end
        % En este layout, col 4 es {1×1 cell -> linIdx}
        linCell = cst_scen{st,4};
        if ~iscell(linCell) || numel(linCell) < 1, continue; end
        linIdx = linCell{1};

        if isempty(linIdx)
            cst_scen{st,4}{1} = [];  % deja vacío
            continue;
        elseif islogical(linIdx)
            if ~isequal(size(linIdx), [Ny0 Nx0 Nz0])
                warning('Estructura %d, escenario %d: máscara lógica tamaño inesperado. Se omite.', st, s);
                cst_scen{st,4}{1} = [];
                continue;
            end
            linIdx = find(linIdx);
        elseif ~isvector(linIdx)
            warning('Estructura %d, escenario %d: máscara no soportada. Se omite.', st, s);
            cst_scen{st,4}{1} = [];
            continue;
        end

        [iy, ix, iz] = ind2sub([Ny0,Nx0,Nz0], linIdx(:));
        inCrop = (mapY(iy)>0) & (mapX(ix)>0) & (mapZ(iz)>0);

        if ~any(inCrop)
            cst_scen{st,4}{1} = [];
            continue;
        end

        iyc = double(mapY(iy(inCrop)));
        ixc = double(mapX(ix(inCrop)));
        izc = double(mapZ(iz(inCrop)));
        linNew = sub2ind([Ny1,Nx1,Nz1], iyc, ixc, izc);

        cst_scen{st,4}{1} = linNew(:);
    end

    cst_out{1,s} = cst_scen;
end

end
