function [tmp_ct_out, commonGrid, cropIdx] = matRad_cropCtToOverlap(tmp_ct_in, tol)
% Recorta múltiples CT (uno por escenario) a la intersección 3D de sus ejes.
%
% INPUT
%   tmp_ct_in{ctPhase} : structs ct individuales (como salen de matRad_importDicom)
%   tol                : tolerancia numérica (opcional, p.ej. 1e-5)
%
% OUTPUT
%   tmp_ct_out         : mismos ct pero recortados a la intersección
%   commonGrid         : struct con x,y,z comunes tras el recorte
%   cropIdx            : cell por escenario con campos Ix,Iy,Iz usados para recortar
%
% NOTAS
%  - Requiere que todos los escenarios tengan el mismo paso en cada eje
%    (dentro de la tolerancia). Si no, habría que re-muestrear (fuera de alcance aquí).
%  - Se asume convención matRad: size(cubeHU) = [numY, numX, numZ].
%
if nargin < 2, tol = 1e-5; end

nScen = numel(tmp_ct_in);
assert(nScen >= 1, 'Se requiere al menos un escenario');

% Recolectar ejes y pasos
allX = cellfun(@(c)c.x, tmp_ct_in, 'uni', 0);
allY = cellfun(@(c)c.y, tmp_ct_in, 'uni', 0);
allZ = cellfun(@(c)c.z, tmp_ct_in, 'uni', 0);

dx = cellfun(@(v) median(diff(v)), allX);
dy = cellfun(@(v) median(diff(v)), allY);
dz = cellfun(@(v) median(diff(v)), allZ);

% Validar pasos (mismo spacing por eje)
if max(dx) - min(dx) > tol
    error('Los escenarios tienen pasos diferentes en X. Re-muestreo requerido.');
end
if max(dy) - min(dy) > tol
    error('Los escenarios tienen pasos diferentes en Y. Re-muestreo requerido.');
end
if max(dz) - min(dz) > tol
    error('Los escenarios tienen pasos diferentes en Z. Re-muestreo requerido.');
end

% Calcular intersección de rangos (en coordenadas físicas)
xMinCommon = max(cellfun(@(v) min(v), allX));
xMaxCommon = min(cellfun(@(v) max(v), allX));
yMinCommon = max(cellfun(@(v) min(v), allY));
yMaxCommon = min(cellfun(@(v) max(v), allY));
zMinCommon = max(cellfun(@(v) min(v), allZ));
zMaxCommon = min(cellfun(@(v) max(v), allZ));

if ~(xMinCommon < xMaxCommon && yMinCommon < yMaxCommon && zMinCommon < zMaxCommon)
    error('No hay región de superposición entre escenarios en al menos un eje.');
end

% Para la malla común elegimos los puntos que caen dentro del rango en TODOS,
% pero usando los vectores del primer escenario como referencia "canónica"
% solo para construir commonGrid (los datos se recortan por escenario).
xRef = tmp_ct_in{1}.x; yRef = tmp_ct_in{1}.y; zRef = tmp_ct_in{1}.z;
xCommon = xRef(xRef >= xMinCommon - tol & xRef <= xMaxCommon + tol);
yCommon = yRef(yRef >= yMinCommon - tol & yRef <= yMaxCommon + tol);
zCommon = zRef(zRef >= zMinCommon - tol & zRef <= zMaxCommon + tol);

% Recortar por escenario
tmp_ct_out = tmp_ct_in;
cropIdx = cell(1,nScen);

for s = 1:nScen
    xS = tmp_ct_in{s}.x; yS = tmp_ct_in{s}.y; zS = tmp_ct_in{s}.z;

    Ix = find(xS >= xMinCommon - tol & xS <= xMaxCommon + tol);
    Iy = find(yS >= yMinCommon - tol & yS <= yMaxCommon + tol);
    Iz = find(zS >= zMinCommon - tol & zS <= zMaxCommon + tol);

    if isempty(Ix) || isempty(Iy) || isempty(Iz)
        error('Escenario %d no tiene índices dentro de la intersección.', s);
    end

    % Recortar el volumen HU (orden [Y, X, Z] en matRad)
    HU = tmp_ct_in{s}.cubeHU{1};
    HUc = HU(Iy, Ix, Iz);

    % Actualizar ct
    tmp_ct_out{s}.x = xS(Ix);
    tmp_ct_out{s}.y = yS(Iy);
    tmp_ct_out{s}.z = zS(Iz);
    tmp_ct_out{s}.cubeHU{1} = HUc;
    tmp_ct_out{s}.cubeDim = [numel(Ix),numel(Iy),numel(Iz)];

    cropIdx{s} = struct('Ix', Ix, 'Iy', Iy, 'Iz', Iz, ...
                        'origDim', tmp_ct_in{s}.cubeDim);
end

% Construir una malla común final como intersección estricta de todos,
% validando que todos terminaron con mismas longitudes (debería ocurrir)
nx = cellfun(@(c) numel(c.x), tmp_ct_out);
ny = cellfun(@(c) numel(c.y), tmp_ct_out);
nz = cellfun(@(c) numel(c.z), tmp_ct_out);
if ~(all(nx == nx(1)) && all(ny == ny(1)) && all(nz == nz(1)))
    error('Tras recorte, las dimensiones aún difieren. Revisar tolerancias.');
end

% Elegimos la malla del primer ct recortado como común (todas iguales ya)
commonGrid = struct('x', tmp_ct_out{1}.x, ...
                    'y', tmp_ct_out{1}.y, ...
                    'z', tmp_ct_out{1}.z);
end
