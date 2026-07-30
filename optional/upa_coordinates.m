function positions = upa_coordinates(Nx, Ny, d)
    %UPA_COORDINATES Centered UPA positions, ordered t=(i-1)*Ny+j.
    arguments
    Nx (1,1) double {mustBeInteger, mustBePositive}
    Ny (1,1) double {mustBeInteger, mustBePositive}
    d (1,1) double {mustBePositive}
end
ix = (1:Nx).' - (Nx + 1) / 2;
iy = (1:Ny).' - (Ny + 1) / 2;
% i is the outer loop and j is the inner loop: t=(i-1)*Ny+j.
positions = [repelem(ix * d, Ny), repmat(iy * d, Nx, 1), zeros(Nx * Ny, 1)];
end
