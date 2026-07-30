function partition = partition_upa(Nx, Ny, partsX, partsY)
    %PARTITION_UPA Regular non-overlapping rectangular partition of a UPA.
    % Inputs are positive integer antenna and partition counts.  Each returned
    % index vector is column-shaped and follows global t=(i-1)*Ny+j ordering.
    % Output partition.subarrays is (partsX*partsY)x1, ordered m=(mx-1)*partsY+my.
    arguments
    Nx (1,1) double {mustBeInteger, mustBePositive}
    Ny (1,1) double {mustBeInteger, mustBePositive}
    partsX (1,1) double {mustBeInteger, mustBePositive}
    partsY (1,1) double {mustBeInteger, mustBePositive}
end
if mod(Nx, partsX) ~= 0 || mod(Ny, partsY) ~= 0
    error('partition_upa:NonDivisible', 'Nx and Ny must be divisible by partsX and partsY.');
end
subNx = Nx / partsX; subNy = Ny / partsY;
globalIndices = repelem((0:Nx-1).' * Ny, 1, Ny) + repmat(1:Ny, Nx, 1);
template = struct('index', [], 'iRange', [], 'jRange', [], 'Nx', subNx, 'Ny', subNy, 'mx', [], 'my', []);
subarrays = repmat(template, partsX * partsY, 1);
for mx = 1:partsX
    for my = 1:partsY
        m = (mx - 1) * partsY + my;
        iRange = (mx - 1) * subNx + (1:subNx);
        jRange = (my - 1) * subNy + (1:subNy);
        localIndices = globalIndices(iRange, jRange);
        subarrays(m) = struct('index', reshape(localIndices.', [], 1), ...
            'iRange', iRange, 'jRange', jRange, 'Nx', subNx, 'Ny', subNy, 'mx', mx, 'my', my);
    end
end
partition = struct('Nx', Nx, 'Ny', Ny, 'partsX', partsX, 'partsY', partsY, ...
    'subNx', subNx, 'subNy', subNy, 'subarrays', subarrays);
end
