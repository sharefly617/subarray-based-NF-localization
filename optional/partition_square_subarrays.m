function partition = partition_square_subarrays(Nx, Ny, partsPerAxis)
    %PARTITION_SQUARE_SUBARRAYS Uniform square UPA subarrays, defaulting to nine.
    % Inputs Nx,Ny and partsPerAxis are integer counts. Output is the regular
    % partition from partition_upa; for 60,60,3 it contains nine 20x20 blocks in
    % m=(mx-1)*partsPerAxis+my order and global t=(i-1)*Ny+j index order.
    arguments
    Nx (1,1) double {mustBeInteger, mustBePositive}
    Ny (1,1) double {mustBeInteger, mustBePositive}
    partsPerAxis (1,1) double {mustBeInteger, mustBePositive} = 3
end
if Nx ~= Ny
    error('partition_square_subarrays:NonSquareArray', 'A uniform square partition requires Nx=Ny.');
end
partition = partition_upa(Nx, Ny, partsPerAxis, partsPerAxis);
if partition.subNx ~= partition.subNy
    error('partition_square_subarrays:NonSquareSubarray', 'Subarrays must be square.');
end
end
