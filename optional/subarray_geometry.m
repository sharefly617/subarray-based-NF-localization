function geometry = subarray_geometry(antennaPositions, partition, position)
    %SUBARRAY_GEOMETRY Centers, distances and paper direction cosines per subarray.
    % antennaPositions is (Nx*Ny)x3 metres in global t order; position is 3x1 m.
    % Outputs centers Mx3 (m), distances Mx1 (m), theta Mx2 dimensionless, with
    % theta=(p_U-p_BS,m)/norm(p_U-p_BS,m) in x and y, precisely NF_Loc Eq. (9).
    arguments
    antennaPositions (:,3) double {mustBeFinite}
    partition struct
    position (3,1) double {mustBeFinite}
end
if size(antennaPositions, 1) ~= partition.Nx * partition.Ny
    error('subarray_geometry:LengthMismatch', 'antennaPositions must contain Nx*Ny rows.');
end
M = numel(partition.subarrays); centers = zeros(M, 3);
for m = 1:M, centers(m, :) = mean(antennaPositions(partition.subarrays(m).index, :), 1); end
offset = position.' - centers; distances = sqrt(sum(offset.^2, 2));
if any(distances <= eps(max(1, norm(position))))
    error('subarray_geometry:CoincidentCenter', 'UE coincides with a subarray center.');
end
theta = offset(:, 1:2) ./ distances;
geometry = struct('centers', centers, 'distance', distances, 'theta', theta);
end
