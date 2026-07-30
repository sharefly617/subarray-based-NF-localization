function ySubarrays = extract_subarray_observations(y, partition)
    %EXTRACT_SUBARRAY_OBSERVATIONS Gather a full UPA vector into subarray vectors.
    % y must be (Nx*Ny)x1 in t=(i-1)*Ny+j order; output is Mx1 cell, whose
    % entries are (subNx*subNy)x1 complex vectors in the same local ordering.
    arguments
    y (:,1) double
    partition struct
end
if numel(y) ~= partition.Nx * partition.Ny
    error('extract_subarray_observations:LengthMismatch', 'y must have Nx*Ny elements.');
end
M = numel(partition.subarrays); ySubarrays = cell(M, 1);
for m = 1:M, ySubarrays{m} = y(partition.subarrays(m).index); end
end
