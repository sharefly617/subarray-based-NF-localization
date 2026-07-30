function [y, steering, noise, antennaPositions] = generate_nearfield_snapshot(position, alpha, lambda, sigma2, arrayConfig, includeNoise)
    %GENERATE_NEARFIELD_SNAPSHOT Eq. (5) alpha-only exact complex-Gaussian snapshot.
    % Inputs: position is 3x1 metres; alpha is scalar; lambda is metres;
    % sigma2 is complex-noise variance; arrayConfig has Nx, Ny, d (metres).
    % Outputs y, steering and noise are (Nx*Ny)x1 in t=(i-1)*Ny+j order.
    % The caller supplies alpha directly; this function intentionally performs no
    % path-loss calculation. Use generate_nf_observation for the explicit
    % 'common-beta' or 'per-element-beta' models of NF_Loc Eqs. (1)--(5).
    arguments
    position (3,1) double {mustBeFinite}
    alpha (1,1) double {mustBeFinite}
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBeNonnegative}
    arrayConfig struct
    includeNoise (1,1) logical = true
end
for field = {'Nx','Ny','d'}
    if ~isfield(arrayConfig, field{1})
        error('generate_nearfield_snapshot:InvalidArrayConfig', 'arrayConfig.%s is required.', field{1});
    end
end
antennaPositions = upa_coordinates(arrayConfig.Nx, arrayConfig.Ny, arrayConfig.d);
[steering, ~] = nearfield_steering(antennaPositions, position, lambda);
if includeNoise
    noise = sqrt(sigma2 / 2) * (randn(size(steering)) + 1j * randn(size(steering)));
else
    noise = zeros(size(steering));
end
y = alpha * steering + noise;
end
