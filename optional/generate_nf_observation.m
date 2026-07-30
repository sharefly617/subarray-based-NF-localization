function [y, model] = generate_nf_observation(position, transmitSymbol, lambda, sigma2, arrayConfig, channelConfig, includeNoise)
    %GENERATE_NF_OBSERVATION Exact one-snapshot observation for NF_Loc (1)--(5).
    % Inputs: position is 3x1 metres; transmitSymbol is scalar; lambda and d are
    % metres; sigma2 is complex-noise variance. arrayConfig has Nx,Ny,d. The
    % explicit channelConfig has pathLossMode ('common-beta' or 'per-element-beta')
    % and nonnegative scalar GT, GR. Outputs y is (Nx*Ny)x1 in t=(i-1)*Ny+j
    % order; model contains steering, beta, h and alpha/equivalentAlpha.
    % Exact steering/channel phase is exp(-1j*2*pi*r/lambda).
    arguments
    position (3,1) double {mustBeFinite}
    transmitSymbol (1,1) double {mustBeFinite}
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBeNonnegative}
    arrayConfig struct
    channelConfig struct
    includeNoise (1,1) logical = true
end
for field = {'Nx','Ny','d'}
    if ~isfield(arrayConfig, field{1}), error('generate_nf_observation:InvalidArrayConfig', 'arrayConfig.%s is required.', field{1}); end
end
for field = {'pathLossMode','GT','GR'}
    if ~isfield(channelConfig, field{1}), error('generate_nf_observation:InvalidChannelConfig', 'channelConfig.%s is required.', field{1}); end
end
if ~(isscalar(channelConfig.GT) && isscalar(channelConfig.GR) && isfinite(channelConfig.GT) && isfinite(channelConfig.GR) && channelConfig.GT >= 0 && channelConfig.GR >= 0)
    error('generate_nf_observation:InvalidGains', 'GT and GR must be finite nonnegative scalars.');
end
antennaPositions = upa_coordinates(arrayConfig.Nx, arrayConfig.Ny, arrayConfig.d);
[steering, distance] = nearfield_steering(antennaPositions, position, lambda);
mode = string(channelConfig.pathLossMode);
scale = sqrt(channelConfig.GT * channelConfig.GR) * lambda / (4*pi);
switch mode
case "per-element-beta"
    beta = scale ./ distance;                  % Eq. (2), Nx*Ny x 1
    equivalentAlpha = [];
case "common-beta"
    r = norm(position);                         % paper's neglected-variation model
    if r <= eps, error('generate_nf_observation:ZeroRange', 'Common-beta model requires nonzero UE-to-array-center range.'); end
    beta = repmat(scale / r, numel(distance), 1); % Eq. (3), explicit common beta
    equivalentAlpha = beta(1) * transmitSymbol;
otherwise
    error('generate_nf_observation:InvalidPathLossMode', 'pathLossMode must be ''common-beta'' or ''per-element-beta''.');
end
h = beta .* steering;
if includeNoise
    noise = sqrt(sigma2/2) * (randn(size(h)) + 1j*randn(size(h)));
else
    noise = zeros(size(h));
end
y = h * transmitSymbol + noise;
model = struct('antennaPositions', antennaPositions, 'distance', distance, 'steering', steering, ...
    'beta', beta, 'h', h, 'noise', noise, 'pathLossMode', char(mode), ...
    'transmitSymbol', transmitSymbol, 'equivalentAlpha', equivalentAlpha);
end
