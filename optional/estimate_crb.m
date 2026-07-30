function result = estimate_crb(position, alpha, lambda, sigma2, arrayConfig, options)
    %ESTIMATE_CRB Exact near-field position CRB with unknown complex gain.
    %   The real parameter vector is eta=[x,y,z,angle(alpha),abs(alpha)]^T.
    %   This implements (44) and J=(2/sigma2)*real(D'*D) from the paper.
    %   The returned position CRB eliminates the two gain parameters by a
    %   Schur complement; alpha is deliberately never treated as known.

    arguments
    position (3,1) double {mustBeFinite}
    alpha (1,1) double {mustBeFinite}
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBePositive}
    arrayConfig struct
    options.maxCondition (1,1) double {mustBePositive} = 1e12
    options.minEigenvalue (1,1) double {mustBeNonnegative} = 1e-10
end

if abs(alpha) == 0
    error('estimate_crb:ZeroGain', 'abs(alpha) must be nonzero because its phase is a nuisance parameter.');
end
required = {'Nx', 'Ny', 'd'};
for k = 1:numel(required)
    if ~isfield(arrayConfig, required{k})
        error('estimate_crb:InvalidArrayConfig', 'arrayConfig.%s is required.', required{k});
    end
end

% Paper system model: the UPA is centered at the origin in the xOy plane.
antennaPositions = upa_coordinates(arrayConfig.Nx, arrayConfig.Ny, arrayConfig.d);
delta = antennaPositions - position.';             % p_BS - p_U
[steering, distance] = nearfield_steering(antennaPositions, position, lambda);
mu = alpha * steering;

% Five columns follow the order [x y z angle(alpha) abs(alpha)].
D = zeros(numel(distance), 5);
D(:, 1:3) = (1j * 2 * pi / lambda) .* (delta ./ distance) .* mu;
D(:, 4) = 1j * mu;
D(:, 5) = (alpha / abs(alpha)) * steering;

J = (2 / sigma2) * real(D' * D);
J = (J + J.') / 2;                                 % remove round-off skew
eigJ = eig(J);
conditionNumber = cond(J);
% The full five-parameter condition number is unit-dependent (metres versus
% radians versus gain), so it is recorded but not used as an invertibility
% test.  The nuisance and Schur-complement blocks below are scale-relevant.
if ~isfinite(conditionNumber) || min(eigJ) < -options.minEigenvalue
    error('estimate_crb:IllConditionedFIM', ...
        'FIM is invalid (cond=%g, minEig=%g).', conditionNumber, min(eigJ));
end

Jpp = J(1:3, 1:3);
Jpn = J(1:3, 4:5);
Jnn = J(4:5, 4:5);
conditionNuisance = cond(Jnn);
if ~isfinite(conditionNuisance) || conditionNuisance > options.maxCondition
    error('estimate_crb:IllConditionedNuisance', 'Nuisance FIM is ill-conditioned (cond=%g).', conditionNuisance);
end
% QR projection is algebraically the Schur complement above, but avoids
% subtracting two large nearly equal matrices in this geometry.
Dreal = [real(D); imag(D)];
A = Dreal(:, 1:3); B = Dreal(:, 4:5);
[Q, R] = qr(B, 0);
if min(abs(diag(R))) <= eps(norm(R, 'fro'))
    error('estimate_crb:IllConditionedNuisance', 'Nuisance derivative columns are rank deficient.');
end
Aprojected = A - Q * (Q.' * A);
effectivePositionFIM = (2 / sigma2) * (Aprojected.' * Aprojected);
effectivePositionFIM = (effectivePositionFIM + effectivePositionFIM.') / 2;
conditionPosition = cond(effectivePositionFIM);
eigPosition = eig(effectivePositionFIM);
if ~isfinite(conditionPosition) || conditionPosition > options.maxCondition || min(eigPosition) <= options.minEigenvalue
    error('estimate_crb:IllConditionedPositionFIM', ...
        'Schur-complement position FIM is not safely invertible (cond=%g, minEig=%g).', ...
    conditionPosition, min(eigPosition));
end
positionCRB = effectivePositionFIM \ eye(3);
positionCRB = (positionCRB + positionCRB.') / 2;

% Full CRB is computed after diagonal equilibration.  This avoids falsely
% rejecting a valid bound merely because eta mixes metres, radians and gain.
parameterScale = 1 ./ sqrt(diag(J));
if any(~isfinite(parameterScale))
    error('estimate_crb:IllConditionedFIM', 'FIM has a zero or invalid diagonal.');
end
scaledFIM = (parameterScale .* J) .* parameterScale.';
scaledFIM = (scaledFIM + scaledFIM.') / 2;
scaledConditionNumber = cond(scaledFIM);
if ~isfinite(scaledConditionNumber) || scaledConditionNumber > options.maxCondition
    error('estimate_crb:IllConditionedFIM', 'Equilibrated FIM is ill-conditioned (cond=%g).', scaledConditionNumber);
end
scaledCRB = scaledFIM \ eye(5);
scaledInverseResidual = norm(scaledFIM * scaledCRB - eye(5), 'fro') / sqrt(5);
fullCRB = (parameterScale .* scaledCRB) .* parameterScale.';
fullCRB = (fullCRB + fullCRB.') / 2;
fullCRBResidual = norm(J * fullCRB - eye(5), 'fro') / sqrt(5);
% A residual below 1e-4 is consistent with the permitted 1e12 condition
% limit in IEEE double precision; it is recorded for every call.
if ~isfinite(scaledInverseResidual) || scaledInverseResidual > 1e-4
    error('estimate_crb:IllConditionedFIM', 'Equilibrated full-CRB inverse residual is too large (%g).', scaledInverseResidual);
end

result = struct('D', D, 'mu', mu, 'antennaPositions', antennaPositions, ...
    'distance', distance, 'fim', J, 'fimEigenvalues', eigJ, ...
    'fullCRB', fullCRB, 'fullCRBInverseResidual', fullCRBResidual, ...
    'scaledFullCRBInverseResidual', scaledInverseResidual, ...
    'fimConditionNumber', conditionNumber, 'nuisanceConditionNumber', conditionNuisance, ...
    'scaledFIMConditionNumber', scaledConditionNumber, ...
    'effectivePositionFIM', effectivePositionFIM, ...
    'effectivePositionFIMConditionNumber', conditionPosition, ...
    'positionCRB', positionCRB, 'positionRMSE', sqrt(trace(positionCRB)));
end
