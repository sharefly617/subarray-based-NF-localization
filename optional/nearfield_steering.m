function [steering, distance] = nearfield_steering(antennaPositions, position, lambda)
    %NEARFIELD_STEERING Exact spherical-wave steering a_t=exp(-j*2*pi*r_t/lambda).
    arguments
    antennaPositions (:,3) double {mustBeFinite}
    position (3,1) double {mustBeFinite}
    lambda (1,1) double {mustBePositive}
end
distance = sqrt(sum((antennaPositions - position.').^2, 2));
if any(distance <= eps(max(1, norm(position))))
    error('nearfield_steering:CoincidentAntenna', 'UE location coincides with or is too close to an antenna.');
end
steering = exp(-1j * 2 * pi * distance / lambda);
end
