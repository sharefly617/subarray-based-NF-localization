function steering = subarray_plane_wave_steering(Nx, Ny, dx, dy, lambda, thetaX, thetaY)
    %SUBARRAY_PLANE_WAVE_STEERING SFM relative steering for one subarray.
    % Output is (Nx*Ny)x1, t=(k-1)*Ny+l. Spacings/wavelength are metres and
    % theta values are dimensionless Eq. (9) direction cosines (toward UE).
    % Its phase is POSITIVE, exp(+j*2*pi*(...)/lambda), as in (12), unlike
    % the exact negative-phase spherical steering.
    arguments
    Nx (1,1) double {mustBeInteger, mustBePositive}
    Ny (1,1) double {mustBeInteger, mustBePositive}
    dx (1,1) double {mustBePositive}
    dy (1,1) double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    thetaX (1,1) double {mustBeFinite}
    thetaY (1,1) double {mustBeFinite}
end
if thetaX^2 + thetaY^2 > 1 + 100*eps
    error('subarray_plane_wave_steering:InvalidDirection', 'thetaX^2 + thetaY^2 must not exceed one.');
end
k = (1:Nx).' - (Nx + 1)/2; l = (1:Ny).' - (Ny + 1)/2;
ax = exp(1j * 2*pi * k * dx * thetaX / lambda);
ay = exp(1j * 2*pi * l * dy * thetaY / lambda);
steering = kron(ax, ay);
end
