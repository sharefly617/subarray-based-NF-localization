function position = spherical_to_cartesian(r, omega, phi)
    %SPHERICAL_TO_CARTESIAN Paper convention: omega azimuth, phi polar angle.
    arguments
    r (:,1) double {mustBeNonnegative}
    omega (:,1) double {mustBeFinite}
    phi (:,1) double {mustBeFinite}
end
if ~isequal(size(r), size(omega), size(phi))
    error('spherical_to_cartesian:SizeMismatch', 'r, omega, and phi must have equal sizes.');
end
position = [r .* cos(omega) .* sin(phi), r .* sin(omega) .* sin(phi), r .* cos(phi)];
end
