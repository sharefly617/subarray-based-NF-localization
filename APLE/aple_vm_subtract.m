function [mu, kappa, reason] = aple_vm_subtract(muPosterior, kappaPosterior, muIncoming, kappaIncoming, minKappa)
    %APLE_VM_SUBTRACT Stable VM natural-parameter subtraction for Eq. (23).
    % All mu values are radians for the VM variable pi*theta; kappa values are
    % nonnegative dimensionless concentrations. Output is finite and uses a
    % uniform VM (kappa=0) when cancellation makes its direction undefined.
    arguments
    muPosterior (1,1) double {mustBeFinite}
    kappaPosterior (1,1) double {mustBeNonnegative,mustBeFinite}
    muIncoming (1,1) double {mustBeFinite}
    kappaIncoming (1,1) double {mustBeNonnegative,mustBeFinite}
    minKappa (1,1) double {mustBePositive} = 1e-8
end
natural = kappaPosterior*exp(1j*muPosterior) - kappaIncoming*exp(1j*muIncoming);
kappa = abs(natural);
if ~isfinite(kappa) || kappa < minKappa
    mu = 0; kappa = 0; reason = 'vm_cancellation_or_uniform';
else
    mu = angle(natural); reason = '';
end
end
