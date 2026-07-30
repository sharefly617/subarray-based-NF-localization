function [F, alphaHat, gradient, position] = eaple_profile_ml(params, y, lambda, sigma2, arrayConfig)
    %EAPLE_PROFILE_ML Exact-model profile ML objective and analytic polar gradient.
    % params=[r;omega;phi] has r in metres and angles in radians. y is (Nx*Ny)x1
    % in t=(i-1)*Ny+j order. This implements NF_Loc Eqs. (36)--(39) only with
    % exact negative-phase spherical steering; no SFM quantity is used.
    arguments
    params (3,1) double {mustBeFinite}
    y (:,1) double
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBePositive}
    arrayConfig struct
end
if params(1)<=0, error('eaple_profile_ml:InvalidRange','r must be positive.'); end
for field={'Nx','Ny','d'}
    if ~isfield(arrayConfig,field{1}), error('eaple_profile_ml:InvalidArrayConfig','arrayConfig.%s is required.',field{1}); end
end
if numel(y)~=arrayConfig.Nx*arrayConfig.Ny, error('eaple_profile_ml:LengthMismatch','y must contain Nx*Ny elements.'); end
r=params(1); omega=mod(params(2),2*pi); phi=min(pi/2,max(0,params(3)));
direction=[cos(omega)*sin(phi);sin(omega)*sin(phi);cos(phi)]; position=r*direction;
ants=upa_coordinates(arrayConfig.Nx,arrayConfig.Ny,arrayConfig.d); [a,distance]=nearfield_steering(ants,position,lambda);
N=numel(a); c=a'*y; alphaHat=c/N; residual=y-alphaHat*a; F=-real(residual'*residual)/sigma2;
dp=[direction, r*[-sin(omega)*sin(phi);cos(omega)*sin(phi);0], r*[cos(omega)*cos(phi);sin(omega)*cos(phi);-sin(phi)]];
gradient=zeros(3,1); k=2*pi/lambda;
for q=1:3
    dd=((position.'-ants)*dp(:,q))./distance;
    da=(-1j*k*dd).*a;
    dc=da'*y;
    gradient(q)=2*real(conj(c)*dc)/(sigma2*N);
end
if ~isfinite(F)||any(~isfinite(gradient)), error('eaple_profile_ml:NonFinite','Profile objective or gradient is non-finite.'); end
end
