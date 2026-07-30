function result = aple_plane_wave_surrogate(y, Nx, Ny, d, lambda, options)
    %APLE_PLANE_WAVE_SURROGATE Self-contained 2-D SFM AoA substitute, not VALSE.
    % y is (Nx*Ny)x1 in t=(k-1)*Ny+l order. It alternates 1-D matched
    % periodograms for positive-phase Eq. (12) steering. A coarse-grid global peak
    % is refined by bounded continuous 1-D optimization; this is not original VALSE.
    % It returns VM parameters for pi*theta, a complex LS gain, and diagnostics.
    arguments
    y (:,1) double
    Nx (1,1) double {mustBeInteger,mustBePositive}
    Ny (1,1) double {mustBeInteger,mustBePositive}
    d (1,1) double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    options.gridSize (1,1) double {mustBeInteger,mustBePositive} = 181
    options.alternations (1,1) double {mustBeInteger,mustBePositive} = 3
    options.minKappa (1,1) double {mustBePositive} = 1e-4
    options.maxKappa (1,1) double {mustBePositive} = 1e8
    options.refinementTolerance (1,1) double {mustBePositive} = 1e-10
    options.useVALSE (1,1) logical = true
end
if numel(y) ~= Nx*Ny, error('aple_plane_wave_surrogate:LengthMismatch','y must have Nx*Ny elements.'); end
if options.gridSize < 11, error('aple_plane_wave_surrogate:GridTooSmall','gridSize must be at least 11.'); end
Y = reshape(y, Ny, Nx).'; % rows k, columns l under t=(k-1)*Ny+l
grid = linspace(-1,1,options.gridSize); coarseStep=grid(2)-grid(1);
% A zero-valued alternating start can cancel the orthogonal dimension when
% thetaY (or thetaX) is nonzero.  A rank-one SVD start uses only y and the
% known SFM separability, and avoids that deterministic low-SNR failure mode.
[thetaX,thetaY] = local_svd_initial(Y,Nx,Ny,d,lambda,grid);
valseDiagnostic = struct('available',false,'used',false,'fallbackReason','', ...
    'x',[],'y',[],'freqs',{{}},'amps',{{}},'noiseVar',[],'iterations',[]);
if options.useVALSE && exist('VALSE','file')==2
    [thetaXv,dx] = local_valse_axis(U_or_col(Y,Nx,Ny,1),Nx,d,lambda);
    [thetaYv,dy] = local_valse_axis(U_or_col(Y,Nx,Ny,2),Ny,d,lambda);
    valseDiagnostic.available = true;
    valseDiagnostic.x = dx; valseDiagnostic.y = dy;
    if dx.used, thetaX=thetaXv; end
    if dy.used, thetaY=thetaYv; end
else
    valseDiagnostic.fallbackReason='VALSE.m not found or disabled';
end
for it=1:options.alternations
    ay = local_steer(Ny,d,lambda,thetaY); zx = Y*conj(ay); coarseX = local_peak(zx,Nx,d,lambda,grid); thetaX = local_refine(zx,Nx,d,lambda,coarseX,coarseStep,options.refinementTolerance);
    ax = local_steer(Nx,d,lambda,thetaX); zy = Y.'*conj(ax); coarseY = local_peak(zy,Ny,d,lambda,grid); thetaY = local_refine(zy,Ny,d,lambda,coarseY,coarseStep,options.refinementTolerance);
end
steering = kron(local_steer(Nx,d,lambda,thetaX),local_steer(Ny,d,lambda,thetaY));
gain = (steering'*y)/(steering'*steering);
residual = norm(y-gain*steering)^2 / max(norm(y)^2,eps);
% Keep the previous curvature/SNR-proxy value as a diagnostic, but form the
% returned VM confidence from a normalized matched-spectrum posterior.
zx = Y*conj(local_steer(Ny,d,lambda,thetaY)); zy = Y.'*conj(local_steer(Nx,d,lambda,thetaX));
[~, rawKappaX, curvatureX] = local_kappa(zx,Nx,d,lambda,thetaX,residual,options);
[~, rawKappaY, curvatureY] = local_kappa(zy,Ny,d,lambda,thetaY,residual,options);
confidenceX = local_posterior_confidence(zx,Nx,d,lambda,thetaX,grid,residual,options);
confidenceY = local_posterior_confidence(zy,Ny,d,lambda,thetaY,grid,residual,options);
result = struct('theta',[thetaX thetaY], 'mu',pi*[thetaX thetaY], 'kappa',[confidenceX.kappa confidenceY.kappa], ...
    'gain',gain,'relativeResidual',residual,'isOriginalVALSE',false,'method','svd_initialized_alternating_1d_matched_periodogram_with_continuous_refinement', ...
    'initialization',ternary(valseDiagnostic.available,'rank1_svd_plus_VALSE','rank1_svd'), ...
    'valseDiagnostic',valseDiagnostic, ...
    'coarseTheta',[coarseX coarseY],'coarseGridStep',coarseStep,'refinementTolerance',options.refinementTolerance, ...
    'rawKappa',[rawKappaX rawKappaY],'logScoreCurvature',[curvatureX curvatureY], ...
    'peakQuality',[confidenceX confidenceY], ...
    'fallbackReason',{ {confidenceX.fallbackReason confidenceY.fallbackReason} });
end

function z=U_or_col(Y,Nx,Ny,which)
    [U,~,V]=svd(Y,'econ');
    if which==1, z=U(:,1); else, z=conj(V(:,1)); end
end

function [theta,diag]=local_valse_axis(z,N,d,lambda)
    diag=struct('used',false,'fallbackReason','','freqs',[],'amps',[],'noiseVar',NaN,'iterations',NaN,'selectedIndex',NaN);
    theta=0;
    try
        out=VALSE(z,(0:N-1).',2,z);
        freqs=out.freqs(:); amps=out.amps(:);
        diag.freqs=freqs; diag.amps=amps;
        if isfield(out,'noise_var'),diag.noiseVar=out.noise_var;end
        if isfield(out,'iterations'),diag.iterations=out.iterations;end
        if isempty(freqs) || isempty(amps) || any(~isfinite(freqs)) || any(~isfinite(amps))
            diag.fallbackReason='VALSE returned no finite component'; return;
        end
        [~,ix]=max(abs(amps));
        % VALSE frequencies are radians/sample.  Eq. (12) has phase
        % 2*pi*d*theta/lambda, hence theta=freq*lambda/(2*pi*d).
        theta=wrap_to_interval(freqs(ix)*lambda/(2*pi*d),-1,1);
        diag.selectedIndex=ix; diag.used=isfinite(theta);
        if ~diag.used,diag.fallbackReason='nonfinite selected theta';theta=0;end
    catch ME
        diag.fallbackReason=[ME.identifier ': ' ME.message];
    end
end

function y=wrap_to_interval(x,lo,hi)
    period=2*(hi-lo); y=mod(x-lo,period)+lo; if y>hi,y=period-y+2*lo;end
end

function y=ternary(c,a,b),if c,y=a;else,y=b;end,end

function a = local_steer(N,d,lambda,theta)
    k=(1:N).'-(N+1)/2; a=exp(1j*2*pi*k*d*theta/lambda);
end
function theta = local_peak(z,N,d,lambda,grid)
    scores=arrayfun(@(q)local_score(z,N,d,lambda,q),grid); [~,ix]=max(scores); theta=grid(ix);
end
function [thetaX,thetaY]=local_svd_initial(Y,Nx,Ny,d,lambda,grid)
    [U,~,V]=svd(Y,'econ');
    if isempty(U) || isempty(V) || ~all(isfinite(U(:,1))) || ~all(isfinite(V(:,1)))
        thetaX=0; thetaY=0; return;
    end
    thetaX=local_peak(U(:,1),Nx,d,lambda,grid);
    thetaY=local_peak(conj(V(:,1)),Ny,d,lambda,grid);
end
function theta = local_refine(z,N,d,lambda,coarseTheta,coarseStep,tolerance)
    lower=max(-1,coarseTheta-coarseStep); upper=min(1,coarseTheta+coarseStep);
    if upper<=lower,theta=coarseTheta;return;end
    theta=fminbnd(@(q)-local_score(z,N,d,lambda,q),lower,upper,optimset('TolX',tolerance,'Display','off'));
end
function score = local_score(z,N,d,lambda,theta)
    a=local_steer(N,d,lambda,theta); score=abs(a'*z)^2/max(a'*a,eps);
end
function [kappa,rawKappa,curvature]=local_kappa(z,N,d,lambda,theta,residual,options)
    h=min(1e-4,max(1e-7,options.refinementTolerance*100));
    if theta-h<-1||theta+h>1,h=min(theta+1,1-theta)/2;end
    if h<=eps,curvature=0;rawKappa=0;kappa=options.minKappa;return;end
    f0=log(max(local_score(z,N,d,lambda,theta),eps));fp=log(max(local_score(z,N,d,lambda,theta+h),eps));fm=log(max(local_score(z,N,d,lambda,theta-h),eps));
    curvature=max(0,-(fp-2*f0+fm)/h^2);snrProxy=max(1,(1-residual)/max(residual,1e-6));rawKappa=curvature*snrProxy;kappa=min(options.maxKappa,max(options.minKappa,rawKappa));
end

function confidence=local_posterior_confidence(z,N,d,lambda,theta,grid,residual,options)
    % Form a discrete posterior from the complete matched spectrum.  The spectrum
    % supplies its shape; the fitted explained-to-residual energy ratio supplies an
    % effective likelihood scale.  This avoids treating array gain (peak/floor) as
    % if it were an unlimited SNR.
    scores=arrayfun(@(q)local_score(z,N,d,lambda,q),grid);
    [mainPeak,peakIndex]=max(scores);
    noiseFloor=max(median(scores),eps);
    peakToFloor=max(mainPeak/noiseFloor,1);
    explainedToResidual=max((1-residual)/max(residual,eps),0);
    evidence=log1p(explainedToResidual);
    spectrumShape=(scores-noiseFloor)/max(mainPeak-noiseFloor,eps);
    logWeights=evidence*spectrumShape;
    weights=exp(logWeights-max(logWeights));
    probability=weights/max(sum(weights),eps);
    resultant=abs(sum(probability.*exp(1j*pi*grid)));
    posteriorKappa=local_vm_kappa(resultant);

    outside=abs(grid-grid(peakIndex))>=2/N;
    secondPeak=max(scores(outside));
    psrDb=10*log10(max(mainPeak,eps)/max(secondPeak,eps));
    boundaryLimited=abs(theta)>1-2/N;
    localWindow=abs(grid-theta)<=1/N;
    localScores=scores(localWindow); localGrid=grid(localWindow);
    [~,localIndex]=min(abs(localGrid-theta));
    tolerance=max(mainPeak,eps)*1e-9;
    localUnimodal=all(diff(localScores(1:localIndex))>=-tolerance) && ...
        all(diff(localScores(localIndex:end))<=tolerance);
    localMaxima=(localScores(2:end-1)>=localScores(1:end-2)) & ...
        (localScores(2:end-1)>=localScores(3:end));
    localPeakCount=sum(localMaxima);

    reasons={};
    if ~isfinite(posteriorKappa) || ~isfinite(mainPeak)
        reasons{end+1}='nonfinite_posterior';
    end
    if peakToFloor<2
        reasons{end+1}='insufficient_peak_to_floor_evidence';
    end
    if psrDb<3 && ~boundaryLimited
        reasons{end+1}='non_dominant_peak';
    end
    if ~localUnimodal || localPeakCount>1
        reasons{end+1}='locally_multimodal_peak';
    end
    if explainedToResidual<2
        reasons{end+1}='insufficient_explained_energy_evidence';
    end
    if isempty(reasons)
        kappa=min(options.maxKappa,max(options.minKappa,posteriorKappa));
        fallbackReason='';
    else
        % Preserve the finite posterior concentration for geometric fusion.  The
        % fallback remains explicit, but replacing a positionally accurate,
        % noise-limited peak by exactly minKappa discards all cross-subarray
        % information and can make the APLE Hessian singular.
        kappa=min(options.maxKappa,max(options.minKappa,posteriorKappa));
        fallbackReason=strjoin(reasons,';');
    end
    confidence=struct('kappa',kappa,'posteriorKappa',posteriorKappa, ...
        'mainPeak',mainPeak,'secondPeak',secondPeak,'psrDb',psrDb, ...
        'peakToFloor',peakToFloor,'localUnimodal',localUnimodal, ...
        'localPeakCount',localPeakCount,'boundaryLimited',boundaryLimited, ...
        'posteriorResultant',resultant, ...
        'effectiveLogEvidence',evidence,'explainedToResidual',explainedToResidual, ...
        'relativeResidual',residual,'fallbackReason',fallbackReason);
end

function kappa=local_vm_kappa(resultant)
    % Standard piecewise inverse of A(kappa)=I1(kappa)/I0(kappa).
    R=min(max(real(resultant),0),1-1e-12);
    if R<.53
        kappa=2*R+R^3+5*R^5/6;
    elseif R<.85
        kappa=-.4+1.39*R+.43/(1-R);
    else
        kappa=1/(R^3-4*R^2+3*R);
    end
end
