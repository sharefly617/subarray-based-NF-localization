function result = estimate_aple(y, lambda, sigma2, arrayConfig, partition, options)
    %ESTIMATE_APLE Sec. IV APLE with a self-contained non-original-VALSE surrogate.
    % y is (Nx*Ny)x1 in t=(i-1)*Ny+j order. lambda,d and positions are metres;
    % sigma2 is noise variance. The output uses NF_Loc's positive Eq. (9) theta.
    % This function does not use truth, E-APLE, or the exact near-field model.
    arguments
    y (:,1) double
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBePositive}
    arrayConfig struct
    partition struct
    options.iterations (1,1) double {mustBeInteger,mustBePositive} = 1
    options.gridSize (1,1) double {mustBeInteger,mustBePositive} = 181
    options.alternations (1,1) double {mustBeInteger,mustBePositive} = 3
    options.maxGradientSteps (1,1) double {mustBeInteger,mustBePositive} = 25
    options.initialStep (1,1) double {mustBePositive} = 0.2
    options.hessianStep (1,1) double {mustBePositive} = 1e-3
    options.maxCovCondition (1,1) double {mustBePositive} = 1e10
    options.minKappa (1,1) double {mustBePositive} = 1e-5
    options.useVALSE (1,1) logical = true
end
if numel(y) ~= arrayConfig.Nx*arrayConfig.Ny, error('estimate_aple:LengthMismatch','y must have Nx*Ny elements.'); end
if partition.Nx~=arrayConfig.Nx || partition.Ny~=arrayConfig.Ny, error('estimate_aple:PartitionMismatch','partition dimensions must match arrayConfig.'); end
ants=upa_coordinates(arrayConfig.Nx,arrayConfig.Ny,arrayConfig.d); geo=subarray_geometry(ants,partition,[0;0;1]); centers=geo.centers;
M=numel(partition.subarrays); posteriorMu=zeros(M,2); posteriorKappa=zeros(M,2); gains=zeros(M,1); aoaResidual=zeros(M,1); aoaMethod=cell(M,1); aoaFallbackReasons=cell(M,2); aoaRawKappa=zeros(M,2); aoaPeakQuality=cell(M,2); aoaValseDiagnostics=cell(M,1);
for m=1:M
    sub=partition.subarrays(m); ym=y(sub.index);
    aoa=aple_plane_wave_surrogate(ym,sub.Nx,sub.Ny,arrayConfig.d,lambda,'gridSize',options.gridSize,'alternations',options.alternations,'minKappa',options.minKappa,'useVALSE',options.useVALSE);
    posteriorMu(m,:)=aoa.mu; posteriorKappa(m,:)=aoa.kappa; gains(m)=aoa.gain; aoaResidual(m)=aoa.relativeResidual; aoaMethod{m}=aoa.method;
    aoaFallbackReasons(m,:)=aoa.fallbackReason;
    aoaRawKappa(m,:)=aoa.rawKappa;
    aoaPeakQuality(m,:)=num2cell(aoa.peakQuality);
    if isfield(aoa,'valseDiagnostic'), aoaValseDiagnostics{m}=aoa.valseDiagnostic; end
end
incomingMu=zeros(M,2); incomingKappa=zeros(M,2); outgoingMu=zeros(M,2); outgoingKappa=zeros(M,2);
history=cell(options.iterations,1); pairDiagnostics=cell(options.iterations,1); fallbackReasons=cell(M,2);
for outer=1:options.iterations
    for m=1:M
        for u=1:2
            [outgoingMu(m,u),outgoingKappa(m,u),reason]=aple_vm_subtract(posteriorMu(m,u),posteriorKappa(m,u),incomingMu(m,u),incomingKappa(m,u),options.minKappa);
            if ~isempty(reason), fallbackReasons{m,u}=reason; end
        end
    end
    [seed,seedDiagnostic]=local_geometric_seed(centers,outgoingMu/pi,outgoingKappa);
    pd=repmat(struct('target',[0 0],'objectiveHistory',[],'stepHistory',[],'mean',[],'covariance',[],'covarianceCondition',NaN,'reason',''),M,2);
    newMu=zeros(M,2); newKappa=zeros(M,2);
    for m=1:M
        for u=1:2
            exclude=[m u]; f=@(p)local_objective(p,centers,outgoingMu,outgoingKappa,exclude);
            [mode,oh,sh,reason]=local_gradient_ascent(f,seed,options);
            [C,condC,cReason]=local_laplace(f,mode,options);
            if ~isempty(cReason), reason=strjoin(nonempty({reason,cReason}),';'); end
            [newMu(m,u),newKappa(m,u),vReason]=local_gaussian_to_vm(mode,C,centers(m,:),u,options.minKappa);
            if ~isempty(vReason), reason=strjoin(nonempty({reason,vReason}),';'); end
            pd(m,u)=struct('target',[m u],'objectiveHistory',oh,'stepHistory',sh,'mean',mode,'covariance',C,'covarianceCondition',condC,'reason',reason);
        end
    end
    incomingMu=newMu; incomingKappa=newKappa; history{outer}=struct('outgoingMu',outgoingMu,'outgoingKappa',outgoingKappa,'incomingMu',incomingMu,'incomingKappa',incomingKappa,'geometricSeed',seedDiagnostic); pairDiagnostics{outer}=pd;
end
fFinal=@(p)local_objective(p,centers,outgoingMu,outgoingKappa,[]);
[finalSeed,finalSeedDiagnostic]=local_geometric_seed(centers,outgoingMu/pi,outgoingKappa);
[position,finalHistory,finalSteps,finalReason]=local_gradient_ascent(fFinal,finalSeed,options);
[positionCovariance,positionCovarianceCondition,finalCovReason]=local_laplace(fFinal,position,options);
if ~isempty(finalCovReason), finalReason=strjoin(nonempty({finalReason,finalCovReason}),';'); end
result=struct('position',position,'positionCovariance',positionCovariance,'positionCovarianceCondition',positionCovarianceCondition, ...
    'isOriginalVALSE',false,'aoaMethod','alternating_1d_matched_periodogram','posteriorMu',posteriorMu,'posteriorKappa',posteriorKappa, ...
    'aoaFallbackReasons',{aoaFallbackReasons},'aoaRawKappa',aoaRawKappa,'aoaPeakQuality',{aoaPeakQuality},'aoaValseDiagnostics',{aoaValseDiagnostics}, ...
    'complexGain',gains,'aoaRelativeResidual',aoaResidual,'outgoingMu',outgoingMu,'outgoingKappa',outgoingKappa, ...
    'geometricSeedDiagnostic',finalSeedDiagnostic, ...
    'incomingMu',incomingMu,'incomingKappa',incomingKappa,'iterationHistory',{history},'pairDiagnostics',{pairDiagnostics}, ...
    'fallbackReasons',{fallbackReasons},'finalObjectiveHistory',finalHistory,'finalStepHistory',finalSteps,'finalReason',finalReason,'sigma2',sigma2);
end

function [seed,diagnostic]=local_geometric_seed(centers,theta,kappa)
    % Select a finite geometric-message mode without using a true range bound.
    % The old weighted line-intersection seed becomes ill-conditioned for grazing
    % directions.  We retain it as a scale only, then maximize the same Eq. (25c)
    % message objective along the observed mean direction over an adaptive log
    % range.  The powers-of-two span is centered on the data-derived ray scale and
    % is not tied to the Fig. 6 [30,33] interval.
    weights=max(sum(kappa,2),1e-6); weights=weights/sum(weights);
    meanTheta=sum(theta.*weights,1);
    radialNorm=norm(meanTheta);
    if ~isfinite(radialNorm) || radialNorm<1e-12
        meanTheta=[0 0]; radialNorm=0;
    end
    if radialNorm>=1, meanTheta=meanTheta/max(radialNorm,1+1e-6)*(.999999); end
    tz=sqrt(max(1-sum(meanTheta.^2),1e-12)); direction=[meanTheta(:);tz]; direction=direction/max(norm(direction),eps);
    rayCandidate=local_ray_intersection(centers,theta,kappa);
    aperture=max(max(centers,[],1)-min(centers,[],1));
    if all(isfinite(rayCandidate)) && rayCandidate(3)>0
        rayScale=max([norm(rayCandidate); aperture; 1e-3]);
    else
        rayScale=max([norm(mean(centers,1)); aperture; 1e-3]);
    end
    qGrid=log(rayScale)+(-16:16)*log(2); radii=exp(qGrid);
    values=zeros(size(radii));
    for i=1:numel(radii), values(i)=local_objective(radii(i)*direction,centers,mu_from_theta(theta),kappa,[]); end
    [bestValue,bestIndex]=max(values); qBest=qGrid(bestIndex);
    if numel(qGrid)>1
        lo=qGrid(max(1,bestIndex-1)); hi=qGrid(min(numel(qGrid),bestIndex+1));
        if hi>lo
            qRefined=fminbnd(@(q)-local_objective(exp(q)*direction,centers,mu_from_theta(theta),kappa,[]),lo,hi,optimset('Display','off'));
            refinedValue=local_objective(exp(qRefined)*direction,centers,mu_from_theta(theta),kappa,[]);
            if isfinite(refinedValue) && refinedValue>=bestValue, qBest=qRefined; bestValue=refinedValue; end
        end
    end
    seed=exp(qBest)*direction;
    rayObjective=-Inf;
    if all(isfinite(rayCandidate)) && rayCandidate(3)>0
        rayObjective=local_objective(rayCandidate,centers,mu_from_theta(theta),kappa,[]);
        if isfinite(rayObjective) && rayObjective>bestValue
            seed=rayCandidate; bestValue=rayObjective;
        end
    end
    if ~all(isfinite(seed)) || seed(3)<=0, seed=rayScale*direction; bestValue=local_objective(seed,centers,mu_from_theta(theta),kappa,[]); end
    diagnostic=struct('method','adaptive_log_range_profile_seed','direction',direction,'rayScale',rayScale, ...
        'logRangeGrid',qGrid,'rangeGrid',radii,'objectiveGrid',values,'selectedRange',norm(seed), ...
        'selectedObjective',bestValue,'rayCandidate',rayCandidate,'rayObjective',rayObjective, ...
        'thetaInput',theta,'kappaInput',kappa);
end
function mu=mu_from_theta(theta), mu=pi*theta; end
function seed=local_ray_intersection(centers,theta,kappa)
    M=size(centers,1); A=zeros(3); b=zeros(3,1);
    for m=1:M
        tx=max(-.999,min(.999,theta(m,1))); ty=max(-.999,min(.999,theta(m,2))); tz=sqrt(max(1-tx^2-ty^2,1e-6)); u=[tx;ty;tz]; w=max(sum(kappa(m,:)),1e-3); P=eye(3)-u*u.'; A=A+w*P; b=b+w*P*centers(m,:).';
    end
    if rcond(A)<1e-10, seed=[NaN;NaN;NaN]; else, seed=A\b; end
end
function value=local_objective(p,centers,mu,kappa,exclude)
    value=0; for m=1:size(centers,1)
    d=p.'-centers(m,:); r=norm(d); if r<1e-8, value=-Inf; return; end
    for u=1:2
        if ~isempty(exclude)&&m==exclude(1)&&u==exclude(2), continue; end
        theta=d(u)/r; value=value+kappa(m,u)*cos(pi*theta-mu(m,u));
    end
end
end
function [p,history,steps,reason]=local_gradient_ascent(f,p,op)
    history=f(p); steps=[]; reason='';
    for it=1:op.maxGradientSteps
        g=local_gradient(f,p,op.hessianStep); if ~all(isfinite(g))||norm(g)<1e-9, break; end
        step=op.initialStep; accepted=false; current=history(end);
        while step>=1e-8
            candidate=p+step*g/max(norm(g),1); candidate(3)=max(candidate(3),1e-6); v=f(candidate);
            if isfinite(v)&&v>=current, p=candidate; history(end+1)=v; steps(end+1)=step; accepted=true; break; end
            step=step/2;
        end
        if ~accepted, reason='line_search_failed'; break; end
    end
end
function g=local_gradient(f,p,h)
    g=zeros(3,1); for i=1:3, e=zeros(3,1); e(i)=h; g(i)=(f(p+e)-f(p-e))/(2*h); end
end
function [C,c,reason]=local_laplace(f,p,op)
    h=op.hessianStep; H=zeros(3); f0=f(p);
    for i=1:3
        ei=zeros(3,1);ei(i)=h; H(i,i)=(f(p+ei)-2*f0+f(p-ei))/h^2;
        for j=i+1:3
            ej=zeros(3,1);ej(j)=h; H(i,j)=(f(p+ei+ej)-f(p+ei-ej)-f(p-ei+ej)+f(p-ei-ej))/(4*h^2);H(j,i)=H(i,j);
        end
    end
    Q=(-H+(-H).')/2; ev=eig(Q); reason='';
    if any(~isfinite(ev))||min(ev)<=1e-9, C=eye(3)*1e6; c=Inf; reason='hessian_not_positive_definite'; return; end
    c=cond(Q); if ~isfinite(c)||c>op.maxCovCondition, C=eye(3)*1e6; reason='hessian_ill_conditioned'; else, C=Q\eye(3); C=(C+C.')/2; end
end
function [mu,kappa,reason]=local_gaussian_to_vm(meanP,C,center,u,minK)
    d=meanP.'-center; nr=norm(d); reason=''; if nr<1e-8, mu=0;kappa=0;reason='center_coincidence';return;end
    theta=d(u)/nr; if abs(theta)>=1-1e-6, mu=pi*max(-1,min(1,theta));kappa=0;reason='theta_near_boundary';return;end
    e=zeros(3,1);e(u)=1; crossv=cross(d.',e.'); v=cross(crossv,d.'); nv=norm(v); if nv<1e-10, mu=pi*theta;kappa=0;reason='cross_product_degenerate';return;end
    v=v(:)/nv; denom=pi^2*(1-theta^2)*(v.'*C*v); if ~isfinite(denom)||denom<=0, mu=pi*theta;kappa=0;reason='invalid_vm_variance';else,mu=pi*theta;kappa=max(minK,min(1e6,nr^2/denom));end
end
function x=nonempty(c), x=c(~cellfun(@isempty,c)); end
