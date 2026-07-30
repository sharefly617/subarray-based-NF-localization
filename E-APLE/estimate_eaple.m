function result = estimate_eaple(y,lambda,sigma2,arrayConfig,initial,options)
    %ESTIMATE_EAPLE Exact-spherical E-APLE BCA for NF_Loc Eqs. (35)--(43).
    % initial accepts APLE output via initial.position, or GRID polar data via
    % initial.r, initial.omega, initial.phi. BCA alternates angles then range,
    % accepts only nondecreasing profile-ML objective values, and defaults T2=50.
    arguments
    y (:,1) double
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBePositive}
    arrayConfig struct
    initial struct
    options.maxIterations (1,1) double {mustBeInteger,mustBePositive} = 50
    options.tolerance (1,1) double {mustBePositive} = 1e-5
    options.angleInitialStep (1,1) double {mustBePositive} = 0.1
    options.rangeInitialStep (1,1) double {mustBePositive} = 0.1
end
p = local_initial(initial);
[F,~,g] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
history = F; angleSteps = []; rangeSteps = []; updates = [];
angleGradientHistory = []; rangeGradientHistory = [];
parameterHistory = p; reasons = {}; stopReason = 'max_iterations';
for it=1:options.maxIterations
    [pA,okA,stepA] = local_backtrack( ...
        p,g,F,y,lambda,sigma2,arrayConfig,options.angleInitialStep,[2 3]);
    if okA
        [FA,~,gA] = eaple_profile_ml(pA,y,lambda,sigma2,arrayConfig);
    else
        pA = p; FA = F; gA = g;
        reasons{end+1} = 'angle_line_search_failed_or_stationary';
    end
    [pR,okR,stepR] = local_backtrack( ...
        pA,gA,FA,y,lambda,sigma2,arrayConfig,options.rangeInitialStep,1);
    if okR
        [FR,~,gR] = eaple_profile_ml(pR,y,lambda,sigma2,arrayConfig);
    else
        pR = pA; FR = FA; gR = gA;
        reasons{end+1} = 'range_line_search_failed_or_stationary';
    end
    update = pR-p; p = pR; F = FR; g = gR;
    history(end+1) = F; angleSteps(end+1) = stepA; rangeSteps(end+1) = stepR;
    updates(:,end+1) = update; angleGradientHistory(:,end+1) = gA;
    rangeGradientHistory(:,end+1) = gR; parameterHistory(:,end+1) = p;
    if norm(update) < options.tolerance
        stopReason = 'early_stop_update_norm'; break
    end
    if ~(okA || okR)
        stopReason = 'no_accepted_coordinate_update'; break
    end
end
[F,alphaHat,~,position]=eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
risk={}; if strcmp(stopReason,'max_iterations'),risk{end+1}='local_maximum_risk_max_iterations';end
if any(p(3)<1e-8|p(3)>pi/2-1e-8),risk{end+1}='phi_boundary_risk';end
result=struct('params',p,'position',position,'alphaHat',alphaHat,'objective',F,'objectiveHistory',history, ...
    'angleStepHistory',angleSteps,'rangeStepHistory',rangeSteps,'updateHistory',updates,'angleGradientHistory',angleGradientHistory,'rangeGradientHistory',rangeGradientHistory,'parameterHistory',parameterHistory,'stopReason',stopReason, ...
    'reasons',{reasons},'localMaximumRisk',{risk},'method','bca_exact_profile_ml','initialSource',local_source(initial), ...
    'maxIterations',options.maxIterations,'iterationsPerformed',numel(angleSteps));
end

function p=local_initial(initial)
    if isfield(initial,'position'),p=eaple_cartesian_to_polar(initial.position);elseif all(isfield(initial,{'r','omega','phi'})),p=[initial.r;initial.omega;initial.phi];else,error('estimate_eaple:InvalidInitial','initial requires position or r, omega, phi.');end
    p=local_project(p);
end
function source=local_source(initial),if isfield(initial,'source'),source=initial.source;elseif isfield(initial,'position'),source='APLE_or_cartesian';else,source='GRID_or_polar';end,end
function [candidate,ok,acceptedStep]=local_backtrack(p,g,F,y,lambda,sigma2,cfg,step,indices)
    candidate=p;ok=false;acceptedStep=0;direction=local_preconditioned_direction(p,g,indices);if norm(direction(indices))<1e-12,return;end
    while step>=1e-10
        trial=p;
        if any(indices==1), trial(1)=p(1)*exp(step*direction(1)); end
        angularIndices=intersect(indices,[2 3]);trial(angularIndices)=trial(angularIndices)+step*direction(angularIndices);trial=local_project(trial);Ft=eaple_profile_ml(trial,y,lambda,sigma2,cfg);
        if isfinite(Ft)&&Ft>=F-1e-12,candidate=trial;ok=true;acceptedStep=step;return;end
        step=step/2;
    end
end
function p=local_project(p),p(1)=max(p(1),1e-8);p(2)=mod(p(2),2*pi);p(3)=min(pi/2,max(0,p(3)));end
function direction=local_preconditioned_direction(p,g,indices)
    direction=zeros(3,1);
    if any(indices==1), logRangeGradient=p(1)*g(1); direction(1)=logRangeGradient/max(1,abs(logRangeGradient)); end
    angularIndices=intersect(indices,[2 3]);if ~isempty(angularIndices),angularGradient=g(angularIndices);direction(angularIndices)=angularGradient/max(1,norm(angularGradient));end
end
