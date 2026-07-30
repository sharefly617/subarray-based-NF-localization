function result = eaple_joint_gradient_ascent(y,lambda,sigma2,arrayConfig,initial,options)
    %EAPLE_JOINT_GRADIENT_ASCENT Fig. 6 GA Ini: APLE control (GA means gradient ascent).
    % Uses exact profile ML and jointly updates [r,omega,phi]; it is not a genetic
    % algorithm. initial may contain position (m) or r,omega,phi. Output histories
    % record only accepted nondecreasing profile-objective updates.
    arguments
    y (:,1) double
    lambda (1,1) double {mustBePositive}
    sigma2 (1,1) double {mustBePositive}
    arrayConfig struct
    initial struct
    options.maxIterations (1,1) double {mustBeInteger,mustBePositive} = 500
    options.tolerance (1,1) double {mustBePositive} = 1e-5
    options.initialStep (1,1) double {mustBePositive} = 0.5
end
p = local_initial(initial);
[F,~,g] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
% The raw polar derivatives have very different curvature scales (especially
% omega/phi versus range).  Estimate the diagonal curvature once in
% q=[log(r),omega,phi] coordinates and use it as a deterministic local
% preconditioner.  This remains gradient ascent with backtracking; it is not a
% Newton or genetic-algorithm update.  Keeping the scale fixed also makes the
% trajectory auditable and avoids silently changing the objective.
[curvatureScale,curvatureDiagnostic] = local_initial_curvature(p,F,y,lambda,sigma2,arrayConfig);
history = F; steps = []; updates = [];
gradientHistory = g; parameterHistory = p; scaledGradientHistory = [];
reason = 'max_iterations';
for it=1:options.maxIterations
    [candidate,accepted,step,scaledDirection] = local_backtrack( ...
        p,g,F,y,lambda,sigma2,arrayConfig,options.initialStep,[1 2 3],curvatureScale);
    if ~accepted
        reason = 'line_search_failed_or_stationary';
        break
    end
    update = candidate-p;
    p = candidate;
    [F,~,g] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
    history(end+1) = F; steps(end+1) = step;
    updates(:,end+1) = update; gradientHistory(:,end+1) = g;
    parameterHistory(:,end+1) = p;
    scaledGradientHistory(:,end+1) = scaledDirection;
    if norm(update) < options.tolerance
        reason = 'early_stop_update_norm';
        break
    end
end
[F,alphaHat,~,position] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
result = struct('params',p,'position',position,'alphaHat',alphaHat, ...
    'objective',F,'objectiveHistory',history,'stepHistory',steps, ...
    'updateHistory',updates,'gradientHistory',gradientHistory, ...
    'scaledGradientHistory',scaledGradientHistory, ...
    'parameterHistory',parameterHistory,'curvatureScale',curvatureScale, ...
    'curvatureDiagnostic',curvatureDiagnostic,'stopReason',reason, ...
    'method','joint_curvature_preconditioned_gradient_ascent_not_genetic', ...
    'maxIterations',options.maxIterations,'iterationsPerformed',numel(steps));
end

function p=local_initial(initial)
    if isfield(initial,'position')
        p = eaple_cartesian_to_polar(initial.position);
    elseif all(isfield(initial,{'r','omega','phi'}))
        p = [initial.r;initial.omega;initial.phi];
    else
        error('eaple_joint_gradient_ascent:InvalidInitial', ...
            'initial requires position or r, omega, phi.');
    end
    p = local_project(p);
end
function [candidate,ok,acceptedStep,direction] = local_backtrack( ...
        p,g,F,y,lambda,sigma2,cfg,step,indices,curvatureScale)
    candidate = p; ok = false; acceptedStep = 0;
    direction = local_preconditioned_direction(p,g,indices,curvatureScale);
    if norm(direction(indices)) < 1e-12, return; end
    while step>=1e-10
        trial=p;
        if any(indices==1), trial(1)=p(1)*exp(step*direction(1)); end
        angularIndices = intersect(indices,[2 3]);
        trial(angularIndices) = trial(angularIndices) + ...
            step*direction(angularIndices);
        trial = local_project(trial);
        Ft = eaple_profile_ml(trial,y,lambda,sigma2,cfg);
        if isfinite(Ft) && Ft >= F-1e-12
            candidate = trial; ok = true; acceptedStep = step; return
        end
        step = step/2;
    end
end
function p = local_project(p)
    p(1) = max(p(1),1e-8);
    p(2) = mod(p(2),2*pi);
    p(3) = min(pi/2,max(0,p(3)));
end
function direction=local_preconditioned_direction(p,g,indices,curvatureScale)
    direction=zeros(3,1);
    qGradient=[p(1)*g(1);g(2);g(3)];
    direction=qGradient./curvatureScale;
    direction(setdiff(1:3,indices))=0;
    normDirection=norm(direction(indices));
    if normDirection>1,direction=direction/normDirection;end
end

function [scale,diagnostic]=local_initial_curvature(p,F,y,lambda,sigma2,cfg)
    % Estimate diagonal curvature of F in q=[log(r),omega,phi] at initialization.
    q=[log(p(1));p(2);p(3)]; h=[1e-3;1e-4;1e-4]; diagonal=nan(3,1); values=nan(3,2);
    for k=1:3
        qp=q; qm=q; qp(k)=qp(k)+h(k); qm(k)=qm(k)-h(k);
        pp=[exp(qp(1));qp(2);qp(3)]; pm=[exp(qm(1));qm(2);qm(3)];
        Fp=eaple_profile_ml(pp,y,lambda,sigma2,cfg); Fm=eaple_profile_ml(pm,y,lambda,sigma2,cfg);
        diagonal(k)=(Fp-2*F+Fm)/(h(k)^2); values(k,:)=[Fp,Fm];
    end
    if any(~isfinite(diagonal)),scale=ones(3,1);else,scale=max(abs(diagonal),1);end
    diagnostic=struct('coordinateNames',{{'logRange','omega','phi'}},'finiteDifferenceStep',h,'diagonalCurvature',diagonal,'curvatureScale',scale,'offsetObjectives',values);
end
