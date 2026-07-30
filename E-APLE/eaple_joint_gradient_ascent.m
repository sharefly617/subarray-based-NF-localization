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
    options.maxIterations (1,1) double {mustBeInteger,mustBePositive} = 50
    options.tolerance (1,1) double {mustBePositive} = 1e-5
    options.initialStep (1,1) double {mustBePositive} = 0.5
end
p = local_initial(initial);
[F,~,g] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
history = F; steps = []; updates = [];
gradientHistory = g; parameterHistory = p; directionHistory = [];
reason = 'max_iterations';
for it=1:options.maxIterations
    [candidate,accepted,step,direction] = local_backtrack( ...
        p,g,F,y,lambda,sigma2,arrayConfig,options.initialStep,[1 2 3]);
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
    directionHistory(:,end+1) = direction;
    if norm(update) < options.tolerance
        reason = 'early_stop_update_norm';
        break
    end
end
[F,alphaHat,~,position] = eaple_profile_ml(p,y,lambda,sigma2,arrayConfig);
result = struct('params',p,'position',position,'alphaHat',alphaHat, ...
    'objective',F,'objectiveHistory',history,'stepHistory',steps, ...
    'updateHistory',updates,'gradientHistory',gradientHistory, ...
    'directionHistory',directionHistory, ...
    'parameterHistory',parameterHistory,'stopReason',reason, ...
    'method','joint_gradient_ascent_eq39_not_genetic', ...
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
        p,g,F,y,lambda,sigma2,cfg,step,indices)
    candidate = p; ok = false; acceptedStep = 0;
    direction = zeros(3,1);
    direction(indices) = g(indices);
    directionNorm = norm(direction(indices));
    if directionNorm > 0
        direction = direction/directionNorm;
    end
    if norm(direction(indices)) < 1e-12, return; end
    while step>=1e-10
        trial=p;
        trial(indices) = p(indices) + step*direction(indices);
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
