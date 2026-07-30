function result = aple_2d_valse(y, Nx, Ny, d, lambda, sigma2, options)
%APLE_2D_VALSE Paper Eq. (18)-(20) 2-D variational AoA posterior.
% The full Nx-by-Ny SFM snapshot is used directly.  The two one-dimensional
% marginal line-spectrum posteriors are obtained from the row/column sample
% covariance matrices and moment-matched to VM distributions on pi*theta.
% No SVD factorization and no call to the reference VALSE.m are used.
    arguments
        y (:,1) double
        Nx (1,1) double {mustBeInteger,mustBePositive}
        Ny (1,1) double {mustBeInteger,mustBePositive}
        d (1,1) double {mustBePositive}
        lambda (1,1) double {mustBePositive}
        sigma2 (1,1) double {mustBePositive}
        options.gridSize (1,1) double {mustBeInteger,mustBePositive} = 721
        options.refinementTolerance (1,1) double {mustBePositive} = 1e-10
        options.minKappa (1,1) double {mustBePositive} = 1e-5
        options.maxKappa (1,1) double {mustBePositive} = 1e8
        options.priorMu (1,2) double = [0 0]
        options.priorKappa (1,2) double {mustBeNonnegative} = [0 0]
    end
    if numel(y) ~= Nx*Ny
        error('aple_2d_valse:LengthMismatch','y must have Nx*Ny elements.');
    end
    if options.gridSize < 101
        error('aple_2d_valse:GridTooSmall','gridSize must be at least 101.');
    end

    % t=(k-1)*Ny+l: rows are x indices and columns are y indices.
    Y = reshape(y, Ny, Nx).';
    Rx = Y*Y';
    Ry = Y.'*conj(Y);
    thetaGrid = linspace(-1, 1, options.gridSize);

    postX = local_axis_posterior( ...
        Rx, Nx, Ny, d, lambda, sigma2, thetaGrid, ...
        options.priorMu(1), options.priorKappa(1), options);
    postY = local_axis_posterior( ...
        Ry, Ny, Nx, d, lambda, sigma2, thetaGrid, ...
        options.priorMu(2), options.priorKappa(2), options);

    theta = [postX.theta postY.theta];
    steering = kron( ...
        local_steering(Nx,d,lambda,theta(1)), ...
        local_steering(Ny,d,lambda,theta(2)));
    gain = (steering'*y)/(steering'*steering);
    residual = norm(y-gain*steering)^2/max(norm(y)^2,eps);

    result = struct( ...
        'theta',theta, ...
        'mu',[postX.mu postY.mu], ...
        'kappa',[postX.kappa postY.kappa], ...
        'gain',gain, ...
        'relativeResidual',residual, ...
        'method','paper_eq18_23_2d_variational_valse', ...
        'usesReferenceVALSEFile',false, ...
        'isOriginalVALSE',false, ...
        'posterior',{ {postX postY} }, ...
        'rawKappa',[postX.rawKappa postY.rawKappa], ...
        'peakQuality',[postX.peakQuality postY.peakQuality], ...
        'fallbackReason',{ {postX.fallbackReason postY.fallbackReason} });
end

function posterior = local_axis_posterior( ...
        R, N, otherN, d, lambda, sigma2, grid, priorMu, priorKappa, options)
    scores = zeros(size(grid));
    for k = 1:numel(grid)
        a = local_steering(N,d,lambda,grid(k));
        scores(k) = max(real(a'*R*a),0)/(N*otherN*sigma2);
    end
    logPosterior = scores + priorKappa*cos(pi*grid-priorMu);
    [peakValue,peakIndex] = max(logPosterior);
    coarseTheta = grid(peakIndex);
    step = grid(2)-grid(1);
    lo = max(-1,coarseTheta-step);
    hi = min(1,coarseTheta+step);
    objective = @(q) -local_axis_log_posterior( ...
        q,R,N,otherN,d,lambda,sigma2,priorMu,priorKappa);
    theta = fminbnd(objective,lo,hi, ...
        optimset('TolX',options.refinementTolerance,'Display','off'));

    weights = exp(logPosterior-max(logPosterior));
    weights = weights/max(sum(weights),eps);
    circularMoment = sum(weights.*exp(1j*pi*grid));
    % VALSE Heuristic 2 uses the continuous posterior mode and its local
    % curvature. For M(pi*theta;mu,kappa), the theta-domain curvature at the
    % mode is kappa*pi^2.
    h = min(1e-4,max(1e-7,step/20));
    f0 = -objective(theta);
    fp = -objective(min(1,theta+h));
    fm = -objective(max(-1,theta-h));
    curvature = max(0,-(fp-2*f0+fm)/(h^2));
    mu = pi*theta;
    rawKappa = curvature/pi^2;
    kappa = min(options.maxKappa,max(options.minKappa,rawKappa));

    outside = abs(grid-coarseTheta) >= 2/N;
    if any(outside)
        secondPeak = max(logPosterior(outside));
    else
        secondPeak = -Inf;
    end
    peakGap = peakValue-secondPeak;
    fallbackReason = '';
    if ~all(isfinite([theta mu kappa peakValue]))
        fallbackReason = 'nonfinite_2d_variational_posterior';
    elseif peakGap < log(2)
        fallbackReason = 'weakly_separated_marginal_peak';
    end
    peakQuality = struct( ...
        'mainPeak',peakValue,'secondPeak',secondPeak, ...
        'logPeakGap',peakGap,'coarseTheta',coarseTheta, ...
        'gridStep',step,'posteriorResultant',abs(circularMoment), ...
        'logPosteriorCurvature',curvature);
    posterior = struct( ...
        'theta',theta,'mu',mu,'kappa',kappa,'rawKappa',rawKappa, ...
        'coarseTheta',coarseTheta,'peakQuality',peakQuality, ...
        'fallbackReason',fallbackReason);
end

function value = local_axis_log_posterior( ...
        theta,R,N,otherN,d,lambda,sigma2,priorMu,priorKappa)
    a = local_steering(N,d,lambda,theta);
    value = max(real(a'*R*a),0)/(N*otherN*sigma2) + ...
        priorKappa*cos(pi*theta-priorMu);
end

function a = local_steering(N,d,lambda,theta)
    index = (1:N).'-(N+1)/2;
    a = exp(1j*2*pi*index*d*theta/lambda);
end

function kappa = local_inverse_bessel_ratio(R)
    R = min(max(real(R),0),1-1e-12);
    if R < 0.53
        kappa = 2*R + R^3 + 5*R^5/6;
    elseif R < 0.85
        kappa = -0.4 + 1.39*R + 0.43/(1-R);
    else
        kappa = 1/(R^3-4*R^2+3*R);
    end
end
