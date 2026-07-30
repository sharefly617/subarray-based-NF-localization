%% RMSEvsSNR - standalone Fig. 6 reproduction
% Edit this block to change the experiment. All functions are resolved from
% this release; the original code/ directory is never called.
clc; close all;
MC = 60*10;
if exist('RELEASE_SMOKE_MC', 'var')
    MC = RELEASE_SMOKE_MC;
end
randomSeed = 12345;
SNRdB = 0:5:30;
fHz = 28e9;
lambda = 299792458/fHz;
Nx = 150; Ny = 150; d = lambda/2; sigma2 = 1;
rRange = [30 45]; 
omegaRange = [0.1 2*pi]; 
phiRange = [0.1, pi/2];
T2 = 50; GAIterations = 50;
useParallel = true;
numWorkers = 60;
outputName = 'fig6_release_results.mat';

%% Release paths and experiment initialization
releaseRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(releaseRoot, 'APLE'));
addpath(fullfile(releaseRoot, 'E-APLE'));
addpath(fullfile(releaseRoot, 'optional'));
outputFile = fullfile(releaseRoot, 'results', outputName);
if isfile(outputFile)
    warning('RMSEvsSNR:Overwrite', ...
        'Output exists and will be overwritten: %s', outputFile);
end
if mod(Nx, 3) ~= 0 || mod(Ny, 3) ~= 0
    error('RMSEvsSNR:Partition', ...
        'Nx and Ny must be divisible by 3.');
end
rng(randomSeed, 'twister');
% rng('shuffle', 'twister');
arrayConfig = struct('Nx', Nx, 'Ny', Ny, 'd', d);
partition = partition_square_subarrays(Nx, Ny, 3);
S = numel(SNRdB);
truePolar = [ ...
    rRange(1) + diff(rRange)*rand(MC, 1), ...
    omegaRange(1) + diff(omegaRange)*rand(MC, 1), ...
    phiRange(1) + diff(phiRange)*rand(MC, 1)];

%% Parallel trial execution
% Each parfor iteration owns one full trial. The explicit per-SNR seed makes
% serial and parallel execution use the same observation/noise realization.
trialData = cell(MC, 1);
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(numWorkers)
            parpool('local');
        else
            localCluster = parcluster('local');
            activeWorkers = min(numWorkers, localCluster.NumWorkers);
            if activeWorkers < numWorkers
                warning('RMSEvsSNR:WorkerCap', ...
                    'Requested %d workers; local profile permits %d.', ...
                    numWorkers, activeWorkers);
            end
            parpool('local', activeWorkers);
        end
    end
    parfor t = 1:MC
        trialData{t} = local_run_trial( ...
            t, truePolar(t,:), SNRdB, randomSeed, lambda, sigma2, ...
            arrayConfig, partition, GAIterations, T2);
    end
else
    for t = 1:MC
        trialData{t} = local_run_trial( ...
            t, truePolar(t,:), SNRdB, randomSeed, lambda, sigma2, ...
            arrayConfig, partition, GAIterations, T2);
        fprintf('RMSEvsSNR: %d/%d trials complete.\n', t, MC);
    end
end

%% Assemble trial records and compute curves
truePosition = zeros(MC, 3);
alpha = zeros(MC, S);
apInitial = nan(MC, S, 3);
estimate = nan(MC, S, 2, 3);
squaredError = nan(MC, S, 2);
crbTrace = nan(MC, S);
crbFailure = cell(MC, S);
apDiagnostics = cell(MC, S);
optimizationDiagnostics = cell(MC, S, 2);
success = false(MC, S, 2);
for t = 1:MC
    z = trialData{t};
    truePosition(t,:) = z.truePosition;
    alpha(t,:) = z.alpha;
    apInitial(t,:,:) = z.apInitial;
    estimate(t,:,:,:) = z.estimate;
    squaredError(t,:,:) = z.squaredError;
    crbTrace(t,:) = z.crbTrace;
    crbFailure(t,:) = z.crbFailure;
    apDiagnostics(t,:) = z.apDiagnostics;
    optimizationDiagnostics(t,:,:) = z.optimizationDiagnostics;
    success(t,:,:) = z.success;
end
rmse = [ ...
    sqrt(mean(squaredError(:,:,1), 1, 'omitnan')); ...
    sqrt(mean(squaredError(:,:,2), 1, 'omitnan')); ...
    sqrt(mean(crbTrace, 1, 'omitnan'))];

%% Save and display the result
trueP = reshape(truePosition, MC, 1, 1, 3);
apInitial4 = reshape(apInitial, MC, S, 1, 3);
% The coordinate dimension is the fourth dimension after explicit reshaping.
initialError = sqrt(sum((apInitial4 - trueP).^2, 4));
meanApleInitialRMSE = squeeze(sqrt(mean(initialError.^2, 1, 'omitnan'))).';
meanSquaredError = squeeze(mean(squaredError, 1, 'omitnan'));
successRate = squeeze(mean(success, 1, 'omitnan')).';
crbValidCount = sum(isfinite(crbTrace), 1);
crbFailureCount = MC - crbValidCount;
config = struct( ...
    'MC', MC, 'randomSeed', randomSeed, 'SNRdB', SNRdB, ...
    'fHz', fHz, 'lambda', lambda, 'Nx', Nx, 'Ny', Ny, 'd', d, ...
    'sigma2', sigma2, 'rRange', rRange, 'omegaRange', omegaRange, ...
    'phiRange', phiRange, 'T2', T2, 'GAIterations', GAIterations, ...
    'useParallel', useParallel, 'numWorkers', numWorkers, ...
    'methods', {{'GA Ini: APLE', 'E-APLE', 'CRB'}});
results = struct( ...
    'config', config, 'rmse', rmse, ...
    'meanApleInitialRMSE', meanApleInitialRMSE, ...
    'meanSquaredError', meanSquaredError, ...
    'successRate', successRate, ...
    'crbValidCount', crbValidCount, ...
    'crbFailureCount', crbFailureCount, ...
    'status', 'complete', ...
    'storageNote', 'Only Monte Carlo summaries are saved; per-trial arrays and diagnostics are intentionally omitted.');
save(outputFile, 'results', '-v7.3');

figure('Color', 'w');
semilogy(SNRdB, rmse(1,:), '-o', 'LineWidth', 1.5);
hold on;
semilogy(SNRdB, rmse(2,:), '-s', 'LineWidth', 1.5);
semilogy(SNRdB, rmse(3,:), '-^', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE (m)');
legend(config.methods, 'Location', 'southwest');
title('Fig. 6: RMSE versus SNR');
fprintf('Saved %s\n', outputFile);

function z = local_run_trial(t, polar, SNRdB, randomSeed, lambda, ...
        sigma2, arrayConfig, partition, GAIterations, T2)
    % One worker executes all SNR points for one true position.
    S = numel(SNRdB);
    p = spherical_to_cartesian(polar(1), polar(2), polar(3));
    % Keep the physical position as a column vector. This prevents implicit
    % 3-by-3 expansion when subtracting a row/column estimate below.
    p = p(:);
    z.truePosition = p(:).';
    z.alpha = zeros(1, S);
    z.apInitial = nan(1, S, 3);
    z.estimate = nan(1, S, 2, 3);
    z.squaredError = nan(1, S, 2);
    z.crbTrace = nan(1, S);
    z.crbFailure = cell(1, S);
    z.apDiagnostics = cell(1, S);
    z.optimizationDiagnostics = cell(1, S, 2);
    z.success = false(1, S, 2);
    for s = 1:S
        z.alpha(s) = sqrt(10^(SNRdB(s)/10));
        rng(randomSeed + 100000*t + s, 'twister');
        [y, ~, ~, ~] = generate_nearfield_snapshot( ...
            p, z.alpha(s), lambda, sigma2, arrayConfig, true);
        ap = estimate_aple(y, lambda, sigma2, arrayConfig, partition);
        z.apInitial(1,s,:) = reshape(ap.position, 1, 1, 3);
        z.apDiagnostics{1,s} = ap;
        initial = struct('position', ap.position, 'source', 'APLE');
        try
            ga = eaple_joint_gradient_ascent( ...
                y, lambda, sigma2, arrayConfig, initial, ...
                'maxIterations', GAIterations);
            z.estimate(1,s,1,:) = reshape(ga.position, 1, 1, 1, 3);
            z.squaredError(1,s,1) = sum((p-ga.position).^2);
            z.success(1,s,1) = all(isfinite(ga.position));
            z.optimizationDiagnostics{1,s,1} = ga;
        catch ME
            z.optimizationDiagnostics{1,s,1} = struct( ...
                'error', [ME.identifier ': ' ME.message]);
        end
        try
            ea = estimate_eaple( ...
                y, lambda, sigma2, arrayConfig, initial, ...
                'maxIterations', T2);
            z.estimate(1,s,2,:) = reshape(ea.position, 1, 1, 1, 3);
            z.squaredError(1,s,2) = sum((p-ea.position).^2);
            z.success(1,s,2) = all(isfinite(ea.position));
            z.optimizationDiagnostics{1,s,2} = ea;
        catch ME
            z.optimizationDiagnostics{1,s,2} = struct( ...
                'error', [ME.identifier ': ' ME.message]);
        end
        try
            crb = estimate_crb(p, z.alpha(s), lambda, sigma2, arrayConfig);
            z.crbTrace(s) = trace(crb.positionCRB);
        catch ME
            % Preserve the exact CRB diagnostic. Do not regularize or replace
            % an invalid FIM with a fabricated finite value.
            z.crbFailure{s} = [ME.identifier ': ' ME.message];
        end
    end
end
