%% Select a release MAT result
releaseRoot = fileparts(mfilename('fullpath'));
resultFile = fullfile(releaseRoot, 'results', 'fig6_release_results.mat');

% If the default file is absent, let the user choose another result.
if ~isfile(resultFile)
    [name, path] = uigetfile( ...
        fullfile(releaseRoot, 'results', '*.mat'), ...
        'Select a result MAT file');
    if isequal(name, 0)
        error('draw:NoFile', 'No MAT file selected.');
    end
    resultFile = fullfile(path, name);
end

loaded = load(resultFile, 'results');
results = loaded.results;
fprintf('Loaded %s\n', resultFile);

%% Fig. 6 RMSE curves
SNRdB = results.config.SNRdB;
rmse = results.rmse;

figure('Color', 'w');
semilogy(SNRdB, rmse(1,:), '-o', 'LineWidth', 1.5);
hold on;
semilogy(SNRdB, rmse(2,:), '-s', 'LineWidth', 1.5);
semilogy(SNRdB, rmse(3,:), '-^', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (m)');
legend('GA Ini: APLE', 'E-APLE', 'CRB', 'Location', 'southwest');
title('Fig. 6 reproduction');
% exportgraphics(gcf, ...
%     fullfile(releaseRoot, 'results', 'fig6_release.png'), ...
%     'Resolution', 200);

%% APLE initial position error
trueP = reshape(results.truePosition, size(results.truePosition, 1), 1, 1, 3);
initialError = sqrt(sum( ...
    (results.apleInitialPosition - trueP).^2, 4));
initialRMSE = squeeze(sqrt(mean(initialError.^2, 1)));

figure('Color', 'w');
semilogy(SNRdB, initialRMSE, '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('APLE initial RMSE (m)');
title('APLE initial-position error');
% exportgraphics(gcf, ...
%     fullfile(releaseRoot, 'results', 'aple_initial_error.png'), ...
%     'Resolution', 200);

%% One-trial objective and parameter trajectory
trial = 1;
snrIndex = 1;
methodIndex = 2;
diag = results.optimizationDiagnostics{trial, snrIndex, methodIndex};

if isfield(diag, 'objectiveHistory')
    figure('Color', 'w');
    iterations = 0:numel(diag.objectiveHistory)-1;
    plot(iterations, diag.objectiveHistory, '-o');
    grid on;
    xlabel('Iteration');
    ylabel('Profile objective');
    title(sprintf( ...
        '%s objective, trial %d, %g dB', ...
        results.config.methods{methodIndex}, trial, SNRdB(snrIndex)));
    % exportgraphics(gcf, ...
    %     fullfile(releaseRoot, 'results', 'single_trial_objective.png'), ...
    %     'Resolution', 200);
end

if isfield(diag, 'parameterHistory')
    figure('Color', 'w');
    plot(diag.parameterHistory.');
    grid on;
    xlabel('Iteration');
    ylabel('Parameter value');
    legend('r (m)', 'omega (rad)', 'phi (rad)', 'Location', 'best');
    title('Single-trial parameter trajectory');
    % exportgraphics(gcf, ...
    %     fullfile(releaseRoot, 'results', 'single_trial_parameters.png'), ...
    %     'Resolution', 200);
end
