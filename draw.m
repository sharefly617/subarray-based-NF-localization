%% Select a release MAT result
clc;clear;close all;
newDefaultFont = 'Times';  % 替换为你想要的字体名称
set(0, 'DefaultAxesFontName', newDefaultFont);
set(0, 'defaultTextFontName', newDefaultFont);
set(groot,'DefaultFigureRenderer','painters');
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

% ---------- 统一绘图参数，放在此处开头 ----------
markerList = {'o', 's', '^'};           % 标记符号
markerSize = 8;                         % 标记大小
lineWidth = 1.5;                        % 线宽
colorList = {"#ff0000","#0723fb",'#040b02', '#994487',"#86ec84","b","m"};
% ------------------------------------------------

figure('Color', 'w');
semilogy(SNRdB, rmse(1,:), '-', ...
    'Marker', markerList{1}, 'MarkerSize', markerSize, ...
    'LineWidth', lineWidth, 'Color', colorList{1}, ...
    'MarkerFaceColor', colorList{1});
hold on;
semilogy(SNRdB, rmse(2,:), '-', ...
    'Marker', markerList{2}, 'MarkerSize', markerSize, ...
    'LineWidth', lineWidth, 'Color', colorList{2}, ...
    'MarkerFaceColor', colorList{2});
semilogy(SNRdB, rmse(3,:), '--', ...
    'Marker', markerList{3}, 'MarkerSize', markerSize, ...
    'LineWidth', lineWidth, 'Color', colorList{3}, ...
    'MarkerFaceColor', colorList{3});
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (m)');
legend('GA Ini: APLE', 'E-APLE', 'CRB', 'Location', 'southwest');
% title('Fig. 6 reproduction');
exportgraphics(gcf, ...
    fullfile(releaseRoot, 'results', 'fig6_release.png'), ...
    'Resolution', 800);

% %% APLE initial position error
% initialRMSE = results.meanApleInitialRMSE;

% figure('Color', 'w');
% semilogy(SNRdB, initialRMSE, '-o', 'LineWidth', 1.5);
% grid on;
% xlabel('SNR (dB)');
% ylabel('APLE initial RMSE (m)');
% title('APLE initial-position error');
% % exportgraphics(gcf, ...
% %     fullfile(releaseRoot, 'results', 'aple_initial_error.png'), ...
% %     'Resolution', 200);

% %% Per-trial trajectory note
% % The compact runner intentionally omits per-trial histories. Run a separate
% % diagnostic experiment if an objective or parameter trajectory is needed.
% fprintf(['No single-trial objective/parameter trajectory is available in ', ...
%     'this summary-only MAT file.\n']);