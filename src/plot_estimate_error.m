close all
clear

%% Set Options

source_path = '../result/estimation/synthetic';
estimate_suffix = '_fix';

parameter_disp_idx = [1:17, 20];

n_loc = 17;
n_parameter = 23;
n_test_param = 100;
n_sim = 10;

parameter_name = {};
for i = 1:n_loc
    parameter_name = [parameter_name, sprintf('\\beta_{%d}', i)];
end
parameter_name = [parameter_name, '\mu', 'Z', '\alpha', 'D^r', 'D^u', 'G'];

color = colororder;
fontsize = 20;

result_path = source_path;

%% Compute Errors

error = zeros(n_parameter, n_test_param);
error_relative = zeros(n_parameter, n_test_param);

for i = 1:n_test_param

    data_suffix = sprintf('_%d', i);

    % Load truth
    load(sprintf('../pipeline/generate_synthetic_data/data%s.mat', ...
        data_suffix), '*_truth')

    % Load estimate
    load(sprintf('%s/estimate%s%s.mat', source_path, data_suffix, estimate_suffix), ...
        'x_post_obs', 'x_post_parameter')

    parameter_truth = [beta_truth; mu_truth; Z_truth; alpha_truth; Dr_truth; Du_truth; G_truth];

    x_post_parameter_mean = squeeze(mean(x_post_parameter, [2, 4]));

    for j = 1:n_parameter
        error(j, i) = x_post_parameter_mean(j, end) - parameter_truth(j);
    end

    error_relative(:, i) = error(:, i) ./  parameter_truth;
end

%% Plot Boxchart

close all

figure(1)
hold on
boxchart(error(parameter_disp_idx, :)', 'BoxFaceColor', color(5, :), ...
    'MarkerColor', color(5, :), 'JitterOutliers', 'on', 'MarkerStyle', '.')
yline(0, 'k')
hold off
set(gca, 'XTickLabel', parameter_name(parameter_disp_idx))
title('Error')
set(gca, 'FontSize', fontsize)
pos = get(gcf, 'OuterPosition');
set(gcf, 'OuterPosition', [0, pos(2), pos(3) * 1.5, pos(4)])
exportgraphics(gcf, sprintf('%s/error%s.png', result_path, estimate_suffix))

figure(2)
hold on
boxchart(error_relative(parameter_disp_idx, :)', 'BoxFaceColor', color(6, :), ...
    'MarkerColor', color(6, :), 'JitterOutliers', 'on', 'MarkerStyle', '.')
yline(0, 'k')
hold off
ylim([-1, 1])
set(gca, 'XTickLabel', parameter_name(parameter_disp_idx))
title('Relative error')
set(gca, 'FontSize', fontsize)
pos = get(gcf, 'OuterPosition');
set(gcf, 'OuterPosition', [0, pos(2), pos(3) * 1.5, pos(4)])
exportgraphics(gcf, sprintf('%s/relative_error%s.png', result_path, estimate_suffix))
