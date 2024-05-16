close all
clear

%% Set Options

dates = [datetime('5-May-2020')
    datetime('3-August-2020')
    datetime('11-October-2020')
    datetime('25-February-2021')];

result_path = '../result/estimation';
suffix = '_fix';

beta_islocal = true;

loc2plot = [1, 3, 9, 15];% location to plot Re

ci_case_isPlot = true;
ci_beta_isPlot = true;
ci_global_isPlot = true;
ci_Rt_isPlot = true;
ci_state_isPlot = false;

tab_beta_is_gen = true;
tab_global_is_gen = true;
tab_Rt_is_gen = true;
tab_total_case_is_gen = true;

ci_level = 95;
ci_level_span = 0:5:100;

% Credible interval lower and upper limits
ci_ll_level = (100 - ci_level) / 2;
ci_ul_level = ci_level + (100 - ci_level) / 2;

%% Set Figure Options

set(0, 'DefaultFigureVisible', 'on')

color = colororder;

linewidth = 2;
facealpha = 0.1;
facealpha_p = 0.3;
fontsize = 20;

resolution = 500;

loc_names_kor = {'서울', '부산', '대구', '인천', '광주', '대전', '울산', '세종', '경기', ...
    '강원', '충북', '충남', '전북', '전남', '경북', '경남', '제주'};
loc_names_eng = {'Seoul', 'Busan', 'Daegu', 'Incheon', 'Gwangju', 'Daejeon', ...
    'Ulsan', 'Sejong', 'Gyeonggi', 'Gangwon', 'North Chungcheong', 'South Chungcheong', ...
    'North Jeolla', 'South Jeolla', 'North Gyeongsang', 'South Gyeongsang', 'Jeju'};

state_name = {'S', 'E', 'Ir', 'Iu', 'Q', 'R'};
state_name_latex = {'S', 'E', 'I^r', 'I^u', 'Q', 'R'};
parameter_name = {'beta', 'mu', 'Z', 'alpha', 'Dr', 'Du', 'G'};
parameter_name_latex = {'\beta', '\mu', 'Z', '\alpha', 'D^r', 'D^u', 'G'};

%% Load Data

% Load population data
pop_data = readmatrix('../data/demographic/population_census_2020.xlsx', ...
    'Sheet', 'Si-Do', 'Range', 'B4:B20')';% 2020 population census

%% Load Estimate

load(sprintf('%s/estimate%s.mat', result_path, suffix), 'case_data', 'obs_data', ...
    '*_lb', '*_ub', 'date_span', 'x_post_state', 'x_post_obs', 'x_post_incidence_u', 'x_post_parameter')

%% Generate Secondary Variables

n_date = size(dates, 1);
n_loc2plot = length(loc2plot);

n_loc = size(x_post_obs, 1);
n_ens = size(x_post_obs, 2);
n_obs = size(x_post_obs, 3);
n_real = size(x_post_obs, 4);
n_state_type = size(x_post_state, 1) / n_loc;

parameter_lb = [mu_lb; Z_lb; alpha_lb; Dr_lb; Du_lb; G_lb];
parameter_ub = [mu_ub; Z_ub; alpha_ub; Dr_ub; Du_ub; G_ub];
n_global_parameter = length(parameter_lb);% Number of parameters other than beta

x_post_state_p = reshape(permute(x_post_state, [1, 3, 2, 4]), n_loc * n_state_type, n_obs, n_ens * n_real);
x_post_obs_p = reshape(permute(x_post_obs, [1, 3, 2, 4]), n_loc, n_obs, n_ens * n_real);
x_post_incidence_u_p = reshape(permute(x_post_incidence_u, [1, 3, 2, 4]), n_loc, n_obs, n_ens * n_real);
x_post_parameter_p = reshape(permute(x_post_parameter, [1, 3, 2, 4]), n_loc + n_global_parameter, n_obs, n_ens * n_real);
x_post_incidence_cum_p = cumsum(x_post_obs_p + x_post_incidence_u_p, 2);
x_post_incidence_cum_p(end + 1, :, :) = sum(x_post_incidence_cum_p, 1);% National cases

x_post_R_p = x_post_state_p(5 * n_loc + (1:n_loc), :, :);
x_post_R_p_mean = mean(x_post_R_p, 3);
x_post_R_p_ll = prctile(x_post_R_p, ci_ll_level, 3);
x_post_R_p_ul = prctile(x_post_R_p, ci_ul_level, 3);

x_post_state_mean = squeeze(mean(x_post_state, 2));
x_post_obs_mean = squeeze(mean(x_post_obs, 2));
x_post_parameter_mean = squeeze(mean(x_post_parameter, 2));

x_post_state_p_mean = squeeze(mean(x_post_state_p, 3));
x_post_obs_p_mean = squeeze(mean(x_post_obs_p, 3));
x_post_parameter_p_mean = squeeze(mean(x_post_parameter_p, 3));
x_post_incidence_cum_p_mean = squeeze(mean(x_post_incidence_cum_p, 3));

x_post_state_ll = squeeze(prctile(x_post_state, ci_ll_level, 2));
x_post_obs_ll = squeeze(prctile(x_post_obs, ci_ll_level, 2));
x_post_parameter_ll = squeeze(prctile(x_post_parameter, ci_ll_level, 2));

x_post_state_p_ll = squeeze(prctile(x_post_state_p, ci_ll_level, 3));
x_post_obs_p_ll = squeeze(prctile(x_post_obs_p, ci_ll_level, 3));
x_post_parameter_p_ll = squeeze(prctile(x_post_parameter_p, ci_ll_level, 3));
x_post_incidence_cum_p_ll = squeeze(prctile(x_post_incidence_cum_p, ci_ll_level, 3));

x_post_state_ul = squeeze(prctile(x_post_state, ci_ul_level, 2));
x_post_obs_ul = squeeze(prctile(x_post_obs, ci_ul_level, 2));
x_post_parameter_ul = squeeze(prctile(x_post_parameter, ci_ul_level, 2));

x_post_state_p_ul = squeeze(prctile(x_post_state_p, ci_ul_level, 3));
x_post_obs_p_ul = squeeze(prctile(x_post_obs_p, ci_ul_level, 3));
x_post_parameter_p_ul = squeeze(prctile(x_post_parameter_p, ci_ul_level, 3));
x_post_incidence_cum_p_ul = squeeze(prctile(x_post_incidence_cum_p, ci_ul_level, 3));

date_init = date_span(1);
date_span_obs = date_span(2:end);
dates_idx = daysact(date_init, dates);

str4legend = {};
for i = 1:length(dates)
    str4legend = [str4legend, sprintf('%s', string(dates(i)))];
end

%% Observation Mean and CI

close all

if ci_case_isPlot

    xx = [date_span_obs, fliplr(date_span_obs)];

    for i = 1:n_loc

        case_i = case_data(:, i);
        data_i = obs_data(:, i);

        figure(i)
        hold on
        y_ll = squeeze(x_post_obs_p_ll(i, :));
        y_ul = squeeze(x_post_obs_p_ul(i, :));
        yy = [y_ll, fliplr(y_ul)];
        f = fill(xx, yy, color(3, :), 'FaceAlpha', facealpha_p, 'EdgeColor', 'none');
        plot(date_span_obs, squeeze(x_post_obs_p_mean(i, :)), 'LineWidth', ...
            linewidth, 'Color', color(3, :))
        p2 = plot(date_span_obs, case_i, '.', 'MarkerSize', 10, ...
            'Color', color(1, :));
        p3 = plot(date_span_obs, data_i, '.', 'MarkerSize', 10, ...
            'Color', color(2, :));
        p4 = plot(date_span_obs, nan(n_obs, 1), 'LineWidth', linewidth, ...
            'Color', color(3, :));
        hold off
        box on
        switch i
            case {3, 15}
                legend([p2, p3, p4, f], 'Data', 'Smoothed data', 'Posterior mean', ...
                    sprintf('%d%% CI', ci_level), 'Location', 'northeast')
            otherwise
                legend([p2, p3, p4, f], 'Data', 'Smoothed data', 'Posterior mean', ...
                    sprintf('%d%% CI', ci_level), 'Location', 'northwest')
        end
        xlim([date_span_obs(1), date_span_obs(end)])
        xlabel('Date')
        ylabel('Daily cases')
        title(sprintf('%s (%s)', loc_names_eng{i}, loc_names_kor{i}))
        set(gca, 'FontSize', fontsize)
        pos = get(gcf, 'OuterPosition');
        set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
        exportgraphics(gcf, sprintf('%s/ci%d_case%d%s.png', result_path, ci_level, i, suffix))
    end
end

%% beta Mean and CI

close all

if ci_beta_isPlot

    xx = [date_span(2:end), fliplr(date_span(2:end))];

    for i = 1:n_loc

        figure(i)
        hold on
        y_ll = squeeze(x_post_parameter_p_ll(i, :));
        y_ul = squeeze(x_post_parameter_p_ul(i, :));
        yy = [y_ll, fliplr(y_ul)];
        f = fill(xx, yy, color(3, :), 'FaceAlpha', facealpha_p, 'EdgeColor', 'none');
        plot(date_span(2:end), squeeze(x_post_parameter_p_mean(i, :)), ...
            'LineWidth', linewidth, 'Color', color(3, :))
        p = plot(date_span(2:end), nan(n_obs, 1), 'LineWidth', linewidth, 'Color', color(3, :));
        hold off
        box on
        legend([p, f], 'Posterior mean', sprintf('%d%% CI', ci_level), ...
            'Location', 'northwest')
        xlim([date_span(2), date_span(end)])
        ylim([beta_lb(1), beta_ub(1)])
        xlabel('Date')
        if beta_islocal
            ylabel(sprintf('%s_{%d}', parameter_name_latex{1}, i))
            title(sprintf('%s (%s)', loc_names_eng{i}, loc_names_kor{i}))
        else
            title(sprintf('%s', parameter_name_latex{1}))
        end
        set(gca, 'FontSize', fontsize)
        pos = get(gcf, 'OuterPosition');
        set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
        if beta_islocal
            exportgraphics(gcf, sprintf('%s/ci%d_%s%d%s.png', result_path, ci_level, parameter_name{1}, i, suffix))
        else
            exportgraphics(gcf, sprintf('%s/ci%d_%s%s.png', result_path, ci_level, parameter_name{1}, suffix))
            break
        end
    end
end

%% beta Table

if tab_beta_is_gen

    mean_ci = cell(n_loc, n_date);

    rowname = {};

    for i = 1:n_loc

        for j = 1:n_date

            mean_ci{i, j} = sprintf('%.2f (%.2f, %.2f)', ...
                x_post_parameter_p_mean(i, dates_idx(j)), ...
                x_post_parameter_p_ll(i, dates_idx(j)), ...
                x_post_parameter_p_ul(i, dates_idx(j)));
        end

        rowname = [rowname, sprintf('beta%d', i)];
    end

    mean_ci_table = cell2table(mean_ci, 'RowNames', rowname, 'VariableNames', str4legend);
    writetable(mean_ci_table, sprintf('%s/table_beta%s.csv', result_path, suffix), ...
        'WriteRowNames', true)
end

%% Global Parameters Mean and CI

close all

if ci_global_isPlot

    xx = [date_span(2:end), fliplr(date_span(2:end))];

    for i = 1:n_global_parameter

        figure(i)
        hold on
        y_ll = squeeze(x_post_parameter_p_ll(n_loc + i, :));
        y_ul = squeeze(x_post_parameter_p_ul(n_loc + i, :));
        yy = [y_ll, fliplr(y_ul)];
        f = fill(xx, yy, color(3, :), 'FaceAlpha', facealpha_p, 'EdgeColor', 'none');
        plot(date_span(2:end), squeeze(x_post_parameter_p_mean(n_loc + i, :)), ...
            'LineWidth', linewidth, 'Color', color(3, :))
        p = plot(date_span(2:end), nan(n_obs, 1), 'LineWidth', linewidth, 'Color', color(3, :));
        hold off
        box on
        legend([p, f], 'Posterior mean', sprintf('%d%% CI', ci_level), ...
            'Location', 'northwest')
        xlim([date_span(2), date_span(end)])
        xlabel('Date')
        ylabel(sprintf('%s', parameter_name_latex{i + 1}))
        set(gca, 'FontSize', fontsize)
        pos = get(gcf, 'OuterPosition');
        set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
        exportgraphics(gcf, sprintf('%s/ci%d_%s%s.png', result_path, ...
            ci_level, parameter_name{i + 1}, suffix))
    end
end

%% Global Table

if tab_global_is_gen

    mean_ci = cell(n_global_parameter, n_date);

    rowname = {};

    for i = 1:n_global_parameter

        for j = 1:n_date

            mean_ci{i, j} = sprintf('%.2f (%.2f, %.2f)', ...
                x_post_parameter_p_mean(n_loc + i, dates_idx(j)), ...
                x_post_parameter_p_ll(n_loc + i, dates_idx(j)), ...
                x_post_parameter_p_ul(n_loc + i, dates_idx(j)));
        end

        rowname = [rowname, sprintf('%s', parameter_name{1 + i})];
    end

    mean_ci_table = cell2table(mean_ci, 'RowNames', rowname, 'VariableNames', str4legend);
    writetable(mean_ci_table, sprintf('%s/table_global%s.csv', result_path, suffix), ...
        'WriteRowNames', true)
end

%% Reproductive Numbers

close all

if ci_Rt_isPlot

    figure(1)
    hold on

    figure(2)
    hold on

    p1 = [];
    p2 = [];
    ii = 0;
    for i = loc2plot

        ii = ii + 1;

        beta_i_post = x_post_parameter_p(i, :, :);
        mu_post = x_post_parameter_p(n_loc + 1, :, :);
        alpha_post = x_post_parameter_p(n_loc + 3, :, :);
        Dr_post = x_post_parameter_p(n_loc + 4, :, :);
        Du_post = x_post_parameter_p(n_loc + 5, :, :);

        S_i_post = x_post_state_p(i, :, :);
        N_i = pop_data(i);

        Rt = beta_i_post .* (alpha_post .* Dr_post ...
            + (1 - alpha_post) .* mu_post .* Du_post);
        Re = beta_i_post .* (alpha_post .* Dr_post ...
            + (1 - alpha_post) .* mu_post .* Du_post) .* S_i_post ./ N_i;
        Re_r = beta_i_post .* alpha_post .* Dr_post .* S_i_post ./ N_i;
        Re_u = beta_i_post .* (1 - alpha_post) .* mu_post .* Du_post .* S_i_post ./ N_i;

        figure(1)
        p_temp = plot(date_span_obs, mean(Rt, 3), 'LineWidth', linewidth, ...
            'Color', color(ii + 3, :));
        p1 = [p1, p_temp];

        figure(2)
        p_temp = plot(date_span_obs, mean(Re, 3), 'LineWidth', linewidth, ...
            'Color', color(ii + 3, :));
        p2 = [p2, p_temp];

    end

    figure(1)
    yline(1, 'k')
    hold off
    box on
    legend(p1, loc_names_eng(loc2plot))
    xlabel('Date')
    ylabel('R_t')
    title('Time-varying basic reproduction numbers')
    xlim([date_span_obs(1), date_span_obs(end)])
    set(gca, 'FontSize', fontsize)
    pos = get(gcf, 'OuterPosition');
    set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
    exportgraphics(gcf, sprintf('%s/Rt%s.eps', result_path, suffix))
    exportgraphics(gcf, sprintf('%s/Rt%s.png', result_path, suffix))

    figure(2)
    yline(1, 'k')
    text(datetime('20-Feb-2020'), 3.6, '1st wave', 'FontSize', 13)
    text(datetime('10-May-2020'), 1.6, 'Sporadic outbreaks', 'FontSize', 13)
    text(datetime('5-Aug-2020'), 2.4, '2nd wave', 'FontSize', 13)
    text(datetime('24-Nov-2020'), 1.8, '3rd wave', 'FontSize', 13)
    hold off
    box on
    legend(p2, loc_names_eng(loc2plot))
    xlabel('Date')
    ylabel('R_e')
    title('Effective reproduction numbers')
    xlim([date_span_obs(1), date_span_obs(end)])
    ylim([0, 4])
    set(gca, 'FontSize', fontsize)
    set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
    exportgraphics(gcf, sprintf('%s/Re%s.eps', result_path, suffix))
    exportgraphics(gcf, sprintf('%s/Re%s.png', result_path, suffix))
end

%% Reproductive Number Table

dates_idx_ = [1; dates_idx];
str4legend_ = [string(date_init + 1), str4legend];

if tab_Rt_is_gen

    mean_ci_Rt = cell(n_loc2plot, n_date + 1);
    mean_ci_Re = cell(n_loc2plot, n_date + 1);

    rowname_Rt = {};
    rowname_Re = {};

    for i = 1:n_loc2plot

        beta_i_post = x_post_parameter_p(loc2plot(i), :, :);
        mu_post = x_post_parameter_p(n_loc + 1, :, :);
        alpha_post = x_post_parameter_p(n_loc + 3, :, :);
        Dr_post = x_post_parameter_p(n_loc + 4, :, :);
        Du_post = x_post_parameter_p(n_loc + 5, :, :);

        S_i_post = x_post_state_p(loc2plot(i), :, :);
        N_i = pop_data(loc2plot(i));

        Rt = squeeze(beta_i_post .* (alpha_post .* Dr_post ...
            + (1 - alpha_post) .* mu_post .* Du_post));

        Re = squeeze(beta_i_post .* (alpha_post .* Dr_post ...
            + (1 - alpha_post) .* mu_post .* Du_post) .* S_i_post ./ N_i);

        Rt_mean = mean(Rt, 2);
        Re_mean = mean(Re, 2);

        Rt_ll = prctile(Rt, ci_ll_level, 2);
        Rt_ul = prctile(Rt, ci_ul_level, 2);

        Re_ll = prctile(Re, ci_ll_level, 2);
        Re_ul = prctile(Re, ci_ul_level, 2);

        for j = 1:n_date + 1
            mean_ci_Rt{i, j} = sprintf('%.2f (%.2f, %.2f)', Rt_mean(dates_idx_(j)), ...
                Rt_ll(dates_idx_(j)), Rt_ul(dates_idx_(j)));
            mean_ci_Re{i, j} = sprintf('%.2f (%.2f, %.2f)', Re_mean(dates_idx_(j)), ...
                Re_ll(dates_idx_(j)), Re_ul(dates_idx_(j)));
        end

        rowname_Rt = [rowname_Rt, sprintf('Rt%d', loc2plot(i))];
        rowname_Re = [rowname_Re, sprintf('Re%d', loc2plot(i))];
    end

    mean_ci_Rt_table = cell2table(mean_ci_Rt, 'RowNames', rowname_Rt, ...
        'VariableNames', str4legend_);
    writetable(mean_ci_Rt_table, sprintf('%s/table_Rt%s.csv', result_path, suffix), ...
        'WriteRowNames', true)

    mean_ci_Re_table = cell2table(mean_ci_Re, 'RowNames', rowname_Re, ...
        'VariableNames', str4legend_);
    writetable(mean_ci_Re_table, sprintf('%s/table_Re%s.csv', result_path, suffix), ...
        'WriteRowNames', true)
end

%% Total Case Table

pop_data_sum = [pop_data, sum(pop_data)];

if tab_total_case_is_gen

    mean_ci = cell(n_loc + 1, 2 * n_date + 1);

    for i = 1:n_loc + 1

        mean_ci{i, 1} = sprintf('%.0f', pop_data_sum(i));

        str4legend_ = {'Population'};

        % Total cases computed from the flows
        for j = 1:n_date

            mean_ci{i, 2 * j} = sprintf('%.0f (%.0f, %.0f)', ...
                x_post_incidence_cum_p_mean(i, dates_idx(j)), ...
                x_post_incidence_cum_p_ll(i, dates_idx(j)), ...
                x_post_incidence_cum_p_ul(i, dates_idx(j)));

            mean_ci{i, 2 * j + 1} = sprintf('%.2f (%.2f, %.2f)', ...
                x_post_incidence_cum_p_mean(i, dates_idx(j)) / pop_data_sum(i) * 100, ...
                x_post_incidence_cum_p_ll(i, dates_idx(j)) / pop_data_sum(i) * 100, ...
                x_post_incidence_cum_p_ul(i, dates_idx(j)) / pop_data_sum(i) * 100);

            str4legend_ = [str4legend_, sprintf('%s', str4legend{j}), ...
                sprintf('%s (%%)', str4legend{j})];
        end
    end

    mean_ci_table = cell2table(mean_ci, 'RowNames', [loc_names_eng, 'Nationwide'], ...
        'VariableNames', str4legend_);
    writetable(mean_ci_table, sprintf('%s/table_total_case%s.csv', ...
        result_path, suffix), 'WriteRowNames', true)

end

%% State Variables Mean and CI

close all

y_ulim = [6e+7, 2.5e+5, 2e+4, 1e+5, 1.5e+4, 5e+7];

if ci_state_isPlot

    for i = 1:n_state_type

        close all

        for j = 1:n_loc

            loc_idx = (i - 1) * n_loc + j;

            figure(j)
            hold on
            xx = [date_span(2:end), fliplr(date_span(2:end))];
            y_ll = x_post_state_p_ll(loc_idx, :);
            y_ul = x_post_state_p_ul(loc_idx, :);
            yy = [y_ll, fliplr(y_ul)];
            fill(xx, yy, color(3, :), 'FaceAlpha', facealpha_p, 'EdgeColor', 'none')
            plot(date_span(2:end), squeeze(x_post_state_p_mean(loc_idx, :)), ...
                'LineWidth', linewidth, 'Color', color(3, :))
            p2 = plot(date_span(2:end), nan(n_obs, 1), 'LineWidth', linewidth, 'Color', color(3, :));
            hold off
            legend(p2, 'Posterior', 'Location', 'best')
            xlim([date_span(2), date_span(end)])
            xlabel('Date')
            title(sprintf('%s in %s (%s)', state_name_latex{i}, loc_names_eng{j}, loc_names_kor{j}))
            set(gca, 'FontSize', fontsize)
            pos = get(gcf, 'OuterPosition');
            set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
            exportgraphics(gcf, sprintf('%s/ci%d_%s%d%s.png', result_path, ci_level, state_name{i}, j, suffix))
        end

        figure(j + 1)
        hold on
        xx = [date_span(2:end), fliplr(date_span(2:end))];
        y_ll = sum(x_post_state_p_ll((i - 1) * n_loc + 1:i * n_loc, :), 1);
        y_ul = sum(x_post_state_p_ul((i - 1) * n_loc + 1:i * n_loc, :), 1);
        yy = [y_ll, fliplr(y_ul)];
        fill(xx, yy, color(3, :), 'FaceAlpha', facealpha_p, 'EdgeColor', 'none')
        plot(date_span(2:end), sum(x_post_state_p_mean((i - 1) * n_loc + 1:i * n_loc, :), 1), ...
            'LineWidth', linewidth, 'Color', color(3, :))
        p2 = plot(date_span(2:end), nan(n_obs, 1), 'LineWidth', linewidth, 'Color', color(3, :));
        hold off
        box on
        legend(p2, 'Posterior', 'Location', 'best')
        xlim([date_span(2), date_span(end)])
        ylim([0, y_ulim(i)])
        xlabel('Date')
        title(sprintf('%s', state_name_latex{i}))
        set(gca, 'FontSize', fontsize)
        pos = get(gcf, 'OuterPosition');
        set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
        exportgraphics(gcf, sprintf('%s/ci%d_%s%s.png', result_path, ci_level, state_name{i}, suffix))
    end
end
