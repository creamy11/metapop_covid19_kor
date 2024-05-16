close all
clear

end_date_idx = 500;

color = colororder;

fontsize = 13;
linewidth = 1;

loc_names_kor = {'서울', '부산', '대구', '인천', '광주', '대전', '울산', '세종', '경기', ...
    '강원', '충북', '충남', '전북', '전남', '경북', '경남', '제주'};
loc_names_eng = {'Seoul', 'Busan', 'Daegu', 'Incheon', 'Gwangju', 'Daejeon', ...
    'Ulsan', 'Sejong', 'Gyeonggi', 'Gangwon', 'North Chungcheong', 'South Chungcheong', ...
    'North Jeolla', 'South Jeolla', 'North Gyeongsang', 'South Gyeongsang', 'Jeju'};

result_path = '../result';

%% Incidence Data

data_date = readmatrix('../data/incidence/COVID-19_confirmed_230303.xlsx', ...
    'Sheet', '시도별(17개시도+검역)', 'Range', 'A7:A1145', 'OutputType', 'string');

data_case = readmatrix('../data/incidence/COVID-19_confirmed_230303.xlsx', ...
    'Sheet', '시도별(17개시도+검역)', 'Range', 'C7:S1145');

%% Process

date_span = datetime(data_date);

day_num = weekday(date_span);
mon_bool = day_num == 2;

data_case(isnan(data_case)) = 0;
n_loc = size(data_case, 2);
data_case(:, end + 1) = sum(data_case, 2);

%% Plot

close all

figure(1)
hold on
plot(date_span(1:end_date_idx), data_case(1:end_date_idx, end), 'LineWidth', 1)
text(datetime('20-Feb-2020'), 1000, '1st wave', 'FontSize', fontsize)
text(datetime('13-May-2020'), 150, 'Sporadic outbreaks', 'FontSize', 10)
text(datetime('15-Aug-2020'), 500, '2nd wave', 'FontSize', fontsize)
text(datetime('28-Nov-2020'), 1300, '3rd wave', 'FontSize', fontsize)
text(datetime('20-Mar-2021'), 850, 'Vaccination', 'FontSize', fontsize)
% text(datetime('7-Jul-2021'), 8, '4th wave', 'FontSize', fontsize)
% text(datetime('24-Jan-2022'), 13.5, 'Omicron', 'FontSize', fontsize)
% xline(datetime('18-Feb-2020'), 'LineWidth', 1)
xline(datetime('5-May-2020'), 'LineWidth', 1)
xline(datetime('3-Aug-2020'), 'LineWidth', 1)
xline(datetime('11-Oct-2020'), 'LineWidth', 1)
xline(datetime('25-Feb-2021'), 'LineWidth', 1)
hold off
set(gca, 'FontSize', fontsize)
xlabel('Date', 'FontSize', 20)
ylabel('Daily cases', 'FontSize', 20)
xlim([date_span(1), date_span(end_date_idx)])
title('COVID-19 Outbreaks in South Korea', 'FontSize', 20)
pos = get(gcf, 'OuterPosition');
set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) + 200, pos(4)])
exportgraphics(gcf, sprintf('%s/case_data.eps', result_path))
exportgraphics(gcf, sprintf('%s/case_data.png', result_path))

end_date_idx = daysact(date_span(1), datetime('25-February-2021'));
mon_idx = find(mon_bool(1:end_date_idx));
n_mon = sum(mon_bool(1:end_date_idx));

for i = 1:n_loc

    figure(i + 1)
    hold on
    b = bar(date_span(1:end_date_idx), data_case(1:end_date_idx, i));
    b.FaceColor = 'flat';
    b.CData(mon_idx, :) = kron(ones(n_mon, 1), color(2, :));
    hold off
    legend('Monday')
    xlabel('Date')
    ylabel('Cases')
    title(sprintf('Daily new cases in %s (%s)', loc_names_eng{i}, loc_names_kor{i}))
    set(gca, 'FontSize', fontsize)
    pos = get(gcf, 'OuterPosition');
    set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) + 300, pos(4)])
    exportgraphics(gcf, sprintf('%s/case_data_loc%d.eps', result_path, i))

end

%% Mobility Matrix

data_mobility1 = readmatrix('../data/mobility/2020-OD-PSN-OBJ-00 전국지역간 목적 OD(250존)(2019-2050)_221121/5. OD/2019년~2050년 목적별OD(250).xlsx', ...
    'Sheet', '대존_2019', 'Range', 'C5:S182');

n_loc = 17;
n_type = 8;

M = zeros(n_loc, n_loc, n_type);% 17 locations 8 types

for i = 1:n_type

    M_i = data_mobility1(23 * (i - 1) + 1:23 * (i - 1) + n_loc, :);

    % Rearrange the order (due to Sejong)
    M_temp = [M_i(1:7, 1:7), M_i(1:7, 17), M_i(1:7, 8:16)
        M_i(17, 1:7), M_i(17, 17), M_i(17, 8:16)
        M_i(8:16, 1:7), M_i(8:16, 17), M_i(8:16, 8:16)];

    M(:, :, i) = M_temp;

end

close all

loc_names_eng = {'Seoul', 'Busan', 'Daegu', 'Incheon', 'Gwangju', 'Daejeon', ...
    'Ulsan', 'Sejong', 'Gyeonggi', 'Gangwon', 'North Chungcheong', 'South Chungcheong', ...
    'North Jeolla', 'South Jeolla', 'North Gyeongsang', 'South Gyeongsang', 'Jeju'};

title_str = {'Going to work', 'Going to school', 'Business', 'Shopping', ...
    'Going back home', 'Leisure', 'ETC', 'Total'};

for i = 1:n_type

    figure(i)
    h = heatmap(M(:, :, i));
    h.XData = loc_names_eng;
    h.YData = loc_names_eng;
    xlabel('Destination')
    ylabel('Origin')
    title(title_str{i})
    set(gca, 'FontSize', fontsize)
    exportgraphics(gcf, sprintf('%s/od_data_type%d.eps', result_path, i))

end

%% Mobility Change

intra_data = readmatrix('../data/mobility/통신모바일 인구이동량_시도별_관내.xlsx', ...
    'Range', 'B10:R63');% 2020 Feb 3rd week
inter_data = readmatrix('../data/mobility/통신모바일 인구이동량_시도별_관외.xlsx', ...
    'Range', 'B10:R63');% 2020 Feb 3rd week

week_name = readmatrix('../data/mobility/통신모바일 인구이동량_시도별_관내.xlsx', ...
    'Range', 'A10:A63', 'OutputType', 'string');

n_week = size(intra_data, 1);
n_loc = size(intra_data, 2);

week_str = {};
for i = 1:n_week

    week_current = week_name(i);
    week_current = strrep(week_current, '''', '');
    week_current = strrep(week_current, '주차', '');
    week_current = strrep(week_current, ' ', ' week');

    week_str = [week_str, week_current];

end

close all

for i = 1:n_loc

    figure(i)
    hold on
    plot(intra_data(:, i) / max(intra_data(:, i)), 'LineWidth', linewidth)
    plot(inter_data(:, i) / max(inter_data(:, i)), 'LineWidth', linewidth)
    hold off
    box on
    legend('Intra-regional movements', 'Inter-regional movements', 'Location', 'southeast')
    ylabel('Average daily movements')
    xlim([1, n_week])
    xtick_idx = get(gca, 'Xtick') + 1;
    set(gca, 'XTickLabel', week_str(xtick_idx))
    title(sprintf('%s(%s)', loc_names_eng{i}, loc_names_kor{i}))
    pos = get(gcf, 'OuterPosition');
    set(gcf, 'OuterPosition', [pos(1), pos(2), pos(3) * 1.5, pos(4)])
    set(gca, 'FontSize', fontsize)

end