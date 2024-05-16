close all
clear

%% Set Paths

restoredefaultpath

%% Load Data

% Load mobility data
M_data = readmatrix('../data/mobility/2020-OD-PSN-OBJ-00 전국지역간 목적 OD(250존)(2019-2050)_230608/5. OD/2019년~2050년 목적별OD(250).xlsx', ...
    'Sheet', '대존_2019', 'Range', 'C5:S182');

%% Process Mobility Data

n_loc = 17;% Number of locations
n_type = 8;% Work, school, business, shopping, ...

M = zeros(n_loc, n_loc, n_type);% 17 locations 8 types
M_og = zeros(n_loc, n_loc, n_type);% Original form with intra mobility

for i = 1:n_type

    M_i = M_data(23 * (i - 1) + 1:23 * (i - 1) + n_loc, :);

    % Rearrange the order (due to Sejong)
    M_temp = [M_i(1:7, 1:7), M_i(1:7, 17), M_i(1:7, 8:16)
        M_i(17, 1:7), M_i(17, 17), M_i(17, 8:16)
        M_i(8:16, 1:7), M_i(8:16, 17), M_i(8:16, 8:16)];

    % Transpose the O\D matrix to make D\O matrix
    M_temp = M_temp';

    M_og(:, :, i) = M_temp;

    % Remove diagonal term
    for j = 1:n_loc
        M_temp(j, j) = 0;
    end

    M(:, :, i) = M_temp;

end

% Check the data
M_sum = sum(M(:, :, 1:end - 1), 3);
bool = sum(M_sum == M(:, :, end), 'all');
if bool == n_loc ^ 2
    disp('Total matrix = Sum of the other matrices')
else
    disp('Error')
end

%% Plot Heatmap

close all

fontsize = 12;

loc_names_eng = {'Seoul', 'Busan', 'Daegu', 'Incheon', 'Gwangju', 'Daejeon', ...
    'Ulsan', 'Sejong', 'Gyeonggi', 'Gangwon', 'North Chungcheong', 'South Chungcheong', ...
    'North Jeolla', 'South Jeolla', 'North Gyeongsang', 'South Gyeongsang', 'Jeju'};

figure(1)
heatmap(loc_names_eng, loc_names_eng, M_og(:, :, 5))
title('')
set(gca, 'FontSize', fontsize)
exportgraphics(gcf, '../result/heatmap_mobility.png')

%% Save as mat

save(sprintf('../pipeline/process_mobility_data/mobility.mat'), 'M', 'M_og')
