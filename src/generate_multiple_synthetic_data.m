close all
clear

%% Set Paths

addpath 'model' 'solver'

%% Set Simulation Options

rng_seed = 0;
rng(rng_seed)

n_test_param = 150;

oev_scaler = 5;

result_path = '../pipeline/generate_synthetic_data';

%% Load Data

% Load population data
pop_data = readmatrix('../data/demographic/population_census_2020.xlsx', ...
    'Sheet', 'Si-Do', 'Range', 'B4:B20')';% 2020 population census

% Load mobility data
load('../pipeline/process_mobility_data/mobility.mat', 'M')

% Load mobility change data
M_change_data = readmatrix('../data/mobility/통신모바일 인구이동량_시도별_합계.xlsx', ...
    'Range', 'B8:R63');% 2020 Feb 1st week

%% Process Mobility Matrix

M = round(M(:, :, 5)');% D\O data going back home matrix

% Get relative mobility change
M_prop = M_change_data(2:end, :) ./ M_change_data(1:end - 1, :);

% Replace diagonal elements
n_loc = length(pop_data);
N0 = M;
for i = 1:n_loc
    N0(i, i) = pop_data(i) - sum(M(:, i));
end

%% Set the Range of State Variables and Parameters

% Initial conditions
E0_truth = zeros(n_loc);
E0_truth(3, 3) = 5;

Ir0_truth = zeros(n_loc);

Iu0_truth = zeros(n_loc);
Iu0_truth(3, 3) = 5;

Q0_truth = zeros(n_loc);

R0_truth = zeros(n_loc);

S0_truth = N0 - E0_truth - Ir0_truth - Iu0_truth;% Omitted Q0_truth and R0_truth;

n_state_type = 6;% Number of state variables (S, E, Ir, Iu, Q, R)

% Parameters
beta_lb = 0.3 * ones(n_loc, 1);% Transmission rate
beta_ub = 1.6 * ones(n_loc, 1);

mu_lb = 0.1;% Relative infectivity of unreported cases
mu_ub = 0.7;

Z_lb = 2;% Infection to infectious
Z_ub = 5;

alpha_lb = 0.1;% Reporting proportion
alpha_ub = 0.75;

Dr_lb = 1;% Infectious to quarantine
Dr_ub = 3;

Du_lb = 2;% Infectious to recovery
Du_ub = 6;

G_lb = 12;% Quarantine to recovery
G_ub = 15;

%% Set Time

t_init = 0;
t_termi = 365;
dt_day = 1 / 3;

%% Generate Test Sets and Time Span

param_lb = [beta_lb; mu_lb; Z_lb; alpha_lb; Dr_lb; Du_lb; G_lb];
param_ub = [beta_ub; mu_ub; Z_ub; alpha_ub; Dr_ub; Du_ub; G_ub];

param_sample = kron(param_lb, ones(1, n_test_param)) ...
    + rand(length(param_lb), n_test_param) .* kron(param_ub - param_lb, ones(1, n_test_param));

beta_grid = ones(n_loc, 1) * param_sample(1, :);
mu_grid = param_sample(n_loc + 1, :);
Z_grid = param_sample(n_loc + 2, :);
alpha_grid = param_sample(n_loc + 3, :);
Dr_grid = param_sample(n_loc + 4, :);
Du_grid = param_sample(n_loc + 5, :);
G_grid = param_sample(n_loc + 6, :);

% Find indices that have reasonable R0
R_0_grid = beta_grid .* (alpha_grid .* Dr_grid + (1 - alpha_grid) .* mu_grid .* Du_grid);
col_bool = (R_0_grid(1, :) > 1) + (R_0_grid(1, :) < 5);
col_idx = find(col_bool == 2);
n_test_param_filtered = length(col_idx);

% Generate tspan
tspan_day = t_init:t_termi;
n_obs = length(tspan_day) - 1;

% Generate tspan including day and night steps
tspan_all = zeros(1, 2 * n_obs + 1);
tspan_all(1) = t_init;
for i = 1:n_obs
    tspan_all(2 * i) = t_init + (i - 1) + dt_day;
    tspan_all(2 * i + 1) = t_init + i;
end

%% Generate Synthetic Data

for i = 1:100

    % Print
    fprintf('Generating %dth data...\n', i)

    beta_truth = beta_grid(1:n_loc, col_idx(i));
    mu_truth = mu_grid(col_idx(i));
    Z_truth = Z_grid(col_idx(i));
    alpha_truth = alpha_grid(col_idx(i));
    Dr_truth = Dr_grid(col_idx(i));
    Du_truth = Du_grid(col_idx(i));
    G_truth = G_grid(col_idx(i));
    R_0_truth = beta_truth * (alpha_truth * Dr_truth + (1 - alpha_truth) * mu_truth * Du_truth);

    % Initialize state variable values
    S0 = S0_truth;
    E0 = E0_truth;
    Ir0 = Ir0_truth;
    Iu0 = Iu0_truth;
    Q0 = Q0_truth;
    R0 = R0_truth;
    N = N0;

    % Initialize variables
    sol_truth = zeros(n_loc, n_loc, n_state_type, n_obs + 1);
    sol_truth(:, :, 1, 1) = S0;
    sol_truth(:, :, 2, 1) = E0;
    sol_truth(:, :, 3, 1) = Ir0;
    sol_truth(:, :, 4, 1) = Iu0;
    sol_truth(:, :, 5, 1) = Q0;
    sol_truth(:, :, 6, 1) = R0;
    truth_data = zeros(n_loc, n_obs);

    for j = 1:n_obs

        %% Adjust Mobility

        if mod(j, 7) == 0

            x = [S0(:); E0(:); Ir0(:); Iu0(:); Q0(:); R0(:)];

            [x_adjusted, N] = adjust_mobility(x, N, n_state_type, M_prop(j / 7, :));

            X = reshape(x_adjusted, n_loc, n_loc, n_state_type);
            S0 = X(:, :, 1);
            E0 = X(:, :, 2);
            Ir0 = X(:, :, 3);
            Iu0 = X(:, :, 4);
            Q0 = X(:, :, 5);
            R0 = X(:, :, 6);
        end

        %% Step Forward

        % Set params
        params_cell = {'S0', S0, false, '$S(0)$'
            'E0', E0, false, '$E(0)$'
            'Ir0', Ir0, false, '$I^r(0)$'
            'Iu0', Iu0, false, '$I^u(0)$'
            'Q0', Q0, false, '$Q(0)$'
            'R0', R0, false, '$R(0)$'
            'n_state_type', n_state_type, false, 'Number of states'

            'N', N, false, '$N$'
            'beta', beta_truth, false, '$\beta$'
            'mu', mu_truth, false, '$\mu$'
            'Z', Z_truth, false, '$Z$'
            'alpha', alpha_truth, false, '$\alpha$'
            'Dr', Dr_truth, false, '$D^r$'
            'Du', Du_truth, false, '$D^u$'
            'G', G_truth, false, '$G$'

            't_init', t_init, false, 'Initial time'
            't_termi', t_termi, false, 'Terminal time'};

        params = cell2struct(params_cell(:, 2:end), params_cell(:, 1), 1);

        tspan = tspan_all(2 * (j - 1) + 1:2 * j + 1);

        % Concatenate initial conditions
        incidence_r = zeros(n_loc, n_loc);% Dummy
        incidence_u = zeros(n_loc, n_loc);% Dummy

        y0 = [S0(:); E0(:); Ir0(:); Iu0(:); Q0(:); R0(:); incidence_r(:); incidence_u(:)];

        % Solve the model
        n_total_state = n_loc * n_loc * n_state_type;
        [sol, incidence_r, incidence_u] = euler_stochastic(@(t, y)seiqr_dayNnight_ens(t, y, params), tspan, y0, n_total_state);

        % Reshape sol
        nt = length(tspan);
        sol_ = reshape(sol, n_loc, n_loc, n_state_type + 2, nt);

        % Compute incidence
        incidence_r_reshaped = squeeze(sum(reshape(incidence_r, n_loc, n_loc, nt - 1), 1));

        n_day = (length(tspan) + 1) / 2;
        sum_operator = kron(eye(n_day - 1), ones(1, 2));
        daily_incidence_r = (sum_operator * incidence_r_reshaped')';        

        truth_data(:, j) = round(daily_incidence_r);

        % Update initial conditions for the next step
        S0 = sol_(:, :, 1, nt);
        E0 = sol_(:, :, 2, nt);
        Ir0 = sol_(:, :, 3, nt);
        Iu0 = sol_(:, :, 4, nt);
        Q0 = sol_(:, :, 5, nt);
        R0 = sol_(:, :, 6, nt);

        sol_truth(:, :, :, j + 1) = sol_(:, :, 1:n_state_type, nt);

    end

    %% Add Observation Error

    obs_error_sample = randn(n_loc, n_obs);% Standard normal

    % Initialize variables
    obs_error_var = zeros(n_loc, n_obs);
    obs_data = zeros(n_loc, n_obs);

    for j = 1:n_obs
        for k = 1:n_loc

            % Prevent zero observation variance
            obs_error_var(k, j) = max(5, (truth_data(k, j) ^ 2) / oev_scaler);

            % Add observation error
            obs_data(k, j) = truth_data(k, j) + obs_error_sample(k, j) * sqrt(obs_error_var(k, j));

            % Remove negative values
            if obs_data(k, j) < 0
                obs_data(k, j) = truth_data(k, j);
            end
        end
    end

    %% Save Results

    save(sprintf('%s/data_%d.mat', result_path, i), ...
        'truth_data', 'obs_data', 'obs_error_var', '*_truth')

end
