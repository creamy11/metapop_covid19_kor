close all
clear

%% Set Paths

addpath 'model' 'solver'

%% Load Data

% Population data
pop_data = readmatrix('../data/demographic/population_census_2020.xlsx', ...
    'Sheet', 'Si-Do', 'Range', 'B4:B20')';% 2020 Census population

% Mobility data
load('../pipeline/process_mobility_data/mobility.mat', 'M')

% Mobility change data
M_change_data = readmatrix('../data/mobility/통신모바일 인구이동량_시도별_합계.xlsx', ...
    'Range', 'B8:R63');% 2020 Feb 1st week

% Incidence data
case_data = readmatrix('../data/incidence/COVID-19_confirmed_230303.xlsx', ...
    'Sheet', '시도별(17개시도+검역)', 'Range', 'C26:S409');% From Feb 8th

case_data(isnan(case_data)) = 0;

%% Process Mobility Matrix

M = round(M(:, :, 5)');% D\O data going back home matrix

% Get relative mobility change
M_prop = M_change_data(2:end, :) ./ M_change_data(1:end - 1, :);

% Replace diagonal elements
n_loc = length(pop_data);% Number of locations
N0 = M;
for i = 1:n_loc
    N0(i, i) = pop_data(i) - sum(M(:, i));
end

%% Smooth the Data

obs_data = movmean(case_data, 7, 1);% 7 day moving average

%% Set State Variable Options

E0_lb = zeros(n_loc, n_loc);
E0_ub = diag(5 * ones(n_loc, 1));
E_is_inflated = false(n_loc);
E_is_reprobed = false(n_loc);

Ir0_lb = zeros(n_loc, n_loc);
Ir0_ub = zeros(n_loc, n_loc);
Ir_is_inflated = false(n_loc);
Ir_is_reprobed = false(n_loc);

Iu0_lb = zeros(n_loc, n_loc);
Iu0_ub = diag(5 * ones(n_loc, 1));
Iu_is_inflated = false(n_loc);
Iu_is_reprobed = false(n_loc);

S0_lb = round(0.9 * N0);
S0_ub = N0 - E0_ub - Ir0_ub - Iu0_ub;% Make sure the sum does not exceed N
S_is_inflated = false(n_loc);
S_is_reprobed = false(n_loc);

Q0_lb = zeros(n_loc, n_loc);
Q0_ub = zeros(n_loc, n_loc);
Q_is_inflated = false(n_loc);
Q_is_reprobed = false(n_loc);

R0_lb = zeros(n_loc, n_loc);% Will be determined by the other state variables
R0_ub = zeros(n_loc, n_loc);
R_is_inflated = false(n_loc);% Always false
R_is_reprobed = false(n_loc);

n_state_type = 6;% Number of state variables (S, E, Ir, Iu, Q, R)

%% Observation

obs0_lb = zeros(n_loc, 1);% Also, this will be determined by state variables
obs0_ub = pop_data';
obs_is_inflated = true(n_loc, 1);
obs_is_reprobed = false(n_loc, 1);

%% Unreported Cases

incidence_u0_lb = zeros(n_loc, 1);
incidence_u0_ub = pop_data';
incidence_u_is_inflated = true(n_loc, 1);
incidence_u_is_reprobed = false(n_loc, 1);

%% Set Parameter Options

beta_lb = zeros(n_loc, 1);% Transmission rate
beta_ub = 5 * ones(n_loc, 1);
beta0_lb = 0.3 * ones(n_loc, 1);% Prior range
beta0_ub = 1.6 * ones(n_loc, 1);
beta_is_local = true;
if beta_is_local
    beta_is_inflated = true(n_loc, 1);
    beta_is_reprobed = true(n_loc, 1);
else
    beta_is_inflated = [true; false(n_loc - 1, 1)];
    beta_is_reprobed = [true; false(n_loc - 1, 1)];
end

mu_lb = 0;% Relative infectivity of unreported case
mu_ub = 1;
mu0_lb = 0.1;
mu0_ub = 0.7;
mu_is_inflated = false;
mu_is_reprobed = false;

Z_lb = 1;% Infection to infectious
Z_ub = 6;
Z0_lb = 2;
Z0_ub = 5;
Z_is_inflated = false;
Z_is_reprobed = false;

alpha_lb = eps;% Reporting proportion
alpha_ub = 1;
alpha0_lb = 0.1;
alpha0_ub = 0.75;
alpha_is_inflated = false;
alpha_is_reprobed = true;

Dr_lb = 0.01;% Infectious to quarantine
Dr_ub = 4;
Dr0_lb = 1;
Dr0_ub = 3;
Dr_is_inflated = false;
Dr_is_reprobed = false;

Du_lb = 1;% Infectious to recovery
Du_ub = 7;
Du0_lb = 2;
Du0_ub = 6;
Du_is_inflated = false;
Du_is_reprobed = false;

G_lb = 10;% Quarantine to recovery
G_ub = 17;
G0_lb = 12;
G0_ub = 15;
G_is_inflated = false;
G_is_reprobed = false;

n_parameter = n_loc + 6;

%% Set EAKF Options

n_ens = 300;
n_real = 10;
rng(0)

% Time
date_init = datetime('7-Feb-2020');% 1 day before initial observation
date_termi = datetime('25-Feb-2021');% Before the vaccination
dt_obs = 1;
dt_day = 1 / 3;

% Inflation
is_inflation = true;
prop_inflation = 0.01;

% Re-probing
is_reprobing = true;
prop_reprobing = 0.01;

result_path = '../result/estimation';
suffix = '_fix';

%% Generate Secondary Variables

n_state_total = n_state_type * n_loc * n_loc;
n_var = n_state_total + 2 * n_loc + n_parameter;

% Concatenate lower and upper bounds
x_lb = [zeros(n_state_total + 2 * n_loc, 1); beta_lb; mu_lb; Z_lb; alpha_lb;
    Dr_lb; Du_lb; G_lb];
x0_lb = [S0_lb(:); E0_lb(:); Ir0_lb(:); Iu0_lb(:); Q0_lb(:); R0_lb(:); obs0_lb;
    incidence_u0_lb; beta0_lb; mu0_lb; Z0_lb; alpha0_lb; Dr0_lb; Du0_lb; G0_lb];
x0_ub = [S0_ub(:); E0_ub(:); Ir0_ub(:); Iu0_ub(:); Q0_ub(:); R0_ub(:); obs0_ub;
    incidence_u0_ub; beta0_ub; mu0_ub; Z0_ub; alpha0_ub; Dr0_ub; Du0_ub; G0_ub];

% Variable indices
S_idx = 1:n_loc * n_loc;
E_idx = n_loc * n_loc + 1:2 * n_loc * n_loc;
Ir_idx = 2 * n_loc * n_loc + 1:3 * n_loc * n_loc;
Iu_idx = 3 * n_loc * n_loc + 1:4 * n_loc * n_loc;
Q_idx = 4 * n_loc * n_loc + 1:5 * n_loc * n_loc;
R_idx = 5 * n_loc * n_loc + 1:6 * n_loc * n_loc;

obs_idx = 6 * n_loc * n_loc + (1:n_loc);
incidence_r_idx = 6 * n_loc * n_loc + n_loc + (1:n_loc);

beta_idx = 6 * n_loc * n_loc + 2 * n_loc + (1:n_loc);
mu_idx = n_var - 5;
Z_idx = n_var - 4;
alpha_idx = n_var - 3;
Dr_idx = n_var - 2;
Du_idx = n_var - 1;
G_idx = n_var;

% Generate tspan
date_span = date_init:caldays(dt_obs):date_termi;
n_obs = length(date_span) - 1;
tspan_obs = 0:dt_obs:n_obs;

tspan4integ = cell(n_obs, 1);
for i = 1:n_obs
    dt_init = tspan_obs(i);
    dt_termi = tspan_obs(i + 1);
    nt_day = length(dt_init:dt_termi);
    tspan = zeros(1, 2 * nt_day - 1);
    tspan(1) = dt_init;
    for j = 1:nt_day - 1
        tspan(2 * j) = dt_init + (j - 1) + dt_day;
        tspan(2 * j + 1) = dt_init + j;
    end
    tspan4integ{i, 1} = tspan;
end

% For averaging Kalman gain
if beta_is_local
    n_loc_var = (n_state_type - 1) * n_loc + 3;% Number of local variables for each location j
    global_idx = alpha_idx;
    %     global_idx = [mu_idx, Z_idx, alpha_idx, Dr_idx, Du_idx, G_idx];
else
    n_loc_var = (n_state_type - 1) * n_loc + 2;
    global_idx = [beta_idx(1), mu_idx, Z_idx, alpha_idx, Dr_idx, Du_idx, G_idx];
end
n_global_parameter = length(global_idx);

% Concatenate is_inflated
is_inflated = [S_is_inflated(:); E_is_inflated(:); Ir_is_inflated(:)
    Iu_is_inflated(:); Q_is_inflated(:); R_is_inflated(:); obs_is_inflated;
    incidence_u_is_inflated; beta_is_inflated; mu_is_inflated; Z_is_inflated;
    alpha_is_inflated; Dr_is_inflated; Du_is_inflated; G_is_inflated];

inflation_idx = find(is_inflated);

% Concatenate is_reprobed
is_reprobed = [S_is_reprobed(:); E_is_reprobed(:); Ir_is_reprobed(:)
    Iu_is_reprobed(:); Q_is_reprobed(:); R_is_reprobed(:); obs_is_reprobed;
    incidence_u_is_reprobed; beta_is_reprobed; mu_is_reprobed; Z_is_reprobed;
    alpha_is_reprobed; Dr_is_reprobed; Du_is_reprobed; G_is_reprobed];

reprobing_idx = find(is_reprobed);

n_reprobing = round(n_ens * prop_reprobing);% Number of ensembles re-probed
n_var_reprobed = sum(is_reprobed);% Number of variables re-probed

%% Main Loop

% Initialize
x_post_state = zeros(n_loc * n_state_type, n_ens, n_obs, n_real);
x_post_obs = zeros(n_loc, n_ens, n_obs, n_real);
x_post_incidence_u = zeros(n_loc, n_ens, n_obs, n_real);
x_post_parameter = zeros(n_parameter, n_ens, n_obs, n_real);

t_start = tic;

% Loop over realizations
for ii = 1:n_real

    % Print realization number
    fprintf('Simulation %d is running...\n', ii)

    %% Generate Initial Ensemble Members

    % Sample initial ensemble memebers from Uniform distribution
    x = kron(x0_lb, ones(1, n_ens)) + rand(n_var, n_ens) .* kron(x0_ub - x0_lb, ones(1, n_ens));
    x(1:n_state_total, :) = round(x(1:n_state_total, :));% Rounding off state variable values
    x(R_idx, :) = kron(N0(:), ones(1, n_ens)) - x(S_idx, :) - x(E_idx, :) ...
        - x(Ir_idx, :) - x(Iu_idx, :) - x(Q_idx, :);% Adjust R size

    % If beta is a global parameter
    if ~beta_is_local
        x(beta_idx, :) = ones(n_loc, 1) * x(beta_idx(1), :);
    end

    %% EAKF

    N = N0;
    obs_model = zeros(n_loc, n_ens, n_obs);
    x_prior = zeros(n_var, n_ens, n_obs);
    x_post = zeros(n_var, n_ens, n_obs);

    % Loop over observations
    for i = 1:n_obs

        %% Adjust Mobility

        % Re-distribute population weekly
        if mod(i, 7) == 0% Every Monday
            [x(1:n_state_total, :), N] = adjust_mobility(x(1:n_state_total, :), N, n_state_type, M_prop(i / 7, :));
        end

        %% Step Forward

        S0 = reshape(x(S_idx, :), n_loc, n_loc, n_ens);
        E0 = reshape(x(E_idx, :), n_loc, n_loc, n_ens);
        Ir0 = reshape(x(Ir_idx, :), n_loc, n_loc, n_ens);
        Iu0 = reshape(x(Iu_idx, :), n_loc, n_loc, n_ens);
        Q0 = reshape(x(Q_idx, :), n_loc, n_loc, n_ens);
        R0 = reshape(x(R_idx, :), n_loc, n_loc, n_ens);

        beta = x(beta_idx, :);
        mu = x(mu_idx, :);
        Z = x(Z_idx, :);
        alpha = x(alpha_idx, :);
        Dr = x(Dr_idx, :);
        Du = x(Du_idx, :);
        G = x(G_idx, :);

        x_ub = [kron(ones(n_state_type, 1), N(:)); obs0_ub; incidence_u0_ub;
            beta_ub; mu_ub; Z_ub; alpha_ub; Dr_ub; Du_ub; G_ub];% Depend on N

        params_cell = {'S0', S0, false, '$S(0)$'
            'E0', E0, false, '$E(0)$'
            'Ir0', Ir0, false, '$I^r(0)$'
            'Iu0', Iu0, false, '$I^u(0)$'
            'Q0', Q0, false, '$Q(0)$'
            'R0', R0, false, '$R(0)$'
            'n_state_type', n_state_type, false, 'Number of states'

            'N', N, false, '$N$'
            'beta', beta, false, '$\beta$'
            'mu', mu, false, '$\mu$'
            'Z', Z, false, '$Z$'
            'alpha', alpha, false, '$\alpha$'
            'Dr', Dr, false, '$D^r$'
            'Du', Du, false, '$D^u$'
            'G', G, false, '$G$'

            'x_lb', x_lb, false, 'Lower bound'
            'x_ub', x_ub, false, 'Upper bound'};

        params = cell2struct(params_cell(:, 2:end), params_cell(:, 1), 1);

        tspan = tspan4integ{i, 1};

        % Concatenate initial conditions
        y0 = zeros(n_loc * n_loc * (n_state_type + 2), n_ens);% Including space for observation
        y0(1:n_state_total, :) = x(1:n_state_total, :);

        % Solve the model
        [sol, incidence_r, incidence_u] = euler_stochastic(@(t, y)seiqr_dayNnight_ens(t, y, params), tspan, y0, n_state_total);

        % Reshape sol
        nt = length(tspan);
        sol_ = reshape(sol, n_loc, n_loc, n_state_type + 2, n_ens, nt);

        % Compute incidence
        incidence_r_reshaped = squeeze(sum(reshape(incidence_r, n_loc, n_loc, n_ens, nt - 1), 1));
        incidence_u_reshaped = squeeze(sum(reshape(incidence_u, n_loc, n_loc, n_ens, nt - 1), 1));

        n_day = (length(tspan) + 1) / 2;
        sum_operator = kron(eye(n_day - 1), ones(2, 1));
        daily_incidence_r = zeros(n_loc, n_ens, n_day - 1);
        daily_incidence_u = zeros(n_loc, n_ens, n_day - 1);
        for j = 1:n_ens
            daily_incidence_r(:, j, :) = squeeze(incidence_r_reshaped(:, j, :)) * sum_operator;
            daily_incidence_u(:, j, :) = squeeze(incidence_u_reshaped(:, j, :)) * sum_operator;
        end

        obs_model(:, :, i) = round(daily_incidence_r(:, :, n_day - 1));% Leave this for the deterministic model

        %% Update Initial Conditions for the Next Step

        sol_end = sol_(:, :, 1:n_state_type, :, end);
        x_next = x;
        x_next(1:n_state_total, :) = reshape(sol_end, n_state_total, n_ens);
        x_next(n_state_total + (1:n_loc), :) = obs_model(:, :, i);
        x_next(n_state_total + n_loc + (1:n_loc), :) = round(daily_incidence_u(:, :, n_day - 1));

        %% Inflation

        % Inflation
        if is_inflation
            x_next(inflation_idx, :) = mean(x_next(inflation_idx, :), 2) ...
                + (1 + prop_inflation) .* (x_next(inflation_idx, :) - mean(x_next(inflation_idx, :), 2));
            x_next = check_bound(params, x_next);
        end

        %% Re-probing

        if is_reprobing
            if mod(i, 7) == 0% Weekly re-probing

                % Sample ensemble indices for re-probing
                reprobing_ens_idx = datasample(1:n_ens, n_reprobing);

                % Sample new ensemble members
                x_next(reprobing_idx, reprobing_ens_idx) = kron(x_lb(reprobing_idx), ones(1, n_reprobing)) ...
                    + rand(n_var_reprobed, n_reprobing) .* kron(x_ub(reprobing_idx) - x_lb(reprobing_idx), ones(1, n_reprobing));
                x_next = check_bound(params, x_next);
            end
        end

        %% Ensemble Adjustment Kalman Filter

        x_prior(:, :, i) = x_next;
        obs_i = obs_model(:, :, i);

        dx_global = zeros(n_global_parameter, n_ens, n_loc);

        % Loop over location
        for j = 1:n_loc

            obs_var = max(5, (obs_data(i, j) ^ 2) / 5);

            prior_var = var(obs_i(j, :));

            post_var = prior_var * obs_var / (prior_var + obs_var);

            if prior_var == 0
                prior_var = 1e-3;
                post_var = 1e-3;
            end

            prior_mean = mean(obs_i(j, :));
            post_mean = post_var * (prior_mean / prior_var + obs_data(i, j) / obs_var);

            % Adjustment
            alpha = sqrt(obs_var / (obs_var + prior_var));
            dy = post_mean + alpha * (obs_i(j, :) - prior_mean) - obs_i(j, :);

            % Indices for location j's state variables and parameters and
            % global parameters
            if beta_is_local
                j_loc_idx = [(j - 1) * n_loc + 1:j * n_loc, ...% Local S
                    n_loc * n_loc + (j - 1) * n_loc + 1:n_loc * n_loc + j * n_loc, ...% Local E
                    2 * n_loc * n_loc + (j - 1) * n_loc + 1:2 * n_loc * n_loc + j * n_loc, ...% Local Ir
                    3 * n_loc * n_loc + (j - 1) * n_loc + 1:3 * n_loc * n_loc + j * n_loc, ...% Local Iu
                    4 * n_loc * n_loc + (j - 1) * n_loc + 1:4 * n_loc * n_loc + j * n_loc, ...% Local Q
                    6 * n_loc * n_loc + j, ...% Local daily new reported cases
                    6 * n_loc * n_loc + n_loc + j, ...% Local daily new unreported cases
                    beta_idx(j)];% Parameters
            else
                j_loc_idx = [(j - 1) * n_loc + 1:j * n_loc, ...% Local S
                    n_loc * n_loc + (j - 1) * n_loc + 1:n_loc * n_loc + j * n_loc, ...% Local E
                    2 * n_loc * n_loc + (j - 1) * n_loc + 1:2 * n_loc * n_loc + j * n_loc, ...% Local Ir
                    3 * n_loc * n_loc + (j - 1) * n_loc + 1:3 * n_loc * n_loc + j * n_loc, ...% Local Iu
                    4 * n_loc * n_loc + (j - 1) * n_loc + 1:4 * n_loc * n_loc + j * n_loc, ...% Local Q
                    6 * n_loc * n_loc + j, ...% Local daily new reported cases
                    6 * n_loc * n_loc + n_loc + j];% Local daily new unreported cases
            end

            % Compute covariance
            rr_local = zeros(n_loc_var, 1);
            for k = 1:n_loc_var
                A = cov(x_next(j_loc_idx(k), :), obs_i(j, :));
                rr_local(k) = A(2, 1) / prior_var;
            end

            rr_global = zeros(n_global_parameter, 1);
            for k = 1:n_global_parameter
                A = cov(x_next(global_idx(k), :), obs_i(j, :));
                rr_global(k) = A(2, 1) / prior_var;
            end
            dx_global(:, :, j) = rr_global * dy;

            % Update x (local variables)
            dx = rr_local * dy;
            x_next(j_loc_idx, :) = x_next(j_loc_idx, :) + dx;

            % Remove aphysicality
            x_next = check_bound(params, x_next);

        end

        % Update x (global variables)
        x_next(global_idx, :) = x_next(global_idx, :) + mean(dx_global, 3);
        x_next = check_bound(params, x_next);

        if ~beta_is_local
            x_next(beta_idx, :) = ones(n_loc, 1) * x_next(beta_idx(1), :);
        end

        x_post(:, :, i) = x_next;
        x = x_next;

    end

    % Save posterior distribution of state variables
    X_post_state = reshape(x_post(1:n_state_total, :, :), n_loc, n_loc, n_state_type, n_ens, n_obs);
    x_post_state(:, :, :, ii) = reshape(sum(X_post_state, 1), n_loc * n_state_type, n_ens, n_obs);

    % Save posterior distribution of observations
    x_post_obs(:, :, :, ii) = x_post(n_state_total + (1:n_loc), :, :);

    % Save posterior distribution of unreported daily cases
    x_post_incidence_u(:, :, :, ii) = x_post(n_state_total + n_loc + (1:n_loc), :, :);

    % Save posterior distribution of parameters
    x_post_parameter(:, :, :, ii) = x_post(n_state_total + 2 * n_loc + (1:n_parameter), :, :);

end

fprintf('Saving results with suffix name: %s\n', suffix)
save(sprintf('%s/estimate%s.mat', result_path, suffix), 'obs_data', 'case_data', ...
    '*_lb', '*_ub', 'date_span', 'n_state_total', 'x_post_state', 'x_post_obs', ...
    'x_post_incidence_u', 'x_post_parameter', '-v7.3')

t_end = toc(t_start);
fprintf('Elapsed time is %.2f seconds.\n', t_end)