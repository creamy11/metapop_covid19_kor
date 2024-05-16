function x_ = check_bound(params, x)
%% Assign Parameters

n_state_type = params(1).n_state_type;
N_ = params(1).N;
x_lb_ = params(1).x_lb;
x_ub_ = params(1).x_ub;

n_loc = size(N_, 1);
n_ens = size(x, 2);
n_var = size(x, 1);

n_state_total = n_loc * n_loc * n_state_type;

%% Remove Negativity

x_ = x;

% Replace negative values to zero
x_(x_ < 0) = 0;

%% Check State Variables

% Adjusting R population size
N_copied = N_(:) * ones(1, n_ens);
R_idx = (n_state_type - 1) * n_loc * n_loc + 1:n_state_type * n_loc * n_loc;
X = reshape(x_(1:n_state_total, :), n_loc * n_loc, n_state_type, n_ens);
x_(R_idx, :) = N_copied - squeeze(sum(X(:, 1:n_state_type - 1, :), 2));

% If S+E+Ir+Iu is greater than N, replace R with 0 and scale down the values
% of the other state variables proportionally to their magnitude.
x_(x_ < 0) = 0;

X = reshape(x_(1:n_state_total, :), n_loc * n_loc, n_state_type, n_ens);

X_sum = squeeze(sum(X, 2));

if sum(X_sum - N_copied, 'all')
    for i = 1:n_state_type
        X(:, i, :) = squeeze(X(:, i, :)) ./ X_sum .* N_copied;
    end
end
x_(1:n_state_total, :) = reshape(X, n_loc * n_loc * n_state_type, n_ens);

x_(isnan(x_)) = 0;% To avoid nan

%% Check Observations

incidence_r_idx = n_state_total + (1:n_loc);
total_pop = sum(N_, 1);
if sum(x_(incidence_r_idx, :) > total_pop', 'all')
    [row, col] = find(x_(incidence_r_idx, :) > total_pop');
    for i = 1:length(row)
        x_(incidence_r_idx(row(i)), col(i)) = total_pop(row(i));
    end
end

%% Check Unreported Cases

incidence_u_idx = n_state_total + n_loc + (1:n_loc);
total_pop = sum(N_, 1);
if sum(x_(incidence_u_idx, :) > total_pop', 'all')
    [row, col] = find(x_(incidence_u_idx, :) > total_pop');
    for i = 1:length(row)
        x_(incidence_u_idx(row(i)), col(i)) = total_pop(row(i));
    end
end

%% Check Parameters

% Boundary conditions
parameter_idx = n_state_total + 2 * n_loc + 1:n_var;

% Resample parameters
for i = parameter_idx
    out_index = find(x_(i, :) < x_lb_(i) | x_(i, :) > x_ub_(i));
    in_index = setdiff(1:n_ens, out_index);
    if isempty(in_index)% If all ensemble members are out of the boundary
        x_(i, :) = x_lb_(i) + rand(1, n_ens) * (x_ub_(i) - x_lb_(i));% re-initialize
    else
        x_(i, out_index) = datasample(x_(i, in_index), length(out_index));
    end
end

end