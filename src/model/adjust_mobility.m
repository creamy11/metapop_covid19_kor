function [x_adjusted, N_adjusted] = adjust_mobility(x_state, N, n_state_type, M_prop)

n_loc = size(N, 1);
n_ens = size(x_state, 2);

M_prop_ext = ones(n_loc, 1) * M_prop;
M_prop_ext = M_prop_ext - diag(diag(M_prop_ext)) + diag(ones(n_loc, 1));

% Adjust state variables
X = reshape(x_state, n_loc, n_loc, n_state_type, n_ens);
for i = 1:n_state_type

    X_i = round(squeeze(X(:, :, i, :)));% Rounding off here is to reduce floating point errors
    X_adjusted = round(M_prop_ext .* X_i);% It does not necessarily have to be integer.
    X_diff = X_i - X_adjusted;

    for j = 1:n_loc
        X_adjusted(j, j, :) = X_i(j, j, :) + sum(X_diff(:, j, :), 1);
    end

    X(:, :, i, :) = X_adjusted;
end

% Adjust N
N_adjusted = round(M_prop_ext .* N);
N_diff = N - N_adjusted;
for j = 1:n_loc
    N_adjusted(j, j) = N(j, j) + sum(N_diff(:, j));
end

% Process rounding off errors (re-adjusted R and S compartments)
X(:, :, n_state_type, :) = X(:, :, n_state_type, :) + N_adjusted - sum(X, 3);
R = squeeze(X(:, :, n_state_type, :));

R_neg = R;
R_neg(R_neg > 0) = 0;
X(:, :, 1, :) = squeeze(X(:, :, 1, :)) + R_neg;

R_pos = R;
R_pos(R_pos < 0) = 0;
X(:, :, n_state_type, :) = R_pos;

% Check errors
if sum(sum(N, 1) == sum(N_adjusted, 1)) < n_loc
    pause
end

if sum(sum(X, 3) == N_adjusted, 'all') < n_loc * n_loc * n_ens
    pause
end

x_adjusted = reshape(X, n_loc * n_loc * n_state_type, n_ens);

end