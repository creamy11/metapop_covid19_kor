function [sol, incidence_r, incidence_u] = euler_stochastic(odefun, tspan, y0, n_state_total)

nt = length(tspan);
n_y0 = size(y0, 1);
n_ens = size(y0, 2);

n_loc_squared = (n_y0 - n_state_total) / 2;% 17 X 17 = 289

incidence_r_idx = n_state_total + (1:n_loc_squared);
incidence_u_idx = n_state_total + n_loc_squared + (1:n_loc_squared);

sol = zeros(n_y0, n_ens, nt);
sol(1:n_state_total, :, 1) = y0(1:n_state_total, :);

incidence_r = zeros(n_loc_squared, n_ens, nt - 1);
incidence_u = zeros(n_loc_squared, n_ens, nt - 1);

for i = 1:nt - 1

    t_current = tspan(i);

    dydt = odefun(t_current, sol(:, :, i));

    sol(:, :, i + 1) = max(sol(:, :, i) + dydt, 0);

    incidence_r(:, :, i) = max(dydt(incidence_r_idx, :), 0);
    incidence_u(:, :, i) = max(dydt(incidence_u_idx, :), 0);

end

end