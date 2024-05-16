function dydt = seiqr_dayNnight_ens(t, y, params)
%% Assign Parameters

n_state_type = params(1).n_state_type;

N_ = params(1).N;
beta_ = params(1).beta;
mu_ = params(1).mu;
Z_ = params(1).Z;
alpha_ = params(1).alpha;
Dr_ = params(1).Dr;
Du_ = params(1).Du;
G_ = params(1).G;

%% Reshape Variables

n_loc = size(N_, 1);
n_ens = size(y, 2);

yy = reshape(y, n_loc * n_loc, n_state_type + 2, n_ens);
S = squeeze(yy(:, 1, :));
E = squeeze(yy(:, 2, :));
Ir = squeeze(yy(:, 3, :));
Iu = squeeze(yy(:, 4, :));
Q = squeeze(yy(:, 5, :));

Ir_ = reshape(Ir, n_loc, n_loc, n_ens);
Iu_ = reshape(Iu, n_loc, n_loc, n_ens);
Q_ = reshape(Q, n_loc, n_loc, n_ens);

Ir_org = squeeze(sum(Ir_, 1));% Ir origin

beta_day = kron(ones(n_loc, 1), beta_);
beta_night = kron(beta_, ones(n_loc, 1));

%%

dyydt = zeros(size(yy));

if mod(t, 1) == 0

    %% Day 

    dt = 1 / 3;

    Ir_des = squeeze(sum(Ir_, 2));% Ir destination
    Ir_des_ext = kron(ones(n_loc, 1), Ir_des);

    Iu_des = squeeze(sum(Iu_, 2));% Iu destination
    Iu_des_ext = kron(ones(n_loc, 1), Iu_des);

    if n_ens == 1
        N_day = sum(N_, 2) + sum(Q_, 1)' - sum(Q_, 2);
        N_day_ext = kron(ones(n_loc, 1), N_day);
    else
        N_day = sum(N_, 2) * ones(1, n_ens) + squeeze(sum(Q_, 1)) - squeeze(sum(Q_, 2));
        N_day_ext = kron(ones(n_loc, 1), N_day);
    end

    S2E = poissrnd(beta_day .* S .* Ir_des_ext ./ N_day_ext * dt) ...
        + poissrnd(mu_ .* beta_day .* S .* Iu_des_ext ./ N_day_ext * dt);

else

    %% Night

    dt = 2 / 3;

    if n_ens == 1
        Ir_org_ext = kron(Ir_org', ones(n_loc, 1));

        Iu_org = sum(Iu_, 1);% Iu origin
        Iu_org_ext = kron(Iu_org', ones(n_loc, 1));

        N_night = sum(N_, 1);
        N_night_ext = kron(N_night', ones(n_loc, 1));
    else
        Ir_org_ext = kron(Ir_org, ones(n_loc, 1));

        Iu_org = squeeze(sum(Iu_, 1));% Iu origin
        Iu_org_ext = kron(Iu_org, ones(n_loc, 1));

        N_night = sum(N_, 1)' * ones(1, n_ens);
        N_night_ext = kron(N_night, ones(n_loc, 1));
    end

    S2E = poissrnd(beta_night .* S .* Ir_org_ext ./ N_night_ext * dt) ...
        + poissrnd(mu_ .* beta_night .* S .* Iu_org_ext ./ N_night_ext * dt);

end

E2Ir = poissrnd(alpha_ .* E ./ Z_ * dt);
E2Iu = poissrnd((1 - alpha_) .* E ./ Z_ * dt);
Ir2Q = poissrnd(Ir ./ Dr_ * dt);
Q2R = poissrnd(Q ./ G_ * dt);
Iu2R = poissrnd(Iu ./ Du_ * dt);

dyydt(:, 1, :) = -S2E;% S
dyydt(:, 2, :) = S2E - E2Ir - E2Iu;% E
dyydt(:, 3, :) = E2Ir - Ir2Q;% Ir
dyydt(:, 4, :) = E2Iu - Iu2R;% Iu
dyydt(:, 5, :) = Ir2Q - Q2R;% Q
dyydt(:, 6, :) = Q2R + Iu2R;% R
dyydt(:, 7, :) = Ir2Q;
dyydt(:, 8, :) = E2Iu;

dydt = reshape(dyydt, n_loc * n_loc * (n_state_type + 2), n_ens);

end