function [freq_SPS, freq_GF, freq_OF] = coverage(Params)
% COVERAGE   Compute coverage frequency SPS, GF, OF

belong_SPS = zeros(1,Params.N_check);
belong_GF = zeros(1,Params.N_check);
belong_OF = zeros(1,Params.N_check);

F = randn(Params.n); % state matrix
F = F/(norm(F)+Params.stab); % stabilized F (for numerical reasons)
% F used in Fig. 1 & Table 1
% F = [0.1540   -0.1998   -0.2220   -0.2658
%      -0.0513    0.1222   -0.0549    0.2303
%       0.7746    0.4134    0.3559   -0.0177
%       0.1532   -0.0254    0.0071   -0.5559];

if Params.IV_case == 2 % IV case
    Params.T = Params.T+Params.T_est;
end

for s = 1:Params.N_check

    x = zeros(Params.n,Params.T+2); % state sequence
    u = randn(Params.n,Params.T+1); % input sequence

    switch Params.noise_case % noise case

        case 1
            w = Params.sigma_nom.*randn(Params.n,Params.T+1);

        case 2
            sigma = zeros(Params.n,Params.T+1);
            for ii = 1:Params.T+1
                for jj = 1:Params.n
                    if rand < Params.prob_mix2
                        sigma(jj,ii) = sqrt((Params.sigma_nom^2-(1-Params.prob_mix2)*Params.sigma_mix2^2)/Params.prob_mix2);
                    else
                        sigma(jj,ii) = Params.sigma_mix2;
                    end
                end
            end
            w = sigma.*randn(Params.n,Params.T+1);

        case 3
            for ii = 1:Params.T+1
                if rand < Params.prob_mix3
                    sigma(ii) = Params.sigma_mix3;
                else
                    sigma(ii) = Params.sigma_nom;
                end
            end
            w = sigma.*randn(Params.n,Params.T+1);

        otherwise
            error('noise_case value not accepted')
    end

    y = zeros(Params.n,Params.T+1);
    % collect data
    Omega = [];
    for i = 1:Params.T+1
        x(:,i+1) = F*x(:,i) + u(:,i) + w(:,i);
        y(:,i) = x(:,i+1)-u(:,i);
        Omega = [Omega; kron(x(:,i)',eye(Params.n))];
    end

    switch Params.IV_case % IV case
        case 1
            % linear system: y = Omega*vec(F) + w
            Omega = Omega(Params.n+1:end,:);
            y = vec(y);
            y = y(Params.n+1:end);
            theta_hat_LS = pinv(Omega)*y;
            F_hat_LS = reshape(theta_hat_LS,[Params.n,Params.n]);

            Omega_IV = [];
            for i = 1:Params.T
                Omega_IV = [Omega_IV; kron(u(:,i)',eye(Params.n))];
            end
        case 2

            % linear system to estimate F for IV
            Omega_est = Omega(1:Params.n*Params.T_est,:);
            y_est = vec(y(:,1:Params.T_est));
            theta_hat_LS_est = pinv(Omega_est)*y_est;
            F_hat_LS_est = reshape(theta_hat_LS_est,[Params.n,Params.n]);

            % linear system to build confidence regions
            Omega = Omega(Params.n*Params.T_est+1:end-Params.n,:);
            y = vec(y(:,Params.T_est+1:end-1));
            theta_hat_LS = pinv(Omega)*y;
            F_hat_LS = reshape(theta_hat_LS,[Params.n,Params.n]);

            % collect data noise-free to define IV
            x_free = zeros(Params.n,Params.T-Params.T_est+1);
            x_free(:,1) = y_est((end-Params.n+1):end)+u(:,Params.T_est);
            Omega_IV = [];
            for i = 1:Params.T-Params.T_est
                x_free(:,i+1) = F_hat_LS_est*x_free(:,i) + u(:,i+Params.T_est);
                zz(:,i) = x_free(:,i+1)-u(:,i);
                Omega_IV = [Omega_IV; kron(x_free(:,i)',eye(Params.n))];
            end
            Omega_IV = Omega_IV(1:end,:);

        otherwise
            error('IV_case value not accepted')
    end

    alpha = zeros(Params.r-1,Params.N);
    D = zeros(Params.N,Params.N,Params.r-1);
    Q = zeros(Params.n^2,Params.n^2,Params.r-1);
    P = zeros(Params.n^2,Params.r-1);
    A = zeros(Params.n^2,Params.n^2,Params.r-1);
    b = zeros(Params.n^2,Params.r-1);
    b_bar = zeros(Params.n^2,Params.r-1);
    c = zeros(1,Params.r-1);
    c_bar = zeros(1,Params.r-1);
    H = (Omega_IV'*Omega_IV)/Params.N;
    V = (Omega_IV'*Omega)/Params.N;

    theta_hat_SPS = (inv(V)*Omega_IV'*y)/Params.N;

    belong_SPS_tmp = 0;

    % SPS confidence region
    for i = 1:Params.r-1

        alpha(i,:) = (rand(1,Params.N)<.5)*2 - 1; % vector of random signs
        D(:,:,i) = diag(alpha(i,:));
        Q(:,:,i) = (Omega_IV'*D(:,:,i)*Omega)/Params.N;
        P(:,i) = Omega_IV'*D(:,:,i)*y/Params.N;
        A(:,:,i) = V'*inv(H)*V-Q(:,:,i)'*inv(H)*Q(:,:,i);
        b(:,i) = Q(:,:,i)'*inv(H)*P(:,i)-V'*inv(H)*V*theta_hat_SPS;
        b_bar(:,i) = -inv(A(:,:,i))*b(:,i);
        c(i) = theta_hat_SPS'*V'*inv(H)*V*theta_hat_SPS-P(:,i)'*inv(H)*P(:,i);
        c_bar(i) = -c(i)+b(:,i)'*inv(A(:,:,i))*b(:,i);

        if ((vec(F)-b_bar(:,i))'*A(:,:,i)*(vec(F)-b_bar(:,i))<=c_bar(i))
            belong_SPS_tmp = 1;
        end

    end

    belong_SPS(s) = belong_SPS_tmp;
    belong_GF(s) = ((vec(F)-theta_hat_LS)'*(Omega'*Omega)*(vec(F)-theta_hat_LS)<=Params.sigma_nom^2*chi2inv(1-Params.q/Params.r,Params.N)-norm(y-Omega*theta_hat_LS)^2);
    belong_OF(s) = ((vec(F)-theta_hat_LS)'*(Omega'*Omega)*(vec(F)-theta_hat_LS)<=Params.sigma_nom^2*chi2inv(1-Params.q/Params.r,Params.n^2));

end

freq_SPS = length(belong_SPS(belong_SPS~=0))/Params.N_check;
freq_GF = length(belong_GF(belong_GF~=0))/Params.N_check;
freq_OF = length(belong_OF(belong_OF~=0))/Params.N_check;

end