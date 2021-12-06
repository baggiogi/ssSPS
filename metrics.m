function [flag, theta, theta_hat_SPS, theta_hat_LS, box_bounds_SPS, box_bounds_GF, box_bounds_OF, dist_SPS, dist_GF, dist_OF, radius_SPS, radius_GF, radius_OF] = metrics(Params)
% METRICS   Compute performance metrics SPS, GF, OF

F = randn(Params.n); % state matrix
F = F/(norm(F)+Params.stab); % stabilized F (for numerical reasons)
% F used in Fig. 1 & Table 1
% F = [0.1540   -0.1998   -0.2220   -0.2658
%      -0.0513    0.1222   -0.0549    0.2303
%       0.7746    0.4134    0.3559   -0.0177
%       0.1532   -0.0254    0.0071   -0.5559];

theta = vec(F);

if Params.IV_case == 2 % IV case
    Params.T = Params.T+Params.T_est;
end

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
theta_hat_LS = pinv(Omega)*y;

Tsel=[];

upper_bound = zeros(Params.n^2,Params.r-1);
lower_bound = zeros(Params.n^2,Params.r-1);

flag = 0;
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

    if any(eig(A(:,:,i))<=0)
        flag = 1;
        % bounds SPS
        upper_bound(:,i) = NaN;
        lower_bound(:,i) = NaN;

    else
        % bounds SPS
        upper_bound(:,i) = b_bar(:,i)+sqrt(c_bar(:,i)*diag(inv(A(:,:,i))));
        lower_bound(:,i) = b_bar(:,i)-sqrt(c_bar(:,i)*diag(inv(A(:,:,i))));

        % sampling from ellipsoid
        X_Cnz=randn(Params.n^2,Params.Runs);
        X_Cnz=X_Cnz./kron(ones(Params.n^2,1),sqrt(sum(X_Cnz.^2))); % Points uniformly distributed on hypersphere
        R=ones(Params.n^2,1)*(rand(1,Params.Runs).^(1/Params.n^2));
        unif_sph=R.*X_Cnz; % Runs points in the hypersphere
        TT=chol(inv(A(:,:,i))); % Cholesky factorization
        unif_ell =TT'*unif_sph ; % Hypersphere to hyperellipsoid mapping
        T2=(unif_ell*sqrt(c_bar(i))+(b_bar(:,i)*ones(1,Params.Runs)))'; % Translation around gate center
        Tsel=[Tsel;T2]; % samples in original coordinates
        logV(i)=log(((pi^(Params.n^2/2))/(gamma(Params.n^2/2+1))))+0.5*sum(log(eig(inv(A(:,:,i))*c_bar(i))));
    end
end

if flag == 0
    % uniform resampling from SPS region
    V=exp(logV);
    Prob=V/sum(V);
    Xsel=Tsel(1,:)';
    selX=1;

    for i=1:Params.r-1
        BelongX(i)=((Tsel(selX,:)'-b_bar(:,i))'*A(:,:,i)*(Tsel(selX,:)'-b_bar(:,i))<=c_bar(i));%ellipsoids flag
    end

    for k=1:Params.nMCMC
        q=randsample(Params.r-1,1,true,Prob);%randomly chosen ellipsoid
        selY=(q-1)*Params.Runs+ceil(Params.Runs*rand);
        Ysel=Tsel(selY,:)';
        for i=1:Params.r-1
            BelongY(i)=((Tsel(selY,:)'-b_bar(:,i))'*A(:,:,i)*(Tsel(selY,:)'-b_bar(:,i))<=c_bar(i)); %ellipsoids flag
        end
        U=rand;
        if U<(sum(BelongX)/sum(BelongY))
            accX(k,:)=Ysel';
            Xsel=Ysel;
            BelongX=BelongY;
        else
            accX(k,:)=Xsel';
        end
    end

    centroid_SPS = mean(accX);
    dist_SPS = norm(accX-mean(accX)); % distance from center (=centroid) SPS
else
    dist_SPS = NaN;
end

% Gaussian Fisher
A_GF = Omega'*Omega;
b_bar_GF = theta_hat_LS;
c_bar_GF = Params.sigma_nom^2*chi2inv(1-Params.q/Params.r,Params.N)-norm(y-Omega*theta_hat_LS)^2;

if c_bar_GF>0
    % uniform sampling from Fisher region
    X_Cnz=randn(Params.n^2,2*Params.Runs);
    X_Cnz=X_Cnz./kron(ones(Params.n^2,1),sqrt(sum(X_Cnz.^2))); % Points uniformly distributed on hypersphere
    R=ones(Params.n^2,1)*(rand(1,2*Params.Runs).^(1/Params.n^2));
    unif_sph=R.*X_Cnz; % Runs points in the hypersphere
    TT=chol(inv(A_GF)); % Cholesky factorization
    unif_ell =TT'*unif_sph ; % Hypersphere to hyperellipsoid mapping
    Tsel=(unif_ell*sqrt(c_bar_GF)+(b_bar_GF*ones(1,2*Params.Runs)))'; % Translation around gate center
    dist_GF = norm(Tsel'-b_bar_GF); % distance from center Fisher
else
    dist_GF = 0;
end
% Observed Fisher
A_OF = Omega'*Omega;
b_bar_OF = theta_hat_LS;
c_bar_OF = Params.sigma_nom^2*chi2inv(1-Params.q/Params.r,Params.n^2);

% uniform sampling from OF region
X_Cnz=randn(Params.n^2,2*Params.Runs);
X_Cnz=X_Cnz./kron(ones(Params.n^2,1),sqrt(sum(X_Cnz.^2))); % Points uniformly distributed on hypersphere
R=ones(Params.n^2,1)*(rand(1,2*Params.Runs).^(1/Params.n^2));
unif_sph=R.*X_Cnz; % Runs points in the hypersphere
TT=chol(inv(A_OF)); % Cholesky factorization
unif_ell =TT'*unif_sph ; % Hypersphere to hyperellipsoid mapping
Tsel=(unif_ell*sqrt(c_bar_OF)+(b_bar_OF*ones(1,2*Params.Runs)))'; % Translation around gate center
dist_OF = norm(Tsel'-b_bar_OF); % distance from center Bayes

% compute l2 ball radius (centered at theta_hat_SPS)
% SPS ball
% if flag == 0
%     radius_SPS_values = [];
%     
%     for i = 1:Params.r-1
%         cvx_begin sdp quiet
%         variable gam
%         variable lambda
%         maximize(gam)
%         lambda>=0
%         P = [-eye(Params.n^2)+lambda*A(:,:,i), theta_hat_SPS+lambda*b(:,i);
%             (theta_hat_SPS+lambda*b(:,i))', -theta_hat_SPS'*theta_hat_SPS+lambda*c(i)-gam];
%         0.5*(P+P')>=0 % symmetrized to prevent roundoff errors
%         cvx_end
%         radius_SPS_values = [radius_SPS_values sqrt(-gam)];
%     end
% 
%     radius_SPS_ord = sort(radius_SPS_values,'descend'); % SPS radius
%     radius_SPS = radius_SPS_ord(Params.q);
% else
%     radius_SPS = Inf;
% end


% compute l2 ball radius (optimizing center)
% SPS ball
if flag == 0
    radius_SPS_values = [];

    cvx_begin sdp quiet
    variable xc(Params.n^2,1)
    variable t
    variable gam
    variable tau(Params.r-1,1)
    minimize(t)
    for i = 1:Params.r-1
        P = [eye(Params.n^2) -xc;
            -xc' gam]-tau(i).*[A(:,:,i) b(:,i); b(:,i)' c(i)];
        0.5*(P+P')<=0 % symmetrized to prevent roundoff errors
        tau(i)>=0
    end
    [eye(Params.n^2) xc; xc' t+gam]>=0;
    cvx_end

    radius_SPS = sqrt(xc'*xc-gam);

else
    radius_SPS = Inf;
end



if c_bar_GF>0
    radius_GF = sqrt(c_bar_GF)*max(sqrt(eig(inv(A_GF)))); % GF radius
else
    radius_GF = 0;
end

radius_OF = sqrt(c_bar_OF)*max(sqrt(eig(inv(A_OF))));  % OF radius

% box bounds SPS
if flag == 1
    box_bounds_SPS = NaN(Params.n^2,2);
else
    upper_bound_SPS = max(upper_bound,[],2);
    lower_bound_SPS = min(lower_bound,[],2);
    box_bounds_SPS = [upper_bound_SPS lower_bound_SPS];
end

% box bounds GF~
upper_bound_GF = b_bar_GF+sqrt(c_bar_GF*diag(inv(A_GF)));
lower_bound_GF = b_bar_GF-sqrt(c_bar_GF*diag(inv(A_GF)));
box_bounds_GF = [upper_bound_GF lower_bound_GF];

% box bounds OF
upper_bound_OF = b_bar_OF+sqrt(c_bar_OF*diag(inv(A_OF)));
lower_bound_OF = b_bar_OF-sqrt(c_bar_OF*diag(inv(A_OF)));
box_bounds_OF = [upper_bound_OF lower_bound_OF];

end