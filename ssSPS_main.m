%% ssSPS_main.m
% compute coverage probability, box bounds, mean distance from centroid, L2 ball radius of ssSPS regions
% include comparison with Gaussian Fisher (GF) and Observed Fisher (OF) regions for three noise scenarios

clc
clear all
close all

%% parameters
Params.n = 4; % state dimension
Params.T = 200; % time horizon
Params.T_est = 40; % time horizon to estimate F for IV_case = 2
Params.q = 1; % define desired coverage probability
Params.r = 20; % define desired coverage probability (q/r=1/20: 95% coverage)
Params.N = Params.n*Params.T; % number of data
Params.stab = 0.1; % define margin of stability state matrix
Params.sigma_nom = 1; % nominal sd noise
Params.Runs = 1000; % number of samples from each ssSPS ellipsoids
Params.nMCMC = 2*Params.Runs; % number of MCMC samples for uniform sampling from union of ellipsoids
Params.N_check = 5000; % trials for checking coverage probability

% Select noise scenario: 
% noise_case = 1: Gaussian with known sd
% noise_case = 2: Mixture of 2 Gaussians with correct var
% noise_case = 3: Mixture of 2 Gaussians with incorrect var
Params.noise_case = 1;
Params.sigma_mix2 = 0.01; % sd of first Gaussian in case 2
Params.prob_mix2 = 0.1; % probability activation first Gaussian in case 2
Params.sigma_mix3 = 2*Params.sigma_nom; % sd of first Gaussian in case 3
Params.prob_mix3 = 0.1; % probability activation first Gaussian in case 3

% Select SPS Instrumental Variable: 
% IV_case = 1: one-step delayed input
% IV_case = 2: noise-free reconstructed output using LS estimate
Params.IV_case = 2;

% Flag variable: 1 to compute coverage probability, 0 otherwise
check_freq = 0;
% Flag variable: 1 to compute performance metrics, 0 otherwise
check_perf = 1;

%% compute coverage probability
if check_freq == 1
    [freq_SPS, freq_GF, freq_OF] = coverage(Params);
    disp('==============')
    freq_SPS, freq_GF, freq_OF
    disp('==============')
end

%% compute box bounds, mean distance from centroid, L2 ball radius
if check_perf == 1
    [flag, theta, theta_hat_SPS, theta_hat_LS, box_bounds_SPS, box_bounds_GF, box_bounds_OF, dist_SPS, dist_GF, dist_OF, radius_SPS, radius_GF, radius_OF] = metrics(Params);
    if flag == 1
        disp('==============')
        disp('WARNING: unbounded ssSPS confidence region')
        disp('==============')
    end
    disp('==============')
    dist_SPS, dist_GF, dist_OF, radius_SPS, radius_GF, radius_OF
    disp('==============')

    % plot box bounds
    h = figure;
    SPS = plot(box_bounds_SPS(:,1),'--r','LineWidth',1.5);
    hold on
    plot(box_bounds_SPS(:,2),'--r','LineWidth',1.5)
    hold on
    GF = plot(box_bounds_GF(:,1),'--b','LineWidth',1.5);
    hold on
    plot(box_bounds_GF(:,2),'--b','LineWidth',1.5)
    hold on
    OF = plot(box_bounds_OF(:,1),'--g','LineWidth',1.5);
    hold on
    plot(box_bounds_OF(:,2),'--g','LineWidth',1.5)
    hold on
    l1 = plot(theta,'-*k','LineWidth',1.5);

    legend([l1 SPS GF OF], 'true', 'ssSPS', 'GF', 'OF','interpreter','latex','fontsize',15);
    xlabel('components of $\mathrm{vec}(F)$','interpreter','latex','fontsize',20)
    xlim([1,Params.n^2])

end

