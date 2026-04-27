%% This testing file is to generate numerical results of iEPPA in Section 4;
%% For 2-marginal capacited constrained optimal transport problems;
%% Two marginals are probability generated from standard uniform distribution on [0,1];
%% The support points of the two marginals are drawn from a Guassian mixture distribution.
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

clear;
clc;
addpath(genpath(pwd));
% add gurobi path (gurobi needs to be installed with a valid license to run gurobi)
addpath('/home/liang/Downloads/gurobi9.5.1_linux64/gurobi951/linux64/matlab');

%% problem size
eg = 1;
if eg == 1
    M = [4000,4000,4000,5000,5000,5000,6000,6000,6000,7000,7000,7000];
    N = [2000,4000,8000,2500,5000,10000,3000,6000,12000,3500,7000,14000];
else
    M = [1000,2000,3000,4000,5000,6000,7000,8000,9000];
    N = [1000,2000,3000,4000,5000,6000,7000,8000,9000];
end

%% specified the stopping tolerance
stoptol = 1e-5;

%% only allow matlab to use 4 threads
maxNumCompThreads(4);

rungurobi      = 0;
runEPPA        = 0;
rundykstra     = 0;
rundykstra_sta = 1;

idx = 1:length(N);
for id = idx(1)
    m = M(id);
    n = N(id);
    ranseed = 480;
    rng(ranseed, 'twister');

    %% generate data
    d = 3;
    a = rand(m,1); a = a/sum(a);
    b = rand(n,1); b = b/sum(b);

    % number of mixtures
    gm_num = 5;

    % mean values for the Gaussians means
    gm_mean = [-20; -10; 0; 10; 20];

    % variances of the Gaussians
    sigma = zeros(1,1,gm_num);
    sigma(1,1,:) = 5*ones(gm_num,1);

    % generate the mixture weights
    gm_ws1 = rand(gm_num,1);
    gm_ws1 = gm_ws1/sum(gm_ws1);
    distrib = gmdistribution(gm_mean, sigma, gm_ws1);

    % generate the support points for distribution 1
    Supp1 = reshape(random(distrib,d*m), d, m);

    % generate the mixture weights
    gm_ws2 = rand(gm_num,1);
    gm_ws2 = gm_ws2/sum(gm_ws2);
    distrib = gmdistribution(gm_mean, sigma, gm_ws2);

    % generate the support points for distribution 2
    Supp2 = reshape(random(distrib,d*n), d, n);

    % compute the cost matrix and normalize
    C = repmat(sum(Supp1.*Supp1,1)',1,n) - 2*Supp1'*Supp2 + repmat(sum(Supp2.*Supp2,1),m,1);
    normC = 1+norm(C(:));
    C = C/max(C(:));

    % generate an upper bound U
    U = 2*a*b';

    %% tesing Gurobi
    if (rungurobi == 1)
        clear model;
        clear params;

        % generate gurobi model
        A = [kron(ones(1,n),speye(m)); kron(speye(n),ones(1,m))];
        bb = [a; b];
        cc = reshape(C, [m*n,1]);
        hh = []; H = [];

        model.obj = cc;
        model.A = sparse(A);
        model.rhs = bb;
        model.sense = '=';
        model.lb = zeros(length(cc),1);
        model.ub = reshape(U,[m*n,1]);

        % set param for gurobi
        params.Presolve =  0;
        params.Crossover = 0;
        params.Method  = 2;
        params.Threads = 4;

        % solve model
        result = gurobi(model,params);

    end

    %% testing iEPPA
    if (runEPPA == 1)
        % set options for iEPPA
        options.stoptol = stoptol;
        options.maxiter = 500;

        % solve
        [X, f, g, W, info] = ieppa_2d(a, b, C, U, options);
    end

    %% testing Dykstra's algorithm
    if (rundykstra == 1)
        Eps = [1e-1,1e-2];
        for iteps = 1:2
            % set options for Dykstra
            epsilon = Eps(iteps);
            options.stoptol = stoptol;
            options.maxiter = 20000;
            options.nrmab = sqrt(norm(a)^2+norm(b)^2);
            options.m = m;
            options.n = n;

            % solve
            [P_dyk,Q_dyk,info_dyk] = dykstra(a,b,C,U,epsilon,options);
        end
    end

    %% testing log-domain and stabilized Dykstra's algorithm for small epsilon
    if (rundykstra_sta == 1)
        Eps = [1e-3,1e-4];
        for iteps = 1:2
            % set options for stablized Dykstra
            epsilon = Eps(iteps);
            options.stoptol = stoptol;
            options.maxiter = 200000;

            % solve
            [P_sta,Q_sta,info_sta] = stabilized_dykstra(a,b,C,U,epsilon,options);
        end
    end


end
%% end of testings