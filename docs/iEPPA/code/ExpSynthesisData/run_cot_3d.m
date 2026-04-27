%% This testing file is to generate numerical results of iEPPA in Section 4;
%% For 3-marginal capacited constrained optimal transport problems;
%% Three marginals are probability generated from standard uniform distribution on [0,1];
%% The support points of the three marginals are drawn from a Guassian mixture distribution.
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang
clear;
clc;

addpath(genpath(pwd));
% need to add path for gurobi (with valid license)
addpath("/home/liang/Downloads/gurobi9.5.1_linux64/gurobi951/linux64/matlab");

%% problem size 
Nf = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
Ng = Nf;
Nh = Nf;

%% stopping tolerance
stoptol = 1e-5;

%% run testings
rungurobi = 0;
runEPPA = 1;
maxNumCompThreads(4);
idx = 1:length(Nf);

for id = idx(1)
    nf = Nf(id);
    ng = Ng(id);
    nh = Nh(id);
    [ii,jj,kk] = ndgrid(1:nf,1:ng,1:nh);
    
    ranseed = 623;
    rng(ranseed, 'twister');
    
    %% generate data
    d = 3;
    gm_num = 5;                          % number of mixtures
    gm_mean = [-20; -10; 0; 10; 20];     % mean values for the Gaussians means
    sigma = zeros(1,1,gm_num);
    sigma(1,1,:) = 5*ones(gm_num,1);     % variances of the Gaussians
    
    gm_weights = rand(gm_num,1);
    gm_weights = gm_weights/sum(gm_weights);  % generate the mixture weights
    distrib = gmdistribution(gm_mean, sigma, gm_weights);
    Supp1 = reshape(random(distrib,d*nf), d, nf);
    
    gm_weights = rand(gm_num,1);
    gm_weights = gm_weights/sum(gm_weights);  % generate the mixture weights
    distrib = gmdistribution(gm_mean, sigma, gm_weights);
    Supp2 = reshape(random(distrib,d*ng), d, ng);
    
    gm_weights = rand(gm_num,1);
    gm_weights = gm_weights/sum(gm_weights);  % generate the mixture weights
    distrib = gmdistribution(gm_mean, sigma, gm_weights);
    Supp3 = reshape(random(distrib,d*nh), d, nh);
    
    CC = zeros(nf,ng,nh);
    for i = 1:nf
        for j = 1:ng
            for k = 1:nh
                CC(i,j,k) = norm(Supp1(:,i)-Supp2(:,j))^2+norm(Supp1(:,i)-Supp3(:,k))^2+norm(Supp3(:,k)-Supp2(:,j))^2;
            end
        end
    end
    
    C = CC/max(CC(:));
    a = rand(nf,1);
    b = rand(ng,1);
    c = rand(nh,1);
    a = a/sum(a);
    b = b/sum(b);
    c = c/sum(c);
    U = 2*(a(ii).*b(jj).*c(kk)); % upper bounded
    
    
    %% Gurobi
    if (rungurobi == 1)
        clear model;
        clear params;        
        
        % model with parameters
        A = [repmat(speye(nf),1,ng*nh); ...
            repmat(kron(speye(ng),ones(1,nf)),1,nh);...
            kron(speye(nh),ones(1,nf*ng))];
        bb = [a; b; c];
        cc = reshape(C, [nf*ng*nh,1]);
        hh = []; H = [];
        
        model.obj = cc;
        model.A = sparse(A);
        model.rhs = bb;
        model.sense = '=';
        model.lb = zeros(length(cc),1);
        model.ub = reshape(U,[nf*ng*nh,1]);        
        params.Crossover = 0;
        params.Method  = 2;
        params.Threads = 4;
        params.Presolve = 0;
        
        % solve 
        result = gurobi(model,params);    
    end
    
    %% iEPPA
    if runEPPA == 1
        % options
        options.stoptol = stoptol;
        options.maxiter = 500;
        
        % solve 
        [X, f, g, h, W, info] = ieppa_3d(a, b, c, C, U, options);        
    end
end