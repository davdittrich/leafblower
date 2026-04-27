%% inexact entropic proximal point algorithm (ieppa) for solving 2-marginal capacited OT problems
%% Input: a, b marginals
%%        C cost matrix
%%        U capacity bound
%%        options struct storing parameters needed for ieppa
%%                options.printyes = 1 (whether to print iterational information)
%%                options.maxiter = 500 (maximum number of iterations)
%%                options.stoptol = 1e-4 (stopping tolerance)
%%                options.epsilon = 5e-2 (entropy parameter for ieppa)
%% Copyright December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function [X, f, g, W, info] = ieppa_2d(a, b, C, U, options)
tstart = clock;
totaliterbcd = 0;
totalbdist = 0;
printyes = 1;
maxiter = 500;
stoptol = 1e-4;
epsilon = 5e-2;
%% options
if isfield(options, 'printyes'); printyes = options.printyes; end
if isfield(options, 'maxiter'); maxiter = options.maxiter; end
if isfield(options, 'stoptol'); stoptol = options.stoptol; end
if isfield(options, 'epsilon'); epsilon = options.epsilon; end

%% norms
norma = norm(a);
normb = norm(b);
normab = 1 + sqrt(norma^2+normb^2);
normC = 1 + norm(C(:));
normU = 1 + norm(U(:));

%% sizes
m = length(a);
n = length(b);

%% initial points
X = a*b'; %% P = ones(m,n);
xi = X(:);
f = zeros(m,1);
g = zeros(n,1);
W = zeros(m,n);

%% initial KKT residual
Rp1 = a - sum(X, 2);
Rp2 = b - sum(X, 1)';
errRp = sqrt(norm(Rp1)^2+norm(Rp2)^2) / normab;
Aty = f*ones(1,n) + ones(m,1)*g';
Z = C - Aty - W;
errRd = norm(min(0, Z(:))) / normC;
errRc = abs(sum(sum(Z.*X))) / normC;
errRu = abs(sum(sum(W .* (U - X)))) / normU;
% maxerr = max([errRp, errRd, errRc, errRu]);
pobj = sum(sum(C.*X));

%% print header
if printyes 
    fprintf("\n ---------------------------------------------------------------------------------------");
    fprintf("\n An inexact EPPA for 2D COT problems");
    fprintf("\n m = %d, n = %d", m, n);
    fprintf("\n ---------------------------------------------------------------------------------------");
    fprintf("\n iter|       objective|    errp    errd    errc    erru    robj|     time| itbcd itbrg|");
    fprintf("\n %4d| %- 9.8e| %2.1e %2.1e %2.1e %2.1e %2.1e| %8.1f| %5d %5d",...
        0, pobj, errRp, errRd, errRc, errRu, 0, etime(clock, tstart), totaliterbcd, totalbdist);
end

%% main loop
breakyes = 0;
msg = [];
opts.compbd = 1;
opts.normab = normab;
opts.normU = normU;
opts.xi = xi;
opts.maxiter = 500;
for iter = 1:maxiter    
    %% BCD method for solving subproblems 
    if iter == 1
        opts.tolRp = 1e0/21^3;       
    else
        opts.tolRp = max(1e-6, errRp/1.5);
    end    
    opts.tolRb = 1e0/iter^1.1;
    [X, f, g, W, infosub] = bcd_2d(a, b, C, U, epsilon, X, g, W, opts);
    totaliterbcd = totaliterbcd + infosub.iter;
    totalbdist = totalbdist + infosub.totalbdist;
    
    %% update kkt residual
    Aty = f*ones(1,n) + ones(m,1)*g';
    Z = C - Aty - W;    
    errRp = infosub.errRp;
    errRd = norm(min(0, Z(:))) / normC;
    errRc = abs(sum(sum(Z.*X))) / normC;
    errRu = abs(sum(sum(W .* (U - X)))) / normU;
    maxerr = max([errRp, errRd, errRc, errRu]);
    pobjold = pobj;
    pobj = sum(sum(C.*X));
    relobj = abs(pobjold-pobj)/(1+abs(pobjold));
    info.pobj(iter) = pobj;
    info.pinf(iter) = errRp;
    info.bdist(iter) = infosub.errRb;
    
    %% check for termination
    if iter > 1 && (maxerr < stoptol || (errRp < stoptol && relobj < stoptol))
        breakyes = 1;
        msg = "converged";
    end
    if iter == maxiter
        breakyes = 2;
        msg = "maxiter reached";
    end
    
    %% print iter
    if printyes
        fprintf("\n %4d| %- 9.8e| %2.1e %2.1e %2.1e %2.1e %2.1e| %8.1f| %5d %5d|",...
            iter, pobj, errRp, errRd, errRc, errRu, relobj, etime(clock, tstart), ...
            totaliterbcd, totalbdist);
        %fprintf(" %2.1e %2.1e %2.1e %2.1e", normC, norm(f), norm(g), norm(W(:)));
    end
    
    %% break
    if breakyes
        if printyes
            fprintf("\n %s\n", msg); 
            fprintf(" ----------------------------------------------------------------------------------------------");
            fprintf("\n");
        end
        break;
    end
end

%% info
info.primobj = pobj;
info.errRp = errRp;
info.errRd = errRd;
info.errRc = errRc;
info.errRu = errRu;
info.totalbdist = totalbdist;
info.totaliterbcd = totaliterbcd;
info.iter = iter;
info.breakyes = breakyes;
info.msg = msg;
info.cputime = etime(clock, tstart);
end