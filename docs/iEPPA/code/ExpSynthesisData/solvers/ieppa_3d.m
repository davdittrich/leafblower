%% inexact entropic proximal point algorithm (ieppa) for solving 3-marginal capacited OT problems
%% Input: a, b, c marginals
%%        C cost matrix
%%        U capacity bound
%%        options struct storing parameters needed for ieppa
%%                options.printyes = 1 (whether to print iterational information)
%%                options.maxiter = 500 (maximum number of iterations)
%%                options.stoptol = 1e-4 (stopping tolerance)
%%                options.epsilon = 5e-2 (entropy parameter for ieppa)
%% Copyright December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function [X, f, g, h, W, info] = ieppa_3d(a, b, c, C, U, options)
tstart = clock;
totaliterbcd = 0;
totalbdist = 0;
printyes = 1;
maxiter = 500;
stoptol = 1e-3;
epsilon = 5e-2;
%% options
if isfield(options, 'printyes'); printyes = options.printyes; end
if isfield(options, 'maxiter'); maxiter = options.maxiter; end
if isfield(options, 'stoptol'); stoptol = options.stoptol; end
if isfield(options, 'epsilon'); epsilon = options.epsilon; end
%% sizes
n1 = length(a);
n2 = length(b);
n3 = length(c);
Ones = ones(n1, n2, n3);
%% norms
norma = norm(a);
normb = norm(b);
normc = norm(c);
normabc = 1 + sqrt(norma^2 + normb^2 + normc^2);
cc = C(:);
normC = 1 + norm(cc); c_ratio = max(cc)/min(cc);
u = U(:);
normU = 1 + norm(u);

%% initial points
aff = bsxfun(@times, a, Ones);
bgg = bsxfun(@times, b', Ones);
chh = bsxfun(@times, reshape(c, 1, 1, n3), Ones);
X = aff .* bgg .* chh;
xi = X(:);
f = zeros(n1, 1);
g = zeros(n2, 1);
h = zeros(n3, 1);
W = zeros(n1, n2, n3);

%% initial KKT residual
errRp1 = norm(sum(sum(X, 2), 3) - a);
errRp2 = norm(sum(sum(X, 1), 3)' - b);
sumX3 = sum(sum(X,1),2);
errRp3 = norm(sumX3(:) - c);
errRp = sqrt(errRp1^2+errRp2^2+errRp3^2)/normabc;
fff = bsxfun(@times, f, Ones);
ggg = bsxfun(@times, g', Ones);
hhh = bsxfun(@times, reshape(h, 1, 1, n3), Ones);
Aty = fff + ggg + hhh;
Z = C - Aty - W;
z = Z(:);
x = X(:);
w = W(:);
errRd = norm(min(0, z))/normC;
errRc = abs(sum(x .* z))/normC;
errRu = abs(sum(w .*(u-x)))/normU;
pobj = sum(cc .* x);

%% print header
if printyes 
    fprintf("\n ----------------------------------------------------------------------------------------------");
    fprintf("\n An inexact EPPA for 3D COT problems");
    fprintf("\n n1 = %d, n2 = %d, n3 = %d", n1, n2, n3);
    fprintf("\n ----------------------------------------------------------------------------------------------");
    fprintf("\n iter|       objective|    errp    errd    errc    erru    robj|     time| itbcd itbrg|");
    fprintf("\n %4d| %- 9.8e| %2.1e %2.1e %2.1e %2.1e %2.1e| %8.1f| %5d %5d|",...
        0, pobj, errRp, errRd, errRc, errRu, 0, etime(clock, tstart), totaliterbcd, totalbdist);
end
%% main loop
breakyes = 0;
msg = [];
opts.compbd = 1;
opts.normabc = normabc;
opts.normU = normU;
opts.xi = xi;
opts.maxiter = 500;
opts.Ones = Ones;
for iter = 1:maxiter
    %% BCD for solving subproblems
    if iter == 1
        opts.tolRp = 1e0/21^3;
    else
        opts.tolRp = max(errRp/1.5, 1e-6);
    end
    opts.tolRb = 1e0/iter^1.1;
    [X, f, g, h, W, infosub] = bcd_3d(a, b, c, C, U, epsilon, X, g, h, W, opts);
    totaliterbcd = totaliterbcd + infosub.iter;
    totalbdist = totalbdist + infosub.totalbdist;
    
    %% update kkt residual
    fff = bsxfun(@times, f, Ones);
    ggg = bsxfun(@times, g', Ones);
    hhh = bsxfun(@times, reshape(h, 1, 1, n3), Ones);
    Aty = fff + ggg + hhh;
    Z = C - Aty - W;
    z = Z(:);
    x = X(:);
    w = W(:);
    errRp = infosub.errRp;
    errRd = norm(min(0, z))/c_ratio;
    errRc = abs(sum(x .* z))/normC;
    errRu = abs(sum(w .*(u-x)))/normU;
    maxerr = max([errRp, errRd, errRc, errRu]);
    pobjold = pobj;
    pobj = sum(cc .* x);
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