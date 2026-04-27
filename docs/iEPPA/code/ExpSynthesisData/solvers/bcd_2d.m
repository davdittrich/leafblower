%% Block corrdinate descent method for solving dual ieppa subproblem for 2d problem
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function [X, f, g, W, info] = bcd_2d(a, b, C, U, epsilon, X, g, W, options)
%% options
tolRp = options.tolRp;
tolRb = options.tolRb;
maxiter = options.maxiter;
normab = options.normab;
normU = options.normU;
xi = options.xi;
compbd = options.compbd;
%% initial points
M0 = exp(-C/epsilon);
Mtld = M0 .* X;
gtld = exp(g/epsilon);
Wtld = exp(W/epsilon);
WMtld = Wtld .* Mtld;
%% main loop
totalbdist = 0;
for iter = 1:maxiter
    %% BCD
    ftld = a ./ (WMtld*gtld);
    gtld = b ./ ((ftld'*WMtld)');
    fgt = ftld*gtld';
    Wtld = min(1, U./(Mtld.*fgt));
    WMtld = Wtld .* Mtld;
    X = max(eps(0), fgt.*WMtld);
    Rp1 = a - sum(X, 2);
    Rp2 = b - sum(X, 1)';
    errRp = sqrt(norm(Rp1)^2+norm(Rp2)^2) / normab;
    bdist = 0;
    
    %% check for termination
    if errRp < tolRp
        if compbd == 1
            bdist = bregdist_2d(X, xi, a, b, U) / normU;
            totalbdist = totalbdist + 1;
            if bdist < tolRb
                break;
            end
        else
            break;
        end
        
    end
end
%% output
f = epsilon*log(ftld);
g = epsilon*log(gtld);
W = epsilon*log(Wtld);
%% info
info.iter = iter;
info.errRp = errRp;
info.errRb = bdist;
info.totalbdist = totalbdist;
end