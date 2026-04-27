%% Block corrdinate descent method for solving dual ieppa subproblem for 3d problem
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function [X, f, g, h, W, info] = bcd_3d(a, b, c, C, U, epsilon, X, g, h, W, options)
tolRp = options.tolRp;
tolRb = options.tolRb;
maxiter = options.maxiter;
normabc = options.normabc;
normU = options.normU;
xi = options.xi;
compbd = options.compbd;
Ones = options.Ones;
%% initial points
n3 = length(c);
M0 = exp(-C/epsilon);
Mtld = M0 .* X;
gtld = exp(g/epsilon);
ggtld = bsxfun(@times, gtld', Ones);
htld = exp(h/epsilon);
hhtld = bsxfun(@times, reshape(htld, 1, 1, n3), Ones);
Wtld = exp(W/epsilon);
WMtld = Wtld .* Mtld;
%% main loop
totalbdist = 0;
for iter = 1:maxiter
    %% BCD    
    ftmp = (ggtld.*hhtld).*WMtld;
    ftld = a ./ sum(sum(ftmp,2),3);
    fftld = bsxfun(@times, ftld, Ones);
    gtmp = (fftld.*hhtld).*WMtld;
    gtld = b ./ (sum(sum(gtmp,1),3)');
    ggtld = bsxfun(@times, gtld', Ones);
    htmp = (fftld.*ggtld).*WMtld;
    sumhtmp = sum(sum(htmp,1),2);    
    htld = c ./ sumhtmp(:);
    hhtld = bsxfun(@times, reshape(htld, 1, 1, n3), Ones);
    fght = fftld .* ggtld .* hhtld;
    Wtld = min(1, U./(Mtld.*fght));
    WMtld = Wtld .* Mtld;
    X = max(eps(0), fght.*WMtld);
    errRp1 = norm(sum(sum(X, 2), 3) - a);
    errRp2 = norm(sum(sum(X, 1), 3)' - b);
    sumX3 = sum(sum(X,1),2);
    errRp3 = norm(sumX3(:) - c);
    errRp = sqrt(errRp1^2+errRp2^2+errRp3^2)/normabc;
    bdist = 0;
    %% check for termination
    if errRp < tolRp
        if compbd == 1
            bdist = bregdist_3d(X, xi, a, b, c, U)/normU;
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
h = epsilon*log(htld);
W = epsilon*log(Wtld);
%% info
info.iter = iter;
info.totalbdist = totalbdist;
info.errRp = errRp;
info.errRb = bdist;
end