%% Compute Bregman distance between X and \tilde{X} for 2d problems
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function bd = bregdist_3d(X, xi, a, b, c, U)
x = X(:);
Z = ot_rounding_3d(X, a, b, c);
z = Z(:);
u = U(:);
%% compute lambda
idx1 = find(z > u);
z1 = z(idx1);
u1 = u(idx1);
xi1 = xi(idx1);
lambda = max((z1-u1)./(z1-xi1));
if lambda < 0 || lambda > 1
    fprintf("\n lambda = %- 2.1", lambda);
    error("lambda out of range");
end
xtld = z + lambda*(xi-z);
%% compute bregman distance
idxp = find(xtld > 1e-12);
xtldtmp = xtld(idxp);
xtmp = x(idxp);
bd = sum(x - xtld) + sum(xtldtmp .* log(xtldtmp./xtmp));
end