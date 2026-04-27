%% round a 3rd-order tensor P onto the OT constraint set
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function Q = ot_rounding_3d(P, a, b, c)
%% check nonnegative of P
if min(P(:)) < 0
    error('Some entries in P are negative!');
end

%% check dimension
if size(a,2) ~= 1; a = a'; end
if size(b,2) ~= 1; b = b'; end
if size(c,2) ~= 1; c = c'; end

%% rounding
r1 = sum(P,[2,3]);
P1 = bsxfun(@times, min(a./r1,1), P);

r2 = sum(P1,[1,3]);
bb = b';
P2 = bsxfun(@times, min(bb./r2,1), P1);

r3 = sum(P2,[1,2]);
n3 = length(c);
cc = reshape(c,1,1,n3);
P3 = bsxfun(@times, min(cc./r3,1), P2);

ea = a  - sum(P3,[2,3]);
eb = bb - sum(P3,[1,3]);
ec = cc - sum(P3,[1,2]);

Q = P3 + bsxfun(@times, ec/sum(abs(ea))^2, repmat(ea*eb,[1,1,n3]));
Q = max(Q, 0);

end