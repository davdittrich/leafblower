%% round a matrix P onto the OT constraint set
%% Copyright: December 2022
%% Hong Chu, Ling Liang, Kim-Chuan Toh, and Lei Yang

function Q = ot_rounding_2d(P, a, b) 

% check nonnegativity of P
if min(P(:)) < 0
    error('Some entries in P are negative!');
end

% check dimension
if size(a,2)~=1, a = a'; end
if size(b,2)~=1, b = b'; end

% rounding 
P_rowsum = sum(P,2);
P1 = bsxfun(@times, min(a./P_rowsum,1), P);    %P1 = repmat(min(a./P_rowsum,1), 1, n) .* P;
P1_colsum = sum(P1,1);
P2 = bsxfun(@times, min(b'./P1_colsum,1), P1); %P2 = repmat(min(b'./P1_colsum,1), m, 1) .* P1;
Da = a - sum(P2,2);
Db = b - sum(P2,1)';
Q = P2 + (Da/sum(Da))*Db';
Q = max(Q, 0);
end