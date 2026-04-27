function [dphi] = Gprocedure(X,Opers,OpersAdj,Proj,U)
%GPROCEDURE Summary of this function goes here
%   Detailed explanation goes here
% 
[N1,N2] = size(X);


if length(Proj) == 2
    alpha = Proj{1}; beta = Proj{2};
    rhs = [sum(X,2)-alpha;sum(X,1)'-beta];
    rhs1 = rhs(1:N1); rhs2 = rhs(N1+1:end);
    sumall = sum(rhs1);
    u = rhs1/N2 - sumall/(2*N1*N2);
    v = rhs2/N1 - sumall/(2*N1*N2);
    Z = X - u - v';
    Xri = alpha*beta';
    J1 = find(Z>U);
    J2 = find(Z<0);

    if isempty(J1) & isempty(J2)
        Xtilde = Z; 
    else
        Term1 = (Z-U)./(Z-Xri);
        Term2 = (-Z)./(Xri - Z);
        lambda1 = max(Term1(J1));
        lambda2 = max(Term2(J2));
        lambda = max([lambda1,lambda2]);
        Xtilde = lambda*Xri + (1-lambda)*Z;
    end
    

else % Hong
    ii = 1; Xtilde = X;
    
    residue = sumcell(cellfun(@(oper,proj) norm(oper(Xtilde)-proj,'fro'),Opers,Proj,'Uni',false))...
        + norm(Xtilde - min(max(Xtilde,0),U),'fro');
    
    while ii < 1e5 & residue > 1e-7
       for nn = 1:length(Proj)
          Xtilde = Xtilde.*OpersAdj{nn}(Proj{nn}./Opers{nn}(Xtilde)); % alternatively project
       end
       Xtilde = min(max(Xtilde,0),U); % project to box contraints
       residue = sumcell(cellfun(@(oper,proj) norm(oper(Xtilde)-proj,'fro'),Opers,Proj,'Uni',false))...
        + norm(Xtilde - min(max(Xtilde,0),U),'fro');
        ii = ii + 1;
    end
    
    

end



dphi = KL(Xtilde,X);






end


function [val] = KL(P,Q)
    p = reshape(P,[],1); p = max(p,eps(0));
    q = reshape(Q,[],1); q = max(q,eps(0));
    p_over_q = p./q; p_over_q(p_over_q<0) = eps(0);
    log_of_p_over_q = log(p_over_q);
    log_of_p_over_q(abs(log_of_p_over_q)<1e-9) = 0;
    val = sum(p.*log_of_p_over_q) - sum(p) + sum(q);
    val = max(val,0);
end










