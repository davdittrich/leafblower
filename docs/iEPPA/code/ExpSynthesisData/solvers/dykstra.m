function [P,Q,info] = dykstra(a,b,C,U,epsilon,options)

printyes = 1;
maxiter  = 100000;
stoptol  = 1e-6;
tstart   = clock;

Fnorm = @(x) norm(x, 'fro');

if isfield(options,'maxiter'); maxiter = options.maxiter; end
if isfield(options,'stoptol'); stoptol = options.stoptol; end
if isfield(options,'nrmab'); nrmab= options.nrmab; end

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
fprintf('\n solving capacitied OT by Dykstra algorithm');
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');

% K0 = options.K0;
K0 = exp(-C/epsilon);
K0(K0<=1e-100) = 1e-100;
m = options.m; n = options.n; em = ones(m,1); en = ones(n,1);
Q0 =  em*en'; Qm1 = Q0; Qm2 = Q0;

P = K0;

err = zeros(maxiter,1); obj = zeros(maxiter,1);
breakyes = 0;

for iter = 1:maxiter
    
    PQ = P.*Qm2;
    P0 = P;
    %     if rem(iter,3) == 1
    %
    %         P = (a./sum(PQ,2)).*PQ;
    %         Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    %
    %     elseif rem(iter,3) == 2
    %
    %         P = (em*(b./sum(PQ,1)')').*PQ;
    %         Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    %
    %     elseif rem(iter,3) == 0
    %
    %         P = min(PQ,U);
    %         Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    %
    %     end
    
    
    
    P = (a./sum(PQ,2)).*PQ;
    Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    
    PQ = P.*Qm2;
    P0 = P;
    P = (em*(b./sum(PQ,1)')').*PQ;
    Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    
    PQ = P.*Qm2;
    P0 = P;
    P = min(PQ,U);
    Qtmp = Q0; Q0 = Qm2.*(P0./P); Qm2 = Qm1; Qm1 = Qtmp;
    
    
    
    nrm1 = Fnorm(sum(P,2)-a);
    nrm2 = Fnorm(sum(P,1)'-b);
    err(iter+1) = sqrt(nrm1^2+nrm2^2)/(1+nrmab);
    obj(iter) = sum(P(:).*C(:));
    
    if err(iter+1) < stoptol && iter > 5 ...
            %&& abs(obj(iter) - obj(iter-1))/abs(obj(iter-1)) < stoptol
        breakyes = 1;
        msg = 'Dykstra accuracy reached';
    end
    
    if printyes == 1 && rem(iter,10) == 1
        ttime = etime(clock,tstart);
        fprintf('\n     %3d| %5.4e %5.4e| %3.1f',iter,err(iter+1),obj(iter),ttime);
    end
    
    if breakyes
        break;
    end
    
end

Q = Q0;
if iter == maxiter
    msg = 'maxiter reached';
end

info.msg   = msg;
info.iter  = iter;
info.ttime = etime(tstart,clock);
info.err   = err(1:iter+1);
info.obj   = obj(1:iter);

fprintf('\n     %s, obj = %5.4e\n',msg,info.obj(end));
end