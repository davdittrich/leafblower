function [P,Q,info] = stabilized_dykstra(a,b,C,U,epsilon,options)

printyes = 1;
maxiter  = 100000;
stoptol  = 1e-6;
tstart   = clock;

Fnorm = @(x) norm(x, 'fro');

if isfield(options,'maxiter'); maxiter = options.maxiter; end
if isfield(options,'stoptol'); stoptol = options.stoptol; end

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
fprintf('\n solving capacitied OT by stabilized-Dykstra algorithm');
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');

m = length(a); n = length(b);
nrmab = sqrt(Fnorm(a)^2+Fnorm(b)^2);

p  = -C; 
u  = epsilon*log(U);
q0 = zeros(m,n); qm1 = q0; qm2 = q0;
aa = epsilon*log(a);
bb = epsilon*log(b);

err = zeros(maxiter,1); obj = zeros(maxiter,1); Time = zeros(maxiter,1);
breakyes = 0;
J = ones(m,n); em = ones(m,1); en = ones(n,1);

for iter = 1:maxiter
    
    p0 = p;
    pq = p+qm2;
    
    %if rem(iter,3)==0
        p = min(pq,u);
        qtmp = q0; q0 = qm2+p0-p; qm2 = qm1; qm1 = qtmp;
        
    %elseif rem(iter,3) == 1
    
        p0 = p;
        pq = p+qm2;
        
        stb = 1;
        if stb == 1
            % choose maximal entry in each row
            D1 = -max(pq,[],2); 
        elseif stb == 2
            % choose minimal entry in each row
            D1 = -min(pq,[],2);   
        elseif stb == 3
            % choose mean value in each row
            D1 = -mean(pq,2);           
        elseif stb == 4
            % choose (max-min)/2
            D1 = -(max(pq,[],2)-min(pq,[],2))/2;  
        else
            D1 = zeros(m,1);
        end
        % D1 is a column vector
        
        pqD1J = pq+D1.*J;
        sum1 = sum(exp(pqD1J/epsilon),2);
        sum1 = max(sum1,eps(0));
        
        % update p
        p = (aa-epsilon*log(sum1)+D1)*en'+pq;
        qtmp = q0; q0 = qm2+p0-p; qm2 = qm1; qm1 = qtmp;
        
    %elseif rem(iter,3) == 2
    
            p0 = p;
        pq = p+qm2;
        
        
        stb = 1;
        if stb == 1
            % choose maximal entry in each colum
            D2 = -max(pq,[],1);          
        elseif stb == 2
            % choose minimal entry in each colum
            D2 = -min(pq,[],1);
        elseif stb == 3
            % choose mean value in each colum
            D2 = -mean(pq,1);     
        elseif stb == 4
            % choose (max-min)/2
            D2 = -(max(pq,[],1)-min(pq,[],1))/2;         
        else
            D2 = zeros(1,n);
        end
        % D2 is a row vector
        
        pqJD2 = pq + J.*D2;
        sum2 = sum(exp(pqJD2/epsilon),1)';
        sum2 = max(sum2,eps(0));
        
        % update p
        p = pq + em*((bb-epsilon*log(sum2)+D2'))';
        qtmp = q0; q0 = qm2+p0-p; qm2 = qm1; qm1 = qtmp;
    %end
    
  
       
    compP = 1;
    if compP == 1
        P = exp(p/epsilon);
        nrm1 = Fnorm(sum(P,2)-a);
        nrm2 = Fnorm(sum(P,1)'-b);
        err(iter) = sqrt(nrm1^2+nrm2^2)/(1+nrmab);
        obj(iter) = sum(P(:).*C(:));
    else
        err(iter) = Fnorm(p-p0)/Fnorm(p0);
    end
   
    
        
    if err(iter) < stoptol && iter > 5 ...
            %&& abs(obj(iter) - obj(iter-1))/abs(obj(iter-1)) < stoptol
        breakyes = 1;
        msg = 'Dykstra accuracy reached';
    end
    
    ttime = etime(clock,tstart); 
    Time(iter) = ttime;
    
    if printyes == 1 && rem(iter,100) == 1
        if compP == 0
            fprintf('\n     %3d| %5.4e %3.1f',iter,err(iter),ttime);
        elseif compP == 1
            fprintf('\n     %3d| %5.4e %5.4e| %3.1f',iter,err(iter),obj(iter),ttime);
        end
    end
    
    if breakyes 
        break;
    end
    
end

Q = exp(q0/epsilon);
if iter == maxiter
    msg = 'maxiter reached';
end


info.msg   = msg;
info.iter  = iter;
info.ttime = etime(tstart,clock);
info.err   = err(1:iter);
info.Time  = Time(1:iter);

if compP == 1 
    info.obj = obj(1:iter);
else
    info.obj   = sum(P(:).*C(:)); 
end


ploterr = 0;
if ploterr == 1
    plot(info.err);
    xlabel('iter');
    ylabel('relative error');
end

plotobj = 0;
if plotobj == 1 && compP == 1
    plot(info.obj);
    xlabel('iter');
    ylabel('objective function value');
end

fprintf('\n     %s, obj = %5.4e\n',msg,info.obj(end));
end