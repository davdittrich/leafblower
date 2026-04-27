function [P,infoE] = EPPA(C,Opers,OpersAdj,Proj,options)


arguments
    C double
    Opers cell
    OpersAdj cell
    Proj cell
    
    options.U           double = ones(size(C));
    options.Tolerance   double = 1e-5;
    options.epsilon     double = 0.01;
    options.FlagEpsilon logical = false;
    options.MaxIter     double = 500;
    options.MaxTime     double = 1800;
    options.FlagPrint   double = 2;
    options.FlagPlot    double = 1;
    options.FlagDphi    logical = false;
    
end

%% Options for EPPA, BCD
U               = options.U;
epsilon         = options.epsilon;
FlagEpsilon     = options.FlagEpsilon;
maxiterEPPA     = options.MaxIter;
MaxTime         = options.MaxTime;
% Options for EPPA
stoptolEPPA     = options.Tolerance;  
printEPPAyes    = options.FlagPrint;
dphiyes         = options.FlagDphi;
FlagPlot        = options.FlagPlot;
% Options BCD
stoptolBCD      = 1e-4; % warm up
minstoptolBCD   = 1e-8;
maxiterBCD      = 50; % warm up
printBCDyes     = printEPPAyes-1;
rangeiter       = 5;
tic
%% Init
Size = size(C);
N = length(Proj); % number of projections

nrmProj = sqrt(sum(cellfun(@(xx) norm(xx,'fro')^2,Proj)));
nrmC    = norm(C,'fro');
nrmU    = norm(U,'fro');
K                   = exp(-C/epsilon);
M                   = exp(-C/epsilon);
P                   = M; % primal solution
Pold                = ones(Size);  % primal solution
Gamma               = ones(Size); % stablized dual variable for upper bound 
Xi                  = cellOnes1D(cellfun(@length,Proj));  % stablized dual variable for constraints

status              = 'NOT OPTIMAL';
%% Recorders
DATAfield          = {'EPPAiter','BCDiter','kkt','pobj','dobj','gap','Dphi'... 
                    'delta1','delta2','delta3','delta4','delta5','delta6','delta7',... 
                    'epsilon','BCDtol','time'}; 
DATAarray            = zeros(maxiterEPPA,length(DATAfield));
BCDiters        = zeros(maxiterEPPA,1); 
%% Print header
if printEPPAyes > 0
    fprintf('\n<strong>Solving by EPPA+BCD </strong>\n');
    fprintf('Parameters: epsilon: %5.4e, tolerance: %5.4e, max iteration: %4d\n\n',epsilon,stoptolEPPA,maxiterEPPA);
end
EPPAmsg = 'EPPA terminates by maxiter';

%% EPPA 
for iterEPPA = 1:maxiterEPPA  
    %% BCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BCDmsg = 'BCD terminates by maxiter'; delta1 = 1;
    PrintVecBCD = [1,round(maxiterBCD.*(1:10)./10)];
    for iterBCD = 1:maxiterBCD
        %% BCD update
        Temp = M./OpersAdj{1}(Xi{1});
        Xi{1} = Proj{1}./Opers{1}(Temp);
        for nn = 2:N
            Temp    = Temp./OpersAdj{nn}(Xi{nn}).*OpersAdj{nn-1}(Xi{nn-1});
            Xi{nn}  = Proj{nn}./Opers{nn}(Temp);
        end  
        Temp = Temp./Gamma.*OpersAdj{N}(Xi{N});
        Gamma = min(1,U./Temp);
        M = Temp.*Gamma; P = M;
        W = epsilon*log(Gamma);
        %% Compute feasibility      
        % Primal feasibility (BCD) (expensive)
        if delta1 <  stoptolBCD*1.001 | ismember(iterBCD,PrintVecBCD)
            cal_delta1_yes = 1;
            delta1 = sqrt(sumcell(cellfun(@(x1,x2) norm(x1(P)-x2,'fro')^2,Opers,Proj,...
                'Uni',false)))/(nrmProj+1); 
        else 
            cal_delta1_yes = 0;
        end
        delta3 = norm(min(P,0),'fro')/(1+nrmU);
        % Dual feasibility (BCD) 
        delta5 = norm(max(W,0),'fro')/(1+norm(W,'fro'));  
        % KKT (BCD)    
        delta6 = abs(sum((U(:)-P(:)).*W(:)))/(1+nrmU);
        deltaBCD = max([delta1,delta3,delta5,delta6]);
        
        ttime = toc;
        %% Record results and print BCD
        if printBCDyes > 0 & cal_delta1_yes & ismember(iterBCD,PrintVecBCD)
            fprintf('\n                 |del1 %5.2e  |del3 %5.2e  |del5 %5.2e  |del6 %5.2e  (%.0f)',delta1,delta3,delta5,delta6,iterBCD);
        end       
        %% Check BCD convergence 
        if deltaBCD < stoptolBCD 
            BCDmsg = 'BCD converges'; break
        elseif ttime > MaxTime
            BCDmsg = 'BCD reaches max time'; break
        end
    end
    % Finish BCD
    if printBCDyes > 0  
        fprintf('\n                 |deltaBCD %5.2e  tolBCD %1.2e                             (%.0f)%s',deltaBCD,stoptolBCD,iterBCD,BCDmsg);
    end
    BCDiters(iterEPPA) = iterBCD;      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %% Recover dual variables 
    Y = cellfun(@(x) epsilon*log(x),Xi,'Uni',false);
    AtYCell = cellfun(@(operadj,Yn) operadj(Yn),OpersAdj,Y,'Uni',false);
    sumAtY = sumcell(AtYCell);
    W = epsilon*log(Gamma);
    %% Compute residues
    % KKT
    SS     = sumAtY + W - C;
    delta2 = norm(max(SS,0),'fro')/(1+nrmC);
    delta3 = norm(min(P,0),'fro')/(1+nrmU);
    delta4 = norm(min(U-P,0),'fro')/(1+nrmU);
    delta5 = norm(max(W,0),'fro')/(1+norm(W,'fro'));  
    delta6 = abs(sum((U(:)-P(:)).*W(:)))/(1+nrmU);
    delta7 = abs(sum(P(:).*SS(:)))/(1+nrmC);
    deltaKKT = max([delta1,delta2,delta3,delta4,delta5,delta6,delta7]);
    % GAP
    pobj    = sum(P(:).*C(:));
    dobj = sum(cellfun(@(cel1,cel2) cel1'*cel2,Proj,Y)) + W(:)'*U(:);
    deltaGAP = abs(pobj-dobj)/(1+abs(pobj)+abs(dobj));
    
    %% KL (optional, very slow)
    if dphiyes 
        dphi = Gprocedure(P,Opers,OpersAdj,Proj,U);
    else
        dphi = NaN;
    end
    %% Print
    if printEPPAyes > 0
        fprintf('\n<strong>                 |del1 %5.2e |del2 %5.2e |del3 %5.2e |del4 %5.2e |del5 %5.2e |del6 %5.2e |del7 %5.2e </strong>',...
                delta1,delta2,delta3,delta4,delta5,delta6,delta7);
    end
    %% Record BCD history
    DATAarray(iterEPPA,:)      = [iterEPPA,iterBCD,deltaKKT,pobj,dobj,deltaGAP,dphi,...
                            delta1,delta2,delta3,delta4,delta5,delta6,delta7,...
                            epsilon,stoptolBCD,ttime];
    %% Check for terminations
    error = max([deltaKKT,deltaGAP]);
    if sum(isnan([deltaKKT,deltaGAP]))
        status = 'NUMERICAL ERROR';
        EPPAmsg = 'EPPA is terminated';
        error = NaN;
        P = Pold;
        break
    elseif error < stoptolEPPA
        status = 'OPTIMAL';
        EPPAmsg = 'EPPA converges';
        break
    elseif ttime > MaxTime
        status = 'NOT OPTIMAL';
        EPPAmsg = 'EPPA reaches max time';
        break
    end  
    %% Update EPPA
    M = M./Pold.*P;
    Pold = P;
    %% Print results
    if printEPPAyes > 0 
        fprintf('\n<strong>EPPA %3d|BCD %3d |KKT  %5.4e |GAP %5.4e |EPS %5.4e |Time %2.1f</strong>\n',iterEPPA,iterBCD,deltaKKT,deltaGAP,epsilon,ttime);
    end
    if FlagPlot > 1
        figure(112);
        imagesc(P); pbaspect([1,1,1]); axis off;
    end
    %% Adjust tol, maxiter for BCD
    decreasing_rate = [0.8;0.85;0.9;0.95;1.05];
    subitervec = max(BCDiters(max(1,iterEPPA-rangeiter):iterEPPA));
    subscore = find(subitervec <= maxiterBCD*[0.25;0.5;0.75;0.9;1],1);
    stoptolBCD = max(minstoptolBCD,stoptolBCD*decreasing_rate(subscore));
    if subscore >= 4
        maxiterBCD = min(500,round(maxiterBCD*1.2));
    else        
        max_val_range = max(cellfun(@(cel) max(cel(:)), Xi)./cellfun(@(cel) min(cel(:)), Xi));
        if FlagEpsilon & max_val_range < 1e128 & min(P(:)) > 1e-128
            epsilon = max(1e-3,epsilon*decreasing_rate(subscore)^0.5);
            Kold = K; K = exp(-C/epsilon);
            M = M./Kold.*K;
        end
    end
    
    
    
    
end

%% Output
DATAarray = DATAarray(1:iterEPPA,:);
DATA = splitvars(table(DATAarray));
DATA.Properties.VariableNames = DATAfield;
infoE.DATA = DATA;
infoE.BCDiters = sum(BCDiters);
infoE.status  = status;
infoE.iter = iterEPPA;
infoE.time = ttime;
infoE.error = error;

%% Plot
if FlagPlot > 0  
    warning('off');
    close(figure(111)); figure(111);
    tiledlayout(2,2,'TileSpacing','tight'); sgtitle('EPPA');
    
    nexttile;
    tableplot(DATA(:,8:14),'semilogy');    hold on    
    semilogy(DATA{:,16},'k-.','DisplayName',DATAfield{16}); hold off
    title('delta');
    
    nexttile;
    tableplot(DATA(:,[3,6,7]),'semilogy');
    title('stopping conditions');
    
    nexttile; 
    tableplot(DATA(:,[4,5]),'plot');
    title('objective values');
    
    nexttile; 
    yyaxis left; bar(DATA{:,2}); hold on; 
    yyaxis right; semilogy(DATA{:,15}); hold off; 
    legend({'BCD iter','epsilon'});
end
% Print
if printEPPAyes > 0
    fprintf(['\n',EPPAmsg,'\n\n']);
    fprintf(1,'\n\n<strong>EPPA+BCD output</strong>');
    fprintf('\nEPPA iters %3d, BCD iters  %3d, running time %.1f ',iterEPPA,sum(BCDiters),ttime);
    fprintf('\nStatus : %s.  ',status);
    fprintf('\nObjective %5.8e ',pobj);
    fprintf('\n<strong>-------------------End EPPA+BCD-------------------</strong>\n\n\n\n');
end




end

