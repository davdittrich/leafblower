function [Opers,OpersAdj,MatAlins] = genCollectionA(Size,Theta,options)
%GENCOLLECTIONA Summary of this function goes here
%   Detailed explanation goes here

arguments
    Size    (1,2)  double   = 1;
    Theta   double          = [0,90];
    options.foldername char = 'databank';
    options.Input char      = 'degree';
    options.Level double    = [32,32];
end

foldername      = options.foldername;
Input           = options.Input;
Level           = options.Level;
size1           = Size(1); size2 = Size(2); 
size12          = size1*size2;
%% Stored data
filename = ['/data_',num2str(size1),'_',num2str(size2),'_',...
    num2str(Level(1)),'_',num2str(Level(2)),'.mat'];
savename = [foldername,filename];
if isfile(savename)
    load(savename,'HalfTans','HalfAngles','HalfUnits',...
        'AllOpers','AllOpersAdj','AllMatAlins'); 
else
    [Ismall,Jsmall]     = meshgrid(1:Level(1),1:Level(2)); % max(p,q)
    Ismall = Ismall(:); Jsmall = Jsmall(:); IJsmall = [Ismall,Jsmall];
    FirstQuarter = IJsmall - [1,1].*ones(size(IJsmall)); 
    FirstRatio = FirstQuarter(:,1)./FirstQuarter(:,2);
    [Num,Dem] = rat(FirstRatio,0); FirstQuarter = [Num,Dem]; 
    [FirstQuarter,FirstIdx] = unique(FirstQuarter,'rows'); FirstRatio = FirstRatio(FirstIdx);
    FirstQuarter(1,:) = []; FirstRatio(1) = [];% remove (0,0)
    FourthQuarter = [FirstQuarter(:,1),-FirstQuarter(:,2)]; FourthRatio = -FirstRatio;
    FourthQuarter(1:2,:) = []; FourthRatio(1:2) = [];  %remove (0,-1) and (1,-0)
    HalfTans = [FirstRatio;FourthRatio];
    HalfAngles = atan(HalfTans);
    HalfUnits = [FirstQuarter;FourthQuarter];
    NumHalf = size(HalfUnits,1);
    AllOpers = cell(NumHalf,1);
    AllOpersAdj = cell(NumHalf,1);
    AllMatAlins = cell(NumHalf,1);
    save(savename,'HalfTans','HalfAngles','HalfUnits',...
        'AllOpers','AllOpersAdj','AllMatAlins', '-v7.3');    
end


%
if strcmp(Input,'degree') % directions are given by degree
    Theta = Theta(:);
    Theta(Theta == -90) = 90;
    Theta = deg2rad(Theta); Theta = atan(tan(Theta));
    [Grid1,Grid2] = meshgrid(HalfAngles(:),Theta(:));
    Grid = abs(Grid1-Grid2);
    [ErrorTheta,ThetaIdx] = min(Grid,[],2);
    fprintf('\nAngles: gap %1.2e, rounding error %1.2e, ',...
        min(abs(diff(Theta))),max(ErrorTheta)); 
    if length(ThetaIdx) ~= length(unique(ThetaIdx))
        [~,newidx] = unique(ThetaIdx);
       Theta = Theta(sort(newidx));
       ThetaIdx = ThetaIdx(sort(newidx));
    end
elseif strcmp(Input,'Ratio') % directions given by units
    if size(Theta,2) ~=2, Dir = Theta'; else, Dir = Theta; end
    ThetaIdx = find(ismember(HalfUnits,Dir,'rows'));  
end

fprintf('total %d angles.\n', length(ThetaIdx));  
Tans     = HalfTans(ThetaIdx);
Angles = HalfAngles(ThetaIdx);
Units = HalfUnits(ThetaIdx,:);
N = length(ThetaIdx);
%
figure(1001); sgtitle('selected angles (in red) and available angles (in yellow)');
subplot(2,1,1); polarplot(HalfAngles,ones(size(HalfAngles)),'y.');
hold on; polarplot(Angles,ones(size(Angles)),'r.'); hold off
subplot(2,1,2); plot(HalfAngles,ones(size(HalfAngles)),'y.');
hold on; plot(Angles,ones(size(Angles)),'r.'); hold off




%% Loop each directions
ReWriteFlag = false;
[II,JJ]     = meshgrid(1:size1,1:size2);
II = II(:); JJ = JJ(:);
for nn = 1:N
    if isempty(AllMatAlins{ThetaIdx(nn)})
        ReWriteFlag = true;
        ibar = Units(nn,2); jbar = Units(nn,1);
        if ibar == 0
            oper = @(X) sum(X,2);
            operAdj = @(y) repmat(y,[1,size2]);
            matAlin = repmat(speye(size1),[1,size2]);
        elseif jbar == 0
            oper = @(X) sum(X,1).';
            operAdj = @(y) repmat(y.',[size1,1]);
            matAlin = sparse(fix((0:size1*size2-1)'/size1) +1,(1:size1*size2)',ones(size1*size2,1));
        else
            GroupVector = NaN(size12,1);
            counter = 1; 
            while sum(isnan(GroupVector))
                start_pt = find(isnan(GroupVector),1);
                diff_J = JJ - JJ(start_pt);
                diff_I = II - II(start_pt);    
                rat_I = diff_I./ibar; rat_J = diff_J./jbar;
                rem_I = rem(diff_I,ibar);
                who_belongs = ( (rat_I-rat_J)== 0  &  rem_I == 0); 
                GroupVector(who_belongs) = counter;
                counter = counter+1;
            end
            GroupVector = findgroups(GroupVector);
            matAlin = sparse(GroupVector,(1:size12)',ones(size12,1));
            oper = @(X) matAlin*reshape(X,[],1); 
            matAlinT = matAlin';
            operAdj = @(y) reshape(matAlinT*y,Size); 

        end

        AllOpers{ThetaIdx(nn)}    = oper; % this is collection of operators A
        AllOpersAdj{ThetaIdx(nn)} = operAdj; % this is collection of operators A^T
        AllMatAlins{ThetaIdx(nn)} = matAlin; 
    end
end

Opers = AllOpers(ThetaIdx);
OpersAdj = AllOpersAdj(ThetaIdx);
MatAlins = AllMatAlins(ThetaIdx);

if ReWriteFlag
save(savename,'HalfTans','HalfAngles','HalfUnits',...
        'AllOpers','AllOpersAdj','AllMatAlins', '-v7.3');  
end
end



