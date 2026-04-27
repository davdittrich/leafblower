%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-D tomography Matlab toolbox
% Chu, H., Liang, L., Toh, K. C., & Yang, L. (2020). 
% An efficient implementable inexact entropic proximal point algorithm 
% for a class of linear programming problems. arXiv preprint arXiv:2011.14312.
% Please refer to readme.pdf
% Press CTRL+Enter to run each section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;
addpath(genpath(pwd));
filename         = './figures/2D_true3.png'; 
Size            = [256,256];
[P0,U0,Const]   = genData2D(filename,...
    'UpperBound',[],'Size',Size,'FlagNormalize',true);  
Size1 = size(P0,1); Size2 = size(P0,2);

%% Cost matrix
vec1        = linspace(-Size(1)/2,Size(1)/2,Size(1));
vec2        = linspace(-Size(2)/2,Size(2)/2,Size(2));
[II,JJ]     = meshgrid(vec1,vec2);
C           = abs(II-JJ).^2;  
C = C./max(abs(C(:)));

%% Get projection (linear contraints)
N = 10; % number of projections
source = 1;
if source == 1
    % either read directions from file (default in our experiments)
    input_file = './databank/Dir.xlsx';
    DIR = readtable(input_file,detectImportOptions(input_file));
    DIR = DIR{:,:};
    Dir = DIR(:,1:N);
    [Opers,OpersAdj,MatAlins]  = genCollectionA(Size,Dir,'Input','Ratio');
else
    % or get directions from degree (user selections)
    Angles = linspace(-90,90,N+1)'; Angles(1) = [];
    [Opers,OpersAdj,MatAlins]  = genCollectionA(Size,Angles,'Input','degree');
end
%% get projections
Proj = cellfun(@(cel) cel(P0),Opers,'Uni',false);

%% run
stoptol = 1e-5;
epsilon = 0.05;
epsilon_adjust_yes = true;
[P_EPPA,infoE]  = EPPA(C,Opers,OpersAdj,Proj,...
    'Tolerance',stoptol,'MaxIter',1e3,'MaxTime',3600,...
    'epsilon',epsilon,'FlagEpsilon',epsilon_adjust_yes,...
    'FlagPrint',2,'FlagPlot',1,'FlagDphi',false);

%% radon
Angles = linspace(-90,90,N+1)'; Angles(1) = [];
[R,xp] = radon(P0,Angles);
JJ = iradon(R,Angles); 
diff1 = floor((size(JJ,2)-Size1)/2); diff2 = floor((size(JJ,2) - Size2)/2);
P_radon = JJ(diff1+1:diff1+Size1, diff2+1:diff2+Size2); 
infoR.time = toc;
%% plot solution
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex')
P_plot = P_EPPA; % P_plot = P_radon;
figure(1); tl = tiledlayout(1,2,'TileSpacing','tight');
color_limit = [min(P0(:)), max(P0(:))];
nexttile; imagesc(P0); caxis(color_limit);
pbaspect([1,1,1]); axis off; title('ground-truth');
nexttile; imagesc(P_plot);  caxis(color_limit);
pbaspect([1,1,1]); axis off;
title(['N = ',num2str(N),', PSNR = ',...
    num2str(psnr2(P_plot,P0,max(P0(:))))]); 
cb = colorbar; cb.Layout.Tile = 'south';

