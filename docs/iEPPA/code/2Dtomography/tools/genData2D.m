function [P0,U0,Const] = genData2D(filename,options)
%   genData function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    filename
    options.UpperBound double = 1;
    options.Size (1,2) double = [128,128]
    options.FlagNormalize logical = true;
end

Size            = options.Size;
UpperBound      = options.UpperBound;
FlagNormalize    = options.FlagNormalize;

if isnumeric(filename)
    P0 = filename;
elseif ischar(filename) & isfile(filename)
    P0 = im2double(imread(filename));
    if size(P0,3) == 3
        P0 = rgb2gray(P0); 
    end
else
    P0 = randn(Size); P0 = P0 + max(abs(P0(:))) + 0.5;
end
size1 = size(P0,1); size2 = size(P0,2); 

% Prior matrix (upper bound)
if isnumeric(UpperBound) & numel(UpperBound) == numel(P0)
    U0 = reshape(UpperBound,size(P0));
elseif ischar(UpperBound) & isfile(UpperBound)
    U0 = imread(UpperBound); 
    if size(U0,3) == 3
        U0 = rgb2gray(U0);
    end
    U0 = double(U0); 
else
    U0 = 2*max(P0(:))*ones(size(P0));
end

%% Reduce dimensions to 'Size'
reduce1 = round(size1/Size(1)); reduce2 = round(size1/Size(2));
P0 = P0(1:reduce1:size1,1:reduce2:size2);
U0 = U0(1:reduce1:size1,1:reduce2:size2);

%% Normalize data
if FlagNormalize
    Const = sum(P0(:));
    P0 = P0./Const;
    U0 = U0./Const;
else
    Const = 1;
end

end



























