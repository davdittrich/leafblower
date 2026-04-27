function [Cell] = cellOnes1D(vec)
%   cellOnes generates cell of ones matrices of size 
%   Detailed explanation goes here
N = length(vec);
Cell = cell(N,1);
for nn = 1:N
   Cell{nn} = ones(vec(nn),1); 
end
end

