function [M] = sumcell(C)
%SUMCELL Summary of this function goes here
%   Detailed explanation goes here
M = cat(3,C{:});
M = sum(M,3);
end

