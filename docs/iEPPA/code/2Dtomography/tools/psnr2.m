function [value] = psnr2(X,Xref,gamma)
%PSNR2 Summary of this function goes here
%   Detailed explanation goes here
MSE = norm(X-Xref,'fro')^2/(size(X,1)*size(X,2));
value = 10*log10(gamma^2/MSE);
end

