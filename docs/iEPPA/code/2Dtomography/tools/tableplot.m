function [] = tableplot(T,varargin)
%TABLEPLOT Summary of this function goes here
%   Detailed explanation goes here
% Hong 08/21

if isempty(varargin)
    plot(T{:,:});
else
    fn = varargin{1};
    feval(fn,T{:,:});
end
legend(T.Properties.VariableNames);

end

