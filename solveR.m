function [rateRatio]=solveR(Rmean)
% based on SH's code to solve Rmean->n
rateRatio=nan(size(Rmean));

Ireal=find(Rmean==Rmean);
for i=1:length(Ireal);
    rateRatio(Ireal(i))=fminsearch(@(s) solveS(s,Rmean(Ireal(i))),1);
end
rateRatio=rateRatio-1;

function[residual]=solveS(s,Rmean)
if s==1
    residual=abs(Rmean-0.5);
else
    residual=abs(Rmean-(s./(s-1).^2)*(1./s + log(s)-1));
end