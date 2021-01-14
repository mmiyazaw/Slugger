function [n_mean,n_high,n_low]=CalculateErr(Nbstrp,Rorig,Amp,Rb,iAmpBins,Nampbin0)

Rvect=Rorig(:);
[R_mean,R_low,R_high]=AverageRAmp_err(Rvect,Amp,iAmpBins,Nbstrp,Nampbin0);

n_high=ComputeRates(R_low,Rb);
n_low=ComputeRates(R_high,Rb);
n_mean=ComputeRates(R_mean,Rb);
