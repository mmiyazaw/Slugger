function [n,nraw,nsimul,RMean,RMeanraw,RMeansimul]=ComputeRates(RMeanraw,RMeansimul);

RMean=0.5+RMeanraw-RMeansimul;
n=solveR(RMean);
nraw=solveR(RMeanraw);
nsimul=solveR(RMeansimul);