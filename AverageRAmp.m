function [Rmean,AmpMean,NMeasurements,iAmpBins]=AverageRAmp(R,Amp,iAmpBins,Nampbin)

if (length(iAmpBins)==1)
    NAmpBins=iAmpBins;
else
    NAmpBins=Nampbin;
end

if (nargin<3)|(length(iAmpBins)==1)
    Iexists=isfinite(R);
    Ampexists=(Amp(Iexists));
    sortedAmp=log10(sort(Ampexists(:)));
    
    NRtotal=length(sortedAmp);
    NAmpBins=Nampbin;
    dAmp=floor(NRtotal/(NAmpBins));
    iAmpBins=sortedAmp(1:dAmp:NRtotal);
end
% average over amplitude bin
RaisedAmp=10.^iAmpBins;
   parfor iAmp=1:length(iAmpBins)-1,%Ampmin:dAmp:Ampmax,
     
        Amp1=RaisedAmp(iAmp);
        Amp2=RaisedAmp(iAmp+1);
        IAmpBIN=(Amp(:)>Amp1)&(Amp(:)<Amp2);
        AmpMean(iAmp)=10.^mean(iAmpBins(iAmp:iAmp+1));
    
        Rmean(iAmp)=nanmean(R(IAmpBIN));
        NMeasurements(iAmp)=sum(isfinite(R(IAmpBIN)));
        end
      
end
