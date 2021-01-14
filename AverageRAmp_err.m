function [Rmean,Rlow,Rhigh,AmpMean,NMeasurements,iAmpBins]=AverageRAmp_err(R,Amp,iAmpBins,Nboot,Nampbin0)
% this version does the bootstrap inside the ampbins

if (length(iAmpBins)==1)
    NAmpBins=iAmpBins;
else
    NAmpBins=Nampbin0;
end

if (nargin<3)|(length(iAmpBins)==1)
    Iexists=isfinite(R);
    Ampexists=(Amp(Iexists));
    sortedAmp=log10(sort(Ampexists(:)));
    
    NRtotal=length(sortedAmp);
    NAmpBins=Nampbin0;
    dAmp=floor(NRtotal/(NAmpBins));
    iAmpBins=sortedAmp(1:dAmp:NRtotal);
end
% average over amplitude bin
RaisedAmp=10.^iAmpBins;
Rmean=nan+iAmpBins(1:end-1);
Rlow=Rmean;
Rhigh=Rmean;
   for iAmp=1:length(iAmpBins)-1,%Ampmin:dAmp:Ampmax,
     
        Amp1=RaisedAmp(iAmp);
        Amp2=RaisedAmp(iAmp+1);
        IAmpBIN=(Amp(:)>Amp1)&(Amp(:)<Amp2);
        AmpMean(iAmp)=10.^mean(iAmpBins(iAmp:iAmp+1));
        NN=sum(IAmpBIN);
        NMeasurements(iAmp)=NN;
        Rvect=R(IAmpBIN);
        if (NN>0)
        parfor ib=1:Nboot,
            Ib=randi(length(Rvect),length(Rvect),1);
            Rboot=Rvect(Ib);
            Rerr_Array(iAmp,ib)=nanmean(Rboot);
        end
        Rmean(iAmp)=nanmean(Rerr_Array(iAmp,:));
        Rlow(iAmp)=prctile(Rerr_Array(iAmp,:),5);
        Rhigh(iAmp)=prctile(Rerr_Array(iAmp,:),95);
        end
   end
      
end
