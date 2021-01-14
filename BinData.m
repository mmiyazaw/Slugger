function [IndexBinsStruct, Bins, Nbin, latbin, lonbin, nevbin]=BinData(t, lat, lon, mag, depth, nevleast, Mth, Mc, latbinMc, lonbinMc);
% bin data

dlat=0.1;
dlon=0.1;

Nbin=0;
Bins=lat*0;

minlat=32;
maxlat=42;
minlon=-124;
maxlon=-114;

latloopvect=[minlat:dlat:maxlat];% vectors to avoid rounding boundary problems
lonloopvect=[minlon:dlon:maxlon];
Nlatloopvect=length(latloopvect)-1;
Nlonloopvect=length(lonloopvect)-1;
% following loop MUST be serial
for ilatloop=1:Nlatloopvect,
    for ilonloop=1:Nlonloopvect;
        %Nbin=Nbin+1;
        Ilat=(lat>=latloopvect(ilatloop))&(lat<latloopvect(ilatloop+1));
        Ilon=(lon>=lonloopvect(ilonloop))&(lon<lonloopvect(ilonloop+1));
        Ibin=Ilat&Ilon;
        latbintmp=latloopvect(ilatloop)+dlat/2;
        lonbintmp=lonloopvect(ilonloop)+dlon/2;
        E1 = abs(lonbinMc-lonbintmp) + abs(latbinMc-latbintmp);
        E2=(E1>0.0001);
        E1=0*E1+1;
        E1(E2)=nan;
        binMc = rmmissing(Mc.*E1);
        if(binMc<=Mth) %
        if(sum(Ibin)>=nevleast) % only increment if bin is populated with at least [nevleast] events
            Nbin=Nbin+1;
            Bins=Bins+Ibin*Nbin;
            latbin(Nbin)=latbintmp;
            lonbin(Nbin)=lonbintmp;
            nevbin(Nbin)=sum(Ibin);
        end
        end
    end
end

IndexBinsStruct=struct('I',{});
parfor ibin=1:max(Bins),
    I=(Bins==ibin);
    IndexBinsStruct(ibin).I=I;
    IndexBinsStruct(ibin).magbin=mag(I);

end

