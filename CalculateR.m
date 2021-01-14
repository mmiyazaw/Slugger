function [R,Amp,Mbefore,Mafter,t0Array,t1Array,t2Array]=CalculateR(lattrigger,lontrigger,ttrigger,Ntrig,BinTrgLon,BinTrgLat,BinTrgTm,BinTrgAm,IndexBinsStruct,Bins,latbin,lonbin,t,minsep,Ampmin,Ampmax,dtmax)

R=nan(length(lattrigger),length(latbin));

t0Array=R;
t1Array=R;
t2Array=R;
Mbefore=R; % arrays to store the magnitudes;
Mafter=R;
Amp=R;

     minsep=km2deg(minsep);
     tic

        parfor ibin=1:max(Bins)

            E1 = abs(BinTrgLon-lonbin(ibin))+abs(BinTrgLat-latbin(ibin));
            E2=(E1>0.0001);
            E1=0*E1+1;
            E1(E2)=nan;
            TmpTm = BinTrgTm.*E1;
            TmpTm2 = TmpTm(~isnan(TmpTm))./3600./24.;
            TmpAm = BinTrgAm.*E1;
            TmpAm2 = TmpAm(~isnan(TmpAm));
            D=distance(lattrigger,lontrigger,latbin(ibin),lonbin(ibin)); %,[6378.137 0.081819190]);
            Itooclose=(D<minsep);
            TmpAm2(Itooclose) = nan;
            
            AmpBin=TmpAm2./3.5e3;
            Itoosmall=(AmpBin<Ampmin);
            AmpBin(Itoosmall)=nan;
            Itoolarge=(AmpBin>Ampmax);
            AmpBin(Itoolarge)=nan;

            I=IndexBinsStruct(ibin).I; 
            tbin=t(I);
            dt=zeros(size(tbin));
            Iwindow1=logical(size(dt));
            Iwindow2=Iwindow1;
           
            for itrig=1:Ntrig,
                Amp(itrig,ibin)=AmpBin(itrig);
            if (isfinite(AmpBin(itrig)))
                t0=ttrigger(itrig)+TmpTm2(itrig); % add traveltime
                dt=tbin-t0;   
                Iwindow1=(dt<0)&(dt>-dtmax);

             if (any(Iwindow1))   
                Iwindow2=(dt>0)&(dt<dtmax);

                if(any(Iwindow2)) % there are enough events in the window before and after
                % record the magnitudes
                IndexAfter=find(dt>0,1,'first'); % local subscripts
                IndexBefore=find(dt<0,1,'last'); 
                t1=-dt(IndexBefore);
                t2=dt(IndexAfter);
                Mbeforetemp=IndexBinsStruct(ibin).magbin(IndexBefore);
                Maftertemp=IndexBinsStruct(ibin).magbin(IndexAfter);
              
                    R_bin=t2./(t1+t2);

                t0Array(itrig,ibin)=t0;
                t1Array(itrig,ibin)=t1;
                t2Array(itrig,ibin)=t2;
                R(itrig,ibin)=R_bin;
                Mbefore(itrig,ibin)=Mbeforetemp;
                Mafter(itrig,ibin)=Maftertemp;
                
                end
            
             end
            
            end
            end
        end
     
    toc
  
