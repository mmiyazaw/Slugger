% Master R Script - Slugger v3.0 - definitely beta
    clear;

%   readdata=true; % read dataset
    readdata=false; 

if (readdata) % read dataset

% Mc for each bin
   load BinMc.mat
   latbinMc = latbin;
   lonbinMc = lonbin;
   
% QTM catalog available at https://scedc.caltech.edu/data/qtm-catalog.html
    fn='qtm_final_9.5dev.hypo';
    mth=0.5;
    magmax=8;
    tmin=datenum('2008-1-1');
    tmax=datenum('2017-12-31');
    depthmin=0;
    depthmax=15;

% Read catalog
    [t, lat, lon, mag, depth]=ReadCat(fn,tmin,tmax,mth,depthmin,depthmax,magmax);    

    nevleast = 2; %least number of events for each bin
    [IndexBinsStruct, Bins, Nbin, latbin, lonbin, nevbin]=BinData(t, lat, lon, mag, depth, nevleast, mth, Mc, latbinMc, lonbinMc);
    
%%
    % Load in triggers
    Ampmin=2e-9; % strain
    Ampmax=1e0;  % strain
    minsep=50; % minimum separation in kilometers

    %save TrgEq lattrigger lontrigger ttrigger Ntrig;
    load TrgEq.mat;
    %-- distant earthquake --
    % lattrigger: latitude
    % lontrigger: longitude
    % ttrigger: origin time
    % Ntrig: number
    
    %windows
    dtmax=365*2;

    %save BinTrg BinTrgLon BinTrgLat BinTrgTm BinTrgAm;
    load BinTrg.mat;
    %-- PGVs at bins for each distant earthquake --
    % BinTrgLon: bin longitude
    % BinTrgLat: bin latitude
    % BinTrgTm: time of PGV(m/s)
    % BinTrgAm: PGV(m/s)

    [Rorig,Amp,Mbefore,Mafter,t0Arrayorig,t1Arrayorig,t2Arrayorig]=CalculateR(lattrigger,lontrigger,ttrigger,Ntrig,BinTrgLon,BinTrgLat,BinTrgTm,BinTrgAm,IndexBinsStruct,Bins,latbin,lonbin,t,minsep,Ampmin,Ampmax,dtmax); % traveltime

else

    load R_Amp.mat;

end



%% output %%
%-- R and n distribution (not-smoothed) --
%    Rmed = nanmedian(Rorig);
%    nmed = solveR(Rmed);
%    npair = nansum(0.*Rorig+1.,1).';
%    T = table(lonbin.',latbin.',Rmed.',nmed.',nevbin.',npair); % R and n table
%    fileout_Rnamp = 'R_n_Amp.dat';
%    writetable(T,fileout_Rnamp,'Delimiter','\t');
%-- peak strain --
%    DateString = datestr(ttrigger,'yyyy-mm-ddTHH:MM:SS');
%    T2 = table(DateString,Amp); 
%    fileout_peak_strain = 'peak_strain.dat';
%    writetable(T2,fileout_peak_strain,'Delimiter','\t');
%-- seismicity before/after triggerer --
%    t0Array=reshape(t0Arrayorig,[],1);
%    Ampre=reshape(Amp,[],1);
%    t1Array = [reshape(t1Arrayorig,[],1),reshape(Mbefore,[],1),Ampre,t0Array];
%    t2Array = [reshape(t2Arrayorig,[],1),reshape(Mafter,[],1),Ampre,t0Array];
%    t0t1t2Array = [reshape(t1Arrayorig,[],1),reshape(Mbefore,[],1),reshape(t2Arrayorig,[],1),reshape(Mafter,[],1),Ampre,t0Array];
%-- t1 and t2 --
%    outArray = rmmissing(t0t1t2Array,1);
%    outDate = datestr(outArray(:,6),'yyyy-mm-ddTHH:MM:SS');
%    outTable_t0t1t2 = table(outArray(:,1:5),outDate);
%    fileID_t0t1t2 = 't0t1t2.dat';
%    writetable(outTable_t0t1t2,fileID_t0t1t2,'Delimiter','\t');
%-- sorted t1 and t2--
%    t1tmp = rmmissing(t1Array,1);
%    I1n = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
%    t1tmp = t1tmp*I1n;
%    [~, Indx] = sort(t1tmp(:,1));
%    t1ren = t1tmp(Indx,:);
%    t2tmp = rmmissing(t2Array,1);
%    [~, Indx] = sort(t2tmp(:,1));
%    t2re = t2tmp(Indx,:);
%    outArray = vertcat(t1ren,t2re);
%    outDate = datestr(outArray(:,4),'yyyy-mm-ddTHH:MM:SS');
%    outTable_t1t2_sort = table(outArray(:,1:3),outDate);
%    fileID_t1t2_all = 't1t2_mag.dat';
%    writetable(outTable_t1t2_sort,fileID_t1t2_all,'Delimiter','\t');


%% triggering intensity - strain relationship

    Nampbin0=50;
    iAmpBins=logspace(log10(1.5),log10(5.5),Nampbin0); 
    iAmpBins =log10(2*pi/(20*3.5e9)*10.^iAmpBins);  

    [RMeanraw,AmpMean,NRMeasurements,iAmpBins]=AverageRAmp(Rorig,Amp,iAmpBins,Nampbin0);

    Rb=RMeanraw*0 + 0.5; % no bias corect
    [n,nraw,nsimul,RMean,RMeanraw,RMeansimul]=ComputeRates(RMeanraw,Rb);

%%
% error calculation
    tic
    Nbstrp=1000; % number of R-array resampling for error calculations (need ~1000 for 90% CI)
    [n_mean,n_high,n_low]=CalculateErr(Nbstrp,Rorig,Amp,Rb,iAmpBins,Nampbin0);
    toc
%

%-- n-strain relationship --
%    T3 = table(AmpMean.',n.',(n-n_low).',(n_high-n).');
%    fileout_fig = 'figure_n_strain_2008-2017.dat';
%    writetable(T3,fileout_fig,'Delimiter','\t');

%%
%-- display results --
    figure(1)
    clf;
    errorbar(AmpMean,n,n-n_low,n_high-n,'r.');
    set(gca,'XScale','Log','YScale','Linear');
    xlabel('Strain');
    ylabel('n')
    title('90% Confidence Interval');
    grid on
    hold on;
%% FIT
    loglog(AmpMean,0.385*(1e6*AmpMean).^0.944,'k--')
    text(1e-7,3,'0.39 \epsilon^{0.94} from MBG 2021')
    hold off;
