function [t,M,lat,lon,depth]= ReadQTM(F)

    id=fopen(F)
    A=textscan(id,'%s %s %s %s %s %s %d %f %f %f %f %*[^\n]','headerlines',1);
    fclose(id);
    
    lat=A{8};
    lon=A{9};
    depth=A{10};
    M=A{11};
    datetimestring=strcat(A{1},{'/'},A{2},{'/'},A{3},A{4},{':'},A{5},{':'},A{6});
    anssDateFormat='yyyy/mm/ddHH:MM:SS.FFF';
    t=datenum(datetimestring, anssDateFormat);
    
    
    
