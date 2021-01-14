%clearvars -except tmin tmax mth;
function [t, lat, lon, mag, depth]=ReadCat(fn,tmin,tmax,mth,depthmin,depthmax,magmax);

[t,mag,lat,lon,depth]=ReadQTM(fn);
latmax=42;
latmin=32;
lonmin=-124;
lonmax=-114;

It=(t>=tmin)&(t<=tmax);
Ispace=(lat>=latmin)&(lat<=latmax)&(lon<=lonmax)&(lon>=lonmin);
Icomplete=(mag>=mth)&(depth<=depthmax)&Ispace&It&(depth>depthmin)&(mag<magmax);
mag=mag(Icomplete);
lat=lat(Icomplete);
lon=lon(Icomplete);
t=t(Icomplete);
depth=depth(Icomplete);

[t Isort]=sort(t);
lat=lat(Isort);
lon=lon(Isort);
mag=mag(Isort);
depth=depth(Isort);
