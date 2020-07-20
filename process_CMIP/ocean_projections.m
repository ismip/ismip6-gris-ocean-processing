% script to generate CMIP5 ocean forcings
clear; close all;

% modelid = 1 for MIROC5 RCP8.5
% modelid = 2 for MIROC5 RCP2.6
% modelid = 3 for NorESM RCP8.5
% modelid = 4 for NorESM RCP2.6
% modelid = 5 for HadGEM RCP8.5
% modelid = 6 for HadGEM RCP2.6
% modelid = 7 for IPSLMR RCP8.5
% modelid = 8 for IPSLMR RCP2.6
% modelid = 9 for CSIRO RCP8.5
% modelid = 10 for CSIRO RCP2.6
% modelid = 11 for CNRM RCP8.5
% modelid = 12 for CNRM RCP2.6
% modelid = 13 for ACCESS RCP8.5 (note access has no RCP2.6)
% modelid = 14 for CNRM-CM6-1 ssp585
% modelid = 15 for CNRM-CM6-1 ssp126
% modelid = 16 for CNRM-ESM2-1 ssp585
% modelid = 17 for UKESM1-0-LL ssp585
% modelid = 18 for CESM2 ssp585
modelid = 18;

% freezing point parameters
l1 = -5.73e-2;
l2 = 8.32e-2;
l3 = -7.53e-4;

% regular grid
dd = 50; % grid resolution in km
xreg = [-3e6:dd*10^3:3e6];
yreg = [-4e6:dd*10^3:0e6];
zreg = [200:25:500];
[Xreg,Yreg] = meshgrid(xreg,yreg);

%% MIROC5 RCP8.5

if modelid == 1,

% path to MIROC5 RCP8.5 temperature output
Tfiletoload = '~/Documents/CMIP5/M5histrcp.nc';
% path to MIROC5 RCP8.5 salinity output
Sfiletoload = '~/Documents/CMIP5/M5histrcp_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');
lon_vert=ncread(Tfiletoload,'lon_vertices');
lat_vert=ncread(Tfiletoload,'lat_vertices');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[1850:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin with NaNs set to 0
% change 0s to NaNs
T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity has NaNs set to 0 so change these
S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of MIROC5 points on BedMachine grid
% gridcell vertices
[x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).MIROC5inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).MIROC5inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).MIROC5inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save MIROC5_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save MIROC5_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% MIROC5 RCP2.6

if modelid == 2,

% path to MIROC5 RCP2.6temperature output
Tfiletoload = '~/Desktop/M5rcp26_thetao.nc';
% path to MIROC5 RCP2.6 salinity output
Sfiletoload = '~/Desktop/M5rcp26_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');
lon_vert=ncread(Tfiletoload,'lon_vertices');
lat_vert=ncread(Tfiletoload,'lat_vertices');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[2006:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin with NaNs set to 0
% change 0s to NaNs
T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity has NaNs set to 0 so change these
S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of MIROC5 points on BedMachine grid
% gridcell vertices
[x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).MIROC5inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).MIROC5inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).MIROC5inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% add historical part to front of run
load MIROC5_historical.mat
year = [year_historical,year];
TF_basins = cat(2,TF_basins_historical,TF_basins);
TF = cat(2,TF_historical,TF);

% save
save MIROC5_RCP26.mat x y TF TF_basins year

end

%% NorESM RCP8.5

if modelid == 3,

% path to NorESM RCP8.5 temperature output
Tfiletoload = '~/Desktop/NorMhistrcp85_thetao.nc';
% path to NorESM RCP8.5 salinity output
Sfiletoload = '~/Desktop/NorMhistrcp85_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');
lon_vert=ncread(Tfiletoload,'lon_vertices');
lat_vert=ncread(Tfiletoload,'lat_vertices');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[1850:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin with NaNs set to 0
% change 0s to NaNs
T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity has NaNs set to 0 so change these
S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of NorESM points on BedMachine grid
% gridcell vertices
[x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).NorESMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).NorESMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).NorESMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save NorESM_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save NorESM_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% NorESM RCP2.6

if modelid == 4,

% path to NorESM RCP2.6 temperature output
Tfiletoload = '~/Desktop/NorMrcp26_thetao.nc';
% path to NorESM RCP2.6 salinity output
Sfiletoload = '~/Desktop/NorMrcp26_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');
lon_vert=ncread(Tfiletoload,'lon_vertices');
lat_vert=ncread(Tfiletoload,'lat_vertices');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[2006:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin with NaNs set to 0
% change 0s to NaNs
T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity has NaNs set to 0 so change these
S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of NorESM points on BedMachine grid
% gridcell vertices
[x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).NorESMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).NorESMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).NorESMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% add historical part to front of run
load NorESM_historical.mat
year = [year_historical,year];
TF_basins = cat(2,TF_basins_historical,TF_basins);
TF = cat(2,TF_historical,TF);

% save
save NorESM_RCP26.mat x y TF TF_basins year

end

%% HadGEM RCP8.5
% note processing for HadGEM is a bit different to the others
% because I had to download the output myself

if modelid == 5,

% path to HadGEM temperature output
Tfiletoload1 = '~/Documents/HadGEM-ES2/HadGEM2-ES-hist-thetao.mat';
Tfiletoload2 = '~/Documents/HadGEM-ES2/HadGEM2-ES-RCP85-thetao.mat';

% path to HadGEM salinity output
Sfiletoload1 = '~/Documents/HadGEM-ES2/HadGEM2-ES-hist-so.mat';
Sfiletoload2 = '~/Documents/HadGEM-ES2/HadGEM2-ES-RCP85-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of HadGEM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).HadGEMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).HadGEMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).HadGEMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save HadGEM_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save HadGEM_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% HadGEM RCP2.6
% note processing for HadGEM is a bit different to the others
% because I had to download the output myself

if modelid == 6,

% path to HadGEM temperature output
Tfiletoload1 = '~/Documents/HadGEM-ES2/HadGEM2-ES-hist-thetao.mat';
Tfiletoload2 = '~/Documents/HadGEM-ES2/HadGEM2-ES-RCP26-thetao.mat';

% path to HadGEM salinity output
Sfiletoload1 = '~/Documents/HadGEM-ES2/HadGEM2-ES-hist-so.mat';
Sfiletoload2 = '~/Documents/HadGEM-ES2/HadGEM2-ES-RCP26-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of HadGEM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).HadGEMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).HadGEMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).HadGEMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save HadGEM_RCP26.mat x y TF TF_basins year

end

%% IPSLCM RCP8.5

if modelid == 7,

% path to IPSLMR RCP8.5 temperature output
Tfiletoload = '~/Desktop/IPSLMRhistrcp85_thetao.nc';
% path to IPSLMR RCP8.5 salinity output
Sfiletoload = '~/Desktop/IPSLMRhistrcp85_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[1850:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of IPSLCM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).IPSLCMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).IPSLCMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).IPSLCMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save IPSLCM_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save IPSLCM_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% IPSLCM RCP2.6

if modelid == 8,

% path to IPSLCM RCP2.6temperature output
Tfiletoload = '~/Desktop/IPSLMRrcp26_thetao.nc';
% path to IPSLCM RCP2.6 salinity output
Sfiletoload = '~/Desktop/IPSLMRrcp26_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[2006:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of IPSLCM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).IPSLCMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).IPSLCMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).IPSLCMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% add historical part to front of run
load IPSLCM_historical.mat
year = [year_historical,year];
TF_basins = cat(2,TF_basins_historical,TF_basins);
TF = cat(2,TF_historical,TF);

% save
save IPSLCM_RCP26.mat x y TF TF_basins year

end

%% CSIRO RCP8.5

if modelid == 9,

% path to CSIRO RCP8.5 temperature output
Tfiletoload = '~/Desktop/CSIROhistrcp85_thetao.nc';
% path to CSIRO RCP8.5 salinity output
Sfiletoload = '~/Desktop/CSIROhistrcp85_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[1850:2100];

% meshgrid lon and lat
[lat,lon] = meshgrid(lat,lon);

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CSIRO points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CSIROinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CSIROinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CSIROinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CSIRO_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save CSIRO_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% CSIRO RCP2.6

if modelid == 10,

% path to CSIRO RCP2.6temperature output
Tfiletoload = '~/Desktop/CSIROhistrcp26_thetao.nc';
% path to CSIRO RCP2.6 salinity output
Sfiletoload = '~/Desktop/CSIROhistrcp26_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[2006:2100];

% meshgrid lon and lat
[lat,lon] = meshgrid(lat,lon);

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CSIRO points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CSIROinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CSIROinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CSIROinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% add historical part to front of run
load CSIRO_historical.mat
year = [year_historical,year];
TF_basins = cat(2,TF_basins_historical,TF_basins);
TF = cat(2,TF_historical,TF);

% save
save CSIRO_RCP26.mat x y TF TF_basins year

end

%% CNRM RCP8.5

if modelid == 11,

% path to CNRM RCP8.5 temperature output
Tfiletoload = '~/Desktop/CNRMhistrcp85_thetao.nc';
% path to CNRM RCP8.5 salinity output
Sfiletoload = '~/Desktop/CNRMhistrcp85_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[1850:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CNRMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CNRMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CNRMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CNRM_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save CNRM_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% CNRM RCP2.6

if modelid == 12,

% path to CNRM RCP2.6temperature output
Tfiletoload = '~/Desktop/CNRMrcp26_thetao.nc';
% path to CNRM RCP2.6 salinity output
Sfiletoload = '~/Desktop/CNRMrcp26_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');

% read only TS between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),1],[Inf,Inf,length(depthinds),Inf]);

% time field is a bit wierd - make own array
year=[2006:2100];

% reshape
T = reshape(T,size(T,1)*size(T,2),size(T,3),size(T,4));
S = reshape(S,size(S,1)*size(S,2),size(S,3),size(S,4));
lon = lon(:);
lat = lat(:);

% trim to area of interest
% spatial trim
inds1 = find(lon>270 & lat>55);
inds2 = find(lon<10 & lat>55);
spatialinds = [inds1;inds2];
T = squeeze(T(spatialinds,:,:));
S = squeeze(S(spatialinds,:,:));
lon = lon(spatialinds);
lat = lat(spatialinds);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CNRMinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CNRMinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CNRMinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% add historical part to front of run
load CNRM_historical.mat
year = [year_historical,year];
TF_basins = cat(2,TF_basins_historical,TF_basins);
TF = cat(2,TF_historical,TF);

% save
save CNRM_RCP26.mat x y TF TF_basins year

end

%% ACCESS RCP8.5
% note processing for ACCESS is a bit different to the others
% because I had to download the output myself

if modelid == 13,

% path to ACCESS temperature output
Tfiletoload1 = '~/Documents/ACCESS1-3/ACCESS1-3-hist-thetao.mat';
Tfiletoload2 = '~/Documents/ACCESS1-3/ACCESS1-3-RCP85-thetao.mat';

% path to ACCESS salinity output
Sfiletoload1 = '~/Documents/ACCESS1-3/ACCESS1-3-hist-so.mat';
Sfiletoload2 = '~/Documents/ACCESS1-3/ACCESS1-3-RCP85-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of ACCESS points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).ACCESSinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).ACCESSinds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).ACCESSinds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save ACCESS_RCP85.mat x y TF TF_basins year

% also save historical part of run
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save ACCESS_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% CNRM-CM6-1 ssp585
% note processing for CNRM-CM6-1 is a bit different to the others
% because I had to download the output myself and also because
% this is a CMIP6 model so the historical period extends to 2014

if modelid == 14,

% path to CNRM-CM6-1 temperature output
Tfiletoload1 = '~/Documents/CNRM-CM6/CNRM-CM6-1-hist-thetao.mat';
Tfiletoload2 = '~/Documents/CNRM-CM6/CNRM-CM6-1-ssp585-thetao.mat';

% path to CNRM-CM6-1 salinity output
Sfiletoload1 = '~/Documents/CNRM-CM6/CNRM-CM6-1-hist-so.mat';
Sfiletoload2 = '~/Documents/CNRM-CM6/CNRM-CM6-1-ssp585-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
% T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM-CM6-1 points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CNRMCM6inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CNRMCM6inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CNRMCM6inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CNRM-CM6-1_ssp585.mat x y TF TF_basins year

% also save historical part of run
% note define historical as <2006 for consistency with CMIP5 runs
% in actuality historical extends to 2014 in CMIP6
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save CNRM-CM6-1_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% CNRM-CM6-1 ssp126
% note processing for CNRM-CM6-1 is a bit different to the others
% because I had to download the output myself and also because
% this is a CMIP6 model so the historical period extends to 2014

if modelid == 15,

% path to CNRM-CM6-1 temperature output
Tfiletoload1 = '~/Documents/CNRM-CM6/CNRM-CM6-1-hist-thetao.mat';
Tfiletoload2 = '~/Documents/CNRM-CM6/CNRM-CM6-1-ssp126-thetao.mat';

% path to CNRM-CM6-1 salinity output
Sfiletoload1 = '~/Documents/CNRM-CM6/CNRM-CM6-1-hist-so.mat';
Sfiletoload2 = '~/Documents/CNRM-CM6/CNRM-CM6-1-ssp126-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
% T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM-CM6-1 points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CNRMCM6inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CNRMCM6inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CNRMCM6inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CNRM-CM6-1_ssp126.mat x y TF TF_basins year

end

%% CNRM-ESM2-1 ssp585
% note processing for CNRM-ESM2-1 is a bit different to the others
% because I had to download the output myself and also because
% this is a CMIP6 model so the historical period extends to 2014
if modelid == 16,

% path to CNRM-ESM2-1 temperature output
Tfiletoload1 = '~/Documents/CNRM-ESM2-1/CNRM-ESM2-1-hist-thetao.mat';
Tfiletoload2 = '~/Documents/CNRM-ESM2-1/CNRM-ESM2-1-ssp585-thetao.mat';

% path to CNRM-ESM2-1 salinity output
Sfiletoload1 = '~/Documents/CNRM-ESM2-1/CNRM-ESM2-1-hist-so.mat';
Sfiletoload2 = '~/Documents/CNRM-ESM2-1/CNRM-ESM2-1-ssp585-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
% T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM-CM6-1 points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CNRMESM2inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CNRMESM2inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CNRMESM2inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CNRM-ESM2-1_ssp585.mat x y TF TF_basins year

% also save historical part of run
% note define historical as <2006 for consistency with CMIP5 runs
% in actuality historical extends to 2014 in CMIP6
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save CNRM-ESM2-1_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% UKESM1-0-LL ssp585
% note processing for UKESM1-0-LL is a bit different to the others
% because I had to download the output myself and also because
% this is a CMIP6 model so the historical period extends to 2014
if modelid == 17,

% path to UKESM1-0-LL temperature output
Tfiletoload1 = '~/Documents/UKESM1-0-LL/UKESM1-0-LL-hist-thetao.mat';
Tfiletoload2 = '~/Documents/UKESM1-0-LL/UKESM1-0-LL-ssp585-thetao.mat';

% path to UKESM1-0-LL salinity output
Sfiletoload1 = '~/Documents/UKESM1-0-LL/UKESM1-0-LL-hist-so.mat';
Sfiletoload2 = '~/Documents/UKESM1-0-LL/UKESM1-0-LL-ssp585-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
% T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM-CM6-1 points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).UKESM1inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).UKESM1inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).UKESM1inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save UKESM1-0-LL_ssp585.mat x y TF TF_basins year

% also save historical part of run
% note define historical as <2006 for consistency with CMIP5 runs
% in actuality historical extends to 2014 in CMIP6
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save UKESM1-0-LL_historical.mat x y TF_historical TF_basins_historical year_historical

end

%% CESM2 ssp585
% note processing for CESM2 is a bit different to the others
% because I had to download the output myself and also because
% this is a CMIP6 model so the historical period extends to 2014
if modelid == 18,

% path to CESM2 temperature output
Tfiletoload1 = '~/Documents/CESM2/CESM2-hist-thetao.mat';
Tfiletoload2 = '~/Documents/CESM2/CESM2-ssp585-thetao.mat';

% path to CESM2 salinity output
Sfiletoload1 = '~/Documents/CESM2/CESM2-hist-so.mat';
Sfiletoload2 = '~/Documents/CESM2/CESM2-ssp585-so.mat';

% read fields
T1 = load(Tfiletoload1);
T2 = load(Tfiletoload2);
S1 = load(Sfiletoload1);
S2 = load(Sfiletoload2);

% merge fields
S = cat(3,S1.S,S2.S);
T = cat(3,T1.T,T2.T);
z = S1.z;
year = [S1.yrs,S2.yrs];
lat = S1.lat;
lon = S1.lon;

% restrict to depths between 200 and 500 m and one point either side
depthinds = [max(find(z<200)):min(find(z>500))];
z = z(depthinds);
T = squeeze(T(:,depthinds,:));
S = squeeze(S(:,depthinds,:));

% check all variables are double
T = double(T);
S = double(S);
lon = double(lon);
lat = double(lat);

% there are some weird duplicates - remove these
[~,lon_inds] = unique(lon);
[~,lat_inds] = unique(lat);
lon_dups = setdiff(1:length(lon),lon_inds);
lat_dups = setdiff(1:length(lat),lat_inds);
dups_inds = intersect(lon_dups,lat_dups);
T(dups_inds,:,:) = [];
S(dups_inds,:,:) = [];
lon(dups_inds) = [];
lat(dups_inds) = [];

% temperature is in kelvin
% NaNs are already NaNs
% T(find(T==0)) = NaN;
% change to celsius
% T = T-273.15;
% salinity NaNs are already NaNs
% S(find(S==0)) = NaN;

% calculate thermal forcing
depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
TF = T - (l1*S+l2+l3*depth);

% get coords of CNRM-CM6-1 points on BedMachine grid
% gridcell vertices
% [x_vert,y_vert] = latlon2utm(lat_vert,lon_vert);
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate TF onto regular grid
for jj=1:size(TF,2),
    for kk=1:size(TF,3),
        TF0 = squeeze(TF(:,jj,kk));
        f = scatteredInterpolant(x,y,TF0,'linear','none');
        TFreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../../ismip6/final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CESM2inds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CESM2inds) = k;
end

% get thermal forcing profile per basin per year
TF_basins_z = NaN(7,length(z),length(year));
TFreg_vec = reshape(TFreg,size(TFreg,1)*size(TFreg,2),size(TFreg,3),size(TFreg,4));
for jj=1:7,
    for kk=1:length(z),
        TF_basins_z(jj,kk,:) = nanmean(TFreg_vec(regions(jj).CESM2inds,kk,:),1);
    end
end

% take mean of thermal forcing over 200-500 m
% first interpolate onto regular grid then take mean
for jj=1:7,
    for kk=1:length(year),
        TF_basins(jj,kk) = nanmean(interp1(z,TF_basins_z(jj,:,kk),zreg,'linear',NaN));
    end
end

% create rough spatial map of TF by taking naive depth mean
TF = nanmean(TF,2);
TF = squeeze(TF(:,1,:));

% save
save CESM2_ssp585.mat x y TF TF_basins year

% also save historical part of run
% note define historical as <2006 for consistency with CMIP5 runs
% in actuality historical extends to 2014 in CMIP6
historical_inds = find(year<2006);
year_historical = year(historical_inds);
TF_historical = TF(:,historical_inds);
TF_basins_historical = TF_basins(:,historical_inds);
save CESM2_historical.mat x y TF_historical TF_basins_historical year_historical

end