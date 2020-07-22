%% THIS EXAMPLE CONTAINS ONLY NORESM

% script to generate forcing files for extrapolation into fjords
% v2 does both temperature and salinity
clear; close all;

% freezing point parameters
% l1 = -5.73e-2;
% l2 = 8.32e-2;
% l3 = -7.53e-4;

% regular grid
dd = 50; % grid resolution in km
xreg = [-3e6:dd*10^3:3e6];
yreg = [-4e6:dd*10^3:0e6];
zreg = [0:50:2000];
[Xreg,Yreg] = meshgrid(xreg,yreg);

forcingid = 1; % 1 for RCP8.5

%% NorESM RCP8.5
if forcingid == 1,

% path to NorESM RCP8.5 temperature output
Tfiletoload = '~/Documents/NorESM/NorMhistrcp85_thetao.nc';
% path to NorESM RCP8.5 salinity output
Sfiletoload = '~/Documents/NorESM/NorMhistrcp85_so.nc';

% read fields
lon=ncread(Tfiletoload,'lon');
lat=ncread(Tfiletoload,'lat');
z=ncread(Tfiletoload,'lev');
time=ncread(Tfiletoload,'time');
lon_vert=ncread(Tfiletoload,'lon_vertices');
lat_vert=ncread(Tfiletoload,'lat_vertices');

% read only TS to 2000 m and one point either side
% also read only from 1950 to save space
year=[1850:2100];
timeinds = find(year>=1950);
year = year(timeinds);
depthinds = [1:max(find(z<2000))+1];
z = z(depthinds);
T=ncread(Tfiletoload,'thetao',[1,1,depthinds(1),timeinds(1)],[Inf,Inf,length(depthinds),Inf]);
S=ncread(Sfiletoload,'so',[1,1,depthinds(1),timeinds(1)],[Inf,Inf,length(depthinds),Inf]);

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
% depth = permute(repmat(z,1,size(S,1),length(year)),[2,1,3]);
% TF = T - (l1*S+l2+l3*depth);

% clear irrelevant variables to speed things up
% clearvars T S;

% get coords of NorESM points on BedMachine grid
% gridcell centres
[x,y] = latlon2utm(lat,lon);

% interpolate T/S onto regular grid
for jj=1:size(T,2),
    for kk=1:size(T,3),
        T0 = squeeze(T(:,jj,kk));
        f = scatteredInterpolant(x,y,T0,'linear','none');
        Treg(:,:,jj,kk) = f(Xreg,Yreg);
        S0 = squeeze(S(:,jj,kk));
        f = scatteredInterpolant(x,y,S0,'linear','none');
        Sreg(:,:,jj,kk) = f(Xreg,Yreg);
    end
end

% load ice-ocean basin definitions
load ../final_region_def/ice_ocean_sectors.mat

% find regular grid points inside ice-ocean basins
basin = NaN(size(Xreg,1),size(Xreg,2));
for k=1:7, % since there are 7 basins
    regions(k).CMIPinds = find(inpolygon(Xreg,Yreg,regions(k).ocean.x,regions(k).ocean.y));
    basin(regions(k).CMIPinds) = k;
end

% get T/S profiles per basin per year
T_basins_z = NaN(7,length(z),length(year));
Treg_vec = reshape(Treg,size(Treg,1)*size(Treg,2),size(Treg,3),size(Treg,4));
S_basins_z = NaN(7,length(z),length(year));
Sreg_vec = reshape(Sreg,size(Sreg,1)*size(Sreg,2),size(Sreg,3),size(Sreg,4));
for jj=1:7,
    for kk=1:length(z),
        T_basins_z(jj,kk,:) = nanmean(Treg_vec(regions(jj).CMIPinds,kk,:),1);
        S_basins_z(jj,kk,:) = nanmean(Sreg_vec(regions(jj).CMIPinds,kk,:),1);
    end
end

% interpolate onto regular vertical grid
for jj=1:7,
    for kk=1:length(year),
        T_basins_reg(jj,:,kk) = interp1(z,T_basins_z(jj,:,kk),zreg,'linear','extrap');
        S_basins_reg(jj,:,kk) = interp1(z,S_basins_z(jj,:,kk),zreg,'linear','extrap');
    end
end

% check if any basin has NaNs - fill these with nearest non-NaN properties
for jj=1:7,
    if ~isempty(find(isnan(T_basins_reg(jj,:,:)))),
        for kk=1:length(year),
            T_basins_reg(jj,find(isnan(T_basins_reg(jj,:,kk))),kk) = T_basins_reg(jj,min(find(isnan(T_basins_reg(jj,:,kk))))-1,kk);
        end
    end
    if ~isempty(find(isnan(S_basins_reg(jj,:,:)))),
        for kk=1:length(year),
            S_basins_reg(jj,find(isnan(S_basins_reg(jj,:,kk))),kk) = S_basins_reg(jj,min(find(isnan(S_basins_reg(jj,:,kk))))-1,kk);
        end
    end
end

% bias correction based on mean 200-500m properties over 1995-2014
% load EN4
load('EN4_ISMIP6_highres.mat');
% project EN4 onto regular z grid
for i=1:7,
    Tbaseline(i,:) = interp1(regions(i).depth,regions(i).meansectorT,zreg,'linear','extrap');
    Sbaseline(i,:) = interp1(regions(i).depth,regions(i).meansectorS,zreg,'linear','extrap');
end
inds200500 = find(zreg>=200 & zreg<=500);
inds19952014 = find(year>=1995 & year<=2014);
for i=1:7,
    Tanomaly(i) = nanmean(nanmean(T_basins_reg(i,inds200500,inds19952014))) - nanmean(Tbaseline(i,inds200500));
    Sanomaly(i) = nanmean(nanmean(S_basins_reg(i,inds200500,inds19952014))) - nanmean(Sbaseline(i,inds200500));
end
% do bias correction
for i=1:7,
    T_basins_reg(i,:,:) = T_basins_reg(i,:,:) - Tanomaly(i);
    S_basins_reg(i,:,:) = S_basins_reg(i,:,:) - Sanomaly(i);
end

% rename for saving with same names as before
% and format exactly as before
z = -zreg;
T = permute(T_basins_reg,[1,3,2]);
S = permute(S_basins_reg,[1,3,2]);

% make some plots
% potential T
figure();
cols = jet(length(year));
region_names = {'SE','SW','CE','CW','NE','NW','NO'};
fs = 8;
lw = 0.25;
for k=1:7,
    subplot(2,4,k); hold on;
    for ii=1:length(year),
        plot(squeeze(T(k,ii,:)),z,'color',cols(ii,:),'linewidth',lw);
    end
    title(region_names{k},'fontsize',fs,'interpreter','latex');
    xlabel('potential temperature ($^{\circ}$C)','fontsize',fs,'interpreter','latex');
    ylabel('depth (m)','fontsize',fs,'interpreter','latex');
    set(gca,'box','on');
    ylim([-2000 0]);
end
h = colorbar('position',[0.8,0.1,0.02,0.25]); colormap(jet); caxis([1950 2100]);
% print
fw = 20;
fh = 10;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 fw fh]);
set(gcf,'PaperSize',[fw fh]);
set(gcf,'units','centimeters','position',[0 0 fw fh]);
set(gcf,'color','w','inverthardcopy','off');
print -dpng -r300 potentialT_RCP85.png
close all;

% practical salinity
figure();
cols = jet(length(year));
region_names = {'SE','SW','CE','CW','NE','NW','NO'};
fs = 8;
lw = 0.25;
for k=1:7,
    subplot(2,4,k); hold on;
    for ii=1:length(year),
        plot(squeeze(S(k,ii,:)),z,'color',cols(ii,:),'linewidth',lw);
    end
    title(region_names{k},'fontsize',fs,'interpreter','latex');
    xlabel('practical salinity (psu)','fontsize',fs,'interpreter','latex');
    ylabel('depth (m)','fontsize',fs,'interpreter','latex');
    set(gca,'box','on');
    ylim([-2000 0]);
end
h = colorbar('position',[0.8,0.1,0.02,0.25]); colormap(jet); caxis([1950 2100]);
% print
fw = 20;
fh = 10;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 fw fh]);
set(gcf,'PaperSize',[fw fh]);
set(gcf,'units','centimeters','position',[0 0 fw fh]);
set(gcf,'color','w','inverthardcopy','off');
print -dpng -r300 practicalS_RCP85.png
close all;

% save relevant variables
% rename for consistency with previous file
for i=1:7,
    basins(i).X = regions(i).ice.x;
    basins(i).Y = regions(i).ice.y;
end
save ocean_extrap_NorESM_RCP85.mat z T S basins year

end
