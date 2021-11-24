clear; close all;
% script to bring together all terminus position records
% and add past runoff, ocean thermal forcing and melt
% for the purpose of calculating kappa distribution
%% Helheim data from Andresen et al. 2011, Nature Geoscience

% comes from excel sheet emailed by Camilla
A = xlsread('andresen/Sermilik Fjord data.xlsx','Glacier length data');
glaciers(1).termpos.L = flipud(A(:,5));
% in absence of more info, assume position refers to middle of year
glaciers(1).termpos.t = flipud(A(:,1))+0.5;
glaciers(1).bjorkid = 'GGN0163';
glaciers(1).rignotid = 3;
glaciers(1).termpos.source = 'andresen';
glaciers(1).lon = -38.3067;
glaciers(1).lat = 66.3735;
glaciers(1).x = [];
glaciers(1).y = [];
%% Jakobshavn data from Steiger et al. 2018, Cryosphere

% sent by email from Nadine
aaa = xlsread('steiger/Retreat.xlsx');
glaciers(2).termpos.t = aaa(:,1);
glaciers(2).termpos.L = aaa(:,3);
glaciers(2).bjorkid = 'GGN0235';
glaciers(2).rignotid = 1;
glaciers(2).termpos.source = 'steiger';
glaciers(2).lat = 69.1898;
glaciers(2).lon = -49.4542;
%% KNS data from Lea et al. 2014, Cryosphere

glaciers(3).termpos.L = [22657.89,17542,8563,8563,10188,8434,7874,7299,3420,...
        1487,2143,1978,1485,2243,2331,2358,1733,1959,1975,1401,1334,1319,...
        1218,1187,749,433,381,145,358,416,238,137,0]/1000;
glaciers(3).termpos.t = [1761,1808,1860,1903,1921,1932,1936,1946,1948,1968,1979,...
        1985,1987,1992,1993,1994,1995,1996,1997,1999,2000,2001,2002,2003,2004,...
        2005,2006,2007,2008,2009,2010,2011,2012]+0.5;
glaciers(3).bjorkid = 'GGN0094';
glaciers(3).rignotid = 36;
glaciers(3).termpos.source = 'lea';
glaciers(3).lat = 64.2966;
glaciers(3).lon = -49.6102;
%% Upernavik data from Haubner et al. 2018, Cryosphere

load haubner/upernavik_retreat_processed.mat
% has data from the three streams; originally these were 1 so make them agree at the start
glaciers(4).termpos.t = t_upernavik+0.5;
glaciers(4).termpos.L = p_upernavik1 - p_upernavik1(1) + 20;
glaciers(4).bjorkid = ''; % not in bjork database
glaciers(4).rignotid = 26;
glaciers(4).termpos.source = 'haubner';
glaciers(4).lat = 73.0449;
glaciers(4).lon = -54.2698;

glaciers(5).termpos.t = t_upernavik+0.5;
glaciers(5).termpos.L = p_upernavik2 - p_upernavik2(1) + 20;
glaciers(5).bjorkid = ''; % not in bjork database
glaciers(5).rignotid = 18;
glaciers(5).termpos.source = 'haubner';
glaciers(5).lat = 72.9451;
glaciers(5).lon = -54.1842;

glaciers(6).termpos.t = t_upernavik+0.5;
glaciers(6).termpos.L = p_upernavik3 - p_upernavik3(1) + 20;
glaciers(6).bjorkid = 'GGN0452';
glaciers(6).rignotid = 22;
glaciers(6).termpos.source = 'haubner';
glaciers(6).lat = 72.8409;
glaciers(6).lon = -54.2837;
%% Catania et al 2018, JGR Earth Surface

load catania/glaciertermini_18-Apr-2018.mat
% ava
glaciers(7).termpos.t = ava.dyear;
glaciers(7).termpos.L = ava.T/1000;
glaciers(7).bjorkid = 'GGN0265';
glaciers(7).rignotid = 53;
glaciers(7).termpos.source = 'catania';
[glaciers(7).lat,glaciers(7).lon] = utm2deg(nanmean(ava.x{end}),nanmean(ava.y{end}),'22 W');
% eqi
glaciers(8).termpos.t = eqi.dyear;
glaciers(8).termpos.L = eqi.T/1000;
glaciers(8).bjorkid = 'GGN0250';
glaciers(8).rignotid = 90;
glaciers(8).termpos.source = 'catania';
[glaciers(8).lat,glaciers(8).lon] = utm2deg(nanmean(eqi.x{end}),nanmean(eqi.y{end}),'22 W');
% ing
glaciers(9).termpos.t = ing.dyear;
glaciers(9).termpos.L = ing.T/1000;
glaciers(9).bjorkid = 'GGN0390';
glaciers(9).rignotid = 88;
glaciers(9).termpos.source = 'catania';
[glaciers(9).lat,glaciers(9).lon] = utm2deg(nanmean(ing.x{end}),nanmean(ing.y{end}),'22 W');
% kan
glaciers(10).termpos.t = kan.dyear;
glaciers(10).termpos.L = kan.T/1000;
glaciers(10).bjorkid = 'GGN0253';
glaciers(10).rignotid = 52;
glaciers(10).termpos.source = 'catania';
[glaciers(10).lat,glaciers(10).lon] = utm2deg(nanmean(kan.x{end}),nanmean(kan.y{end}),'22 W');
% kas
glaciers(11).termpos.t = kas.dyear;
glaciers(11).termpos.L = kas.T/1000;
glaciers(11).bjorkid = 'GGN0326';
glaciers(11).rignotid = 47;
glaciers(11).termpos.source = 'catania';
[glaciers(11).lat,glaciers(11).lon] = utm2deg(nanmean(kas.x{end}),nanmean(kas.y{end}),'22 W');
% kng
glaciers(12).termpos.t = kng.dyear;
glaciers(12).termpos.L = kng.T/1000;
glaciers(12).bjorkid = 'GGN0295';
glaciers(12).rignotid = 70;
glaciers(12).termpos.source = 'catania';
[glaciers(12).lat,glaciers(12).lon] = utm2deg(nanmean(kng.x{end}),nanmean(kng.y{end}),'22 W');
% kss
glaciers(13).termpos.t = kss.dyear;
glaciers(13).termpos.L = kss.T/1000;
glaciers(13).bjorkid = 'GGN0313';
glaciers(13).rignotid = 191;
glaciers(13).termpos.source = 'catania';
[glaciers(13).lat,glaciers(13).lon] = utm2deg(nanmean(kss.x{end}),nanmean(kss.y{end}),'22 W');
% kuj
glaciers(14).termpos.t = kuj.dyear;
glaciers(14).termpos.L = kuj.T/1000;
glaciers(14).bjorkid = 'GGN0256';
glaciers(14).rignotid = 13;
glaciers(14).termpos.source = 'catania';
[glaciers(14).lat,glaciers(14).lon] = utm2deg(nanmean(kuj.x{end}),nanmean(kuj.y{end}),'22 W');
% lik
glaciers(15).termpos.t = lik.dyear;
glaciers(15).termpos.L = lik.T/1000;
glaciers(15).bjorkid = 'GGN0290';
glaciers(15).rignotid = 119;
glaciers(15).termpos.source = 'catania';
[glaciers(15).lat,glaciers(15).lon] = utm2deg(nanmean(lik.x{end}),nanmean(lik.y{end}),'22 W');
% lil
glaciers(16).termpos.t = lil.dyear;
glaciers(16).termpos.L = lil.T/1000;
glaciers(16).bjorkid = 'GGN0287';
glaciers(16).rignotid = 134;
glaciers(16).termpos.source = 'catania';
[glaciers(16).lat,glaciers(16).lon] = utm2deg(nanmean(lil.x{end}),nanmean(lil.y{end}),'22 W');
% prd
glaciers(17).termpos.t = prd.dyear;
glaciers(17).termpos.L = prd.T/1000;
glaciers(17).bjorkid = 'GGN0303';
glaciers(17).rignotid = NaN;
glaciers(17).termpos.source = 'catania';
[glaciers(17).lat,glaciers(17).lon] = utm2deg(nanmean(prd.x{end}),nanmean(prd.y{end}),'22 W');
% rnk
glaciers(18).termpos.t = rnk.dyear;
glaciers(18).termpos.L = rnk.T/1000;
glaciers(18).bjorkid = 'GGN0347';
glaciers(18).rignotid = 5;
glaciers(18).termpos.source = 'catania';
[glaciers(18).lat,glaciers(18).lon] = utm2deg(nanmean(rnk.x{end}),nanmean(rnk.y{end}),'22 W');
% sil
glaciers(19).termpos.t = sil.dyear;
glaciers(19).termpos.L = sil.T/1000;
glaciers(19).bjorkid = 'GGN0299';
glaciers(19).rignotid = 25;
glaciers(19).termpos.source = 'catania';
[glaciers(19).lat,glaciers(19).lon] = utm2deg(nanmean(sil.x{end}),nanmean(sil.y{end}),'22 W');
% str
glaciers(20).termpos.t = str.dyear;
glaciers(20).termpos.L = str.T/1000;
glaciers(20).bjorkid = 'GGN0279';
glaciers(20).rignotid = 6;
glaciers(20).termpos.source = 'catania';
[glaciers(20).lat,glaciers(20).lon] = utm2deg(nanmean(str.x{end}),nanmean(str.y{end}),'22 W');
% umi
glaciers(21).termpos.t = umi.dyear;
glaciers(21).termpos.L = umi.T/1000;
glaciers(21).bjorkid = 'GGN0348';
glaciers(21).rignotid = 68;
glaciers(21).termpos.source = 'catania';
[glaciers(21).lat,glaciers(21).lon] = utm2deg(nanmean(umi.x{end}),nanmean(umi.y{end}),'22 W');
%% Cowton et al. 2018, PNAS

ii0 = length(glaciers);
cowtonid = {'M3','T1','AB','HG','KG',...
            'BG','VG','DJ','WG','HK'};
bjorkid = {'','','GGN0089','GGN0163','GGN0218',...
           'GGN0223','GGN0275','GGN0345','GGN0502','GGN0541'};
rignotid = [31,15,12,3,2,91,48,8,215,226];
cowtonlat = [62.52,62.76,63.82,66.35,68.59,...
             68.60,70.38,71.90,73.83,75.16];
cowtonlon = [-43.04,-43.21,-41.69,-38.18,-32.86,...
             -28.05,-29.09,-28.59,-24.30,-22.46];
[A,B] = xlsread('cowton/PNAS retreat data.xlsx');
for ii=1:10,
    glaciers(ii0+ii).termpos.t = A(:,1)+0.5; % mean annual position
    glaciers(ii0+ii).termpos.L = A(:,ii+1);
    glaciers(ii0+ii).termpos.source = 'cowton';
    glaciers(ii0+ii).bjorkid = bjorkid{ii};
    glaciers(ii0+ii).rignotid = rignotid(ii);
    glaciers(ii0+ii).lat = cowtonlat(ii);
    glaciers(ii0+ii).lon = cowtonlon(ii);
end
%% NSIDC (keep only those with a rignot id)

ii0 = length(glaciers);
load nsidc/nsidc_glaciers_bowmethod.mat
for ii=1:length(nsidc_glaciers),
    if ~isempty(nsidc_glaciers(ii).rignotid) & ~isnan(nsidc_glaciers(ii).rignotid),
        ii0=ii0+1;
        glaciers(ii0).termpos.t = nsidc_glaciers(ii).t;
        glaciers(ii0).termpos.L = nsidc_glaciers(ii).L;
        glaciers(ii0).rignotid = nsidc_glaciers(ii).rignotid;
        glaciers(ii0).termpos.source = 'nsidc';
        glaciers(ii0).x = nsidc_glaciers(ii).x;
        glaciers(ii0).y = nsidc_glaciers(ii).y;
    end
end
%% Bunce 2018 - J. Glac. (only those with rignot id)

% SE glaciers
ii0 = length(glaciers);
load bunce/bunce_se_glaciers.mat
for ii=1:length(bunce_se_glaciers),
    if ~isempty(bunce_se_glaciers(ii).rignotid) & ~isnan(bunce_se_glaciers(ii).rignotid),
        ii0=ii0+1;
        glaciers(ii0).termpos.t = bunce_se_glaciers(ii).t;
        glaciers(ii0).termpos.L = bunce_se_glaciers(ii).L;
        glaciers(ii0).rignotid = bunce_se_glaciers(ii).rignotid;
        glaciers(ii0).termpos.source = 'bunce';
        glaciers(ii0).x = bunce_se_glaciers(ii).x;
        glaciers(ii0).y = bunce_se_glaciers(ii).y;
    end
end
% NW glaciers
ii0 = length(glaciers);
load bunce/bunce_nw_glaciers.mat
for ii=1:length(bunce_nw_glaciers),
    if ~isempty(bunce_nw_glaciers(ii).rignotid) & ~isnan(bunce_nw_glaciers(ii).rignotid),
        ii0=ii0+1;
        glaciers(ii0).termpos.t = bunce_nw_glaciers(ii).t;
        glaciers(ii0).termpos.L = bunce_nw_glaciers(ii).L;
        glaciers(ii0).rignotid = bunce_nw_glaciers(ii).rignotid;
        glaciers(ii0).termpos.source = 'bunce';
        glaciers(ii0).x = bunce_nw_glaciers(ii).x;
        glaciers(ii0).y = bunce_nw_glaciers(ii).y;
    end
end
%% Carr 2017, A. Glac. (only those with rignot id)

ii0 = length(glaciers);
load carr/carr_glaciers.mat
for ii=1:length(carr_glaciers),
    if ~isempty(carr_glaciers(ii).rignotid) & ~isnan(carr_glaciers(ii).rignotid),
        ii0=ii0+1;
        glaciers(ii0).termpos.t = carr_glaciers(ii).t;
        glaciers(ii0).termpos.L = carr_glaciers(ii).L;
        glaciers(ii0).rignotid = carr_glaciers(ii).rignotid;
        glaciers(ii0).termpos.source = 'carr';
        glaciers(ii0).x = carr_glaciers(ii).x;
        glaciers(ii0).y = carr_glaciers(ii).y;
    end
end
%% make raw plot

% % raw terminus position
% figure(); hold on;
% for ii=1:length(glaciers),
%     plot(glaciers(ii).termpos.t,glaciers(ii).termpos.L,'.-');
% end
% xlabel('year'); ylabel('terminus position (km)');
%% processing

% remove any glacier without a rignot id
glaciers_to_remove = [];
for ii=1:length(glaciers),
    if isempty(glaciers(ii).rignotid) | isnan(glaciers(ii).rignotid),
        glaciers_to_remove = [glaciers_to_remove,ii];
    end
end
glaciers(glaciers_to_remove) = [];

% remove any NaNs from terminus position records
for ii=1:length(glaciers),
    inds = find(~isnan(glaciers(ii).termpos.L) & ~isnan(glaciers(ii).termpos.t));
    glaciers(ii).termpos.L = glaciers(ii).termpos.L(inds);
    glaciers(ii).termpos.t = glaciers(ii).termpos.t(inds);
end

% remove ice shelf glaciers
ids = [glaciers.rignotid];
glaciers_to_remove = [];
glaciers_to_remove = [glaciers_to_remove,find(ids==78)]; % petermann
glaciers_to_remove = [glaciers_to_remove,find(ids==144)]; % ryder
glaciers_to_remove = [glaciers_to_remove,find(ids==74)]; % 79N
% remove storestrommen as it surges (Mouginot paper)
glaciers_to_remove = [glaciers_to_remove,find(ids==216)];
% remove 81 due to complications with it splitting to form 2 glaciers
glaciers_to_remove = [glaciers_to_remove,find(ids==81)];
% remove 164 as it is essentially land-terminating
glaciers_to_remove = [glaciers_to_remove,find(ids==164)];
glaciers(glaciers_to_remove) = [];
%% deal with duplicates case by case

ids = [glaciers.rignotid];
[aa,bb] = unique(ids);
glaciers_to_remove = [];

% make example plot for Helheim of how to deal with duplicates
helheimid = 3;
helheimid = find(ids==helheimid);
for jj=1:length(helheimid),
    helheimsources{jj} = glaciers(helheimid(jj)).termpos.source;
end
andresen_t = glaciers(helheimid(1)).termpos.t;
andresen_L = glaciers(helheimid(1)).termpos.L;
cowton_t = glaciers(helheimid(2)).termpos.t;
cowton_L = glaciers(helheimid(2)).termpos.L;
nsidc_t = glaciers(helheimid(3)).termpos.t;
nsidc_L = glaciers(helheimid(3)).termpos.L;
bunce_t = glaciers(helheimid(4)).termpos.t;
bunce_L = glaciers(helheimid(4)).termpos.L;
carr_t = glaciers(helheimid(5)).termpos.t;
carr_L = glaciers(helheimid(5)).termpos.L;
% start with cowton as final record
t = cowton_t;
L = cowton_L;
% add nsidc to the end of cowton
offset = L(end) - interp1(nsidc_t,nsidc_L,t(end));
offset_inds = find(nsidc_t>t(end));
t = [t;nsidc_t(offset_inds)'];
L = [L;nsidc_L(offset_inds)'+offset];
t1 = t;
L1 = L;
t = andresen_t;
L = andresen_L;
% replace long term record with cowton/nsidc
offset = L1(1) - interp1(t,L,t1(1));
offset_inds = find(t>t1(1));
t(offset_inds) = [];
L(offset_inds) = [];
t = [t;t1];
L = [L;L1-offset];
% plot
lw = 1;
ss = 7;
fs = 8;
figure(); hold on;
plot(cowton_t,cowton_L+5,'.-','linewidth',lw,'markersize',ss);
plot(bunce_t,bunce_L+4,'.-','linewidth',lw,'markersize',ss);
plot(carr_t,carr_L+3,'.-','linewidth',lw,'markersize',ss);
plot(nsidc_t,nsidc_L+3,'.-','linewidth',lw,'markersize',ss);
plot(andresen_t,andresen_L-18,'.-','linewidth',lw,'markersize',ss);
plot(t,L-20,'k.-','linewidth',lw,'markersize',ss);
ylabel('relative terminus position (km)','fontsize',fs);
set(gca,'fontsize',fs,'box','on');
l = legend('Cowton18','Bunce18','Carr17','NSIDC','Andresen11','Merged');
set(l,'location','northwest','fontsize',fs);
saveplot(12,7,300,'../../ismip6/write_up/plots/helheim_merge.png');
close all;

% deal first with glaciers with a cowton or catania record
for ii=1:length(aa),
    inds = find(ids==aa(ii));
    sources = [];
    catania_ind = [];
    cowton_ind = [];
    nsidc_ind = [];
    bunce_ind = [];
    carr_ind = [];
    for jj=1:length(inds),
        sources{jj} = glaciers(inds(jj)).termpos.source;
    end
    % cowton records
    if any(ismember(sources,'cowton')),
        % remove carr and bunce as they are redundant for these glaciers
        cowton_ind = find(ismember(sources,'cowton'));
        nsidc_ind = find(ismember(sources,'nsidc'));
        carr_ind = find(ismember(sources,'carr'));
        bunce_ind = find(ismember(sources,'bunce'));
        if ~isempty(carr_ind), glaciers_to_remove = [glaciers_to_remove,inds(carr_ind)]; end
        if ~isempty(bunce_ind), glaciers_to_remove = [glaciers_to_remove,inds(bunce_ind)]; end
        % add NSIDC to the end of cowton if NSIDC is more recent
        % then remove NSIDC
        if ~isempty(nsidc_ind) & glaciers(inds(nsidc_ind)).termpos.t(end) > glaciers(inds(cowton_ind)).termpos.t(end),
            offset = glaciers(inds(cowton_ind)).termpos.L(end) - interp1(glaciers(inds(nsidc_ind)).termpos.t,glaciers(inds(nsidc_ind)).termpos.L,glaciers(inds(cowton_ind)).termpos.t(end));
            offset_inds = find(glaciers(inds(nsidc_ind)).termpos.t>glaciers(inds(cowton_ind)).termpos.t(end));
            glaciers(inds(cowton_ind)).termpos.t = [glaciers(inds(cowton_ind)).termpos.t;glaciers(inds(nsidc_ind)).termpos.t(offset_inds)'];
            glaciers(inds(cowton_ind)).termpos.L = [glaciers(inds(cowton_ind)).termpos.L;glaciers(inds(nsidc_ind)).termpos.L(offset_inds)'+offset];
            glaciers_to_remove = [glaciers_to_remove,inds(nsidc_ind)];
            glaciers(inds(cowton_ind)).termpos.source = 'cowton/nsidc';
        end
    end    
    % catania records
    if any(ismember(sources,'catania')),
        % remove carr and bunce as they are redundant for these glaciers
        catania_ind = find(ismember(sources,'catania'));
        nsidc_ind = find(ismember(sources,'nsidc'));
        carr_ind = find(ismember(sources,'carr'));
        bunce_ind = find(ismember(sources,'bunce'));
        if ~isempty(carr_ind), glaciers_to_remove = [glaciers_to_remove,inds(carr_ind)]; end
        if ~isempty(bunce_ind), glaciers_to_remove = [glaciers_to_remove,inds(bunce_ind)]; end
        % add NSIDC to the end of catania if NSIDC is more recent
        % then remove NSIDC
        if ~isempty(nsidc_ind) & glaciers(inds(nsidc_ind)).termpos.t(end) > glaciers(inds(catania_ind)).termpos.t(end),
            offset = glaciers(inds(catania_ind)).termpos.L(end) - interp1(glaciers(inds(nsidc_ind)).termpos.t,glaciers(inds(nsidc_ind)).termpos.L,glaciers(inds(catania_ind)).termpos.t(end));
            offset_inds = find(glaciers(inds(nsidc_ind)).termpos.t>glaciers(inds(catania_ind)).termpos.t(end));
            glaciers(inds(catania_ind)).termpos.t = [glaciers(inds(catania_ind)).termpos.t';glaciers(inds(nsidc_ind)).termpos.t(offset_inds)'];
            glaciers(inds(catania_ind)).termpos.L = [glaciers(inds(catania_ind)).termpos.L';glaciers(inds(nsidc_ind)).termpos.L(offset_inds)'+offset];
            glaciers_to_remove = [glaciers_to_remove,inds(nsidc_ind)];
            glaciers(inds(catania_ind)).termpos.source = 'catania/nsidc';
        end
    end    
end
glaciers(glaciers_to_remove) = [];

% every catania/cowton glacier now has only one record, except helhiem
ids = [glaciers.rignotid];
[aa,bb] = unique(ids);
glaciers_to_remove = [];
inds = find(ids==3);
offset = glaciers(inds(2)).termpos.L(1) - interp1(glaciers(inds(1)).termpos.t,glaciers(inds(1)).termpos.L,glaciers(inds(2)).termpos.t(1));
offset_inds = find(glaciers(inds(1)).termpos.t>glaciers(inds(2)).termpos.t(1));
glaciers(inds(1)).termpos.t(offset_inds) = [];
glaciers(inds(1)).termpos.L(offset_inds) = [];
glaciers(inds(1)).termpos.t = [glaciers(inds(1)).termpos.t;glaciers(inds(2)).termpos.t];
glaciers(inds(1)).termpos.L = [glaciers(inds(1)).termpos.L;glaciers(inds(2)).termpos.L-offset];
glaciers(inds(1)).termpos.source = 'andresen/cowton/nsidc';
glaciers_to_remove = [glaciers_to_remove,inds(2)];
glaciers(glaciers_to_remove) = [];

% deal with remaining 5 glaciers with long records
% all have a long record, nsidc and carr
ids = [glaciers.rignotid];
[aa,bb] = unique(ids);
glaciers_to_remove = [];
id = [36,1,26,18,22];
for kk=1:length(id),
    inds = find(ids==id(kk));
    sources = [];
    nsidc_ind = [];
    for jj=1:length(inds),
        sources{jj} = glaciers(inds(jj)).termpos.source;
    end
    % first remove carr record
    glaciers_to_remove = [glaciers_to_remove,inds(find(ismember(sources,'carr')))];
    % then combine nsidc record with long record if nsidc is more recent
    % then remove nsidc
    nsidc_ind = find(ismember(sources,'nsidc'));
    if ~isempty(nsidc_ind) & glaciers(inds(nsidc_ind)).termpos.t(end) > glaciers(inds(1)).termpos.t(end),
        offset = glaciers(inds(1)).termpos.L(end) - interp1(glaciers(inds(nsidc_ind)).termpos.t,glaciers(inds(nsidc_ind)).termpos.L,glaciers(inds(1)).termpos.t(end));
        offset_inds = find(glaciers(inds(nsidc_ind)).termpos.t>glaciers(inds(1)).termpos.t(end));
        if size(glaciers(inds(1)).termpos.t,1) > size(glaciers(inds(1)).termpos.t,2),
            glaciers(inds(1)).termpos.t = [glaciers(inds(1)).termpos.t;glaciers(inds(nsidc_ind)).termpos.t(offset_inds)'];
            glaciers(inds(1)).termpos.L = [glaciers(inds(1)).termpos.L;glaciers(inds(nsidc_ind)).termpos.L(offset_inds)'+offset];
        else
            glaciers(inds(1)).termpos.t = [glaciers(inds(1)).termpos.t';glaciers(inds(nsidc_ind)).termpos.t(offset_inds)'];
            glaciers(inds(1)).termpos.L = [glaciers(inds(1)).termpos.L';glaciers(inds(nsidc_ind)).termpos.L(offset_inds)'+offset];
        end
        glaciers_to_remove = [glaciers_to_remove,inds(nsidc_ind)];
        glaciers(inds(1)).termpos.source = [glaciers(inds(1)).termpos.source,'/nsidc'];
    end
end
glaciers(glaciers_to_remove) = [];

% everything else has some combination of nsidc/carr/bunce
ids = [glaciers.rignotid];
[aa,bb] = unique(ids);
glaciers_to_remove = [];
for ii=1:length(aa),

    inds = find(ids==aa(ii));
    
    if length(inds)>1,
        sources = [];
        bunce_ind = [];
        carr_ind = [];
        nsidc_ind = [];
        for jj=1:length(inds),
            sources{jj} = glaciers(inds(jj)).termpos.source;
        end
        bunce_ind = find(ismember(sources,'bunce'));
        carr_ind = find(ismember(sources,'carr'));
        nsidc_ind = find(ismember(sources,'nsidc'));
        
        % if bunce exists, start with this as root
        if ~isempty(bunce_ind),
            % if additionally nsidc exists, add more recent positions
            if ~isempty(nsidc_ind),
                offset = glaciers(inds(bunce_ind)).termpos.L(end) - interp1(glaciers(inds(nsidc_ind)).termpos.t,glaciers(inds(nsidc_ind)).termpos.L,glaciers(inds(bunce_ind)).termpos.t(end));
                offset_inds = find(glaciers(inds(nsidc_ind)).termpos.t>glaciers(inds(bunce_ind)).termpos.t(end));
                glaciers(inds(bunce_ind)).termpos.t = [glaciers(inds(bunce_ind)).termpos.t;glaciers(inds(nsidc_ind)).termpos.t(offset_inds)'];
                glaciers(inds(bunce_ind)).termpos.L = [glaciers(inds(bunce_ind)).termpos.L;glaciers(inds(nsidc_ind)).termpos.L(offset_inds)'+offset];
                glaciers_to_remove = [glaciers_to_remove,inds(nsidc_ind)];
                glaciers(inds(bunce_ind)).termpos.source = 'bunce/nsidc';
            end
            % if additionally carr exists, add 1992 position
            if ~isempty(carr_ind),
                offset = glaciers(inds(bunce_ind)).termpos.L(1) - interp1(glaciers(inds(carr_ind)).termpos.t,glaciers(inds(carr_ind)).termpos.L,glaciers(inds(bunce_ind)).termpos.t(1));
                offset_inds = find(glaciers(inds(carr_ind)).termpos.t<glaciers(inds(bunce_ind)).termpos.t(1));
                glaciers(inds(bunce_ind)).termpos.t = [glaciers(inds(carr_ind)).termpos.t(offset_inds)';glaciers(inds(bunce_ind)).termpos.t];
                glaciers(inds(bunce_ind)).termpos.L = [glaciers(inds(carr_ind)).termpos.L(offset_inds)'+offset;glaciers(inds(bunce_ind)).termpos.L];
                glaciers_to_remove = [glaciers_to_remove,inds(carr_ind)];
                glaciers(inds(bunce_ind)).termpos.source = 'carr/bunce/nsidc';
            end
            
        else % remaining glaciers must have nsidc and carr            
            
            if glaciers(inds(nsidc_ind)).termpos.t(1)>glaciers(inds(carr_ind)).termpos.t(end), % if no overlap, keep carr not nsidc
                glaciers_to_remove = [glaciers_to_remove,inds(nsidc_ind)];
            else % add start of carr to nsidc
                offset = glaciers(inds(nsidc_ind)).termpos.L(1) - interp1(glaciers(inds(carr_ind)).termpos.t,glaciers(inds(carr_ind)).termpos.L,glaciers(inds(nsidc_ind)).termpos.t(1));
                offset_inds = find(glaciers(inds(carr_ind)).termpos.t<glaciers(inds(nsidc_ind)).termpos.t(1));
                glaciers(inds(nsidc_ind)).termpos.t = [glaciers(inds(carr_ind)).termpos.t(offset_inds)';glaciers(inds(nsidc_ind)).termpos.t'];
                glaciers(inds(nsidc_ind)).termpos.L = [glaciers(inds(carr_ind)).termpos.L(offset_inds)'+offset;glaciers(inds(nsidc_ind)).termpos.L'];
                glaciers_to_remove = [glaciers_to_remove,inds(carr_ind)];
                glaciers(inds(nsidc_ind)).termpos.source = 'carr/nsidc';            
            end
            
        end
                
    end
end
glaciers(glaciers_to_remove) = [];
%% convert lat-lon to x-y and vice versa

for ii=1:length(glaciers),
    if ~isempty(glaciers(ii).lat) & isempty(glaciers(ii).x),
        [glaciers(ii).x,glaciers(ii).y] = latlon2utm(glaciers(ii).lat,glaciers(ii).lon);
    elseif ~isempty(glaciers(ii).x) & isempty(glaciers(ii).lat),
        [glaciers(ii).lat,glaciers(ii).lon] = polarstereo_inv(glaciers(ii).x,glaciers(ii).y);
    end
end
%% get distance to land terminating

% % mouginot to rignot id conversion
% m2r = xlsread('rignotid2mouginotbasin.xlsx');
% 
% % load mouginot basins from Heiko
% load('GrIS_01000m_REGIONS.mat');
% 
% % load distances from ocean from Heiko
% dist=ncread('dist_01000m.nc','dist');
% dist(find(dist>=9999)) = NaN;
% 
% for ii=1:length(glaciers),
%     
%     % get mouginot basin number for this rignotid
%     mid = m2r(find(m2r(:,2)==glaciers(ii).rignotid,1));
%     % get basin indices
%     if ~isempty(mid),
%         basin_ids = find(id==mid);
%     else basin_ids = [];
%     end
%     
%     if ~isempty(basin_ids) & ~isempty(find(~isnan(dist(basin_ids)))),
%         disti = NaN*dist;
%         disti(basin_ids) = dist(basin_ids);
%         disti(find(~isnan(disti)))=1;
%         disti(find(isnan(disti)))=0;
%         aa = bwconncomp(disti);
%         % find connected group with terminus in it
%         termgroup = NaN;
%         for j=1:aa.NumObjects,
%             if ismember(min(dist(basin_ids)),dist(aa.PixelIdxList{j})),
%                 termgroup = j;
%             end
%         end
%         % manual exceptions...
%         if glaciers(ii).rignotid==71, termgroup=2; end
%         if glaciers(ii).rignotid==109, termgroup=5; end
%         if glaciers(ii).rignotid==169, termgroup=1; end
%         if glaciers(ii).rignotid==228, termgroup=NaN; end
%         if glaciers(ii).rignotid==11, termgroup=1; end
%         if glaciers(ii).rignotid==57, termgroup=1; end
%         if glaciers(ii).rignotid==9, termgroup=1; end
%         if glaciers(ii).rignotid==126, termgroup=11; end
%         if glaciers(ii).rignotid==179, termgroup=4; end
%         if ~isnan(termgroup),
%             glaciers(ii).dist2land = max(dist(aa.PixelIdxList{termgroup}));
%         else
%             glaciers(ii).dist2land = NaN;
%         end
%     else
%         glaciers(ii).dist2land = NaN;
%     end
%     
%     % set 0s to NaNs
%     if glaciers(ii).dist2land == 0, glaciers(ii).dist2land=NaN; end;
% 
% end
%% associate water depth/ice thickness

% A = xlsread('water_depth.xlsx','rignot');
% id = A(:,1);
% waterdepth = A(:,2);
% for ii=1:length(glaciers),
%     if ~isempty(glaciers(ii).rignotid),
%         xref = find(id==glaciers(ii).rignotid);
%         if ~isempty(xref), glaciers(ii).h = waterdepth(xref); end
%     end
% end
%% associate runoff

load ../runoff/RACMO2.3p2_runoff_MM_tidewater_allbasins_rignotid_withGlacierIDs.mat
for ii=1:length(glaciers),
    for jj=1:length(runoff),
        if ismember(glaciers(ii).rignotid,runoff(jj).rignotGlacierID),
            glaciers(ii).rignotname = runoff(jj).rignotGlacierName;
            glaciers(ii).RACMO.t = 1958+(round(runoff(jj).time)+0.5)/12; % decimal years
            yrs = floor(glaciers(ii).RACMO.t);
            mnths = floor(12*(glaciers(ii).RACMO.t-floor(glaciers(ii).RACMO.t)))+1;
            glaciers(ii).RACMO.Q = runoff(jj).runoff./(86400*eomday(yrs,mnths)); % m3/s
        end
    end
end
% remove any glacier which does not have runoff
% rignotid=149 is the only one  - due to confusion about position it didn't
% make it into the drainage basins
glaciers_to_remove = [];
for ii=1:length(glaciers),
    if isempty(glaciers(ii).RACMO), glaciers_to_remove = [glaciers_to_remove,ii]; end
end
glaciers(glaciers_to_remove) = [];
% add JJA runoff
for ii=1:length(glaciers),
    [glaciers(ii).RACMO.tJJA,glaciers(ii).RACMO.QJJA] = time_bin_JJArunoff(glaciers(ii).RACMO.t,glaciers(ii).RACMO.Q,[1958:2018]);
end
%% associate ocean thermal forcing

load ../EN4/EN4_ISMIP6.mat
basins = {'SE','SW','CE','CW','NE','NW','NO'};
% threshold at which we don't believe ocean temperatures
Toi_threshold = 0.1;
for i=1:7,
    regions(i).TF(find(regions(i).Toi<Toi_threshold)) = NaN;
end
for ii=1:length(glaciers),
    % find ocean basin in which glacier is located
    basin = 0; j = 1;
    while basin == 0,
        if inpolygon(glaciers(ii).x,glaciers(ii).y,regions(j).ice.x,regions(j).ice.y),
            basin = j;
        else j = j+1;
        end
    end
    glaciers(ii).EN4.TF = regions(basin).TF;
    glaciers(ii).EN4.t = regions(basin).t+0.5; % for annual means
    glaciers(ii).sector = basins{basin};
    glaciers(ii).sectornum = basin;
end
%% calculate melt rate from Q and TF

baseline = [1995:2014];
for ii=1:length(glaciers),
    glaciers(ii).melt.t = glaciers(ii).RACMO.tJJA;
    glaciers(ii).melt.m = glaciers(ii).RACMO.QJJA.^0.4.*interp1(glaciers(ii).EN4.t,glaciers(ii).EN4.TF,glaciers(ii).melt.t);
end
baselineinds = find(ismember(floor(glaciers(1).melt.t),baseline));
for ii=1:length(glaciers),
    glaciers(ii).melt.meltbaseline = mean(glaciers(ii).melt.m(baselineinds));
end
%% associate 2000-2010 ice flux from Enderlin 2014

[A,B]=xlsread('../iceflux/flux_vs_retreat.xlsx');
rid = A(:,1);
iceflux = A(:,18);
for ii=1:length(glaciers),
    id = find(rid==glaciers(ii).rignotid);
    if ~isempty(id),
        glaciers(ii).iceflux.enderlin = iceflux(id);
    else
        glaciers(ii).iceflux.enderlin = NaN;
    end
end
%% associate 2000-2010 ice flux from King 2018

[A,B]=xlsread('../iceflux/Greenland_Discharge2000-2010.xlsx');
rid = A(:,7);
iceflux = A(:,3);
for ii=1:length(glaciers),
    id = find(rid==glaciers(ii).rignotid);
    if ~isempty(id),
        glaciers(ii).iceflux.king = iceflux(id);
    else
        glaciers(ii).iceflux.king = NaN;
    end
end
%% fill final ice flux based on these datasets and correlation with runoff

% fill based on datasets with King 2018 taking priority
for ii=1:length(glaciers),
    if ~isnan(glaciers(ii).iceflux.king),
        glaciers(ii).iceflux.final = glaciers(ii).iceflux.king;
    elseif ~isnan(glaciers(ii).iceflux.enderlin),
        glaciers(ii).iceflux.final = glaciers(ii).iceflux.enderlin;
    else
        glaciers(ii).iceflux.final = NaN;
    end
end
% fill those still without an iceflux based on correlation with runoff
for ii=1:length(glaciers),
    runoff_all(ii) = nanmean(glaciers(ii).RACMO.QJJA(find(ismember(floor(glaciers(ii).RACMO.tJJA),[2000:2010]))));
    iceflux_all(ii) = glaciers(ii).iceflux.final;
end
inds = find(~isnan(runoff_all)&~isnan(iceflux_all));
X = runoff_all(inds); Y = iceflux_all(inds);
p = polyfit(X,Y,1);
[R,pval] = corr(X',Y');
for ii=1:length(glaciers),
    if isnan(glaciers(ii).iceflux.final),
        glaciers(ii).iceflux.final = p(1)*runoff_all(ii) + p(2);
    end
end
%% save

save glaciers.mat glaciers