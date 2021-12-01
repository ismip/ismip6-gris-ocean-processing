clear; close all;
% script to add future runoff, ocean TF and melt to glaciers structure
load glaciers.mat

%% associate future thermal forcing from MIROC5 RCP8.5

% load MIROC5 RCP8.5 ocean
load ../process_CMIP/MIROC5_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
miroc5_baseline_inds = find(ismember(year,baseline));
TF0_miroc5 = nanmean(squeeze(TF_basins(:,miroc5_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% miroc5 bias
bias = TF0_miroc5 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_miroc5(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).MIROC5.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).MIROC5.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).MIROC5.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end

% add glacier TF baseline
baselineinds = find(ismember(floor(glaciers(1).MIROC5.RCP85.tTF),baseline));
for ii=1:length(glaciers),
    glaciers(ii).EN4.TFbaseline = mean(glaciers(ii).MIROC5.RCP85.TF(baselineinds));
end
%% associate future thermal forcing from MIROC5 RCP2.6

% load MIROC5 RCP2.6 ocean
load ../process_CMIP/MIROC5_RCP26.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
miroc5_baseline_inds = find(ismember(year,baseline));
TF0_miroc5 = nanmean(squeeze(TF_basins(:,miroc5_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% miroc5 bias
bias = TF0_miroc5 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_miroc5(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).MIROC5.RCP26.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).MIROC5.RCP26.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).MIROC5.RCP26.bias_TF = bias(glaciers(ii).sectornum);
end

% add glacier TF baseline
baselineinds = find(ismember(floor(glaciers(1).MIROC5.RCP26.tTF),baseline));
for ii=1:length(glaciers),
    glaciers(ii).EN4.TFbaseline = mean(glaciers(ii).MIROC5.RCP26.TF(baselineinds));
end
%% associate future thermal forcing from NorESM RCP8.5

% load NorESM RCP8.5 ocean
load ../process_CMIP/NorESM_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
noresm_baseline_inds = find(ismember(year,baseline));
TF0_noresm = nanmean(squeeze(TF_basins(:,noresm_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% noresm bias
bias = TF0_noresm - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_noresm(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).NorESM.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from NorESM RCP2.6

% load NorESM RCP2.6 ocean
load ../process_CMIP/NorESM_RCP26.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
noresm_baseline_inds = find(ismember(year,baseline));
TF0_noresm = nanmean(squeeze(TF_basins(:,noresm_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% noresm bias
bias = TF0_noresm - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_noresm(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).NorESM.RCP26.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).NorESM.RCP26.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).NorESM.RCP26.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from UKESM1-0-LL ssp585

% load UKESM1-0-LL ssp585 ocean
load ../process_CMIP/UKESM1-0-LL_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ukesm1_baseline_inds = find(ismember(year,baseline));
TF0_ukesm1 = nanmean(squeeze(TF_basins(:,ukesm1_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% ukesm1 bias
bias = TF0_ukesm1 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ukesm1(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).UKESM1.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).UKESM1.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).UKESM1.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CESM2 ssp585

% load CESM2 ssp585 ocean
load ../process_CMIP/CESM2_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cesm2_baseline_inds = find(ismember(year,baseline));
TF0_cesm2 = nanmean(squeeze(TF_basins(:,cesm2_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cesm2 bias
bias = TF0_cesm2 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cesm2(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CESM2.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CNRM-ESM2-1 ssp585

% load CNRM-ESM2-1 ssp585 ocean
load ../process_CMIP/CNRM-ESM2-1_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmesm2_baseline_inds = find(ismember(year,baseline));
TF0_cnrmesm2 = nanmean(squeeze(TF_basins(:,cnrmesm2_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cnrmesm2 bias
bias = TF0_cnrmesm2 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmesm2(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CNRMESM2.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CNRMESM2.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CNRMESM2.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CNRM-CM6-1 ssp585

% load CNRM-CM6-1 ssp585 ocean
load ../process_CMIP/CNRM-CM6-1_ssp585.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_cnrmcm6 = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cnrmcm6 bias
bias = TF0_cnrmcm6 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmcm6(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CNRMCM6.ssp585.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CNRMCM6.ssp585.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CNRMCM6.ssp585.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CNRM-CM6-1 ssp126

% load CNRM-CM6-1 ssp126 ocean
load ../process_CMIP/CNRM-CM6-1_ssp126.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
cnrmcm6_baseline_inds = find(ismember(year,baseline));
TF0_cnrmcm6 = nanmean(squeeze(TF_basins(:,cnrmcm6_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% cnrmcm6 bias
bias = TF0_cnrmcm6 - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_cnrmcm6(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CNRMCM6.ssp126.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CNRMCM6.ssp126.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CNRMCM6.ssp126.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from HadGEM RCP8.5

% load HadGEM RCP8.5 ocean
load ../process_CMIP/HadGEM_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
hadgem_baseline_inds = find(ismember(year,baseline));
TF0_hadgem = nanmean(squeeze(TF_basins(:,hadgem_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% hadgem bias
bias = TF0_hadgem - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_hadgem(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).HadGEM.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).HadGEM.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).HadGEM.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from HadGEM RCP2.6

% load HadGEM RCP2.6 ocean
load ../process_CMIP/HadGEM_RCP26.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
hadgem_baseline_inds = find(ismember(year,baseline));
TF0_hadgem = nanmean(squeeze(TF_basins(:,hadgem_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% hadgem bias
bias = TF0_hadgem - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_hadgem(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).HadGEM.RCP26.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).HadGEM.RCP26.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).HadGEM.RCP26.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CSIRO RCP8.5

% load CSIRO RCP8.5 ocean
load ../process_CMIP/CSIRO_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
csiro_baseline_inds = find(ismember(year,baseline));
TF0_csiro = nanmean(squeeze(TF_basins(:,csiro_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% csiro bias
bias = TF0_csiro - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_csiro(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CSIRO.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CSIRO.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CSIRO.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from CSIRO RCP2.6

% load CSIRO RCP2.6 ocean
load ../process_CMIP/CSIRO_RCP26.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
csiro_baseline_inds = find(ismember(year,baseline));
TF0_csiro = nanmean(squeeze(TF_basins(:,csiro_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% csiro bias
bias = TF0_csiro - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_csiro(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).CSIRO.RCP26.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).CSIRO.RCP26.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).CSIRO.RCP26.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from IPSLCM RCP8.5

% load IPSLCM RCP8.5 ocean
load ../process_CMIP/IPSLCM_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ipslcm_baseline_inds = find(ismember(year,baseline));
TF0_ipslcm = nanmean(squeeze(TF_basins(:,ipslcm_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% IPSLCM bias
bias = TF0_ipslcm - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ipslcm(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).IPSLCM.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from IPSLCM RCP2.6

% load IPSLCM RCP2.6 ocean
load ../process_CMIP/IPSLCM_RCP26.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ipslcm_baseline_inds = find(ismember(year,baseline));
TF0_ipslcm = nanmean(squeeze(TF_basins(:,ipslcm_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% IPSLCM bias
bias = TF0_ipslcm - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ipslcm(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).IPSLCM.RCP26.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).IPSLCM.RCP26.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).IPSLCM.RCP26.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future thermal forcing from ACCESS RCP8.5 (note there is no ACCESS RCP2.6)

% load ACCESS RCP8.5 ocean
load ../process_CMIP/ACCESS_RCP85.mat
% load EN4
load EN4_ISMIP6.mat

% present day baseline period
baseline = [1995:2014];
ACCESS_baseline_inds = find(ismember(year,baseline));
TF0_ACCESS = nanmean(squeeze(TF_basins(:,ACCESS_baseline_inds)),2);

% en4 baseline
EN4_baseline_inds = find(ismember(regions(1).t,baseline));
for l=1:7,
    TF0_EN4(l) = nanmean(regions(l).TF(EN4_baseline_inds));
end

% ACCESS bias
bias = TF0_ACCESS - TF0_EN4';

% future T for forcing
clearvars Tfuture
for l=1:7,
Tfuture(l,:) = TF0_EN4(l) + (TF_basins(l,:) - TF0_ACCESS(l));
end

% assign to glaciers
for ii=1:length(glaciers),
    glaciers(ii).ACCESS.RCP85.tTF = year+0.5; % follow convention of annual means
    glaciers(ii).ACCESS.RCP85.TF = Tfuture(glaciers(ii).sectornum,:);
    glaciers(ii).ACCESS.RCP85.bias_TF = bias(glaciers(ii).sectornum);
end
%% associate future runoff from MAR/MIROC5

% load MIROC5-MAR RCP2.6
load ../runoff/MAR3.9-MIROC5-rcp26-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_26 = runoff;
% load MIROC5-MAR RCP8.5
load ../runoff/MAR3.9-MIROC5-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_26(ii).time = [1950:2100];
    runoff_85(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp2.6
    for jj=1:length(runoff_26),
        if ismember(glaciers(ii).rignotid,runoff_26(jj).rignotGlacierID),
            glaciers(ii).MIROC5.RCP26.tJJA = runoff_26(jj).time+0.5; % annual mean convention
            glaciers(ii).MIROC5.RCP26.QJJA = []; 
            for kk=1:length(glaciers(ii).MIROC5.RCP26.tJJA),
                yr = floor(glaciers(ii).MIROC5.RCP26.tJJA(kk));
                glaciers(ii).MIROC5.RCP26.QJJA(kk) = 3*runoff_26(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).MIROC5.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            glaciers(ii).MIROC5.RCP85.QJJA = []; 
            for kk=1:length(glaciers(ii).MIROC5.RCP85.tJJA),
                yr = floor(glaciers(ii).MIROC5.RCP85.tJJA(kk));
                glaciers(ii).MIROC5.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).MIROC5.RCP26.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MIROC5.RCP26.bias_QJJA = mean(glaciers(ii).MIROC5.RCP26.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MIROC5.RCP85.bias_QJJA = mean(glaciers(ii).MIROC5.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).MIROC5.RCP26.QJJA = glaciers(ii).MIROC5.RCP26.QJJA - glaciers(ii).MIROC5.RCP26.bias_QJJA;
    glaciers(ii).MIROC5.RCP85.QJJA = glaciers(ii).MIROC5.RCP85.QJJA - glaciers(ii).MIROC5.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).MIROC5.RCP26.QJJA(find(glaciers(ii).MIROC5.RCP26.QJJA<0)) = 0;
    glaciers(ii).MIROC5.RCP85.QJJA(find(glaciers(ii).MIROC5.RCP85.QJJA<0)) = 0;
end
%% associate future runoff from MAR/NorESM

% load NorESM-MAR RCP8.5
load ../runoff/MAR3.9-NorESM1-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_85(ii).time = [1950:2100]; % time variable is off 
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).NorESM.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).NorESM.RCP85.tJJA),
                yr = floor(glaciers(ii).NorESM.RCP85.tJJA(kk));
                glaciers(ii).NorESM.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).NorESM.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM.RCP85.bias_QJJA = mean(glaciers(ii).NorESM.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).NorESM.RCP85.QJJA = glaciers(ii).NorESM.RCP85.QJJA - glaciers(ii).NorESM.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).NorESM.RCP85.QJJA(find(glaciers(ii).NorESM.RCP85.QJJA<0)) = 0;
end
%% associate future runoff from MAR/CNRM-CM6-1
% load CNRM-CM6-1/MAR ssp126
load ../runoff/MAR3.9-CNRM-CM6-ssp126-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_26 = runoff;
% load CNRM-CM6-1/MAR ssp585
load ../runoff/MAR3.9-CNRM-CM6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_26(ii).time = [1950:2100];
    runoff_85(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp126
    for jj=1:length(runoff_26),
        if ismember(glaciers(ii).rignotid,runoff_26(jj).rignotGlacierID),
            glaciers(ii).CNRMCM6.ssp126.tJJA = runoff_26(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CNRMCM6.ssp126.tJJA),
                yr = floor(glaciers(ii).CNRMCM6.ssp126.tJJA(kk));
                glaciers(ii).CNRMCM6.ssp126.QJJA(kk) = 3*runoff_26(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
    % ssp585
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).CNRMCM6.ssp585.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CNRMCM6.ssp585.tJJA),
                yr = floor(glaciers(ii).CNRMCM6.ssp585.tJJA(kk));
                glaciers(ii).CNRMCM6.ssp585.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CNRMCM6.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMCM6.ssp585.bias_QJJA = mean(glaciers(ii).CNRMCM6.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMCM6.ssp585.QJJA = glaciers(ii).CNRMCM6.ssp585.QJJA - glaciers(ii).CNRMCM6.ssp585.bias_QJJA;
    glaciers(ii).CNRMCM6.ssp126.bias_QJJA = mean(glaciers(ii).CNRMCM6.ssp126.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMCM6.ssp126.QJJA = glaciers(ii).CNRMCM6.ssp126.QJJA - glaciers(ii).CNRMCM6.ssp126.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CNRMCM6.ssp585.QJJA(find(glaciers(ii).CNRMCM6.ssp585.QJJA<0)) = 0;
    glaciers(ii).CNRMCM6.ssp126.QJJA(find(glaciers(ii).CNRMCM6.ssp126.QJJA<0)) = 0;
end
%% associate future runoff from MAR/CNRM-ESM2-1 (ssp 585 only)

% load CNRM-ESM2-1/MAR ssp585
load ../runoff/MAR3.9-CNRM-ESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_85(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).CNRMESM2.ssp585.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CNRMESM2.ssp585.tJJA),
                yr = floor(glaciers(ii).CNRMESM2.ssp585.tJJA(kk));
                glaciers(ii).CNRMESM2.ssp585.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CNRMESM2.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMESM2.ssp585.bias_QJJA = mean(glaciers(ii).CNRMESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CNRMESM2.ssp585.QJJA = glaciers(ii).CNRMESM2.ssp585.QJJA - glaciers(ii).CNRMESM2.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CNRMESM2.ssp585.QJJA(find(glaciers(ii).CNRMESM2.ssp585.QJJA<0)) = 0;
end
%% associate future runoff from MAR/UKESM1-0-LL (ssp 585 only)

% load UKESM1-0-LL/MAR ssp585
load ../runoff/MAR3.9-UKESM1-CM6-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_85(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).UKESM1.ssp585.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).UKESM1.ssp585.tJJA),
                yr = floor(glaciers(ii).UKESM1.ssp585.tJJA(kk));
                glaciers(ii).UKESM1.ssp585.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).UKESM1.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1.ssp585.bias_QJJA = mean(glaciers(ii).UKESM1.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).UKESM1.ssp585.QJJA = glaciers(ii).UKESM1.ssp585.QJJA - glaciers(ii).UKESM1.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).UKESM1.ssp585.QJJA(find(glaciers(ii).UKESM1.ssp585.QJJA<0)) = 0;
end
%% associate future runoff from MAR/CESM2 (ssp 585 only)

% load CESM2/MAR ssp585
load ../runoff/MAR3.9-CESM2-ssp585-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat
runoff_85 = runoff;

% sort time variable
for ii=1:length(runoff_85),
    runoff_85(ii).time = [1950:2100];
end

% assign to glaciers
for ii=1:length(glaciers),
    % ssp585
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).CESM2.ssp585.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CESM2.ssp585.tJJA),
                yr = floor(glaciers(ii).CESM2.ssp585.tJJA(kk));
                glaciers(ii).CESM2.ssp585.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CESM2.ssp585.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2.ssp585.bias_QJJA = mean(glaciers(ii).CESM2.ssp585.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CESM2.ssp585.QJJA = glaciers(ii).CESM2.ssp585.QJJA - glaciers(ii).CESM2.ssp585.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CESM2.ssp585.QJJA(find(glaciers(ii).CESM2.ssp585.QJJA<0)) = 0;
end
%% associate future runoff from MAR/HadGEM

% load HadGEM-MAR RCP8.5 run
load ../runoff/MAR3.9-HadGEM2-ES-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat

% combine into continuous time series for rcp8.5
for ii=1:length(runoff),
    runoff_85(ii).time = runoff(ii).time;
    runoff_85(ii).runoff = runoff(ii).runoff;
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).HadGEM.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).HadGEM.RCP85.tJJA),
                yr = floor(glaciers(ii).HadGEM.RCP85.tJJA(kk));
                glaciers(ii).HadGEM.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).HadGEM.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).HadGEM.RCP85.bias_QJJA = mean(glaciers(ii).HadGEM.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).HadGEM.RCP85.QJJA = glaciers(ii).HadGEM.RCP85.QJJA - glaciers(ii).HadGEM.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).HadGEM.RCP85.QJJA(find(glaciers(ii).HadGEM.RCP85.QJJA<0)) = 0;
end
%% associate future runoff from MAR/CSIRO

% load CSIRO-MAR RCP8.5 run
load ../runoff/MAR3.9-CSIRO-Mk3.6-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat

% combine into continuous time series for rcp8.5
for ii=1:length(runoff),
    runoff_85(ii).time = runoff(ii).time;
    runoff_85(ii).runoff = runoff(ii).runoff;
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).CSIRO.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).CSIRO.RCP85.tJJA),
                yr = floor(glaciers(ii).CSIRO.RCP85.tJJA(kk));
                glaciers(ii).CSIRO.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).CSIRO.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CSIRO.RCP85.bias_QJJA = mean(glaciers(ii).CSIRO.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).CSIRO.RCP85.QJJA = glaciers(ii).CSIRO.RCP85.QJJA - glaciers(ii).CSIRO.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).CSIRO.RCP85.QJJA(find(glaciers(ii).CSIRO.RCP85.QJJA<0)) = 0;
end
%% associate future runoff from MAR/IPSLCM

% load IPSLCM-MAR RCP8.5 run
load ../runoff/MAR3.9-IPSL-CM5-MR-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat

% combine into continuous time series for rcp8.5
for ii=1:length(runoff),
    runoff_85(ii).time = runoff(ii).time;
    runoff_85(ii).runoff = runoff(ii).runoff;
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).IPSLCM.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).IPSLCM.RCP85.tJJA),
                yr = floor(glaciers(ii).IPSLCM.RCP85.tJJA(kk));
                glaciers(ii).IPSLCM.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).IPSLCM.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM.RCP85.bias_QJJA = mean(glaciers(ii).IPSLCM.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).IPSLCM.RCP85.QJJA = glaciers(ii).IPSLCM.RCP85.QJJA - glaciers(ii).IPSLCM.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).IPSLCM.RCP85.QJJA(find(glaciers(ii).IPSLCM.RCP85.QJJA<0)) = 0;
end
%% associate future runoff from MAR/ACCESS

% load ACCESS-1-3 RCP8.5 run
load ../runoff/MAR3.9-ACCESS1.3-rcp85-JJA_mean-tidewaterbasins_rignotid_withGlacierIDs.mat

% combine into continuous time series for rcp8.5
for ii=1:length(runoff),
    runoff_85(ii).time = runoff(ii).time;
    runoff_85(ii).runoff = runoff(ii).runoff;
end

% assign to glaciers
for ii=1:length(glaciers),
    % rcp8.5
    for jj=1:length(runoff_85),
        if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
            glaciers(ii).ACCESS.RCP85.tJJA = runoff_85(jj).time+0.5; % annual mean convention
            for kk=1:length(glaciers(ii).ACCESS.RCP85.tJJA),
                yr = floor(glaciers(ii).ACCESS.RCP85.tJJA(kk));
                glaciers(ii).ACCESS.RCP85.QJJA(kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
            end
        end
    end
end

% present day baseline period
baseline = [1995:2014];
MARbaselineinds = find(ismember(floor(glaciers(1).ACCESS.RCP85.tJJA),baseline));
RACMObaselineinds = find(ismember(floor(glaciers(1).RACMO.tJJA),baseline));

% calculate and account for bias over baseline and add Q baseline
for ii=1:length(glaciers),
    glaciers(ii).RACMO.Qbaseline = mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).ACCESS.RCP85.bias_QJJA = mean(glaciers(ii).ACCESS.RCP85.QJJA(MARbaselineinds)) - mean(glaciers(ii).RACMO.QJJA(RACMObaselineinds));
    glaciers(ii).ACCESS.RCP85.QJJA = glaciers(ii).ACCESS.RCP85.QJJA - glaciers(ii).ACCESS.RCP85.bias_QJJA;
    % make sure no runoff values less than 0
    glaciers(ii).ACCESS.RCP85.QJJA(find(glaciers(ii).ACCESS.RCP85.QJJA<0)) = 0;
end
%% calculate future MIROC5 RCP8.5 melt - old method

% % NOTE RELOAD RUNOFF AND TF FROM SCRATCH AS WANT TO DO BIAS FROM SCRATCH
% % load MIROC5-MAR historical run
% load ../runoff/MAR3.9MIROC5histo19502005tidewaterbasinsrignotidwithGlacierIDs.mat
% runoff_histo = runoff;
% % load MIROC5-MAR RCP8.5
% load ../runoff/MAR3.9MIROC5rcp8520062100tidewaterbasinsrignotidwithGlacierIDs.mat
% runoff_85 = runoff;
% 
% % combine into continuous time series for rcp2.6 and rcp8.5
% for ii=1:length(runoff_histo),
%     runoff_85(ii).time = [runoff_histo(ii).time,runoff_85(ii).time];
%     runoff_85(ii).runoff = [runoff_histo(ii).runoff,runoff_85(ii).runoff];
% end
% 
% % get raw runoff array
% for ii=1:length(glaciers),
%     % rcp8.5
%     for jj=1:length(runoff_85),
%         if ismember(glaciers(ii).rignotid,runoff_85(jj).rignotGlacierID),
%             t_Q(ii,:) = runoff_85(jj).time+0.5; % annual mean convention
%             for kk=1:length(t_Q(ii,:)),
%                 yr = floor(t_Q(ii,kk));
%                 Q(ii,kk) = 3*runoff_85(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
%             end
%         end
%     end
% end
% 
% % get raw TF array
% % load MIROC5
% load ../CMIP5/MIROC5_RCP85.mat
% for ii=1:length(glaciers),
%     t_TF(ii,:) = year+0.5; % follow convention of annual means
%     TF(ii,:) = TF_basins(glaciers(ii).sectornum,:);
% end
% 
% % make raw melt
% for ii=1:length(glaciers),
%     t_melt(ii,:) = t_Q(ii,:);
%     melt(ii,:) = Q(ii,:).^0.4.*interp1(t_TF(ii,:),TF(ii,:),t_melt(ii,:));
% end
% 
% % baseline bias
% baseline = [1995:2014];
% MIROC5baselineinds = find(ismember(floor(t_melt(1,:)),baseline));
% for ii=1:length(glaciers),
%     MIROC5meltbaseline(ii) = mean(melt(ii,MIROC5baselineinds));
%     bias(ii) = MIROC5meltbaseline(ii) - glaciers(ii).melt.meltbaseline;
% end
% 
% % assign to glaciers
% for ii=1:length(glaciers),
%     glaciers(ii).MIROC5.RCP85.tmelt = t_melt(ii,:);
%     glaciers(ii).MIROC5.RCP85.bias_melt = bias(ii);
%     glaciers(ii).MIROC5.RCP85.melt = melt(ii,:) - glaciers(ii).MIROC5.RCP85.bias_melt;
%     % make sure no melt values less than 0
%     glaciers(ii).MIROC5.RCP85.melt(find(glaciers(ii).MIROC5.RCP85.melt<0)) = 0;
% end
%% calculate future MIROC5 melt - new method

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).MIROC5.RCP85.tmelt = glaciers(ii).MIROC5.RCP85.tJJA;
    glaciers(ii).MIROC5.RCP85.melt = (glaciers(ii).MIROC5.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).MIROC5.RCP85.tTF,glaciers(ii).MIROC5.RCP85.TF,glaciers(ii).MIROC5.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).MIROC5.RCP85.melt(find(glaciers(ii).MIROC5.RCP85.melt<0)) = 0;
    % RCP2.6
    glaciers(ii).MIROC5.RCP26.tmelt = glaciers(ii).MIROC5.RCP26.tJJA;
    glaciers(ii).MIROC5.RCP26.melt = (glaciers(ii).MIROC5.RCP26.QJJA.^0.4).*...
        interp1(glaciers(ii).MIROC5.RCP26.tTF,glaciers(ii).MIROC5.RCP26.TF,glaciers(ii).MIROC5.RCP26.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).MIROC5.RCP26.melt(find(glaciers(ii).MIROC5.RCP26.melt<0)) = 0;
end
%% calculate future MIROC5 RCP2.6 melt - old method

% % NOTE RELOAD RUNOFF AND TF FROM SCRATCH AS WANT TO DO BIAS FROM SCRATCH
% % load MIROC5-MAR historical run
% load ../runoff/MAR3.9MIROC5histo19502005tidewaterbasinsrignotidwithGlacierIDs.mat
% runoff_histo = runoff;
% % load MIROC5-MAR RCP2.6
% load ../runoff/MAR3.9MIROC5rcp2620062100tidewaterbasinsrignotidwithGlacierIDs.mat
% runoff_26 = runoff;
% 
% % combine into continuous time series for rcp2.6 and rcp8.5
% for ii=1:length(runoff_histo),
%     runoff_26(ii).time = [runoff_histo(ii).time,runoff_26(ii).time];
%     runoff_26(ii).runoff = [runoff_histo(ii).runoff,runoff_26(ii).runoff];
% end
% 
% % get raw runoff array
% for ii=1:length(glaciers),
%     % rcp2.6
%     for jj=1:length(runoff_26),
%         if ismember(glaciers(ii).rignotid,runoff_26(jj).rignotGlacierID),
%             t_Q(ii,:) = runoff_26(jj).time+0.5; % annual mean convention
%             for kk=1:length(t_Q(ii,:)),
%                 yr = floor(t_Q(ii,kk));
%                 Q(ii,kk) = 3*runoff_26(jj).runoff(kk)/(86400*sum(eomday(yr,[6:8])));
%             end
%         end
%     end
% end
% 
% % get raw TF array
% % load MIROC5
% load ../CMIP5/MIROC5_RCP26.mat
% for ii=1:length(glaciers),
%     t_TF(ii,:) = year+0.5; % follow convention of annual means
%     TF(ii,:) = TF_basins(glaciers(ii).sectornum,:);
% end
% 
% % make raw melt
% for ii=1:length(glaciers),
%     t_melt(ii,:) = t_Q(ii,:);
%     melt(ii,:) = Q(ii,:).^0.4.*interp1(t_TF(ii,:),TF(ii,:),t_melt(ii,:));
% end
% 
% % baseline bias
% baseline = [1995:2014];
% MIROC5baselineinds = find(ismember(floor(t_melt(1,:)),baseline));
% for ii=1:length(glaciers),
%     MIROC5meltbaseline(ii) = mean(melt(ii,MIROC5baselineinds));
%     bias(ii) = MIROC5meltbaseline(ii) - glaciers(ii).melt.meltbaseline;
% end
% 
% % assign to glaciers
% for ii=1:length(glaciers),
%     glaciers(ii).MIROC5.RCP26.tmelt = t_melt(ii,:);
%     glaciers(ii).MIROC5.RCP26.bias_melt = bias(ii);
%     glaciers(ii).MIROC5.RCP26.melt = melt(ii,:) - glaciers(ii).MIROC5.RCP26.bias_melt;
%     % make sure no melt values less than 0
%     glaciers(ii).MIROC5.RCP26.melt(find(glaciers(ii).MIROC5.RCP26.melt<0)) = 0;
% end
%% calculate future NorESM melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).NorESM.RCP85.tmelt = glaciers(ii).NorESM.RCP85.tJJA;
    glaciers(ii).NorESM.RCP85.melt = (glaciers(ii).NorESM.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).NorESM.RCP85.tTF,glaciers(ii).NorESM.RCP85.TF,glaciers(ii).NorESM.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).NorESM.RCP85.melt(find(glaciers(ii).NorESM.RCP85.melt<0)) = 0;
%     % RCP2.6
%     glaciers(ii).NorESM.RCP26.tmelt = glaciers(ii).NorESM.RCP26.tJJA;
%     glaciers(ii).NorESM.RCP26.melt = (glaciers(ii).NorESM.RCP26.QJJA.^0.4).*...
%         interp1(glaciers(ii).NorESM.RCP26.tTF,glaciers(ii).NorESM.RCP26.TF,glaciers(ii).NorESM.RCP26.tJJA);
%     % make sure no melt values less than 0
%     glaciers(ii).NorESM.RCP26.melt(find(glaciers(ii).NorESM.RCP26.melt<0)) = 0;
end
%% calculate future HadGEM melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).HadGEM.RCP85.tmelt = glaciers(ii).HadGEM.RCP85.tJJA;
    glaciers(ii).HadGEM.RCP85.melt = (glaciers(ii).HadGEM.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).HadGEM.RCP85.tTF,glaciers(ii).HadGEM.RCP85.TF,glaciers(ii).HadGEM.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).HadGEM.RCP85.melt(find(glaciers(ii).HadGEM.RCP85.melt<0)) = 0;
end
%% calculate future CNRM-CM6-1 melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CNRMCM6.ssp585.tmelt = glaciers(ii).CNRMCM6.ssp585.tJJA;
    glaciers(ii).CNRMCM6.ssp585.melt = (glaciers(ii).CNRMCM6.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CNRMCM6.ssp585.tTF,glaciers(ii).CNRMCM6.ssp585.TF,glaciers(ii).CNRMCM6.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CNRMCM6.ssp585.melt(find(glaciers(ii).CNRMCM6.ssp585.melt<0)) = 0;
    % ssp126
    glaciers(ii).CNRMCM6.ssp126.tmelt = glaciers(ii).CNRMCM6.ssp126.tJJA;
    glaciers(ii).CNRMCM6.ssp126.melt = (glaciers(ii).CNRMCM6.ssp126.QJJA.^0.4).*...
        interp1(glaciers(ii).CNRMCM6.ssp126.tTF,glaciers(ii).CNRMCM6.ssp126.TF,glaciers(ii).CNRMCM6.ssp126.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CNRMCM6.ssp126.melt(find(glaciers(ii).CNRMCM6.ssp126.melt<0)) = 0;
end
%% calculate future UKESM1-0-LL melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).UKESM1.ssp585.tmelt = glaciers(ii).UKESM1.ssp585.tJJA;
    glaciers(ii).UKESM1.ssp585.melt = (glaciers(ii).UKESM1.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).UKESM1.ssp585.tTF,glaciers(ii).UKESM1.ssp585.TF,glaciers(ii).UKESM1.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).UKESM1.ssp585.melt(find(glaciers(ii).UKESM1.ssp585.melt<0)) = 0;
end
%% calculate future CESM2 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CESM2.ssp585.tmelt = glaciers(ii).CESM2.ssp585.tJJA;
    glaciers(ii).CESM2.ssp585.melt = (glaciers(ii).CESM2.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CESM2.ssp585.tTF,glaciers(ii).CESM2.ssp585.TF,glaciers(ii).CESM2.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CESM2.ssp585.melt(find(glaciers(ii).CESM2.ssp585.melt<0)) = 0;
end
%% calculate future CNRM-ESM2-1 melt (ssp585 only)

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % ssp585
    glaciers(ii).CNRMESM2.ssp585.tmelt = glaciers(ii).CNRMESM2.ssp585.tJJA;
    glaciers(ii).CNRMESM2.ssp585.melt = (glaciers(ii).CNRMESM2.ssp585.QJJA.^0.4).*...
        interp1(glaciers(ii).CNRMESM2.ssp585.tTF,glaciers(ii).CNRMESM2.ssp585.TF,glaciers(ii).CNRMESM2.ssp585.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CNRMESM2.ssp585.melt(find(glaciers(ii).CNRMESM2.ssp585.melt<0)) = 0;
end
%% calculate future CSIRO melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).CSIRO.RCP85.tmelt = glaciers(ii).CSIRO.RCP85.tJJA;
    glaciers(ii).CSIRO.RCP85.melt = (glaciers(ii).CSIRO.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).CSIRO.RCP85.tTF,glaciers(ii).CSIRO.RCP85.TF,glaciers(ii).CSIRO.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).CSIRO.RCP85.melt(find(glaciers(ii).CSIRO.RCP85.melt<0)) = 0;
end
%% calculate future IPSLCM melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).IPSLCM.RCP85.tmelt = glaciers(ii).IPSLCM.RCP85.tJJA;
    glaciers(ii).IPSLCM.RCP85.melt = (glaciers(ii).IPSLCM.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).IPSLCM.RCP85.tTF,glaciers(ii).IPSLCM.RCP85.TF,glaciers(ii).IPSLCM.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).IPSLCM.RCP85.melt(find(glaciers(ii).IPSLCM.RCP85.melt<0)) = 0;
end
%% calculate future ACCESS melt

% just have to combine runoff and thermal forcing
for ii=1:length(glaciers),
    % RCP8.5
    glaciers(ii).ACCESS.RCP85.tmelt = glaciers(ii).ACCESS.RCP85.tJJA;
    glaciers(ii).ACCESS.RCP85.melt = (glaciers(ii).ACCESS.RCP85.QJJA.^0.4).*...
        interp1(glaciers(ii).ACCESS.RCP85.tTF,glaciers(ii).ACCESS.RCP85.TF,glaciers(ii).ACCESS.RCP85.tJJA);
    % make sure no melt values less than 0
    glaciers(ii).ACCESS.RCP85.melt(find(glaciers(ii).ACCESS.RCP85.melt<0)) = 0;
end
%% save

save glaciers.mat glaciers
