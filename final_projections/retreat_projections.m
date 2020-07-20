% script to do final ISMIP6 retreat projections
clear; close all;

% smoothing - choose backwards = 10 and forwards = 9
% for 20 year centred mean
% backwards
sb = 10;
% forwards
sf = 9;

% zero year
yr0 = 2014;

%% calculate retreat

% load glaciers
load ../glaciers/glaciers.mat
for ii=1:length(glaciers),
    iceflux(ii) = glaciers(ii).iceflux.final;
end
sectors = [glaciers.sectornum];
% sector indices
for l=1:7,
    ids(l).inds = find(sectors==l);
end
% make the 8th sector the largest glacier by ice flux in each region
ids(8).inds = [];
for l=1:7,
    ids(8).inds = [ids(8).inds,find(iceflux==max(iceflux(find(sectors==l))))];
end

% load K samples
% this comes from create_K_samples.m
load Ksamples.mat
N = length(k)/length(glaciers);

% loop over RCPs and model scenarios
% mm = 1 for MIROC5 RCP2.6
% mm = 2 for MIROC5 RCP8.5
% mm = 3 for NorESM RCP8.5
% mm = 4 for HadGEM RCP8.5
% mm = 5 for CSIRO RCP8.5
% mm = 6 for IPSLCM RCP8.5
% mm = 7 for ACCESS RCP8.5
% mm = 8 for CNRM-CM6-1 ssp585
% mm = 9 for CNRM-CM6-1 ssp126
% mm = 10 for CNRM-ESM2-1 ssp585
% mm = 11 for UKESM1-0-LL ssp585
% mm = 12 for CESM2
modelscenario = {'MIROC5_RCP2.6','MIROC5_RCP8.5','NorESM_RCP8.5','HadGEM_RCP8.5','CSIRO_RCP8.5','IPSLCM_RCP8.5','ACCESS_RCP8.5'...
                 'CNRM-CM6-1_ssp585','CNRM-CM6-1_ssp126','CNRM-ESM2-1_ssp585','UKESM1-0-LL_ssp585','CESM2_ssp585'};
sectorname = {'SE','SW','CE','CW','NE','NW','NO'};
for mm = 1:12,

    % make array of projected sample retreat
    if mm == 1,
        t = glaciers(1).MIROC5.RCP26.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP26.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).MIROC5.RCP26.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 2,
        t = glaciers(1).MIROC5.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).MIROC5.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 3,
        t = glaciers(1).NorESM.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).NorESM.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 4,
        t = glaciers(1).HadGEM.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).HadGEM.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 5,
        t = glaciers(1).CSIRO.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).CSIRO.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 6,
        t = glaciers(1).IPSLCM.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).IPSLCM.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 7,
        t = glaciers(1).ACCESS.RCP85.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).ACCESS.RCP85.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 8,
        t = glaciers(1).CNRMCM6.ssp585.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).CNRMCM6.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 9,
        t = glaciers(1).CNRMCM6.ssp126.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).CNRMCM6.ssp126.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 10,
        t = glaciers(1).CNRMESM2.ssp585.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).CNRMESM2.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 11,
        t = glaciers(1).UKESM1.ssp585.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).UKESM1.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    elseif mm == 12,
        t = glaciers(1).CESM2.ssp585.tmelt;
        for ii=1:length(glaciers),
%             MA(ii,:) = glaciers(ii).MIROC5.RCP85.melt-glaciers(ii).melt.meltbaseline;
            MA(ii,:) = movmean(glaciers(ii).CESM2.ssp585.melt,[sb,sf]);
            MA(ii,:) = MA(ii,:) - MA(ii,find(floor(t)==yr0));
        end
    end
    
    MA_N = repmat(MA,N,1);
    sampleretreat = k.*MA_N;

    % take flux-weighted means
    for l=1:8,
        for j=1:N,
            fluxweightedmean(l,j,:) = nansum(iceflux(ids(l).inds)'.*sampleretreat(ids(l).inds+(j-1)*length(glaciers),:))/nansum(iceflux(ids(l).inds));
        end
    end

    % identify the ensemble members which have quartile retreats at 2100
    % and extract retreat timeseries for these members
    for l=1:8,
        fwm0 = fluxweightedmean(l,:,end);
        fwm0_sorted = sort(fwm0);
        id25(l) = min(find(fwm0==fwm0_sorted(floor(0.25*length(fwm0_sorted)))));
        id50(l) = min(find(fwm0==fwm0_sorted(floor(0.50*length(fwm0_sorted)))));
        id75(l) = min(find(fwm0==fwm0_sorted(floor(0.75*length(fwm0_sorted)))));
        fw50(l,:) = fluxweightedmean(l,id50(l),:);
        if fw50(l,end)<=0,
            fw25(l,:) = fluxweightedmean(l,id25(l),:);
            fw75(l,:) = fluxweightedmean(l,id75(l),:);
        elseif fw50(l,end)>0,
            fw75(l,:) = fluxweightedmean(l,id25(l),:);
            fw25(l,:) = fluxweightedmean(l,id75(l),:);
        end
    end

    % assign to output arrays
    Rhigh(mm,:,:) = fw25';
    Rmed(mm,:,:) = fw50';
    Rlow(mm,:,:) = fw75';

end

% permute into more logical order
Rhigh = permute(Rhigh,[1,3,2]);
Rmed = permute(Rmed,[1,3,2]);
Rlow = permute(Rlow,[1,3,2]);

% assign to output arrays
retreat.time = t;
% MIROC5 RCP2.6
retreat.MIROC5.RCP26.low = squeeze(Rlow(1,:,:));
retreat.MIROC5.RCP26.med = squeeze(Rmed(1,:,:));
retreat.MIROC5.RCP26.high = squeeze(Rhigh(1,:,:));
% MIROC5 RCP8.5
retreat.MIROC5.RCP85.low = squeeze(Rlow(2,:,:));
retreat.MIROC5.RCP85.med = squeeze(Rmed(2,:,:));
retreat.MIROC5.RCP85.high = squeeze(Rhigh(2,:,:));
% NorESM RCP8.5
retreat.NorESM.RCP85.low = squeeze(Rlow(3,:,:));
retreat.NorESM.RCP85.med = squeeze(Rmed(3,:,:));
retreat.NorESM.RCP85.high = squeeze(Rhigh(3,:,:));
% HadGEM RCP8.5
retreat.HadGEM.RCP85.low = squeeze(Rlow(4,:,:));
retreat.HadGEM.RCP85.med = squeeze(Rmed(4,:,:));
retreat.HadGEM.RCP85.high = squeeze(Rhigh(4,:,:));
% CSIRO RCP8.5
retreat.CSIRO.RCP85.low = squeeze(Rlow(5,:,:));
retreat.CSIRO.RCP85.med = squeeze(Rmed(5,:,:));
retreat.CSIRO.RCP85.high = squeeze(Rhigh(5,:,:));
% IPSLCM RCP8.5
retreat.IPSLCM.RCP85.low = squeeze(Rlow(6,:,:));
retreat.IPSLCM.RCP85.med = squeeze(Rmed(6,:,:));
retreat.IPSLCM.RCP85.high = squeeze(Rhigh(6,:,:));
% ACCESS RCP8.5
retreat.ACCESS.RCP85.low = squeeze(Rlow(7,:,:));
retreat.ACCESS.RCP85.med = squeeze(Rmed(7,:,:));
retreat.ACCESS.RCP85.high = squeeze(Rhigh(7,:,:));
% CNRM-CM6-1 ssp585
retreat.CNRMCM6.ssp585.low = squeeze(Rlow(8,:,:));
retreat.CNRMCM6.ssp585.med = squeeze(Rmed(8,:,:));
retreat.CNRMCM6.ssp585.high = squeeze(Rhigh(8,:,:));
% CNRM-CM6-1 ssp126
retreat.CNRMCM6.ssp126.low = squeeze(Rlow(9,:,:));
retreat.CNRMCM6.ssp126.med = squeeze(Rmed(9,:,:));
retreat.CNRMCM6.ssp126.high = squeeze(Rhigh(9,:,:));
% CNRM-ESM2-1 ssp585
retreat.CNRMESM2.ssp585.low = squeeze(Rlow(10,:,:));
retreat.CNRMESM2.ssp585.med = squeeze(Rmed(10,:,:));
retreat.CNRMESM2.ssp585.high = squeeze(Rhigh(10,:,:));
% UKESM1-0-LL ssp585
retreat.UKESM1.ssp585.low = squeeze(Rlow(11,:,:));
retreat.UKESM1.ssp585.med = squeeze(Rmed(11,:,:));
retreat.UKESM1.ssp585.high = squeeze(Rhigh(11,:,:));
% CESM2 ssp585
retreat.CESM2.ssp585.low = squeeze(Rlow(12,:,:));
retreat.CESM2.ssp585.med = squeeze(Rmed(12,:,:));
retreat.CESM2.ssp585.high = squeeze(Rhigh(12,:,:));

%% add sector TF for models with ice shelves
for l=1:7,    
    retreat.MIROC5.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MIROC5.RCP85.tTF,glaciers(ids(l).inds(1)).MIROC5.RCP85.TF,retreat.time);
    retreat.MIROC5.RCP26.TF(l,:) = interp1(glaciers(ids(l).inds(1)).MIROC5.RCP26.tTF,glaciers(ids(l).inds(1)).MIROC5.RCP26.TF,retreat.time);
    retreat.NorESM.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).NorESM.RCP85.tTF,glaciers(ids(l).inds(1)).NorESM.RCP85.TF,retreat.time);
    retreat.HadGEM.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).HadGEM.RCP85.tTF,glaciers(ids(l).inds(1)).HadGEM.RCP85.TF,retreat.time);
    retreat.CSIRO.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CSIRO.RCP85.tTF,glaciers(ids(l).inds(1)).CSIRO.RCP85.TF,retreat.time);
    retreat.IPSLCM.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).IPSLCM.RCP85.tTF,glaciers(ids(l).inds(1)).IPSLCM.RCP85.TF,retreat.time);
    retreat.ACCESS.RCP85.TF(l,:) = interp1(glaciers(ids(l).inds(1)).ACCESS.RCP85.tTF,glaciers(ids(l).inds(1)).ACCESS.RCP85.TF,retreat.time);
    retreat.CNRMCM6.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMCM6.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMCM6.ssp585.TF,retreat.time);
    retreat.CNRMCM6.ssp126.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMCM6.ssp126.tTF,glaciers(ids(l).inds(1)).CNRMCM6.ssp126.TF,retreat.time);
    retreat.CNRMESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CNRMESM2.ssp585.tTF,glaciers(ids(l).inds(1)).CNRMESM2.ssp585.TF,retreat.time);
    retreat.UKESM1.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).UKESM1.ssp585.tTF,glaciers(ids(l).inds(1)).UKESM1.ssp585.TF,retreat.time);
    retreat.CESM2.ssp585.TF(l,:) = interp1(glaciers(ids(l).inds(1)).CESM2.ssp585.tTF,glaciers(ids(l).inds(1)).CESM2.ssp585.TF,retreat.time);
end

%% plotting

% retreat time series plot
lspace = 0.06;
rspace = 0.025;
bspace = 0.055;
tspace = 0.03;
hspace = 0.03;
vspace = 0.06;
pw = (1-lspace-rspace-3*hspace)/4;
ph = (1-bspace-tspace-vspace)/2;
lwthick = 1;
lwthin = 0.5;
transp = 0.3;
fs = 8;

% colours
cols = [0.000,0.447,0.741;...
        0.850,0.325,0.098;...
        0.929,0.694,0.125;...
        0.494,0.184,0.556;...
        0.466,0.674,0.188;...
        0.301,0.745,0.933;...
        0.635,0.078,0.184;...
        0.000,0.000,0.000;...
        0.500,0.500,0.500;...
        1,0,0;...
        0,1,0;...
        0,0,1];

figure();
a(1) = axes('position',[lspace+0*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(2) = axes('position',[lspace+1*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(3) = axes('position',[lspace+2*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(4) = axes('position',[lspace+3*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(5) = axes('position',[lspace+0*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(6) = axes('position',[lspace+1*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(7) = axes('position',[lspace+2*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);

for l=1:7,
    
    axes(a(l)); hold on;
    plot(t,0*t,'--','linewidth',lwthin,'color',0.5*[1,1,1]);
    % MIROC5 RCP2.6
    plot(t,retreat.MIROC5.RCP26.high(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.MIROC5.RCP26.med(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.MIROC5.RCP26.low(l,:),'color',cols(1,:),'linewidth',lwthick);
    % MIROC5 RCP8.5
    plot(t,retreat.MIROC5.RCP85.high(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.MIROC5.RCP85.med(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.MIROC5.RCP85.low(l,:),'color',cols(2,:),'linewidth',lwthick);
    % NorESM RCP8.5
    plot(t,retreat.NorESM.RCP85.high(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.NorESM.RCP85.med(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.NorESM.RCP85.low(l,:),'color',cols(3,:),'linewidth',lwthick);
    % HadGEM RCP8.5
    plot(t,retreat.HadGEM.RCP85.high(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.HadGEM.RCP85.med(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.HadGEM.RCP85.low(l,:),'color',cols(4,:),'linewidth',lwthick);    
    % CSIRO RCP8.5
    plot(t,retreat.CSIRO.RCP85.high(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.CSIRO.RCP85.med(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.CSIRO.RCP85.low(l,:),'color',cols(5,:),'linewidth',lwthick);    
    % IPSLCM RCP8.5
    plot(t,retreat.IPSLCM.RCP85.high(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.IPSLCM.RCP85.med(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.IPSLCM.RCP85.low(l,:),'color',cols(6,:),'linewidth',lwthick);
    % ACCESS RCP8.5
    plot(t,retreat.ACCESS.RCP85.high(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.ACCESS.RCP85.med(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.ACCESS.RCP85.low(l,:),'color',cols(7,:),'linewidth',lwthick); 
    % CNRM-CM6-1 ssp585
    plot(t,retreat.CNRMCM6.ssp585.high(l,:),'color',cols(8,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.med(l,:),'color',cols(8,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.low(l,:),'color',cols(8,:),'linewidth',lwthick); 
    % CNRM-CM6-1 ssp126
    plot(t,retreat.CNRMCM6.ssp126.high(l,:),'color',cols(9,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp126.med(l,:),'color',cols(9,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp126.low(l,:),'color',cols(9,:),'linewidth',lwthick);
    % CNRM-ESM2-1 ssp585
    plot(t,retreat.CNRMESM2.ssp585.high(l,:),'color',cols(10,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.med(l,:),'color',cols(10,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.low(l,:),'color',cols(10,:),'linewidth',lwthick); 
    % UKESM1-0-LL ssp585
    plot(t,retreat.UKESM1.ssp585.high(l,:),'color',cols(11,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.med(l,:),'color',cols(11,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.low(l,:),'color',cols(11,:),'linewidth',lwthick); 
    % CESM2 ssp585
    plot(t,retreat.CESM2.ssp585.high(l,:),'color',cols(12,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.med(l,:),'color',cols(12,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.low(l,:),'color',cols(12,:),'linewidth',lwthick); 
    xlim([2015 2100]); ylim([-45 5]);
    set(gca,'fontsize',fs,'box','on');
    if l==1 | l==5, ylabel('proj. retreat $\Delta L$ (km)','fontsize',fs); end
    if l<5, set(gca,'xticklabel',[]); end
    if l~=1 & l~=5, set(gca,'yticklabel',[]); end
    text(2020,-23,sectorname{l},'fontsize',fs,'color','k');
    
end

tx = 0.035;
ty = 0.51;
dx = -0.02;
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+1*(ph+vspace),0,0],'string','a','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+1*(ph+vspace),0,0],'string','b','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+1*(ph+vspace),0,0],'string','c','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+3*(pw+hspace),ty+1*(ph+vspace),0,0],'string','d','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+0*(ph+vspace),0,0],'string','e','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+0*(ph+vspace),0,0],'string','f','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+0*(ph+vspace),0,0],'string','g','fontsize',fs,'fontweight','bold','edgecolor','none');

% manual legend
tx = 0.78;
ty = 0.15;
dy = 0.04;
fs0 = 4;
annotation('textbox','position',[tx,ty+0*dy,0,0],'string','MIROC5-RCP2.6','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(1,:));
annotation('textbox','position',[tx,ty+1*dy,0,0],'string','MIROC5-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(2,:));
annotation('textbox','position',[tx,ty+2*dy,0,0],'string','NorESM-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(3,:));
annotation('textbox','position',[tx,ty+3*dy,0,0],'string','HadGEM-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(4,:));
annotation('textbox','position',[tx,ty+4*dy,0,0],'string','CSIRO-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(5,:));
annotation('textbox','position',[tx,ty+5*dy,0,0],'string','IPSLCM-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(6,:));
annotation('textbox','position',[tx,ty+6*dy,0,0],'string','ACCESS-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(7,:));
annotation('textbox','position',[tx,ty+7*dy,0,0],'string','CNRM-CM6-1-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(8,:));
annotation('textbox','position',[tx,ty+8*dy,0,0],'string','CNRM-CM6-1-ssp126','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(9,:));
annotation('textbox','position',[tx,ty+9*dy,0,0],'string','CNRM-ESM2-1-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(10,:));
annotation('textbox','position',[tx,ty+10*dy,0,0],'string','UKESM1-0-LL-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(11,:));
annotation('textbox','position',[tx,ty+11*dy,0,0],'string','CESM2-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(12,:));

saveplot(17,8,300,'retreat_projections.png');
close all;

% TF timeseries plot
figure();
a(1) = axes('position',[lspace+0*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(2) = axes('position',[lspace+1*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(3) = axes('position',[lspace+2*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(4) = axes('position',[lspace+3*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(5) = axes('position',[lspace+0*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(6) = axes('position',[lspace+1*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(7) = axes('position',[lspace+2*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);

for l=1:7,
    
    axes(a(l)); hold on;
    plot(t,retreat.MIROC5.RCP26.TF(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.MIROC5.RCP85.TF(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.NorESM.RCP85.TF(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.HadGEM.RCP85.TF(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.CSIRO.RCP85.TF(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.IPSLCM.RCP85.TF(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.ACCESS.RCP85.TF(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.TF(l,:),'color',cols(8,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp126.TF(l,:),'color',cols(9,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.TF(l,:),'color',cols(10,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.TF(l,:),'color',cols(11,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.TF(l,:),'color',cols(12,:),'linewidth',lwthick);
    xlim([2015 2100]); ylim([1 11]);
    set(gca,'fontsize',fs,'box','on');
    if l==1 | l==5, ylabel('proj. TF ($^{\circ}$C)','fontsize',fs); end
    if l<5, set(gca,'xticklabel',[]); end
    if l~=1 & l~=5, set(gca,'yticklabel',[]); end
    text(2020,9.8,sectorname{l},'fontsize',fs,'color','k');
    
end

tx = 0.035;
ty = 0.51;
dx = -0.02;
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+1*(ph+vspace),0,0],'string','a','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+1*(ph+vspace),0,0],'string','b','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+1*(ph+vspace),0,0],'string','c','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+3*(pw+hspace),ty+1*(ph+vspace),0,0],'string','d','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+0*(ph+vspace),0,0],'string','e','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+0*(ph+vspace),0,0],'string','f','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+0*(ph+vspace),0,0],'string','g','fontsize',fs,'fontweight','bold','edgecolor','none');

% manual legend
tx = 0.78;
ty = 0.5;
dy = 0.03;
fs = 5;
annotation('textbox','position',[tx,ty-0*dy,0,0],'string','MIROC5-RCP2.6','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(1,:));
annotation('textbox','position',[tx,ty-1*dy,0,0],'string','MIROC5-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(2,:));
annotation('textbox','position',[tx,ty-2*dy,0,0],'string','NorESM-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(3,:));
annotation('textbox','position',[tx,ty-3*dy,0,0],'string','HadGEM-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(4,:));
annotation('textbox','position',[tx,ty-4*dy,0,0],'string','CSIRO-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(5,:));
annotation('textbox','position',[tx,ty-5*dy,0,0],'string','IPSLCM-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(6,:));
annotation('textbox','position',[tx,ty-6*dy,0,0],'string','ACCESS-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(7,:));
annotation('textbox','position',[tx,ty-7*dy,0,0],'string','CNRM-CM6-1-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(8,:));
annotation('textbox','position',[tx,ty-8*dy,0,0],'string','CNRM-CM6-1-ssp126','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(9,:));
annotation('textbox','position',[tx,ty-9*dy,0,0],'string','CNRM-ESM2-1-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(10,:));
annotation('textbox','position',[tx,ty-10*dy,0,0],'string','UKESM1-0-LL-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(11,:));
annotation('textbox','position',[tx,ty-11*dy,0,0],'string','CESM2-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(12,:));

saveplot(17,8,300,'TF_projections.png');
close all;

%% save outputs

load ../final_region_def/ice_ocean_sectors.mat

retreat.regions = regions;
for l=1:7,
    retreat.regions(l).name = sectorname{l};
end

save projected_retreat.mat retreat
