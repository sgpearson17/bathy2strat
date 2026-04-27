%% bathyChangeEnvelope_Weighted.m
% determines the age of sediment in different deposits
% v001 - SGP/2020-10-19 - Original
% v002 - SGP/2021-06-10 - Updated for paper to include weighted means
% v003 - SGP/2021-06-16 - Now using 2020 bathy
% v004 - SGP/2021-07-21 - Changed weighting from midpoint to preceding
%                         period as per Helena's suggestion
% v005 - SGP/2021-07-26 - Added contour plots as per Ad's suggestion
% v006 - SGP/2021-08-08 - Updated with 1979, 1982, and 2021
% v008 - SGP/2021-09-29 - NOT with Edwin's new dataset, but with polar vol
% computation

clear all
close all
format compact
clc

% load bathymetric data
% load('..\data\ameland_1926-2019.mat')
% load('..\data\ame_filled-1926-2020.mat')
load('..\data\ame_filled-1926-2021.mat')
outDir = 'D:\PhD\500_Analysis\519_HovmollerAndStratigraphy\results\';

MLW = -1.4; % mean low water (estimated visually based on Terschelling Noordzee station) 
MHW = 1.2; % mean high water (estimated visually based on Terschelling Noordzee station)
% D:\PhD\500_Analysis\512_Hydrodynamics\tidalStations\tidalStationAnalysisRWS_v002_TNZ.m

if exist('bathyChangeEnvelope_Ameland.mat')==2
    load('bathyChangeEnvelope_Ameland.mat')
else
    % create minimum/maximum surface
    X = grd.x./1000; % convert to km
    Y = grd.y./1000;
    Z = -grd.dp; % convert from depth to surface elevation
    T = grd.year;
    % Z(:,:,[18 20]) = []; % get rid of half-year surveys since they are smaller
    % T([18 20]) = [];
    dx = 20;
    iT0=2;
    medianSurf = nanmedian(Z(:,:,iT0:end),3);
    % stdSurf = nanstd(Z(:,:,iT0:end),3); %% NOTE!!!! THIS SHOULD BE WEIGHTED AS NOT ALL TIMESTEPS ARE EQUAL!
    maxSurf = max(Z(:,:,iT0:end),[],3);
    clear grd
    minSurf = min(Z(:,:,iT0:end),[],3);
    meanSurf = nanmean(Z(:,:,iT0:end),3); %% NOTE!!!! THIS SHOULD BE WEIGHTED AS NOT ALL TIMESTEPS ARE EQUAL!
    
    % perform weighted mean surface calc
    weightMidpt = 0;
    if weightMidpt
        Tmidpt(1)=T(iT0);
        for tt = iT0:length(T)-1
            Tmidpt(tt) = (T(tt)+T(tt+1))/2;
        end
        Tmidpt(length(T))=(T(end)+T(end)+1)/2; % assume current year has weight of one year, given annual surveys
        Tweight = diff(Tmidpt);
    else % otherwise just take diff between surveys, which Helena thinks is more correct because bathy represents cumulative erosion/deposition up until that point, and not into the future
        Tweight = [diff(T(iT0:end)) 1];
    end
    
    weightedMeanSurf=zeros(size(X));
    for tt = iT0:length(T)
        weightedMeanSurf = weightedMeanSurf + squeeze(Z(:,:,tt)).*Tweight(tt-1);
    end
    weightedMeanSurf =  weightedMeanSurf./(T(end)-T(iT0));
    
    % estimate volume anomaly
    volumeAnomaly = nansum(nansum((squeeze(Z(:,:,end))-weightedMeanSurf)))*dx^2;
    
end


%% 2020 Bathymetry
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,squeeze(Z(:,:,end)),[-30:0.2:5],'LineStyle','none')
contour(X,Y,Z(:,:,end),[MLW MLW],'Color', [0.5 0.5 0.5],'LineWidth',1.0)
contour(X,Y,Z(:,:,end),[-6 -6],':k','LineWidth',0.5)
% xl=[min(min(X(~isnan(Z(:,:,end))))) max(max(X(~isnan(Z(:,:,end)))))];
% yl=[min(min(Y(~isnan(Z(:,:,end))))) max(max(Y(~isnan(Z(:,:,end)))))];
set(gca,'XTick',[160:177]);
set(gca,'YTick',[605:615]);
axis equal
% set(gca,'XTick',[xLim(1):xLim(2)]);
% set(gca,'YTick',[yLim(1):yLim(2)]);
% xtl = {''}; ytl = {''};
% for ii = 1:length([xLim(1):xLim(2)])
%     xtl{ii,1}=num2str(xLim(1)+ii-1);
%     if mod(ii,2)==0; xtl{ii}=''; end
% end
% for ii = 1:length([yLim(1):yLim(2)])
%     ytl{ii,1}=num2str(yLim(1)+ii-1);
%     if mod(ii,2)==0; ytl{ii}=''; end
% end
% set(gca,'XTickLabel',xtl);
% set(gca,'YTickLabel',ytl);
% xlim([xLim(1) xLim(2)]);
% ylim([yLim(1) yLim(2)]);

set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
% [cmap] = kg2Colormap;
[cmap] = oldWaddenSeaColormap;
colormap(cmap);
clim([-20 10])
axis equal
xlim([160 175]);
ylim([605 615]);
cb=colorbar;
ylabel(cb,'Depth [m]','FontName','Myriad Pro')
title(['Ameland Ebb-Tidal Delta Bathymetry [' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'FontName','Myriad Pro')

figName = ['Bathy' num2str(T(end))];
export_fig([outDir figName],'-png','-r300');

%% plot mean surface
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,meanSurf,[-30:5],'LineStyle','none')
contour(X,Y,meanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,meanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation [m NAP]','FontName','Myriad Pro');
title(['Mean Bathymetric Elevation [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
set(gca,'FontName','Myriad Pro')
figName = 'MeanBathyChange';
export_fig([outDir figName],'-png','-r300');

%% plot weighted mean surface
figure(2)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,weightedMeanSurf,[-30:5],'LineStyle','none')
contour(X,Y,weightedMeanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,weightedMeanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation [m NAP]','FontName','Myriad Pro');
title(['Weighted Mean Bathymetric Elevation [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
set(gca,'FontName','Myriad Pro')
figName = 'WeightedMeanBathyChange';
export_fig([outDir figName],'-png','-r300');

%% plot difference weighted mean surface
figure(3)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,weightedMeanSurf-meanSurf,[-3:0.1:3],'LineStyle','none')
contour(X,Y,weightedMeanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,weightedMeanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation Difference [m]','FontName','Myriad Pro');
title(['Weighted Mean - Unweighted Mean  [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
colormap(parula);
clim([-3 3])
set(gca,'FontName','Myriad Pro')
figName = 'WeightedMeanDiffBathyChange';
export_fig([outDir figName],'-png','-r300');

%% median
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,medianSurf,[-30:5],'LineStyle','none')
contour(X,Y,medianSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,medianSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation [m NAP]','FontName','Myriad Pro');
title(['Median Bathymetric Elevation  [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
set(gca,'FontName','Myriad Pro')
figName = 'MedianBathyChange';
export_fig([outDir figName],'-png','-r300');

%% min
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,minSurf,[-30:5],'LineStyle','none')
contour(X,Y,minSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,minSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation [m NAP]','FontName','Myriad Pro');
title(['Min Bathymetric Elevation  [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
set(gca,'FontName','Myriad Pro')
figName = 'MinBathyChange';
export_fig([outDir figName],'-png','-r300');

%% max
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,maxSurf,[-30:5],'LineStyle','none')
contour(X,Y,maxSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,maxSurf,[-6 -6],':k','LineWidth',0.5)
% contour(X,Y,maxSurf,[-6 -3],':r','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation [m NAP]','FontName','Myriad Pro');
title(['Max Bathymetric Elevation  [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
set(gca,'FontName','Myriad Pro')
figName = 'MaxBathyChange';
export_fig([outDir figName],'-png','-r300');


%% Envelope of Bathymetric Change [1975-2020]
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,maxSurf-minSurf,[0:0.1:10],'LineStyle','none')
contour(X,Y,weightedMeanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,weightedMeanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Elevation Difference [m]','FontName','Myriad Pro');
title(['Envelope of Bathymetric Change [1975-' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
clim([0 10])
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
colormap(parula);
set(gca,'FontName','Myriad Pro')
figName = 'EnvelopeOfBathyChange';
export_fig([outDir figName],'-png','-r300');

envelopeVol = nansum(nansum((maxSurf-minSurf)))*dx^2;

%% active delta thickness (2020)
figure(1)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,squeeze(Z(:,:,end))-minSurf,[0:0.1:10],'LineStyle','none')
contour(X,Y,weightedMeanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,weightedMeanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Thickness [m]','FontName','Myriad Pro');
title(['Active Delta Thickness  [' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
clim([0 10])
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
colormap(parula);
set(gca,'FontName','Myriad Pro')
figName = 'ActiveDeltaThickness2020';
export_fig([outDir figName],'-png','-r300');


activeDeltaVol = nansum(nansum((squeeze(Z(:,:,end))-minSurf)))*dx^2;

%% surface anomaly (2020)
surfaceAnomaly = squeeze(Z(:,:,end))-weightedMeanSurf;

figure(5)
clf
set(gcf,'Color','w');
hold on; box on; grid on;
contourf(X,Y,surfaceAnomaly,[-20:0.5:20],'LineStyle','none')
contour(X,Y,weightedMeanSurf,[MLW MLW],'LineColor','k','LineWidth',1.5)
contour(X,Y,weightedMeanSurf,[-6 -6],':k','LineWidth',0.5)
axis equal
cb=colorbar('southoutside');
ylabel(cb,'Surface Anomaly [m]','FontName','Myriad Pro');
title(['Surface Anomaly [' num2str(T(end)) ']'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'XTick',[160:175]);
set(gca,'YTick',[605:615]);
xlim([160 175]);
ylim([605 615]);
clim([-5 5])
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.3);
colormap(parula);
set(gca,'FontName','Myriad Pro')
figName = 'SurfaceAnomaly2020';
export_fig([outDir figName],'-png','-r300');

% Save everything
save('bathyChangeEnvelope_Ameland.mat');

% %% plot contours for Ad
% 
% plotContours = [-24:2:2];
% figure(999)
% set(gcf,'Color','w');
% % loop through all years
% for ii = 1:length(plotContours)
%     clf
%     hold on; box on; grid on;
%     contourf(X,Y,weightedMeanSurf,[-30:0.2:5],'LineStyle','none')
%     for tt=iT0:length(T)
%         contour(X,Y,squeeze(Z(:,:,tt)),[plotContours(ii) plotContours(ii)], 'Color', [0.5 0.5 0.5],'LineWidth',0.5)
%     end
%     contour(X,Y,maxSurf,[plotContours(ii) plotContours(ii)],':k','LineWidth',1.0)
%     contour(X,Y,minSurf,[plotContours(ii) plotContours(ii)],'-k','LineWidth',1.0)
%     
%     set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
%     [cmap] = kg2Colormap;
%     colormap(cmap);
%     clim([-20 10])
%     axis equal
%     set(gca,'XTick',[160:175]);
%     set(gca,'YTick',[605:615]);
%     xlim([160 175]);
%     ylim([605 615]);
%     set(gca,'XTickLabelMode','manual','XTickMode','manual','YTickLabelMode','manual','YTickMode','manual');
%     cb=colorbar;
%     ylabel(cb,'Depth [m]','FontName','Myriad Pro')
%     title([num2str(plotContours(ii)) ' m Contour'])
%     xlabel('Easting [km]');
%     ylabel('Northing [km]');
%     set(gca,'FontName','Myriad Pro')
%     
%     % Export figure to png file
%     pngName = ['Contours_' num2str(plotContours(ii),'%03.1f')];
%     dimensions = [29.7 21]./2; % [width height]
%     printFigs(dimensions,[outDir filesep pngName],0);
% end

%% VOL CALC RESTRICTED TO POLAR GRID

% load polar data
load('polarData.mat')
thetaDeg = rad2deg(theta);

% coordinate correction for Polar system
x0=169.5;
y0=607;
XP=X-x0;
YP=Y-y0;
dx=20;

plot([xPolar(1,:) xPolar(1,1)],[yPolar(1,:) yPolar(1,1)],'-k')

% check which XY bathy coordinates lie inside a given grid cell
[in] = inpolygon(XP,YP,[xPolar(end,:) xPolar(end,1)],[yPolar(end,:) yPolar(end,1)]);

activeDeltaThk_presentDay = Z(:,:,end)-minSurf;
activeDeltaVol_Polar = nansum(nansum(squeeze(activeDeltaThk_presentDay(in))))*dx^2

envelopeVol_presentDay = maxSurf-minSurf;
envelopeVol_Polar = nansum(nansum(envelopeVol_presentDay(in)))*dx^2

netVol_presentDay = Z(:,:,end)-Z(:,:,9);
netVol_Polar = nansum(nansum(netVol_presentDay(in)))*dx^2

grossPos = nansum(netVol_presentDay(netVol_presentDay>0))*dx^2;
grossNeg = nansum(netVol_presentDay(netVol_presentDay<0))*dx^2;
grossTotal = grossPos + abs(grossNeg)

close all
figure
subplot(2,1,1)
scatter(XP(in),YP(in),10,envelopeVol_presentDay(in),'filled')
axis equal
colorbar
subplot(2,1,2)
scatter(XP(in),YP(in),10,activeDeltaThk_presentDay(in),'filled')
axis equal
colorbar
