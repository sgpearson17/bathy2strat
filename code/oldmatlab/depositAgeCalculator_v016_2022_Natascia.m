%% depositAgeCalculator.m
% determines the age of sediment in different deposits
% v001 - SGP/2020-10-16 - Original
% v004 - SGP/2020-10-21 - Tidied Up
% v005 - SGP/2020-10-22 - Cleaned up into functions
% v006 - SGP/2020-10-22 - Added interpolation function
% v007 - SGP/2020-10-22 - Added GUI input
% v008 - SGP/2021-06-08 - Updates based on Helena's comments
% v009 - SGP/2021-06-10 - Revised transect plotting
% v010 - SGP/2021-06-11 - Made layers of zero thickness for intermediate years
% v011 - SGP/2021-06-11 - Fixed scaling of printed plots
% v012 - SGP/2021-06-15 - Added little black lines to stratigraphy XSs
% v013 - SGP/2021-06-16 - Now loading 2020 data
% v014 - SGP/2021-08-06 - With Frame 2 in official code and with parula as standard
% v016 - SGP/2021-08-08 - Updated with 1979, 1982, and 2021
%        SGP/2023-11-01 - Updated for Natascia Pannozzo with 2022 bathy

clear all
close all
format compact
clc

addpath(genpath('C:\oet\')); % open earth tools

% load bathymetric data
% load('..\data\ameland_1926-2019.mat')
% load('..\data\ame_filled-1926-2020.mat')
% load('..\data\ame_filled-1926-2021.mat')
load('..\data\ame_filled-1926-2022_stuart.mat')

MLW = -1.4; % mean low water (estimated visually based on Terschelling Noordzee station) 
MHW = 1.2; % mean high water (estimated visually based on Terschelling Noordzee station)
% D:\PhD\500_Analysis\512_Hydrodynamics\tidalStations\tidalStationAnalysisRWS_v002_TNZ.m
cm=parula(5); deepBlue = cm(1,:);

%% Natascia's transect locations
% AM23-1: 53’24.641; 005’40.380
% AM23-3: 53%24.498; 005 40.853

transectXY_raw = [5+40.380/60	53+24.641/60  5+40.853/60	53+24.498/60;...
    5.655843	53.392967	5.695531	53.424770];

[xStart,yStart,~]=convertCoordinates(transectXY_raw(:,1),transectXY_raw(:,2),'CS1.code',4326,'CS2.code',28992);
[xEnd,yEnd,~]=convertCoordinates(transectXY_raw(:,3),transectXY_raw(:,4),'CS1.code',4326,'CS2.code',28992);

%% create initial variables
initialDate = 2; % 13=2010
X = grd.x./1000; % convert to km
Y = grd.y./1000;
Z = -grd.dp; % convert from depth to surface elevation

% % optional: remove specific years
% removeTimesteps = logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0]);
Z(:,:,3)=[];
T = grd.year;
T([3])=[];
clear grd

% set any NaN cell from one year equal to NaN in all years
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if sum(isnan(Z(ii,jj,:)))>0
            Z(ii,jj,:) = NaN;
        end
    end
end

% calculate the maximum and minimum surfaces
minSurf = min(Z(:,:,initialDate:end),[],3);
maxSurf = max(Z(:,:,initialDate:end),[],3);


%% calculate deposition time series
surfaceDate = ones(size(Z(:,:,1))) .* T(initialDate);
depositElev = nan(size(Z));
depositElev(:,:,initialDate) = Z(:,:,initialDate);
legendEntries{1} = num2str(T(initialDate));

% loop through each timestep and determine elevation of each deposit
for tt = initialDate+1:size(Z,3)
    
    % check for erosion or accretion of the surface
    dz= Z(:,:,tt)-Z(:,:,tt-1);
    
    % loop through array and find cells with erosion
    for ii = 1:size(X,1)
        for jj = 1:size(X,2)
            
            % assume accreting
            depositElev(ii,jj,tt) = Z(ii,jj,tt);
            
            if dz(ii,jj) < 0 % EROSION
                erosion = dz(ii,jj); %Z(ii,jj,tt)-Z(ii,jj,tt-1)
                qq = tt-1; % timestep for marching backward in time

                for qq = 1:((tt-1)-initialDate)+1
                    if ~isnan(depositElev(ii,jj,tt-qq))
                        if depositElev(ii,jj,tt) < depositElev(ii,jj,tt-qq)
                            depositElev(ii,jj,tt-qq) = depositElev(ii,jj,tt);
                        end
                    end
                end
                depositElev(ii,jj,:);
            end
            
        end
    end
    
    % assign legend entries
    legendEntries{tt-initialDate+1} = num2str(T(tt));
end

%% convert deposit elevations to deposit thicknesses
depositThk = diff(depositElev,1,3);
depositThk(:,:,1:initialDate-1) = [];
depositThk = cat(3,[minSurf],depositThk);

% add layer of zero thickness for each intermediate year
T_full = T(2):T(end);
yearCount = initialDate;
depositThkFull =  zeros([size(X) length(T_full)]);
for tt = 1:length(T_full)
    if T_full(tt) == T(yearCount)
        depositThkFull(:,:,tt) = depositThk(:,:,yearCount-1);
        yearCount = yearCount + 1;
    end
end
depositThkOriginal = depositThk; 
depositThk = depositThkFull; 

%% MEAN SURFACE AGE
% calculate weighted mean age of top 30 cm of sediment (e.g. in a box core)

coreDepth = 1.0; % depth of sample/core [m]
[meanSampleAge] = calcMeanSurfaceAge(depositElev,Z(:,:,end),T,coreDepth);

% SURFACE DATE MAP
figure(3)
set(gcf, 'Color','w');
clf
cmap=parula(length(T_full));
% cmap(1:8,:)=[];
colormap(flipud(cmap))

hold on; box on; grid on;
% contourf(X,Y,meanSampleAge,[0:40],'LineStyle','none')
contourf(X,Y,meanSampleAge,[0:40],'LineStyle','none')
contour(X,Y,Z(:,:,end),[MLW MLW],'k','LineWidth',1.5)
% contour(X,Y,Z(:,:,end),[MHW MHW],'w','LineWidth',1.5)
contour(X,Y,Z(:,:,end),[-6 -6],':k','LineWidth',0.5)
% surf(X,Y,meanSampleAge,'FaceColor','interp','LineStyle','none')
% contour3(X,Y,Z(:,:,end)+1000,[MLW MLW],'k','LineWidth',1.0)
% contour3(X,Y,Z(:,:,end)+1000,[-6 -6],':k','LineWidth',0.5)
set(gca,'XTick',[165:185]);
set(gca,'YTick',[600:610]);
axis equal
cb=colorbar;
caxis([0 T(end)-T(initialDate)])
% title(['Mean Surface Age (Top ' num2str(coreDepth) ' m Sample)'])
xlabel('Easting [km]');
ylabel('Northing [km]');
% xlim([min(min(X(~isnan(Z(:,:,end))))) max(max(X(~isnan(Z(:,:,end)))))]);
% ylim([min(min(Y(~isnan(Z(:,:,end))))) max(max(Y(~isnan(Z(:,:,end)))))]);
xlim([160 185]);
ylim([600 610]);
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
ylabel(cb,'Deposit Age [y]','FontName','Myriad Pro')
set(gca,'FontName','Myriad Pro')

figName = ['Ameland_MeanSurfaceAge(Top' num2str(coreDepth) 'mSample)'];
export_fig(['..\results\Natascia\'  figName],'-png','-r300');

%% Correlation between surface age and deposit thickness

minSurf = min(Z(:,:,2:end),[],3);

figure(747)
set(gcf, 'Color','w');
clf
hold on; box on; grid on;
% scatter(reshape((squeeze(Z(:,:,end))-minSurf),[],1),reshape(meanSampleAge,[],1),10,'k')

[N,c]=hist3([reshape((squeeze(Z(:,:,end))-minSurf),[],1),...
    reshape(meanSampleAge,[],1)],'Centres',{0:2:20 0:5:45},'CdataMode','auto');
dx=20; %dx=20 m
pcolor(c{1,1},c{1,2},log10(N*dx^2))
clim([0 8])
cb=colorbar;
ylabel(cb,'log_{10}(Area) [log_{10}(m^2)]','FontName','Myriad Pro')
view(2)
set(gca,'FontName','Myriad Pro')
xlabel('Active Delta Thickness (2021) [m]');
ylabel('Mean Surface Sample Age (0.9 m) [years]');

figName = ['Surface Age vs Active Deposit Thickness'];
export_fig(['..\results\Natascia\'  figName],'-png','-r300');

% %% load polar data
% load('polarData.mat')
% thetaDeg = rad2deg(theta);

%% coordinate correction for Polar system
x0=0;%169.5;
y0=0;%607;
XP=X-x0;
YP=Y-y0;

%% Make big plot that shows XS locations
figure(3333)
set(gcf, 'Color','w');
clf

% colormap(parula)
hold on; box on; grid on;
% contourf(X,Y,Z(:,:,end),[-30:5],'LineStyle','none')
pc=pcolor(X,Y,Z(:,:,end));
pc.FaceColor = 'flat';
pc.EdgeColor = 'none';
contour(X,Y,Z(:,:,end),[0 0],'--k','LineWidth',1.0)
contour(X,Y,Z(:,:,end),[-6 -6],':k','LineWidth',0.5)
% xl=[min(min(X(~isnan(Z(:,:,end))))) max(max(X(~isnan(Z(:,:,end)))))];
% yl=[min(min(Y(~isnan(Z(:,:,end))))) max(max(Y(~isnan(Z(:,:,end)))))];
set(gca,'XTick',[160:190]);
set(gca,'YTick',[600:610]);

set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
axis equal
xlim([164 179]);
ylim([600 610]);
cb=colorbar;
ylabel(cb,'Depth [m]','FontName','Myriad Pro')
title(['Transect Locations (' num2str(T(end)) ' Bathymetry)'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'FontName','Myriad Pro')

% cross-sections for paper

% section AA: 
transectLetter = 'A';
curvy = 0; % is the transect curved (1) or straight (0)
transectFromOrigin=0; % is this a straight transect from the origin?
if transectFromOrigin
    transectLetter1 = 'O';
    transectLetter2 = transectLetter;
else
    transectLetter1 = transectLetter;
    transectLetter2 = transectLetter;
end

% pick out transects
n=1;
xx = [xStart(n) xEnd(n)]./1000-x0;
yy = [yStart(n) yEnd(n)]./1000-y0;

% resample along cross-section
q=1; % resampling factor
xq = resample(xx,q,1); % note: this command requires signal toolbox
yq = resample(yy,q,1);

% distance from origin along cross-section
if transectFromOrigin
    % % interpolate along cross-section
    np=500; % number of points along cross-section
    xq=linspace(xx(1),xx(end),np); % x-coordinates along cross-section
    yq=linspace(yy(1),yy(end),np); % y-coordinates along cross-section
    dq=sqrt((xq-xx(1)).^2+(yq-yy(1)).^2); % distance from origin along cross-section
else
    if curvy
        dq = nan(size(xq));
        dq(1)=0;
        for ii = 2:length(dq)
            dq_int=sqrt((xq(ii)-xq(ii-1)).^2+(yq(ii)-yq(ii-1)).^2);
            dq(ii) = dq(ii-1)+dq_int;
        end
    else
        np=500; % number of points along cross-section
        xq=linspace(xx(1),xx(end),np); % x-coordinates along cross-section
        yq=linspace(yy(1),yy(end),np); % y-coordinates along cross-section
        dq=sqrt((xq-xx(1)).^2+(yq-yy(1)).^2); % distance from origin along cross-section
        dq=dq-dq(1);
    end
end

Zq = interp2(XP(:,1),YP(1,:),Z(:,:,end)',xq,yq); % surface elevation along cross-section 
Zq_max = interp2(XP(:,1),YP(1,:),maxSurf',xq,yq); % surface elevation along cross-section 

% interpolate deposit thicknesses
clear depThkXS
for qq = 1:size(depositThk,3)
    depThkXS(qq,:) = interp2(XP(:,1),YP(1,:),depositThk(:,:,qq)',xq,yq); % deposit thickness along cross-section 
end
%--------------------------------------------
% plot on XS location map
figure(3333)
hold on;

plot(xq+x0,yq+y0,'-k')
scatter(xq(1)+x0,yq(1)+y0,10,'w','filled','MarkerEdgeColor','k');
for ii = 2:length(dq)
    if floor(dq(ii))>floor(dq(ii-1))
        scatter(xq(ii)+x0,yq(ii)+y0,5,'w','filled','MarkerEdgeColor','k');
    end
end
scatter(xq(end)+x0,yq(end)+y0,10,'w','filled','MarkerEdgeColor','k');
if ~transectFromOrigin
    text(xq(1)+x0-0.4,yq(1)+y0+0.3,transectLetter,'FontWeight','bold','FontName','Myriad Pro','Color','k');
end
text(xq(end)+x0+0.15,yq(end)+y0-0.3,[transectLetter ''''],'FontWeight','bold','FontName','Myriad Pro','Color','k');
%------------------------------------------------------------

% stratigraphic cross-section
figure(30)
set(gcf, 'Color','w');
clf
cmap=parula(length(T_full));
% cmap(1:8,:)=[];
colormap(cmap)
hold on; box on; grid on;
area(dq,depThkXS',-35,'LineWidth',0.2,'FaceColor','flat'); % plot stratigraphy
plot(dq,depThkXS(1,:)','-k', 'LineWidth',1.0); % min surface elev
plot(dq,Zq_max,':k', 'LineWidth',1.0); % max surface elev
plot(dq,Zq,'-k', 'LineWidth',1.5); % most recent surface elev
XT=get(gca,'XTick');
set(gca,'XTick',[XT(1):0.1:XT(end)]);
set(gca,'YTick',[-30:2:20]);
cb=colorbar('southoutside');
title(['Cross-section ' transectLetter1 '-' transectLetter2 '''']);
xlabel('Distance [km]')
ylabel('Elevation [m NAP]')
ylim([min(min(depThkXS))-1 max([max(Zq)+1 3])]);
ylabel(cb,'Deposit Age [y]','FontName','Myriad Pro')
set(cb,'YTick',0:5:size(depositThk,3)+1,'TickLabels',[size(depositThk,3)-(0:5:size(depositThk,3))]);
set(gca,'FontName','Myriad Pro')
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
MLW_line=plot([XT(1) XT(end)],[MLW MLW],'--','Color',deepBlue,'LineWidth',0.5); % plot MLW
MHW_line=plot([XT(1) XT(end)],[MHW MHW],'--','Color',deepBlue,'LineWidth',0.5); % plot MHW
uistack(MHW_line,'bottom'); uistack(MLW_line,'bottom');
text(6.1,MHW+0.5,'MHW','FontName','Myriad Pro','Color',deepBlue);
text(6.1,MLW+0.5,'MLW','FontName','Myriad Pro','Color',deepBlue);
xlim([0 dq(find(~isnan(Zq),1,'last'))]);

% Export figure to png file
pngName = ['Cross-section ' transectLetter1 '-' transectLetter2 ''''];
maxTransectLength=1; maxDepthRange=6; 
maxFigWidth=25; maxFigHeight=5;
dimensions = [diff(get(gca,'XLim'))/maxTransectLength*maxFigWidth...
    diff(get(gca,'YLim'))/maxDepthRange*maxFigHeight].*2; % [width height]
printFigs(dimensions,['..\results\Natascia\' pngName],0);

clear xq yq dq
%
% =======================================================================
% section BB': 
transectLetter = 'B';
curvy = 0; % is the transect curved (1) or straight (0)
transectFromOrigin=0; % is this a straight transect from the origin?
if transectFromOrigin
    transectLetter1 = 'O';
    transectLetter2 = transectLetter;
else
    transectLetter1 = transectLetter;
    transectLetter2 = transectLetter;
end

% pick out transects
n=2;
xx = [xStart(n) xEnd(n)]./1000-x0;
yy = [yStart(n) yEnd(n)]./1000-y0;

% resample along cross-section
q=1; % resampling factor
xq = resample(xx,q,1); % note: this command requires signal toolbox
yq = resample(yy,q,1);

% distance from origin along cross-section
if transectFromOrigin
    % % interpolate along cross-section
    np=500; % number of points along cross-section
    xq=linspace(xx(1),xx(end),np); % x-coordinates along cross-section
    yq=linspace(yy(1),yy(end),np); % y-coordinates along cross-section
    dq=sqrt((xq-xx(1)).^2+(yq-yy(1)).^2); % distance from origin along cross-section
else
    if curvy
        dq = nan(size(xq));
        dq(1)=0;
        for ii = 2:length(dq)
            dq_int=sqrt((xq(ii)-xq(ii-1)).^2+(yq(ii)-yq(ii-1)).^2);
            dq(ii) = dq(ii-1)+dq_int;
        end
    else
        np=500; % number of points along cross-section
        xq=linspace(xx(1),xx(end),np); % x-coordinates along cross-section
        yq=linspace(yy(1),yy(end),np); % y-coordinates along cross-section
        dq=sqrt((xq-xx(1)).^2+(yq-yy(1)).^2); % distance from origin along cross-section
        dq=dq-dq(1);
    end
end

Zq = interp2(XP(:,1),YP(1,:),Z(:,:,end)',xq,yq); % surface elevation along cross-section 
Zq_max = interp2(XP(:,1),YP(1,:),maxSurf',xq,yq); % surface elevation along cross-section 

% interpolate deposit thicknesses
clear depThkXS
for qq = 1:size(depositThk,3)
    depThkXS(qq,:) = interp2(XP(:,1),YP(1,:),depositThk(:,:,qq)',xq,yq); % deposit thickness along cross-section 
end
%--------------------------------------------
% plot on XS location map
figure(3333)
hold on;

plot(xq+x0,yq+y0,'-k')
scatter(xq(1)+x0,yq(1)+y0,10,'w','filled','MarkerEdgeColor','k');
for ii = 2:length(dq)
    if floor(dq(ii))>floor(dq(ii-1))
        scatter(xq(ii)+x0,yq(ii)+y0,5,'w','filled','MarkerEdgeColor','k');
    end
end
scatter(xq(end)+x0,yq(end)+y0,10,'w','filled','MarkerEdgeColor','k');
if ~transectFromOrigin
    text(xq(1)+x0-0.4,yq(1)+y0+0.3,transectLetter,'FontWeight','bold','FontName','Myriad Pro','Color','k');
end
text(xq(end)+x0+0.15,yq(end)+y0-0.3,[transectLetter ''''],'FontWeight','bold','FontName','Myriad Pro','Color','k');
%------------------------------------------------------------

% stratigraphic cross-section
figure(30)
set(gcf, 'Color','w');
clf
cmap=parula(length(T_full));
% cmap(1:8,:)=[];
colormap(cmap)
hold on; box on; grid on;
area(dq,depThkXS',-35,'LineWidth',0.2,'FaceColor','flat'); % plot stratigraphy
plot(dq,depThkXS(1,:)','-k', 'LineWidth',1.0); % min surface elev
plot(dq,Zq_max,':k', 'LineWidth',1.0); % max surface elev
plot(dq,Zq,'-k', 'LineWidth',1.5); % most recent surface elev
XT=get(gca,'XTick');
set(gca,'XTick',[XT(1):0.5:XT(end)]);
set(gca,'YTick',[-30:2:20]);
% cb=colorbar('southoutside');
title(['Cross-section ' transectLetter1 '-' transectLetter2 '''']);
xlabel('Distance [km]')
ylabel('Elevation [m NAP]')
ylim([min(min(depThkXS))-1 max([max(Zq)+1 3])]);
% ylabel(cb,'Deposit Age [y]','FontName','Myriad Pro')
% set(cb,'YTick',0:5:size(depositThk,3)+1,'TickLabels',[size(depositThk,3)-(0:5:size(depositThk,3))]);
set(gca,'FontName','Myriad Pro')
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
MLW_line=plot([XT(1) XT(end)],[MLW MLW],'--','Color',deepBlue,'LineWidth',0.5); % plot MLW
MHW_line=plot([XT(1) XT(end)],[MHW MHW],'--','Color',deepBlue,'LineWidth',0.5); % plot MHW
uistack(MHW_line,'bottom'); uistack(MLW_line,'bottom');
text(6.1,MHW+0.5,'MHW','FontName','Myriad Pro','Color',deepBlue);
text(6.1,MLW+0.5,'MLW','FontName','Myriad Pro','Color',deepBlue);
xlim([0 dq(find(~isnan(Zq),1,'last'))]);

% Export figure to png file
pngName = ['Cross-section ' transectLetter1 '-' transectLetter2 ''''];
maxTransectLength=10; maxDepthRange=30; 
maxFigWidth=25; maxFigHeight=5;
dimensions = [diff(get(gca,'XLim'))/maxTransectLength*maxFigWidth...
    diff(get(gca,'YLim'))/maxDepthRange*maxFigHeight].*2; % [width height]
printFigs(dimensions,['..\results\Natascia\' pngName],0);

clear xq yq dq
%
% =======================================================================

%==========================================================================
% print map
figure(3333)

% % plot origin
% transectFromOrigin=0; % is this a straight transect from the origin?
% if transectFromOrigin
%     scatter(x0,y0,20,'w','filled','MarkerEdgeColor','k');
%     text(x0-0.4,y0-0.4,'O','FontWeight','bold','FontName','Myriad Pro','Color','w');
% end

% % plot frame 2!
% scatter(xyFrame2(1),xyFrame2(2),20,'^y','filled','MarkerEdgeColor','k')

% Export figure to png file
pngName = ['Stratigraphic Transect Map - Natascia'];
dimensions = [29.7 21]./1.2; % [width height]
printFigs(dimensions,['..\results\Natascia\' pngName],0);

