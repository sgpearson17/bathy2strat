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
% v018 - SGP/2022-11-02 - Bogue Inlet data and plotting
% v019- SGP/2023-02-23 - Now layers of zero thickness for every month

clear all
close all
format compact
clc

% load bathymetric data
dataDir = '..\data\';
plotDir = '..\plots\stratigraphy\';
mkdir(plotDir)
load([dataDir filesep 'ame_filled-1926-2021_paperVersion.mat'])

addpath(genpath('C:\oet\')); % open earth tools

MLW = -1.4; % mean low water (estimated visually based on Terschelling Noordzee station) 
MHW = 1.2; % mean high water (estimated visually based on Terschelling Noordzee station)
% D:\PhD\500_Analysis\512_Hydrodynamics\tidalStations\tidalStationAnalysisRWS_v002_TNZ.m
cm=parula(5); deepBlue = cm(1,:);

%% Tjitske's transect locations
transectXY_raw = [5.5559645	53.4370181	5.6096702	53.4450156;...
5.5640548	53.4176334	5.5757599	53.4193774;...
5.609333	53.424667	5.640333	53.424667;...
5.6113916	53.4230702	5.6339411	53.4275833;...
5.6196325	53.4155944	5.6615901	53.4154341;...
5.6131129	53.3982399	5.6676794	53.4122985;...
5.6248333	53.3975	5.7646667	53.3975;...
5.6707347	53.3917028	5.6844516	53.4085406];

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

% add layer of zero thickness for each intermediate month

% get january of second surveyed year
T_startMonth = datenum(year(T(1)),1,1); % second year because first year serves as baseline for comparison
% get december of last surveyed year
T_endMonth = datenum(year(T(end)),12,31);

% loop through and create time series
monthCount = 1;
for yy = year(T(1)):year(T(end))
    for mm = 1:12
        T_full(monthCount) = datenum(yy,mm,1);
        monthCount = monthCount + 1;
    end
end

% get total number of months and make array that big
depositThkFull =  zeros([size(X) length(T_full)]); 

surveyCount = 1; % start with first survey
for tt = 1:length(T_full)-1 % loop through ALL months in entire surveyed period
    if T_full(tt) <= T(surveyCount) ...
            && T(surveyCount) <= T_full(tt+1) % if a survey is between this month and next month, add its deposit thickness
        depositThkFull(:,:,tt) = depositThk(:,:,surveyCount);

        if surveyCount >= length(T); break; end % exit loop once last survey is finished

        surveyCount = surveyCount + 1; % go to next month
    end
end
depositThkOriginal = depositThk; 
depositThk = depositThkFull; 

% %% load polar data
% load('polarData.mat')
% thetaDeg = rad2deg(theta);
% 
% coordinate correction for Polar system
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
set(gca,'XTick',[160:180]);
set(gca,'YTick',[600:610]);

set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
[cmap] = kg2Colormap;
colormap(cmap);
clim([-20 10])
axis equal
xlim([160 180]);
ylim([600 610]);
cb=colorbar;
ylabel(cb,'Depth [m]','FontName','Myriad Pro')
title(['Transect Locations (' num2str(T(end)) ' Bathymetry)'])
xlabel('Easting [km]');
ylabel('Northing [km]');
set(gca,'FontName','Myriad Pro')

% plot settings for transect dots and labels
dotSpacing = 0.1; % intermediate dot spacing along transects in location plan [km]
endDotSize = 40;
midDotSize = 15;
endLetterSize = 16;

% ===================================================================
% cross-sections for paper

% section AA: 
transectLetter = 'GUI';
curvy = 1; % is the transect curved (1) or straight (0)
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

% choose points and cut cross-section
disp('Left-click points on the map to define a transect and right-click to finish')
pp = 1;
button = 1;
while button == 1 % stop if you click anything other than the left mouse button
    [xx(pp),yy(pp),button] = ginput(1); % select point on the map
    plot(xx(pp),yy(pp),'xk'); % plot that point
    pp = pp + 1; % next point
end

xx = xx-x0;
yy = yy-y0;
pt_raw = [xx' yy'];

% remove duplicate points since interparc doesn't like that
% (i.e., if you click the same point twice by accident)
pt_clean = unique(pt_raw,'rows','stable'); % preserve the originally clicked order of points via 'stable'

% interpolate a smooth spline along the transect
numSplinePoints = 1000;
pt_spline = interparc(numSplinePoints,pt_clean(:,1),pt_clean(:,2),'spline');
xq = pt_spline(:,1);
yq = pt_spline(:,2);

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

Zq = interp2(XP,YP,Z(:,:,end),xq,yq); % surface elevation along cross-section 
Zq_max = interp2(XP,YP,maxSurf,xq,yq); % surface elevation along cross-section 
Zq_min = interp2(XP,YP,minSurf,xq,yq); % surface elevation along cross-section 

% interpolate deposit thicknesses
clear depThkXS 
for qq = 1:size(depositThkFull,3)
    depThkXS(qq,:) = interp2(XP,YP,depositThkFull(:,:,qq),xq,yq); % deposit thickness along cross-section 
end

%--------------------------------------------
% plot on XS location map
figure(3333)
hold on;

plot(xq+x0,yq+y0,'-k')
scatter(xq(1)+x0,yq(1)+y0,endDotSize,'w','filled','MarkerEdgeColor','k');
for ii = 2:length(dq)
    if floor(dq(ii)./dotSpacing)>floor(dq(ii-1)./dotSpacing)
        scatter(xq(ii)+x0,yq(ii)+y0,midDotSize,'w','filled','MarkerEdgeColor','k');
    end
end
scatter(xq(end)+x0,yq(end)+y0,endDotSize,'w','filled','MarkerEdgeColor','k');
if ~transectFromOrigin
    text(xq(1)+x0-0.08,yq(1)+y0+0.03,transectLetter,'FontWeight','bold','FontAngle','italic','Color','k','FontSize',endLetterSize);
end
text(xq(end)+x0-0.08,yq(end)+y0+0.03,[transectLetter ''''],'FontWeight','bold','FontAngle','italic','Color','k','FontSize',endLetterSize);
%------------------------------------------------------------

% stratigraphic cross-section
figure(30)
set(gcf, 'Color','w');
clf
cmap=parula(length(T_full)); %%% MODIFY COLORMAP HERE? OR YTICK./12???
% cmap(1:8,:)=[];
colormap(cmap)
hold on; box on; grid on;
area(dq,depThkXS',-35,'LineWidth',0.2,'FaceColor','flat'); % plot stratigraphy
% plot(dq,depThkXS(1,:)','-k', 'LineWidth',1.0); % min surface elev
plot(dq,Zq_min,'-k', 'LineWidth',1.0); % min surface elev
plot(dq,Zq_max,':k', 'LineWidth',1.0); % max surface elev
plot(dq,Zq,'-k', 'LineWidth',1.5); % most recent surface elev
XT=get(gca,'XTick');
set(gca,'XTick',[XT(1):0.1:XT(end)]);
set(gca,'YTick',[-30:2:20]);

title(['Cross-section ' transectLetter1 '-' transectLetter2 '''']);
xlabel('Distance [km]')
ylabel('Elevation [m NAVD88]')
ylim([min(min(Zq_min))-1 max([max(Zq)+1 3])]);

% % COLOUR BAR!!!
% cb=colorbar('southoutside');
% ylabel(cb,'Deposit Age [y]','FontWeight','bold','FontAngle','italic')
% set(cb,'YTick',0:12:size(depositThk,3)+1./12,...
%     'TickLabels',[size(depositThk,3)./12-round((0:12:size(depositThk,3)+1)./12)]);

set(gca,'FontWeight','bold','FontAngle','italic')
set(gca,'Layer','top','GridColor',[0.5 0.5 0.5],'GridAlpha',0.4);
MLW_line=plot([XT(1) XT(end)],[MLW MLW],'--','Color',deepBlue,'LineWidth',0.5); % plot MLW
MHW_line=plot([XT(1) XT(end)],[MHW MHW],'--','Color',deepBlue,'LineWidth',0.5); % plot MHW
uistack(MHW_line,'bottom'); uistack(MLW_line,'bottom');
text(0.02,MHW+0.3,'MHW','FontWeight','bold','FontAngle','italic','Color',deepBlue);
text(0.02,MLW+0.3,'MLW','FontWeight','bold','FontAngle','italic','Color',deepBlue);
xlim([0 dq(find(~isnan(Zq),1,'last'))]);

% Export figure to png file
pngName = ['Cross-section ' transectLetter1 '-' transectLetter2 ''''];
maxTransectLength=1.4; maxDepthRange=14; 
maxFigWidth=20; maxFigHeight=5;
dimensions = [diff(get(gca,'XLim'))/maxTransectLength*maxFigWidth...
    diff(get(gca,'YLim'))/maxDepthRange*maxFigHeight].*2; % [width height]
printFigs(dimensions,[plotDir filesep pngName],0);

clear xq yq dq


%==========================================================================
% print map
figure(3333)

xlim([307 309]);
ylim([3834.3 3836.3]);
set(gca,'XTick',[307:0.2:309]);
set(gca,'YTick',[3834.3:0.2:3836.1]);

% % plot origin
% transectFromOrigin=0; % is this a straight transect from the origin?
% if transectFromOrigin
%     scatter(x0,y0,20,'w','filled','MarkerEdgeColor','k');
%     text(x0-0.4,y0-0.4,'O','FontWeight','bold','FontName','Myriad Pro','Color','w');
% end

% % plot frame 2!
% scatter(xyFrame2(1),xyFrame2(2),20,'^y','filled','MarkerEdgeColor','k')

% Export figure to png file
pngName = ['Stratigraphic Transect Map (GUI-Defined) - Ameland Inlet'];
dimensions = [29.7 21]./1.2; % [width height]
printFigs(dimensions,[plotDir filesep pngName],0);


