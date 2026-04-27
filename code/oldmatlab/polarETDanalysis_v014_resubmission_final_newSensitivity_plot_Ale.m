%% polarETDanalysis_v013_resubmission.m
% track volumetric changes to ebb-tidal delta in polar coordinates from the
% centre of the inlet
% based on bathyChangeEnvelope.m
% v001 - SGP/2021-02-04 - Original in cartesian coordinates
% v002 - SGP/2021-02-08 - Updated to run in polar coordinates.
% v003 - SGP/2021-02-19 - Now with volume calcs in polar coords too
% v004 - SGP/2021-02-22 - Updated after talking to Edwin (incl 2020)
% v005 - SGP/2021-02-22 - With corrected volume calculation and weighted mean
% v006 - SGP/2021-06-18 - New version for paper with better weighting calc, also V_anomaly instead of V_active
% v007 - SGP/2021-08-08 - Updated with 1979, 1982, and 2021
% v008 - SGP/2021-09-21 - Updated with Edwin's new 1926 dataset
% v009 - SGP/2022-01-17 - Updated version for paper resubmission with peak tracking
% v010 - SGP/2022-01-19 - Updated cross-correlation and peak-tracking -> used in resubmission
% v013 - SGP/2022-02-07 - Directly descended from v010, skip v011/012.
%                         NEW SENSITIVITY: just go directly from PolarData.mat to be more
%                         comparable.
% v014 - SGP/2022-02-09 - Now with Alejandra's check for velocity at centre

clear all
close all
format compact
clc

versionCode = 'v013';

% load bathymetric data
% load('..\data\ame_filled-1926-2020.mat')
% load('..\data\ameland_1926-2021_merged.mat')

% set output path
outPathRoot = '..\results\newBathy\resubmission_v003\';

makeMiddlePlots = 0; % 1 = make plots of hAnomaly etc for all surveys
trajectoryPlots = 0; % 1 = make plots of vol anomaly trajectories

o0 = [169.5,607];
origins = [o0(1),o0(2);...
    o0(1),o0(2)+1;...
    o0(1)+1,o0(2);...
    o0(1),o0(2)-1;...
    o0(1)-1,o0(2)];
originName = {'O';'N1';'E1';'S1';'W1'};
createNew = 0;

for oo = 2:length(origins)
    outPath = [outPathRoot filesep originName{oo}];
    mkdir(outPath)
    
    if ~createNew
        
        load(['polarData_' versionCode '_' originName{oo} '.mat'])
        
        % set origin of polar grid
        x0=origins(oo,1);
        y0=origins(oo,2);
        Xp=X-x0;
        Yp=Y-y0;
        XLim = [162 176];
        YLim = [605 615];
        XLim = XLim-x0;
        YLim = YLim-y0;
        TLim = [1 length(T)];
    else
        
        
        % if exist('polarData.mat')
        %     load('polarData.mat')
        
        % set origin of polar grid
        x0=169.5;%x0=169.5;
        y0=607;
        %     Xp=X-x0;
        %     Yp=Y-y0;
        
        XLim = [162 176];
        YLim = [605 615];
        XLim = XLim-x0;
        YLim = YLim-y0;
        TLim = [1 length(T)];
        
        clear hAnomalyPolar maxSurf minSurf theta vAnomalyPolar vAnomalyRho ...
            vAnomalyRhoEast vAnomalyRhoWest vAnomalyTheta weightedMeanSurf ...
            xPolar yPolar zPolar
        
        % else
        %     %% create minimum surface
        %     X = grd.x./1000; % convert to km
        %     Y = grd.y./1000;
        %     Z = -grd.dp; % convert from depth to surface elevation
        %     % Z(:,:,20)=[];
        % %     Z(:,:,3)=[];
        %     T = grd.year;
        % %     T([3])=[];
        %     clear grd
        %
        %     % duplicate first entry of array
        %     Z = cat(3,Z(:,:,1),Z);
        %     T = cat(2,T(1),T);
        %
        %     % Z(:,:,[18 20]) = []; % get rid of half-year surveys since they are smaller
        %     % T([18 20]) = [];
        dx = 20;
        iT0=2;
        
        TLim = [1 length(T)];
        
        % perform weighted mean surface calc
        Tmidpt(1)=T(iT0);
        for tt = iT0:length(T)-1
            Tmidpt(tt) = (T(tt)+T(tt+1))/2;
        end
        Tmidpt(length(T))=(T(end)+T(end)+1)/2; % assume current year has weight of one year, given annual surveys
        Tweight = diff(Tmidpt);
        
        weightedMeanSurf=zeros(size(X));
        for tt = iT0:length(T)
            weightedMeanSurf = weightedMeanSurf + squeeze(Z(:,:,tt)).*Tweight(tt-1);
        end
        weightedMeanSurf =  weightedMeanSurf./(T(end)-T(iT0));
        
        % perform calculations
        minSurf = min(Z(:,:,TLim(1):end),[],3);
        unweightedMeanSurf = nanmean(Z(:,:,TLim(1):end),3);
        % medianSurf = nanmedian(Z(:,:,TLim(1):end),3);
        % stdSurf = nanstd(Z(:,:,TLim(1):end),3); %% NOTE!!!! THIS SHOULD BE WEIGHTED AS NOT ALL TIMESTEPS ARE EQUAL!
        maxSurf = max(Z(:,:,TLim(1):end),[],3);
        
        o0 = [169.5,607];
        origins = [o0(1),o0(2)+1;...
            o0(1)+1,o0(2);...
            o0(1),o0(2)-1;...
            o0(1)-1,o0(2)];
        
        originName = {'N1';'E1';'S1';'W1'};
        for oo = 1:length(origins)
            % set origin of polar grid
            x0=origins(oo,1);
            y0=origins(oo,2);
            Xp=X-x0;
            Yp=Y-y0;
            
            XLim = [162 176];
            YLim = [605 615];
            XLim = XLim-x0;
            YLim = YLim-y0;
            
            % limit calculations to area of interest
            tPlot=T(TLim(1):TLim(2));
            CLim = [-0.15 0.15];
            mask = false(size(X));
            mask(or(X<XLim(1),X>XLim(2))) = true; % focus within X limits
            mask(or(Y<YLim(1),Y>YLim(2))) = true; % focus within Y limits
            elevMask = mask;
            elevMask(Z(:,:,end)>0) = false; % remove all areas >0 m in elev
            mask = repmat(mask,1,1,TLim(2)-TLim(1)+1);
            
            % calculate grid cell surface area
            dx=unique(diff(X(:,1)))*1000; % here, grid is all uniform
            dy=unique(diff(Y(1,:)))*1000; % here, grid is all uniform
            areaGridCell = abs(dx(1)*dy(1));
            
            % calculate the height above the minimum surface (anomaly height)
            hAnomaly = nan(size(Z));
            hAnomaly(:,:,TLim(1):TLim(2)) = Z(:,:,TLim(1):TLim(2))-weightedMeanSurf;
            hActive = nan(size(Z));
            hActive(:,:,TLim(1):TLim(2)) = Z(:,:,TLim(1):TLim(2))-Z(:,:,TLim(1));
            vActive = hActive .* areaGridCell ./10^6; % in Mm3
            % hActive(~mask)=NaN;
            
            % calculate volume above minimum surface (Anomaly volume)
            vAnomaly = hAnomaly .* areaGridCell ./10^6; % in Mm3
            
            
            %% set up grid in polar coordinates
            
            % distribution of angles
            dTheta = -1; % grid angle spacing (in negative deg)
            dRho = 40; % grid spacing in radial direction (in m)
            theta = [190:dTheta:10].*(2*pi/360); % where 0deg =  East and direction is CCW
            rho = [1000:dRho:7000]; % distance from the origin [x0,y0]
            [xPolar,yPolar] = pol2cart(repmat(theta,length(rho),1),repmat(rho,length(theta),1)');
            
            % convert polar coordinates to km
            xPolar = xPolar./1000;
            yPolar = yPolar./1000;
            
            %% Calculate bathymetric properties of polar grid
            
            % XYZ coordinates of bathymetry data (query points)
            xq = reshape(X,[],1);
            yq = reshape(Y,[],1);
            
            % initialize arrays for speed (allegedly...)
            zPolar = nan(length(rho)-1,length(theta)-1,length(T));
            hAnomalyPolar = nan(length(rho)-1,length(theta)-1,length(T));
            vAnomalyPolar = nan(length(rho)-1,length(theta)-1,length(T));
            vAnomalyPolarAlt = nan(length(rho)-1,length(theta)-1,length(T));
            
            tic
            for tt = TLim(1):TLim(2)
                disp(['tt=' num2str(T(tt))]);
                % reshape list of Z/hAnomaly/vAnomaly for a given year
                zq = reshape(Z(:,:,tt),[],1);
                hq = reshape(hAnomaly(:,:,tt),[],1);
                vq = reshape(vAnomaly(:,:,tt),[],1);
                % Loop through and tabulate contents of each grid cell at each timestep
                for jj = 1:length(theta)-1
                    if mod(theta(jj)./(2*pi/360),10)==0
                        disp(['T=' num2str(T(tt)) ', theta=' num2str(theta(jj)./(2*pi/360))]);
                    end
                    for ii = 1:length(rho)-1
                        
                        % identify coordinates of a given polar grid cell (4 corners and repeat first point)
                        xv = [xPolar(ii,jj) xPolar(ii+1,jj) xPolar(ii+1,jj+1) xPolar(ii,jj+1) xPolar(ii,jj)]+x0;
                        yv = [yPolar(ii,jj) yPolar(ii+1,jj) yPolar(ii+1,jj+1) yPolar(ii,jj+1) yPolar(ii,jj)]+y0;
                        
                        % check which XY bathy coordinates lie inside a given grid cell
                        [in] = inpolygon(xq,yq,xv,yv);
                        
                        % multiply by grid cell area and don't forget unit stuff
                        zPolar(ii,jj,tt) = nanmean(zq(in));
                        hAnomalyPolar(ii,jj,tt) = nanmean(hq(in));
                        vAnomalyPolar(ii,jj,tt) = nansum(vq(in));
                        vAnomalyPolarAlt(ii,jj,tt) = nanmean(hq(in)).*areaGridCell.*sum(in)./10^6;
                        % note: the latter is a check to make sure the volume calc is
                        % correct, and they both match as of 2021-02-22
                    end
                end
            end
            toc
            
            % calculate dynamic volumes
            % collapse along Rho dimension
            vAnomalyRho = squeeze(nansum(vAnomalyPolar,2));
            vAnomalyRhoWest = squeeze(nansum(vAnomalyPolar(:,1:33,:),2));
            vAnomalyRhoNorthwest = squeeze(nansum(vAnomalyPolar(:,34:78,:),2));
            vAnomalyRhoNortheast = squeeze(nansum(vAnomalyPolar(:,79:124,:),2));
            vAnomalyRhoEast = squeeze(nansum(vAnomalyPolar(:,125:180,:),2));
            
            % collapse along Theta dimension
            vAnomalyTheta = squeeze(nansum(vAnomalyPolar,1));
            
            hActivePolar = nan(size(zPolar));
            hActivePolar(:,:,TLim(1):TLim(2)) = zPolar(:,:,TLim(1):TLim(2))-min(zPolar,[],3);
            vActivePolar = hActivePolar .* areaGridCell ./10^6; % in Mm3
            
            
            save(['polarData_' versionCode '_' originName{oo} '.mat'],...
                'xPolar','yPolar','theta','rho',...
                'zPolar','hAnomalyPolar','vAnomalyPolar','vAnomalyRho','vAnomalyTheta','vAnomalyRhoWest','vAnomalyRhoEast',...
                'X','Y','Z','T','weightedMeanSurf','minSurf','maxSurf');
        end
    end
    
    %% PLOT Z in POLAR COORDINATES
    
    fontSize = 14;
    
    % ---------------------------------------------------------------------
    % PLOT RESAMPLED DATA IN REAL XY SPACE
    
    for tt = 1:length(T)
        figure(44)
        clf
        set(gcf,'Color','w')
        subplot(1,2,1)
        hold on; grid on; box on;
        pc2=pcolor(xPolar(1:end-1,1:end-1),yPolar(1:end-1,1:end-1),zPolar(:,:,tt));
        pc2.FaceColor = 'flat';
        pc2.EdgeColor = 'none';
        axis equal
        xlabel('Easting [km]');
        ylabel('Northing [km]');
        xlim([XLim(1) XLim(2)]);
        ylim([YLim(1) YLim(2)]);
        [cmap] = kg2Colormap;
        colormap(cmap);
        cb=colorbar('location','southoutside');
        set(cb,'YLim',[-20 10],'FontName','Myriad Pro','FontSize',fontSize);
        caxis([-20 10])
        
        if T(tt) == 2021
            title(['(b) X-Y Space [' num2str(T(tt)) ']']);
        elseif T(tt) == 1975
            title(['(a) X-Y Space [' num2str(T(tt)) ']']);
        else
            title(['X-Y Space [' num2str(T(tt)) ']']);
        end
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        
        % plot radial grid on top
        dThetaPlotGrid = -22.5; % grid angle spacing (in negative deg)
        dRhoPlotGrid = 1000; % grid spacing in radial direction (in m)
        thetaPlotGrid = [360:dThetaPlotGrid:0].*(2*pi/360); % where 0deg =  East and direction is CCW
        thetaPlotGridFine = [360:-1:0].*(2*pi/360); % where 0deg =  East and direction is CCW
        rhoPlotGrid = [0:dRhoPlotGrid:7000]; % distance from the origin [x0,y0]
        [xPolarPlotGrid,yPolarPlotGrid] = pol2cart(repmat(thetaPlotGrid,length(rhoPlotGrid),1),...
            repmat(rhoPlotGrid,length(thetaPlotGrid),1)');
        xPolarPlotGrid = xPolarPlotGrid./1000;
        yPolarPlotGrid = yPolarPlotGrid./1000;
        [xPolarPlotGridFine,yPolarPlotGridFine] = pol2cart(repmat(thetaPlotGridFine,length(rhoPlotGrid),1),...
            repmat(rhoPlotGrid,length(thetaPlotGridFine),1)');
        xPolarPlotGridFine = xPolarPlotGridFine./1000;
        yPolarPlotGridFine = yPolarPlotGridFine./1000;
        
        for ii = 1:length(thetaPlotGrid)
            plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); % plot radial lines
            
        end
        for ii = 1:length(rhoPlotGrid)
            plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); % plot circles
        end
        
        
        % ---------------------------------------------------------------------
        % PLOT IN THETA-RHO SPACE
        thetaDegNorth = theta./(2*pi/360); % convert back from radians to degrees
        thetaDegNorth = 90-thetaDegNorth; % correct to 0=north, CW pos (nautical)
        % thetaDegNorth(thetaDegNorth<0)=thetaDegNorth(thetaDegNorth<0)+360; %correct negative numbers
        [thetaGrid,rhoGrid] = meshgrid(thetaDegNorth(1:end-1),rho(1:end-1)./1000);
        
        subplot(1,2,2)
        hold on; grid on; box on;
        pc2=pcolor(thetaGrid,rhoGrid,zPolar(:,:,tt));
        pc2.FaceColor = 'flat';
        pc2.EdgeColor = 'none';
        xlabel('\theta [{\circ}]');
        ylabel('\rho [km]');
        set(gca,'XTick',[-180:22.5:180])
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4,'FontName','Myriad Pro','FontSize',fontSize);
        xlim([min(thetaDegNorth) 79]);
        ylim([min(rho./1000) max(rho./1000)]);
        [cmap] = kg2Colormap;
        colormap(cmap);
        caxis([-20 10])
        
        if T(tt) == 2021
            title(['(e) \rho-\theta Space [' num2str(T(tt)) ']'])
        elseif T(tt) == 1975
            title(['(d) \rho-\theta Space [' num2str(T(tt)) ']'])
        else
            title(['\rho-\theta Space [' num2str(T(tt)) ']'])
        end
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        
        % Export figure to png file
        pngName = ['polarMap_Overview_' num2str(T(tt))];
        dimensions = [20 10].*1.5; % [width height]
        printFigs(dimensions,[outPath filesep pngName],0);
        
    end
    
    
    if makeMiddlePlots
        % ---------------------------------------------------------------------
        % PLOT RESAMPLED hAnomaly DATA IN REAL XY SPACE
        
        for tt = TLim(1):TLim(2)%1:TLim(2)-TLim(1)+1
            figure(55)
            clf
            set(gcf,'Color','w')
            subplot(1,2,1)
            hold on; grid on; box on;
            %     pc2=pcolor(xPolar(1:end-1,1:end-1),yPolar(1:end-1,1:end-1),hAnomaly(:,:,tt-TLim(1)+1));
            pc2=pcolor(xPolar(1:end-1,1:end-1),yPolar(1:end-1,1:end-1),hAnomalyPolar(:,:,tt));
            pc2.FaceColor = 'flat';
            pc2.EdgeColor = 'none';
            axis equal
            xlabel('Easting [km]');
            ylabel('Northing [km]');
            xlim([XLim(1) XLim(2)]);
            ylim([YLim(1) YLim(2)]);
            [cmap] = parula(100);
            colormap(cmap);
            caxis([-5 5])
            title(['X-Y Space [' num2str(T(tt)) ']']);
            cb=colorbar('southoutside');
            ylabel(cb,'Deviation from Z_{mean}','FontName','Myriad Pro','FontSize',fontSize);
            set(gca,'FontName','Myriad Pro','FontSize',fontSize)
            
            % plot radial grid on top
            dThetaPlotGrid = -22.5; % grid angle spacing (in negative deg)
            dRhoPlotGrid = 1000; % grid spacing in radial direction (in m)
            thetaPlotGrid = [360:dThetaPlotGrid:0].*(2*pi/360); % where 0deg =  East and direction is CCW
            thetaPlotGridFine = [360:-1:0].*(2*pi/360); % where 0deg =  East and direction is CCW
            rhoPlotGrid = [0:dRhoPlotGrid:7000]; % distance from the origin [x0,y0]
            [xPolarPlotGrid,yPolarPlotGrid] = pol2cart(repmat(thetaPlotGrid,length(rhoPlotGrid),1),...
                repmat(rhoPlotGrid,length(thetaPlotGrid),1)');
            xPolarPlotGrid = xPolarPlotGrid./1000;
            yPolarPlotGrid = yPolarPlotGrid./1000;
            [xPolarPlotGridFine,yPolarPlotGridFine] = pol2cart(repmat(thetaPlotGridFine,length(rhoPlotGrid),1),...
                repmat(rhoPlotGrid,length(thetaPlotGridFine),1)');
            xPolarPlotGridFine = xPolarPlotGridFine./1000;
            yPolarPlotGridFine = yPolarPlotGridFine./1000;
            
            for ii = 1:length(thetaPlotGrid)
                plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); % plot radial lines
                
            end
            for ii = 1:length(rhoPlotGrid)
                plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); % plot circles
            end
            
            
            % ---------------------------------------------------------------------
            % PLOT IN THETA-RHO SPACE
            thetaDegNorth = theta./(2*pi/360); % convert back from radians to degrees
            thetaDegNorth = 90-thetaDegNorth; % correct to 0=north, CW pos (nautical)
            % thetaDegNorth(thetaDegNorth<0)=thetaDegNorth(thetaDegNorth<0)+360; %correct negative numbers
            [thetaGrid,rhoGrid] = meshgrid(thetaDegNorth(1:end-1),rho(1:end-1)./1000);
            
            subplot(1,2,2)
            hold on; grid on; box on;
            %     pc2=pcolor(thetaGrid,rhoGrid,zPolar(:,:,tt-TLim(1)+1));
            pc2=pcolor(thetaGrid,rhoGrid,hAnomalyPolar(:,:,tt));
            pc2.FaceColor = 'flat';
            pc2.EdgeColor = 'none';
            xlabel('\theta [{\circ}]');
            ylabel('\rho [km]');
            set(gca,'XTick',[-180:22.5:180])
            set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
            xlim([min(thetaDegNorth) 79]);
            ylim([min(rho./1000) max(rho./1000)]);
            [cmap] = parula(100);
            colormap(cmap);
            caxis([-5 5])
            title(['\rho-\theta Space [' num2str(T(tt)) ']'])
            set(gca,'FontName','Myriad Pro','FontSize',fontSize)
            
            % Export figure to png file
            pngName = ['polarMap_hAnomaly_' num2str(T(tt))];
            dimensions = [20 10].*1.5; % [width height]
            printFigs(dimensions,[outPath filesep pngName],0);
            
        end
        
        % PLOT vAnomaly in POLAR COORDINATES
        
        % ---------------------------------------------------------------------
        % PLOT RESAMPLED vAnomaly DATA IN REAL XY SPACE
        
        for tt = TLim(1):TLim(2)%1:TLim(2)-TLim(1)+1
            figure(66)
            clf
            set(gcf,'Color','w')
            subplot(1,2,1)
            hold on; grid on; box on;
            %     pc2=pcolor(xPolar(1:end-1,1:end-1),yPolar(1:end-1,1:end-1),hAnomaly(:,:,tt-TLim(1)+1));
            pc2=pcolor(xPolar(1:end-1,1:end-1),yPolar(1:end-1,1:end-1),vAnomalyPolar(:,:,tt).*1000);
            pc2.FaceColor = 'flat';
            pc2.EdgeColor = 'none';
            axis equal
            xlabel('Easting [km]');
            ylabel('Northing [km]');
            xlim([XLim(1) XLim(2)]);
            ylim([YLim(1) YLim(2)]);
            [cmap] = parula(100);
            colormap(cmap);
            caxis([-10 10]);
            % 	caxis([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
            if T(tt) == 2021
                title(['(c) X-Y Space [' num2str(T(tt)) ']']);
            else
                title(['X-Y Space [' num2str(T(tt)) ']']);
            end
            cb=colorbar('southoutside');
            ylabel(cb,'V_{Anomaly} [10^3 m^3]','FontName','Myriad Pro','FontSize',fontSize);
            set(gca,'FontName','Myriad Pro','FontSize',fontSize)
            
            % plot radial grid on top
            dThetaPlotGrid = -22.5; % grid angle spacing (in negative deg)
            dRhoPlotGrid = 1000; % grid spacing in radial direction (in m)
            thetaPlotGrid = [360:dThetaPlotGrid:0].*(2*pi/360); % where 0deg =  East and direction is CCW
            thetaPlotGridFine = [360:-1:0].*(2*pi/360); % where 0deg =  East and direction is CCW
            rhoPlotGrid = [0:dRhoPlotGrid:7000]; % distance from the origin [x0,y0]
            [xPolarPlotGrid,yPolarPlotGrid] = pol2cart(repmat(thetaPlotGrid,length(rhoPlotGrid),1),...
                repmat(rhoPlotGrid,length(thetaPlotGrid),1)');
            xPolarPlotGrid = xPolarPlotGrid./1000;
            yPolarPlotGrid = yPolarPlotGrid./1000;
            [xPolarPlotGridFine,yPolarPlotGridFine] = pol2cart(repmat(thetaPlotGridFine,length(rhoPlotGrid),1),...
                repmat(rhoPlotGrid,length(thetaPlotGridFine),1)');
            xPolarPlotGridFine = xPolarPlotGridFine./1000;
            yPolarPlotGridFine = yPolarPlotGridFine./1000;
            
            for ii = 1:length(thetaPlotGrid)
                plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); % plot radial lines
                
            end
            for ii = 1:length(rhoPlotGrid)
                plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); % plot circles
            end
            
            
            % ---------------------------------------------------------------------
            % PLOT IN THETA-RHO SPACE
            thetaDegNorth = theta./(2*pi/360); % convert back from radians to degrees
            thetaDegNorth = 90-thetaDegNorth; % correct to 0=north, CW pos (nautical)
            % thetaDegNorth(thetaDegNorth<0)=thetaDegNorth(thetaDegNorth<0)+360; %correct negative numbers
            [thetaGrid,rhoGrid] = meshgrid(thetaDegNorth(1:end-1),rho(1:end-1)./1000);
            
            subplot(1,2,2)
            hold on; grid on; box on;
            %     pc2=pcolor(thetaGrid,rhoGrid,zPolar(:,:,tt-TLim(1)+1));
            pc2=pcolor(thetaGrid,rhoGrid,vAnomalyPolar(:,:,tt).*1000);
            pc2.FaceColor = 'flat';
            pc2.EdgeColor = 'none';
            xlabel('\theta [{\circ}]');
            ylabel('\rho [km]');
            set(gca,'XTick',[-180:22.5:180])
            set(gca,'Layer','top','GridColor','k','GridAlpha',1);
            xlim([min(thetaDegNorth) 79]);
            ylim([min(rho./1000) max(rho./1000)]);
            [cmap] = parula(100);
            colormap(cmap);
            caxis([-10 10]);
            % 	caxis([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
            if T(tt) == 2021
                title(['(f) \rho-\theta Space [' num2str(T(tt)) ']'])
            else
                title(['\rho-\theta Space [' num2str(T(tt)) ']'])
            end
            %     ylabel(cb,'V_{anomaly} [Mm^3]')
            set(gca,'FontName','Myriad Pro','FontSize',fontSize)
            
            % Export figure to png file
            pngName = ['polarMap_vAnomaly_' num2str(T(tt))];
            dimensions = [20 10].*1.5; % [width height]
            printFigs(dimensions,[outPath filesep pngName],0);
            
        end
        
    end
    
    %% hovmuller diagram (Theta)
    figure(1)
    clf
    set(gcf,'Color','w');
    hold on; box on; grid on;
    pp=pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T)),repmat(T,length(thetaDegNorth(1:end-1)),1),vAnomalyTheta);
    pp.EdgeColor='none';
    % cb=colorbar('Location','southoutside');
    title('(g) Changes in Volume Anomaly [\theta]')
    xlabel('\theta [\circ]');
    ylabel('Time [y]');
    ylim([1975 2021]);
    xlim([-100 79]);
    caxis([-0.51 0.51])
    set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    % ylabel(cb,'V_{anomaly} [Mm^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    % cb.Ticks = [-0.5:0.1:0.5];
    % clim([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
    % xlim([XLim(1) XLim(2)]);
    % clim([CLim(1) CLim(2)]);
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta'];
    dimensions = [20 18]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% hovmuller diagram (Rho)
    figure(2)
    clf
    set(gcf,'Color','w');
    hold on; box on; grid on;
    
    pp=pcolor(repmat(T,length(rho(1:end-1)),1),repmat(rho(1:end-1)./1000,length(T),1)',vAnomalyRho);
    pp.EdgeColor='none';
    % cb=colorbar('Location','southoutside');
    % ylabel(cb,'V_{anomaly} [Mm^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    title('(h) Changes in Volume Anomaly [\rho]')
    ylabel('\rho [km]');
    xlabel('Time [y]');
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    xlim([1975 2021]);
    % ylim([YLim(1) YLim(2)]);
    caxis([-0.31 0.31])
    % clim([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
    
    % Export figure to png file
    pngName = ['HovmullerPolarRho'];
    dimensions = [20 18]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% Joy Division Plot
    
    figure(666);
    clf
    set(gcf,'Color','k');
    
    scaleFactor = 8;
    
    for ii = size(vAnomalyTheta,2):-1:1
        hold on
        patch([thetaDegNorth(1:end-1) thetaDegNorth(end-1) thetaDegNorth(1)],...
            [(vAnomalyTheta(:,ii)-min(min(vAnomalyTheta))).*scaleFactor+T(ii); -100; -100],'k');
        plot(thetaDegNorth(1:end-1),(vAnomalyTheta(:,ii)-min(min(vAnomalyTheta))).*scaleFactor+T(ii),'w')
        set(gca,'Color','k');
    end
    ylim([1975 2045]);
    
    mainTitle='UNKNOWN SEDIMENTS';
    dim = [0.05,0.96,0.9,0.05];
    a = annotation('textbox',dim,'String',mainTitle,'FitBoxToText','on',...
        'Linestyle','none','FontSize',30,'HorizontalAlignment','Center',...
        'Color','k');
    axis off
    
    % Export figure to png file
    pngName = ['JoyDivision_Hovmoeller'];
    dimensions = [30 25]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    
    %% Peakfinding analysis
    
    figure(667)
    set(gcf,'Color','w');
    clf
    
    
    scaleFactor = 2;
    minPeakDist = 8;
    prom = 0.1;
    numPeaks = 6;
    
    for ii=1:22
        hold on
        % plot(vAnomalyTheta(:,ii))
        findpeaks(vAnomalyTheta(:,ii).*scaleFactor+T(ii),thetaDegNorth(1:end-1),...
            'MinPeakDistance',minPeakDist,'MinPeakProminence',prom.*scaleFactor,...
            'Annotate','extents');%,...
        %         'SortStr','descend','NPeaks',numPeaks)
        [peaks(ii).pks,peaks(ii).locs,peaks(ii).width,peaks(ii).prom] = findpeaks(vAnomalyTheta(:,ii),thetaDegNorth(1:end-1),...
            'MinPeakDistance',minPeakDist,'MinPeakProminence',prom);%,...
        %         'SortStr','descend','NPeaks',numPeaks);
        text(peaks(ii).locs+.02,peaks(ii).pks+T(ii),num2str((1:numel(peaks(ii).pks))'))
    end
    title('Changes in Volume Anomaly [\theta]')
    xlabel('\theta [\circ]');
    ylabel('Time [y]');
    % ylim([1975 2021]);
    xlim([-100 79]);
    set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_Peakfinder'];
    dimensions = [30 25]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% Troughfinding analysis
    
    figure(6672)
    set(gcf,'Color','w');
    clf
    
    scaleFactor = 2;
    minPeakDist = 8;
    prom = 0.1;
    numPeaks = 6;
    
    for ii=1:22
        hold on
        % plot(vAnomalyTheta(:,ii))
        findpeaks(-vAnomalyTheta(:,ii).*scaleFactor+T(ii),thetaDegNorth(1:end-1),...
            'MinPeakDistance',minPeakDist,'MinPeakProminence',prom.*scaleFactor,...
            'Annotate','extents');%,...
        %         'SortStr','descend','NPeaks',numPeaks)
        [troughs(ii).pks,troughs(ii).locs,troughs(ii).width,troughs(ii).prom] = findpeaks(-vAnomalyTheta(:,ii),thetaDegNorth(1:end-1),...
            'MinPeakDistance',minPeakDist,'MinPeakProminence',prom);%,...
        %         'SortStr','descend','NPeaks',numPeaks);
        text(troughs(ii).locs+.02,troughs(ii).pks+T(ii),num2str((1:numel(troughs(ii).pks))'))
    end
    title('Changes in Negative Volume Anomaly [\theta]')
    xlabel('\theta [\circ]');
    ylabel('Time [y]');
    % ylim([1975 2021]);
    xlim([-100 79]);
    set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_Troughfinder'];
    dimensions = [30 25]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    
    %% peakplotter
    
    if trajectoryPlots
        trendlabels = 0; % for labelling trendlines
        %     if oo == 12
        
        
        % print trajectory statistics to file
        fid = fopen([outPath filesep 'trajectory_stats.txt'],'w');
        fprintf(fid,[strrep(outPath,'\','/'),'\n']);
        
        figure(668)
        clf
        set(gcf,'Color','w');
        hold on; box on; grid on;
        
        pp=pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T)),repmat(T,length(thetaDegNorth(1:end-1)),1),vAnomalyTheta);
        pp.EdgeColor='none';
        % cb=colorbar('Location','southoutside');
        title('Changes in Volume Anomaly [\theta]')
        xlabel('\theta [\circ]');
        ylabel('Time [y]');
        ylim([1975 2021]);
        xlim([-100 79]);
        caxis([-0.51 0.51])
        set(gca,'XTick',[-180:22.5:180])
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
        % ylabel(cb,'V_{anomaly} [Mm^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (i - boschplaat erosion)
        trajA = [peaks(2).locs(1);...
            peaks(3).locs(1);...
            peaks(4).locs(1);...
            peaks(5).locs(1);...
            peaks(6).locs(1);...
            peaks(7).locs(1);...
            peaks(8).locs(1);...
            peaks(9).locs(1);...
            peaks(10).locs(1);...
            peaks(11).locs(1);...
            peaks(12).locs(1);...
            peaks(13).locs(1);...
            peaks(14).locs(1);...
            peaks(15).locs(1);...
            peaks(16).locs(1);...
            peaks(17).locs(1);...
            peaks(18).locs(1);...
            peaks(19).locs(1);...
            peaks(20).locs(1);...
            peaks(21).locs(1)];
        TsubA = T(2:21);
        plot(trajA,TsubA,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubA'),1) TsubA']\trajA;
        trajA_pred = B(2).*TsubA+B(1);
        R = corrcoef(trajA,trajA_pred);
        Rsq = R(2,1).^2;
        plot(trajA_pred,TsubA,'--r','Linewidth',2.5)
        disp(['Slope i = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope i = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajA_pred(end-5)+3,TsubA(end-5),'i','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (ii - westgat infilling)
        trajB = [peaks(5).locs(2);...
            peaks(6).locs(2);...
            peaks(7).locs(2);...
            peaks(8).locs(2);...
            peaks(9).locs(2);...
            peaks(10).locs(2);...
            peaks(11).locs(2);...
            peaks(12).locs(2);...
            peaks(13).locs(2);...
            peaks(14).locs(2);...
            peaks(15).locs(2);...
            peaks(16).locs(2);...
            peaks(17).locs(2);...
            peaks(18).locs(2);...
            peaks(19).locs(2);...
            peaks(20).locs(2);...
            peaks(21).locs(2);...
            peaks(22).locs(1)];
        TsubB = T(5:end);
        plot(trajB,TsubB,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubB'),1) TsubB']\trajB;
        trajB_pred = B(2).*TsubB+B(1);
        R = corrcoef(trajB,trajB_pred);
        Rsq = R(2,1).^2;
        plot(trajB_pred,TsubB,'--r','Linewidth',2.5)
        disp(['Slope ii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope ii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajB_pred(end-5)+3,TsubB(end-5),'ii','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (iii - boschplaat to ebb lobe 1)
        trajC1 = [peaks(3).locs(2);...
            peaks(4).locs(2);...
            peaks(5).locs(3);...
            peaks(6).locs(3);...
            peaks(7).locs(3);...
            peaks(8).locs(3);...
            peaks(9).locs(3);...
            peaks(10).locs(3);...
            peaks(11).locs(3);...
            peaks(12).locs(3);...
            peaks(13).locs(4);...
            peaks(14).locs(4);...
            peaks(15).locs(4);...
            peaks(16).locs(4);...
            peaks(17).locs(4);...
            peaks(18).locs(4);...
            peaks(19).locs(4);...
            peaks(20).locs(4);...
            peaks(21).locs(4);...
            peaks(22).locs(3)];
        TsubC1 = T(3:end);
        plot(trajC1,TsubC1,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubC1'),1) TsubC1']\trajC1;
        trajC1_pred = B(2).*TsubC1+B(1);
        R = corrcoef(trajC1,trajC1_pred);
        Rsq = R(2,1).^2;
        plot(trajC1_pred,TsubC1,'--r','Linewidth',2.5)
        disp(['Slope iii-1 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope iii-1 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajC1_pred(end-7)+9,TsubC1(end-7)+2,'iii-1','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (iii - boschplaat to ebb lobe 2)
        trajC2 = [peaks(13).locs(3);...
            peaks(14).locs(3);...
            peaks(15).locs(3);...
            peaks(16).locs(3);...
            peaks(17).locs(3);...
            peaks(18).locs(3);...
            peaks(19).locs(3);...
            peaks(20).locs(3);...
            peaks(21).locs(3);...
            peaks(22).locs(2)];
        TsubC2 = T(13:end);
        plot(trajC2,TsubC2,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubC2'),1) TsubC2']\trajC2;
        trajC2_pred = B(2).*TsubC2+B(1);
        R = corrcoef(trajC2,trajC2_pred);
        Rsq = R(2,1).^2;
        plot(trajC2_pred,TsubC2,'--r','Linewidth',2.5)
        disp(['Slope iii-2 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope iii-2 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajC2_pred(end-5)+3,TsubC2(end-5),'iii-2','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following troughs (iv - akkepollegat migration)
        trajT1 = [troughs(2).locs(2);...
            troughs(3).locs(2);...
            troughs(4).locs(2);...
            troughs(5).locs(2);...
            troughs(6).locs(3);...
            troughs(7).locs(3);...
            troughs(8).locs(3);...
            troughs(9).locs(3);...
            troughs(10).locs(3);...
            troughs(11).locs(3);...
            troughs(12).locs(3);...
            troughs(13).locs(4);...
            troughs(14).locs(4);...
            troughs(15).locs(4);...
            troughs(16).locs(4);...
            troughs(17).locs(5);...
            troughs(18).locs(5);...
            troughs(19).locs(4);...
            troughs(20).locs(4);...
            troughs(21).locs(4);...
            troughs(22).locs(5)];
        TsubT1 = T(2:22);
        plot(trajT1,TsubT1,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubT1'),1) TsubT1']\trajT1;
        trajT1_pred = B(2).*TsubT1+B(1);
        R = corrcoef(trajT1,trajT1_pred);
        Rsq = R(2,1).^2;
        plot(trajT1_pred,TsubT1,'--r','Linewidth',2.5)
        disp(['Slope iv = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope iv = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajT1_pred(end-5)+3,TsubT1(end-5),'iv','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (v - bornrif (central))
        trajD = [peaks(2).locs(4);...
            peaks(3).locs(4);...
            peaks(4).locs(4);...
            peaks(5).locs(4);...
            peaks(6).locs(4);...
            peaks(7).locs(5);...
            peaks(8).locs(4);...
            peaks(9).locs(4);...
            peaks(10).locs(4);...
            peaks(11).locs(4);...
            peaks(12).locs(4);...
            peaks(13).locs(5);...
            peaks(14).locs(5);...
            peaks(15).locs(5);...
            peaks(16).locs(5);...
            peaks(17).locs(5);...
            peaks(18).locs(5);...
            peaks(19).locs(5);...
            peaks(20).locs(5);...
            peaks(21).locs(5);...
            peaks(22).locs(4)];
        TsubD = T(2:end);
        plot(trajD,TsubD,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubD'),1) TsubD']\trajD;
        trajD_pred = B(2).*TsubD+B(1);
        R = corrcoef(trajD,trajD_pred);
        Rsq = R(2,1).^2;
        plot(trajD_pred,TsubD,'--r','Linewidth',2.5)
        disp(['Slope v = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope v = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajD_pred(end-3)+3,TsubD(end-3),'v','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following troughs (vi - oostgat migration)
        trajT2 = [troughs(5).locs(3);...
            troughs(6).locs(4);...
            troughs(7).locs(5);...
            troughs(8).locs(4);...
            troughs(9).locs(4);...
            troughs(10).locs(4);...
            troughs(11).locs(4);...
            troughs(12).locs(4);...
            troughs(13).locs(5);...
            troughs(14).locs(5);...
            troughs(15).locs(5);...
            troughs(16).locs(5);...
            troughs(17).locs(6);...
            troughs(18).locs(6);...
            troughs(19).locs(5);...
            troughs(20).locs(5);...
            troughs(21).locs(5);...
            troughs(22).locs(6)];
        TsubT2 = T(5:22);
        plot(trajT2,TsubT2,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubT2'),1) TsubT2']\trajT2;
        trajT2_pred = B(2).*TsubT2+B(1);
        R = corrcoef(trajT2,trajT2_pred);
        Rsq = R(2,1).^2;
        plot(trajT2_pred,TsubT2,'--r','Linewidth',2.5)
        disp(['Slope vi-1 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope vi-1 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajT2_pred(end-5)+3,TsubT2(end-5),'vi-1','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following troughs (vi-b - oostgat migration)
        trajT3 = [troughs(5).locs(3);...
            troughs(6).locs(4);...
            troughs(7).locs(5);...
            troughs(8).locs(5);...
            troughs(9).locs(5);...
            troughs(10).locs(5);...
            troughs(11).locs(5);...
            troughs(12).locs(5);...
            troughs(13).locs(6);...
            troughs(14).locs(6);...
            troughs(15).locs(6);...
            troughs(16).locs(6);...
            troughs(17).locs(7);...
            troughs(18).locs(7);...
            troughs(19).locs(6);...
            troughs(20).locs(6);...
            troughs(21).locs(6);...
            troughs(22).locs(7)];
        TsubT3 = T(5:22);
        plot(trajT3,TsubT3,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubT3'),1) TsubT3']\trajT3;
        trajT3_pred = B(2).*TsubT3+B(1);
        R = corrcoef(trajT2,trajT3_pred);
        Rsq = R(2,1).^2;
        plot(trajT3_pred,TsubT3,'--r','Linewidth',2.5)
        disp(['Slope vi-2 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope vi-2 = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajT3_pred(end-5)+3,TsubT3(end-5),'vi-2','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (vii - bornrif strandhaak)
        trajF = [peaks(4).locs(5);...
            peaks(5).locs(5);...
            peaks(6).locs(5);...
            peaks(7).locs(7);...
            peaks(8).locs(6);...
            peaks(9).locs(7);...
            peaks(10).locs(6);...
            peaks(11).locs(7);...
            peaks(12).locs(6);...
            peaks(13).locs(7);...
            peaks(14).locs(7);...
            peaks(15).locs(7);...
            peaks(16).locs(7);...
            peaks(17).locs(7)];
        TsubF = T(4:17);
        plot(trajF,TsubF,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubF'),1) TsubF']\trajF;
        trajF_pred = B(2).*TsubF+B(1);
        R = corrcoef(trajF,trajF_pred);
        Rsq = R(2,1).^2;
        plot(trajF_pred,TsubF,'--r','Linewidth',2.5)
        disp(['Slope vii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope vii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajF_pred(end-5)+3,TsubF(end-5),'vii','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        % -----------------------------------------------------------------
        % define trajectory following peaks (viii - bornrif bankje)
        trajE = [peaks(13).locs(6);...
            peaks(14).locs(6);...
            peaks(15).locs(6);...
            peaks(16).locs(6);...
            peaks(17).locs(6);...
            peaks(18).locs(6);...
            peaks(19).locs(6);...
            peaks(20).locs(6);...
            peaks(21).locs(6);...
            peaks(22).locs(5)];
        TsubE = T(13:22);
        plot(trajE,TsubE,'.-k')
        
        % linear regression: [B1 B2]=X\Y where B2=slope, B1=intercept
        B=[ones(length(TsubE'),1) TsubE']\trajE;
        trajE_pred = B(2).*TsubE+B(1);
        R = corrcoef(trajE,trajE_pred);
        Rsq = R(2,1).^2;
        plot(trajE_pred,TsubE,'--r','Linewidth',2.5)
        disp(['Slope viii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f')])
        fprintf(fid,['Slope viii = ' num2str(B(2)*10,'%.1f') ' deg/decade, Rsq = '...
            num2str(Rsq,'%.2f') '\n']);
        if trendlabels
            text(trajE_pred(end-7)+3,TsubE(end-7),'viii','Color','w','FontWeight','bold','FontName','Myriad Pro','FontSize',fontSize);
        end
        
        
        % -------------------------------------------------------------------
        
        
        fclose('all');
        
        % Export figure to png file
        pngName = ['HovmullerPolarTheta_Trajectories'];
        dimensions = [20 18]./1.5; % [width height]
        printFigs(dimensions,[outPath filesep pngName],0);
        
        title('(g) Changes in Volume Anomaly [\theta]')
        % Export figure to png file
        pngName = ['HovmullerPolarTheta_Trajectories_fig4g'];
        dimensions = [20 18]./1.5; % [width height]
        printFigs(dimensions,[outPath filesep pngName],0);
    end
    
    %% Cross correlation
    % xcorr
    
    timegap = 1;
    iCutoff = 31; % bin at which to cut off calculations (31 = -70deg)
    for tt=3:length(T)
        
        figure(777)
        set(gcf,'Color','w');
        clf
        
        subplot(2,1,1)
        hold on; box on; grid on;
        plot(thetaDegNorth(iCutoff:end-1),vAnomalyTheta(iCutoff:end,tt-timegap),'color',[0.5 0.5 0.5])
        plot(thetaDegNorth(iCutoff:end-1),vAnomalyTheta(iCutoff:end,tt),'color',[0 0 0])
        title(['(a) Volume Anomaly (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ') [\theta]'])
        xlabel('\theta [\circ]');
        ylabel('Time [y]');
        % ylim([1975 2021]);
        xlim([-100 79]);
        set(gca,'XTick',[-180:22.5:180])
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        legend(num2str(T(tt-timegap)), num2str(T(tt)),'location','southeast');
        
        subplot(2,1,2)
        hold on; box on; grid on;
        
        [xc(tt,:),lags(tt,:)] = xcorr(vAnomalyTheta(iCutoff:end,tt),...
            vAnomalyTheta(iCutoff:end,tt-timegap),25,'normalized');
        stem(lags(tt,:),xc(tt,:));
        peakXC(tt) = max(xc(tt,:));
        peakLag(tt) = lags(tt,find(xc(tt,:)==max(xc(tt,:)),1));
        migRate(tt) = peakLag(tt)./(T(tt)-T(tt-timegap));
        title(['(b) Cross Correlation (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ')'])
        xlabel('\theta [\circ]');
        ylabel('R_{norm} [-]');
        ylim([-1 1]);
        xlim([-25 25]);
        set(gca,'XTick',[-180:5:180])
        text(-19,0.8,['Peak Cross Correlation: ' num2str(peakLag(tt)) '^o'],...
            'FontName','Myriad Pro','FontSize',10);
        text(-19,0.6,['Migration Rate: ' num2str(migRate(tt),'%.1f') '^o/year'],...
            'FontName','Myriad Pro','FontSize',10);
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        
        % Export figure to png file
        pngName = ['HovmullerPolarTheta_CrossCorrelation_' ...
            num2str(T(tt-timegap)) '-' num2str(T(tt))];
        dimensions = [30 25]./1.5; % [width height]
        printFigs(dimensions,[outPath filesep pngName],0);
    end
    
    
    %%
    figure(778)
    set(gcf,'Color','w');
    clf
    
    scaleFactor = 1;
    minPeakDist = 8;
    prom = 0.1;
    numPeaks = 1;
    
    for ii=3:22
        hold on
        % plot(vAnomalyTheta(:,ii))
        findpeaks(xc(ii,:),lags(ii,:),'Annotate','extents',...
            'SortStr','descend','NPeaks',numPeaks)
        [xc_peaks(ii).pks,xc_peaks(ii).locs,xc_peaks(ii).width,xc_peaks(ii).prom] ...
            = findpeaks(xc(ii,:),lags(ii,:),...
            'SortStr','descend','NPeaks',numPeaks);
        text(xc_peaks(ii).locs+.02,xc_peaks(ii).pks+T(ii),num2str((1:numel(xc_peaks(ii).pks))'))
    end
    title('Cross Correlation Through Time [\theta]')
    xlabel('lag [\circ]');
    ylabel('R_{norm} [-]');
    % ylim([1975 2021]);
    % xlim([-100 79]);
    % set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_XCorr_Peakfinder'];
    dimensions = [30 25]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% hovmuller diagram of cross-correlation
    
    weighted_mean_lag = nanmean(peakLag(3:end)./diff(T(2:end)));
    mean(migRate(3:end));
    
    wt_mean_lag = nanmean(peakLag(3:end)./diff(T(2:end)));
    wt_std_lag = nanstd(peakLag(3:end)./diff(T(2:end)));
    wt_mean_peak = nanmean(peakXC(3:end)./diff(T(2:end)));
    
    fid = fopen([outPath filesep 'xcorr_stats.txt'],'w');
    fprintf(fid,[strrep(outPath,'\','/'),'\n']);
    fprintf(fid,['mean mig rate = ' num2str(nanmean(migRate(3:end)).*10,'%.2f') '\n']);
    fprintf(fid,['mean width = ' num2str(nanmean([xc_peaks(3:end).width]),'%.2f') '\n']);
    fprintf(fid,['mean peak = ' num2str(nanmean([xc_peaks(3:end).pks]),'%.2f') '\n']);
    fprintf(fid,['mean prom = ' num2str(nanmean([xc_peaks(3:end).prom]),'%.2f') '\n']);
    fclose('all');
    
    
    nanmean([xc_peaks(3:end).width])
    
    figure(1001)
    clf
    set(gcf,'Color','w');
    
    subplot(2,1,1)
    hold on; box on; grid on;
    pp=pcolor(repmat(lags(end,:)',1,length(T)),repmat(T,length(lags),1),xc');
    pp.EdgeColor='none';
    plot([0 0],[1975 2021],'-k','Linewidth',2.5)
    plot([xc_peaks(3:end).locs],T(3:end),'.-k','Linewidth',1.0)
    cb=colorbar('Location','eastoutside');
    title('Cross-correlation Through Time')
    xlabel('Lag [\circ]');
    ylabel('Time [y]');
    ylim([1975 2021]);
    xlim([-25 25]);
    caxis([-1 1])
    set(gca,'XTick',[-25:5:25])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    ylabel(cb,'R_{norm}','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
    subplot(2,1,2)
    hold on; box on; grid on;
    plot(T(3:end),migRate(3:end).*10,'color',[0 0 0])
    plot([1975 2021],[weighted_mean_lag weighted_mean_lag].*10,':r','Linewidth',2.5)
    title(['Migration Rates (Time-weighted Peak Lags)'])
    xlabel('Time [y]');
    ylabel('d\theta/dt [\circ/yr]');
    % ylim([1975 2021]);
    % xlim([-100 79]);
    % set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    % legend(num2str(T(tt-timegap)), num2str(T(tt)),'location','southeast');
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_XCorr'];
    dimensions = [20 18]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    
    
    %% REVIEWER #2 QUESTION: HOW DOES GRID SIZE VARY AS F'N OF DIST FROM ORIGIN?
    figure(1000)
    set(gcf,'Color','w')
    clf
    
    subplot(1,2,1)
    nearCell = [xPolar(1,1),yPolar(1,1);...
        xPolar(1,2),yPolar(1,2);...
        xPolar(2,2),yPolar(2,2);...
        xPolar(2,1),yPolar(2,1);...
        xPolar(1,1),yPolar(1,1)].*1000;
    plot(nearCell(:,1),nearCell(:,2))
    cellAreaNear = polyarea(nearCell(:,1),nearCell(:,2))
    sqrt(cellAreaNear)
    lNear = sqrt((nearCell(2,1)-nearCell(1,1)).^2+...
        (nearCell(2,2)-nearCell(1,2)).^2);
    wNear = sqrt((nearCell(3,1)-nearCell(2,1)).^2+...
        (nearCell(3,2)-nearCell(2,2)).^2);
    
    axis equal
    title(['Near Cell: A=' num2str(cellAreaNear,'%.0f') 'm^2 LxW='...
        num2str(lNear,'%.0f') 'x' num2str(wNear,'%.0f') 'm']);
    
    subplot(1,2,2)
    farCell = [xPolar(150,1),yPolar(150,1);...
        xPolar(150,2),yPolar(150,2);...
        xPolar(151,2),yPolar(151,2);...
        xPolar(151,1),yPolar(151,1);...
        xPolar(150,1),yPolar(150,1)].*1000;
    plot(farCell(:,1),farCell(:,2))
    cellAreaFar = polyarea(farCell(:,1),farCell(:,2))
    sqrt(cellAreaFar)
    lFar = sqrt((farCell(2,1)-farCell(1,1)).^2+...
        (farCell(2,2)-farCell(1,2)).^2);
    wFar = sqrt((farCell(3,1)-farCell(2,1)).^2+...
        (farCell(3,2)-farCell(2,2)).^2);
    
    axis equal
    title(['Far Cell: A=' num2str(cellAreaFar,'%.0f') 'm^2 LxW='...
        num2str(lFar,'%.0f') 'x' num2str(wFar,'%.0f') 'm']);
    
    
    %% Check velocity of rotation at centre versus outside?
    
    vAnomalyRhoEast = squeeze(nansum(vAnomalyPolar(:,125:180,:),2));
            
    % collapse along Theta dimension
    vAnomalyThetaDistal = squeeze(nansum(vAnomalyPolar(76:end,:,:),1));
    vAnomalyThetaProximal = squeeze(nansum(vAnomalyPolar(1:3,:,:),1));
    
    
    figure(17)
    clf
    set(gcf,'Color','w');
    
    % hovmuller diagram (Theta - proximal)
    subplot(1,2,1)
    hold on; box on; grid on;
    pp=pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T)),repmat(T,length(thetaDegNorth(1:end-1)),1),vAnomalyThetaProximal);
    pp.EdgeColor='none';
    % cb=colorbar('Location','southoutside');
    title('(a) Proximal Rotation')
    xlabel('\theta [\circ]');
    ylabel('Time [y]');
    ylim([1975 2021]);
    xlim([-100 79]);
    caxis([-0.51 0.51])
    set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    % ylabel(cb,'V_{anomaly} [Mm^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    % cb.Ticks = [-0.5:0.1:0.5];
    % clim([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
    % xlim([XLim(1) XLim(2)]);
    % clim([CLim(1) CLim(2)]);
    
    % hovmuller diagram (Theta - distal)
    subplot(1,2,2)
    hold on; box on; grid on;
    pp=pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T)),repmat(T,length(thetaDegNorth(1:end-1)),1),vAnomalyThetaDistal);
    pp.EdgeColor='none';
    % cb=colorbar('Location','southoutside');
    title('(b) Distal Rotation')
    xlabel('\theta [\circ]');
    ylabel('Time [y]');
    ylim([1975 2021]);
    xlim([-100 79]);
    caxis([-0.51 0.51])
    set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    % ylabel(cb,'V_{anomaly} [Mm^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    % cb.Ticks = [-0.5:0.1:0.5];
    % clim([-max(abs(get(gca,'CLim'))) max(abs(get(gca,'CLim')))]);
    % xlim([XLim(1) XLim(2)]);
    % clim([CLim(1) CLim(2)]);
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_InnerVsOuter'];
    dimensions = [20 18]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% Cross correlation - proximal
    
    timegap = 1;
    iCutoff = 31; % bin at which to cut off calculations (31 = -70deg)
    for tt=3:length(T)
        
        figure(17171)
        set(gcf,'Color','w');
        clf
        
        subplot(2,1,1)
        hold on; box on; grid on;
        plot(thetaDegNorth(iCutoff:end-1),vAnomalyThetaProximal(iCutoff:end,tt-timegap),'color',[0.5 0.5 0.5])
        plot(thetaDegNorth(iCutoff:end-1),vAnomalyThetaProximal(iCutoff:end,tt),'color',[0 0 0])
        title(['(a) Proximal Volume Anomaly (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ') [\theta]'])
        xlabel('\theta [\circ]');
        ylabel('Time [y]');
        % ylim([1975 2021]);
        xlim([-100 79]);
        set(gca,'XTick',[-180:22.5:180])
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        legend(num2str(T(tt-timegap)), num2str(T(tt)),'location','southeast');
        
        subplot(2,1,2)
        hold on; box on; grid on;
        
        [proximalXc(tt,:),proximalLags(tt,:)] = xcorr(vAnomalyThetaProximal(iCutoff:end,tt),...
            vAnomalyThetaProximal(iCutoff:end,tt-timegap),25,'normalized');
        stem(proximalLags(tt,:),proximalXc(tt,:));
        proximalPeakXC(tt) = max(proximalXc(tt,:));
        proximalPeakLag(tt) = proximalLags(tt,find(proximalXc(tt,:)==max(proximalXc(tt,:)),1));
        proximalMigRate(tt) = proximalPeakLag(tt)./(T(tt)-T(tt-timegap));
        title(['(b) Proximal Cross Correlation (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ')'])
        xlabel('\theta [\circ]');
        ylabel('R_{norm} [-]');
        ylim([-1 1]);
        xlim([-25 25]);
        set(gca,'XTick',[-180:5:180])
        text(-19,0.8,['Peak Cross Correlation: ' num2str(proximalPeakLag(tt)) '^o'],...
            'FontName','Myriad Pro','FontSize',10);
        text(-19,0.6,['Migration Rate: ' num2str(proximalMigRate(tt),'%.1f') '^o/year'],...
            'FontName','Myriad Pro','FontSize',10);
        set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
        set(gca,'FontName','Myriad Pro','FontSize',fontSize)
        
        % Export figure to png file
        pngName = ['HovmullerPolarTheta_ProximalCrossCorrelation_' ...
            num2str(T(tt-timegap)) '-' num2str(T(tt))];
        dimensions = [30 25]./1.5; % [width height]
        printFigs(dimensions,[outPath filesep pngName],0);
    end
    
    
    %%
    figure(1717)
    set(gcf,'Color','w');
    clf
    
    scaleFactor = 1;
    minPeakDist = 8;
    prom = 0.1;
    numPeaks = 1;
    
    for ii=3:22
        hold on
        % plot(vAnomalyTheta(:,ii))
        findpeaks(proximalXc(ii,:),proximalLags(ii,:),'Annotate','extents',...
            'SortStr','descend','NPeaks',numPeaks)
        [proximalXc_peaks(ii).pks,proximalXc_peaks(ii).locs,proximalXc_peaks(ii).width,proximalXc_peaks(ii).prom] ...
            = findpeaks(proximalXc(ii,:),proximalLags(ii,:),...
            'SortStr','descend','NPeaks',numPeaks);
        text(proximalXc_peaks(ii).locs+.02,proximalXc_peaks(ii).pks+T(ii),num2str((1:numel(proximalXc_peaks(ii).pks))'))
    end
    title('Proximal Cross Correlation Through Time [\theta]')
    xlabel('lag [\circ]');
    ylabel('R_{norm} [-]');
    % ylim([1975 2021]);
    % xlim([-100 79]);
    % set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
     % Export figure to png file
    pngName = ['HovmullerPolarTheta_XCorr_Peakfinder'];
    dimensions = [30 25]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
    %% hovmuller diagram of PROXIMAL cross-correlation
    
    proximalWeighted_mean_lag = nanmean(proximalPeakLag(3:end)./diff(T(2:end)));
    mean(migRate(3:end));
    
    proximalWt_mean_lag = nanmean(proximalPeakLag(3:end)./diff(T(2:end)));
    proximalWt_std_lag = nanstd(proximalPeakLag(3:end)./diff(T(2:end)));
    proximalWt_mean_peak = nanmean(proximalPeakXC(3:end)./diff(T(2:end)));
    
    fid = fopen([outPath filesep 'proximal_xcorr_stats.txt'],'w');
    fprintf(fid,[strrep(outPath,'\','/'),'\n']);
    fprintf(fid,['mean mig rate = ' num2str(nanmean(proximalMigRate(3:end)).*10,'%.2f') '\n']);
    fprintf(fid,['mean width = ' num2str(nanmean([proximalXc_peaks(3:end).width]),'%.2f') '\n']);
    fprintf(fid,['mean peak = ' num2str(nanmean([proximalXc_peaks(3:end).pks]),'%.2f') '\n']);
    fprintf(fid,['mean prom = ' num2str(nanmean([proximalXc_peaks(3:end).prom]),'%.2f') '\n']);
    fclose('all');
    
    
    nanmean([proximalXc_peaks(3:end).width])
    
    figure(1001)
    clf
    set(gcf,'Color','w');
    
    subplot(2,1,1)
    hold on; box on; grid on;
    pp=pcolor(repmat(proximalLags(end,:)',1,length(T)),repmat(T,length(proximalLags),1),xc');
    pp.EdgeColor='none';
    plot([0 0],[1975 2021],'-k','Linewidth',2.5)
    plot([proximalXc_peaks(3:end).locs],T(3:end),'.-k','Linewidth',1.0)
    cb=colorbar('Location','eastoutside');
    title('Cross-correlation Through Time')
    xlabel('Lag [\circ]');
    ylabel('Time [y]');
    ylim([1975 2021]);
    xlim([-25 25]);
    caxis([-1 1])
    set(gca,'XTick',[-25:5:25])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    ylabel(cb,'R_{norm}','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    
    subplot(2,1,2)
    hold on; box on; grid on;
    plot(T(3:end),proximalMigRate(3:end).*10,'color',[0 0 0])
    plot([1975 2021],[proximalWeighted_mean_lag proximalWeighted_mean_lag].*10,':r','Linewidth',2.5)
    title(['Migration Rates (Time-weighted Peak Lags)'])
    xlabel('Time [y]');
    ylabel('d\theta/dt [\circ/yr]');
    % ylim([1975 2021]);
    % xlim([-100 79]);
    % set(gca,'XTick',[-180:22.5:180])
    set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
    set(gca,'FontName','Myriad Pro','FontSize',fontSize)
    % legend(num2str(T(tt-timegap)), num2str(T(tt)),'location','southeast');
    
    % Export figure to png file
    pngName = ['HovmullerPolarTheta_Proximal_XCorr'];
    dimensions = [20 18]./1.5; % [width height]
    printFigs(dimensions,[outPath filesep pngName],0);
    
end


