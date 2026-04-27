function [cmap] = kg2Colormap
    %% kg2Colormap.m Stuart's fun bathymetry colourmap
    %   The fun colourmap I used for the kustgenese bathymetric maps

    dz=0.1; % vertical spacing

    % colours in colormap
    baseColours = ...
        [0 67 143;... % dark blue: -20 m NAP
        13 182 255;... % light blue: -10 m NAP
        255 255 255;... % white: -5 m NAP
        199 181 181;... % light brown: -2 m NAP
        158 144 144;... % muddy brown: 2 m NAP
        29 89 74;... % dark green: 10 m NAP
        29 89 74]./255; % dark green: 10 m NAP

    % intervals at which colour changes
    vertSpacing = [-20; -10; -5; -2; 2; 3; 10]; 

    % loop through and create new colourmap
    cc=1;
    for ii = 1:length(baseColours)-1
        for jj = 1:floor((vertSpacing(ii+1)-vertSpacing(ii))./dz)
            for kk = 1:3
                cmap(cc,kk) = interp1([vertSpacing(ii) vertSpacing(ii+1)],...
                    [baseColours(ii,kk) baseColours(ii+1,kk)],...
                    vertSpacing(ii)+jj.*dz);
            end
            cc=cc+1;
        end
    end

end

