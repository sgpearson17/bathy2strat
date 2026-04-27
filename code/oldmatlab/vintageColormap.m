function [cmap] = vintageColormap(vertSpacing)
    %% vintageColormap.m 
    % Vintage bathymetry colourmap from Edwin/Albert's old 1950s Wadden Sea Maps
    % See Elias et al 2019
    %
    % takes in vertSpacing, a 10-column array with elevations of each
    % contour
        % intervals at which colour changes
%    vertSpacing = [-20; -16; -12; -8; -5; -3; -1.2; 0; 1.4; 10]; 

    dz=0.1; % vertical spacing

    % colours in colormap
    baseColours = ...
        [27 126 129;... % Deepest (darkest green-blue): -20 m NAP
        41 155 151;... % Deeper (dark green-blue): -16 m NAP
        56 170 164;... % Deep (green-blue): -12 m NAP
        130 199 180;... % Subtidal (light green-blue): -8 m NAP
        220 231 194;... % Shallow Subtidal (lighter green blue): -5 m NAP
        255 240 196;... % Shoal (off-white): -3 m NAP
        244 214 176;... % Low Tide (light brown): -1.2 m NAP
        217 188 146;... % High Tide (darker brown): 0 m NAP
        255 221 146;... % Beach (yellowy brown): 1.4 m NAP
        226 129 61]./255; % Dune (orange): 10 m NAP



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

