% GLACIAL_SPLITTING.M - and example script for modelling glacial anisotropy
% 
% This script allows textural data measured from samples recovered
% from an ice sheet core to be used to calculate the expected 
% shear wave splitting parameters for a various possible observations.
% Input data is in the form of measured LPOs of ice and layer information.
% This is converted into the polycrystal elasticity and then into layer
% by layer splitting parameters. These parameters are combined to form
% the effective splitting.
% 
% The script is packaged as a function taking a number of optional
% arguments but it can be run without arguments. The arguments can
% be given in any order and are:
% 
%   glacial_splitting(..., 'phasepole')
%       Generate phase velocity / splitting pole figures
%   glacial_splitting(..., 'splittingplot')
%       Plot the wavelets for each effective splitting calc
%   glacial_splitting(..., 'mode', eff_split_mode)
%       String passed into the mode argument of MS_effective_splitting_N.
%       Note that setting this to 'GaussianWavelet' results in a much 
%       longer execution time.

% (C) Alan Baird and Andrew Walker, 2014, 2015
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [fast_eff,tlag_eff] = glacial_splitting(varargin)

    % Setup defaults for options
    plot_pole = 0; % Phase velocity plot
    plot_waves = 0; % Gaussian wavelemt plot
    eff_split_mode = 's&s';
    min_azi = 0.0;
    max_azi = 180.0;
    del_azi = 5.0;
    quiet = 0;

    % Process those optional arguments
    iarg = 1;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'phasepole'
                plot_pole = 1;
                iarg = iarg + 1;
            case 'splittingplot'
                plot_waves = 1;
                iarg = iarg + 1;
            case 'mode'
                eff_split_mode = lower(varargin{iarg+1});
                iarg = iarg + 2;
            case 'min_azi'
                min_azi = varargin{iarg+1};
                iarg = iarg + 2;
            case 'max_azi'
                max_azi = varargin{iarg+1};
                iarg = iarg + 2;
            case 'del_azi'
                del_azi = varargin{iarg+1};
                iarg = iarg + 2;
            case 'quiet'
                quiet = 1;
                iarg = iarg + 1;
            otherwise
                warning(['Unknown option: ' varargin{iarg}]) ;
                iarg = iarg + 1;
        end
    end

    % Report all run conditions
    tf{1} = 'false';
    tf{2} = 'true';
    fprintf('\n\n');
    fprintf('GLACIAL SPLITTING\n');
    fprintf('+++++++++++++++++\n');
    fprintf('\n');
    fprintf('Run conditions\n');
    fprintf('--------------\n');
    fprintf('Ice texture data from: %s\n', 'internal function');
    fprintf('Will plot pole figures: %s\n', tf{plot_pole+1});
    fprintf(['Effective splitting for azimuths between \n%4f and %4f '...
        'degrees with a step size \nof %4f degrees\n'], ...
        min_azi, max_azi, del_azi);
    fprintf('Source polarization equal to azimuth...\n');
    fprintf('Effective splitting mode: %s\n', eff_split_mode);
    fprintf('Will plot particle motion in splitting calc: %s\n', ...
        tf{plot_waves+1});
    
    
    all_data = get_data(); % This should be optional - argument needed
    
    % Build Cij and thickness arrays - bottom to top
    thickness = zeros(1,length(all_data));
    Cs = zeros(6,6,length(all_data));
    rhos = zeros(1,length(all_data));
    ddts = zeros(1,length(all_data));
    ddbs = zeros(1,length(all_data));
    Vpiso = zeros(1,length(all_data));
    Vsiso = zeros(1,length(all_data));
    AVp = zeros(1,length(all_data));
    AVs_max = zeros(1,length(all_data));
    j = 0;
    fprintf('\nLayer data (bottom to top) \n');
    fprintf('---------------------------\n');
    for i = length(all_data):-1:1
        j = j + 1; % switch array order (read i, write j)
        
        % Thickness of this layer
        if (i==1) 
            thickness(j) = (all_data(i).dtb) / 1000.0; % In km!
            dtt = 0;
        else
            thickness(j) = (all_data(i).dtb - all_data(i-1).dtb) / 1000.0;
            dtt = all_data(i-1).dtb;
        end
        
        dtts(j) = dtt;
        dtbs(j) = all_data(i).dtb;
        temps(j) = all_data(i).Tav;
        
        % Single crystal elasticity of this layer
        [C, rho] = ice_cij(all_data(i).Tav);
        
        % Calculate poly xtal elasticity of this layer 
        [Cs(:,:,j), rhos(j)] = Cs_from_EBSD_file(C, rho, ...
                                             all_data(i).tex_file);
        
        % Calculate average isotropic (VRH) velocity for this layer        
        [ K_vrh, G_vrh ] = MS_polyaverage( C );
        Vpiso(j) = 1e-3.*sqrt(((K_vrh*1e9)+((4.0/3.0)*G_vrh)*1e9)./rhos(j));
        Vsiso(j) = 1e-3.*sqrt((G_vrh*1e9)./rhos(j));
        % We may need the upper and lower bounds (V and R, respectivly)

        % Anisotropy data for plot
        [AVp(j), AVs_max(j)] = anisotropic_data(Cs(:,:,j), rhos(j));
        
        report_layer(i, all_data(i).tex_file, Cs(:,:,j), rhos(j), ...
            all_data(i).dtb, dtt, all_data(i).Tav, plot_pole)
    end

    
    % Plot stratigraphic column
    
    figure();
    subplot(1,6,1);
    bar([fliplr(thickness);nan(size(thickness))],'stack');xlim([0.5,1.5]);ylim([0,4]);
    set(gca,'YDir','Reverse');
    hold on;
    samples=[97 248 399 650 1201 1500 1849 2082 2749 2874 3300 3311 3321 3329 3399 3416]./1000.0;
    %scatter(ones(size(samples)),samples,'MarkerEdgeColor','k','MarkerFaceColor','w')
    scatter(linspace(0.7,1.3,max(size(samples))),samples,'MarkerEdgeColor','k','MarkerFaceColor','w')
    ylabel('Depth (km)')
    hold off;
    
    subplot(1,6,2);
    stairs(temps,dtts./1000);xlim([-25,0]);ylim([0,4]);
    ylabel('Depth (km)');
    xlabel('Temperature (C)');
    set(gca,'YDir','Reverse');
    
    subplot(1,6,3);
    scatter(Vpiso,dtts./1000);ylim([0,4]);
    ylabel('Depth (km)');
    xlabel('Isotropic p-wave velocity (km/s)');
    set(gca,'YDir','Reverse');

    subplot(1,6,4);
    scatter(Vpiso,dtts./1000);ylim([0,4]);
    ylabel('Depth (km)');
    xlabel('Isotropic p-wave velocity (km/s)');
    set(gca,'YDir','Reverse');
    
    subplot(1,6,5);
    scatter(AVp,dtts./1000);ylim([0,4]);
    ylabel('Depth (km)');
    xlabel('P-wave anisotropy (%)');
    set(gca,'YDir','Reverse');

    subplot(1,6,6);
    scatter(AVs_max,dtts./1000);ylim([0,4]);
    ylabel('Depth (km)');
    xlabel('Maximum s-wave anisotropy (%)');
    set(gca,'YDir','Reverse');
    
    % Create model for ATRAK
    create_modsimple_rsp('modsimple.rsp', dtts, dtbs, Cs, rhos)
    

    % The effective splitting calculation for vertical ray path
    fprintf('\nEffective splitting calculation \n');
    fprintf('--------------------------------\n');
    % Set up inclination and azimuth and effective splitting
    % params. These two could be functions of spol, for eg.
    inc = 90.0;
    azi = 0.0;
    % Loop over source polarizations
    % and frequencies
    spol = min_azi:del_azi:max_azi; % Deg
    freq = [3, 30, 3000]; % Hz
    freq_c = ['r', 'g', 'b'];
    freq_names = {'3 Hz', '30 Hz', '3000 Hz'};
    fast_eff = zeros(length(freq),length(spol));
    tlag_eff = zeros(length(freq),length(spol));
    for f = 1:length(freq)
        for s = 1:length(spol)
    
            if ~quiet
                fprintf('\n');
                fprintf('For source polarization: %f (deg)\n', spol(s));
                fprintf('and frequency: %f (Hz)\n', freq(f));
            end
            
            [fast_eff(f,s), tlag_eff(f,s)] = do_effective_splitting(Cs, ...
                     rhos, thickness, inc, azi, freq(f), spol(s), ...
                     plot_waves, eff_split_mode, quiet);           

            if ~quiet
                fprintf('Effective fast direction: %f (deg)\n', fast_eff(f,s));
                fprintf('Effective delay time:     %f (s)\n', tlag_eff(f,s));
            end
        end
    end 

    % Make a figure to plot the results
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;

    % Plot lag times
    subplot(2,1,1)
    for f = 1:length(freq)
        plot(spol,tlag_eff(f,:), freq_c(f))
        hold on
    end
    xlabel('Initial polarisation (deg)')
    ylabel('Lag time (s)')

    % Plot fast directions
    subplot(2,1,2)
    for f = 1:length(freq)
        plot(spol,fast_eff(f,:), freq_c(f))
        hold on
    end
    legend(freq_names)
    xlabel('Initial polarisation (deg)')
    ylabel('Fast directions (deg)')
    
    
    % The effective splitting calculation for inclined rays
    fprintf('\nEffective splitting calculation - inclined rays\n');
    fprintf('------------------------------------------------\n');
    % Set up inclination and azimuth and effective splitting
    % params. These two could be functions of spol, for eg.
    inc = [60 70 80] ;
    azi = [0 45 90] ;
    % Loop over source polarizations
    % and frequencies
    spol = 0:1:180; % Deg
    freq = 30.0; % Hz - no point looking lower for SWS
    inc_c = ['r', 'g', 'b'];
    inc_names = {'60 degrees', '70 degrees', '80 degrees'};
    fast_eff_2 = zeros(length(spol), length(inc), length(azi));
    tlag_eff_2 = zeros(length(spol), length(inc), length(azi));
    for i = 1:length(inc)
        for a = 1:length(azi)
            for s = 1:length(spol)
            
                if ~quiet
                    fprintf('\n');
                    fprintf('For source polarization: %f (deg)\n', spol(s));
                    fprintf('azimuth: %f (deg)\n', azi(a));
                    fprintf('inclination: %f (deg)\n', inc(i));
                end
            
                % Recalculate thickness etc
                
                [fast_eff_2(s,i,a), tlag_eff_2(s,i,a)] = do_effective_splitting(Cs, ...
                        rhos, thickness, inc(i), azi(a), freq, spol(s), ...
                         plot_waves, eff_split_mode, quiet); 
           

                if ~quiet
                    fprintf('Effective fast direction geographical frame: %f (deg)\n', fast_eff_2(s,i,a));
                    fprintf('Effective delay time:     %f (s)\n', tlag_eff_2(s,i,a));
                end
            end
        end
    end
    
    % Make a figure to plot the results
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;
    % Plot lag times
    subplot(2,3,1)
    for i = 1:length(inc)
        plot(spol, tlag_eff_2(:,i,1), inc_c(i))
        hold on
    end
    xlabel('Source polarization (deg)')
    ylabel('Lag time (s)')
    thetitle = sprintf('Splitting as a function of source polarization and inclination, azi = %f (deg)', azi(1));
    title(thetitle);
    
    % Plot fast directions
    subplot(2,3,4)
    for i = 1:length(inc)
        plot(spol ,fast_eff_2(:,i,1), inc_c(i))
        hold on
    end
    legend(inc_names)
    xlabel('Source polarization (deg)')
    ylabel('Fast directions (deg)')

    
    subplot(2,3,2)
    for i = 1:length(inc)
        plot(spol, tlag_eff_2(:,i,2), inc_c(i))
        hold on
    end
    xlabel('Source polarization (deg)')
    ylabel('Lag time (s)')
    thetitle = sprintf('Splitting as a function of source polarization and inclination, azi = %f (deg)', azi(2));
    title(thetitle);

    % Plot fast directions
    subplot(2,3,5)
    for i = 1:length(inc)
        plot(spol ,fast_eff_2(:,i,2), inc_c(i))
        hold on
    end
    legend(inc_names)
    xlabel('Source polarization (deg)')
    ylabel('Fast directions (deg)')
    
    subplot(2,3,3)
    for i = 1:length(inc)
        plot(spol, tlag_eff_2(:,i,3), inc_c(i))
        hold on
    end
    xlabel('Source polarization (deg)')
    ylabel('Lag time (s)')
    thetitle = sprintf('Splitting as a function of source polarization and inclination, azi = %f (deg)', azi(3));
    title(thetitle);

    % Plot fast directions
    subplot(2,3,6)
    for i = 1:length(inc)
        plot(spol ,fast_eff_2(:,i,3), inc_c(i))
        hold on
    end
    legend(inc_names)
    xlabel('Source polarization (deg)')
    ylabel('Fast directions (deg)')
    
end


function [fast_eff, tlag_eff] = do_effective_splitting(Cs, rhos, ...
             thickness, inc, azi, freq, spol, plot_waves, eff_split_mode,...
             quiet)

        % Header line for table...
        if ~quiet
            fprintf('time lag (s)     fast direction in ray frame (deg)\n');
        end
        % Loop over layers and calculate splitting parameters 
        fast = zeros(1,length(rhos));
        tlag = zeros(1,length(rhos));
        for j = 1:length(rhos)

            % Calculate the splitting parameters for this layer
            [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( Cs(:,:,j), ...
                rhos(j), inc, azi );
            fast(j) = pol; 
            % Calculate tlag, ajusting for inclination (vertical ray has
            % inc=90, so this is just length = thickness/1. Gets thicker
            % as ray becomes shallower
            tlag(j) = (thickness(j)/cosd(90-inc))/vs2 - ...
                (thickness(j)/cosd(90-inc))/vs1 ; 
            
            if ~quiet
                fprintf('  %7f        %7f \n', tlag(j), fast(j));
            end
        end
    
        % Calculate effective splitting
        if plot_waves
            [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,...
                tlag,'mode',eff_split_mode,'PlotWavelet');
        else 
            [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,...
                tlag,'mode',eff_split_mode);
        end
        
        % Put fast orentation into geographic reference frame. For 
        % Vertical ray azi = 0, so no change.
        fast_eff = MS_unwind_pm_90(azi+fast_eff);
end


function [C_poly, rho_poly] = Cs_from_EBSD_file(C_single, ...
                                                rho_single, filename)
    % Calculate poly-xtal elasticity of ice 
    % sample from EBSD data and single crystal 
    % elasticity
    
    % Read Eulers from data file
    [eulers, nxtls] = read_EBSD_txt(filename);
          
    % List of elasticities and densities for each crystal 
    Cs = zeros(6,6,nxtls);
    rhos = zeros(1,nxtls);
    for i = 1:nxtls
       Cs(:,:,i) = C_single(:,:);
       rhos(i) = rho_single;
    end

    % Rotate the Cs to give texture
    Cs = MS_rotEuler(Cs, eulers(1,:)', eulers(2,:)', eulers(3,:)', ...
        'sense', 'passive');
  
    % Form VRH mean - each grain has the same volume.
    [C_poly, rho_poly] = MS_VRH(ones(nxtls,1), Cs, rhos);
end


function [C, r] = ice_cij(T)
    % Return the single crystal elasticity and 
    % density as a function of temperature in deg. C
    % results are in GPa and kg/m^3. 
    % We ignore pressure. 3 km of ice is ~0.03 GPa. Void
    % space is a much bigger approximation.
        
    % Get C at -16C
    [C0, ~] = MS_elasticDB('ice');
    
    % Correction to T from e.q. 9 of 
    % Gammon et al 1983.
    C = C0*(1-1.42E-3*(T+16.0));
    
    % Density from e.g. 2 of above paper: Note T in deg C!
    fac = 1 + 1.576E-4*T - 2.778E-7*T^2 + 8.850E-9*T^3 - 1.778E-10*T^4;
    r = (1.0/999.840)*fac; %table 2.3 of physics of ice. In useful units
    r = 1/r;
end


function [data] = get_data()
    % Return default data structure
    % elements run from surface to base
    % (i.e. they need reversing for 
    % MS_effective_splitting_N)
    
    % Deal with file seperators in a cross platform way.
    dat_dir = fullfile('~', 'Dropbox', 'Ice_EBSD_data');

    % EBSD data from Obbard and Baker (2007). Depths to deeper layers from Lipenkov and Barkov (1998)
    data(1) = struct('tex_file', fullfile(dat_dir, 'V97.txt'), ... % R-A
                     'dtb', 200, ... % depth to base of layer m
                     'Tav', -24, ... % Average temperature (deg. C)
                     'azi', 0); % Horizonal rotation for texture
    data(2) = struct('tex_file', fullfile(dat_dir, 'V248.txt'), ... % R-B
                     'dtb', 300, ...
                     'Tav', -23, ...
                     'azi', 0);               
    data(3) = struct('tex_file', fullfile(dat_dir, 'V399-3.txt'), ... % R-C
                     'dtb', 454, ...
                     'Tav', -22, ...
                     'azi', 0);
    data(4) = struct('tex_file', fullfile(dat_dir, 'V650.txt'), ... % A0-A
                     'dtb', 900, ... 
                     'Tav', -19, ... 
                     'azi', 0);    
    data(5) = struct('tex_file', fullfile(dat_dir, 'V1201.txt'), ... % A0-B
                     'dtb', 1300, ... 
                     'Tav', -16, ... 
                     'azi', 0);    
    data(6) = struct('tex_file', fullfile(dat_dir, 'V1500.txt'), ... % A0-C
                     'dtb', 1700, ... 
                     'Tav', -14, ... 
                     'azi', 0);    
    data(7) = struct('tex_file', fullfile(dat_dir, 'V1849.txt'), ... % A0-D
                     'dtb', 1900, ... 
                     'Tav', -12, ... 
                     'azi', 0);    
    data(8) = struct('tex_file', fullfile(dat_dir, 'V2082-M.txt'), ... % A0-E
                     'dtb', 2700, ... 
                     'Tav', -7, ... 
                     'azi', 0);    
    data(9) = struct('tex_file', fullfile(dat_dir, 'V2749-132.txt'), ... % A1
                     'dtb', 2835, ... 
                     'Tav', -6, ... 
                     'azi', 0); 
    data(10) = struct('tex_file', fullfile(dat_dir, 'V2874-P9.txt'), ... % B1
                     'dtb', 2905, ... 
                     'Tav', -5, ... 
                     'azi', 0);    
    data(11) = struct('tex_file', fullfile(dat_dir, 'V2082-M.txt'), ... % A2 ------ girdle
                     'dtb', 3140, ... 
                     'Tav', -4, ... 
                     'azi', 0); 
    data(12) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ... % B2 ------- cluster
                     'dtb', 3220, ... 
                     'Tav', -3.5, ... 
                     'azi', 0);    
    data(13) = struct('tex_file', fullfile(dat_dir, 'V3300_M.txt'), ... % A3
                     'dtb', 3310, ... 
                     'Tav', -3, ... 
                     'azi', 0); 
    data(14) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ... % CA  
                     'dtb', 3330, ... 
                     'Tav', -2.8, ... 
                     'azi', 0);  
    data(15) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ... % CB
                     'dtb', 3350, ... 
                     'Tav', -2.6, ... 
                     'azi', 0);  
    data(16) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ... % CA
                     'dtb', 3370, ... 
                     'Tav', -2.5, ... 
                     'azi', 0); 
    data(17) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ... % A4
                     'dtb', 3460, ... 
                     'Tav', -1.9, ... 
                     'azi', 0); 
    data(18) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ... % DA+B
                     'dtb', 3540, ... 
                     'Tav', -1.4, ... 
                     'azi', 0); 
    data(19) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ... % E
                     'dtb', 3605, ... 
                     'Tav', -1, ... 
                     'azi', 0);     
    data(20) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ... % ?
                     'dtb', 3750, ... 
                     'Tav', 0, ... 
                     'azi', 0); 
end


function [eulers, nxtl] = read_EBSD_txt(filename)
    % Given a file name, read the Euler angles 
    % given in the format used by the ice data
    % files. Output is a (3,nxtl) array of 
    % Euler angles (hopefully in Bunge convention
    % and degrees) and the number of crystals, nxtl.

    fid = fopen(filename); % Read - the default
    assert((fid~=-1), 'Could not open file %s', filename);
    head = fgetl(fid); % Header line - use to distiguish the two 
                       % kinds of file (ID as the first word or Index)
    if all(head(1:6) == '    ID')
        % First data format - from Geoff                       
        data = fscanf(fid, '%f', [12 inf]);
        nxtl = length(data(1,:));
        eulers = zeros(3,nxtl);
        eulers(1,:) = data(3,:); % phi1
        eulers(2,:) = data(4,:); % phi2
        eulers(3,:) = data(5,:); % phi3
    elseif all(head(1:5) == 'Index')
        % New data format - from Rachel                       
        data = fscanf(fid, '%f', [12 inf]);
        nxtl = length(data(1,:));
        eulers = zeros(3,nxtl);
        eulers(1,:) = data(5,:); % euler1
        eulers(2,:) = data(6,:); % euler2
        eulers(3,:) = data(7,:); % euler3
    else
        % Don't know what to do
        error('Could not identify the EBSD data format from the 1st line')
    end
    fclose(fid);
end


function report_layer(layernum, filename, Cvrh, rho, dtb, dtt, tav, ...
    plot_pole)

    fprintf('\nLayer: %i\n', layernum);
    fprintf(['data file: %s, \ntop: %f m, base: %f m, \n'...
        'thickness: %f m, temperature: %f C\n'], ...
        filename, dtt, dtb, dtb-dtt, tav);
   
    if plot_pole
        MS_plot(Cvrh, rho, 'wtitle', filename, 'fontsize', 11, ...
            'avscontours', 0:0.2:10.16, 'pcontours', 3.60:0.01:3.80, ...
            'polsize', 0.18, 0.16, 2.0, 1.0, 'limitsonpol', 'quiet');
    end
end


function create_modsimple_rsp(filename, dtts, dtbs, Cs, rhos)

    model_name = 'icesheet.mod';
    assert(length(dtts)==length(dtbs), 'Tops and bottoms must agree')
    assert(length(dtts)==length(rhos), 'Tops and density must agree')
    assert(length(dtts)==length(Cs), 'Tops and Cijs must agree')

    num_layers = length(dtts);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', model_name);
    fprintf(fid, '%d\n', num_layers);

    % Hard code X, Y and Z nodes here. We could make this
    % an argument, but I don't know what the rules are.
    % (10 is the discritisation of the model)
    fprintf(fid, '%d %f %f\n', 10, 0.0, 10000.0);
    fprintf(fid, '%d %f %f\n', 10, -5000.0, -5000.0);
    fprintf(fid, '%d %f\n', 10, -dtbs(num_layers)); % Base of hole
    
    % Now add each layer in turn (remembering we need to go backwards)
    j = 0;
    for i = num_layers:-1:1
        j = j + 1;
        fprintf(fid, '%d %s\n', j, 'flat');
        fprintf(fid, '%f %f\n', -dtts(i), -dtts(i));
        fprintf(fid, '%s\n', 'ani');
        fprintf(fid, '%d\n', 0); 
        fprintf(fid, '%d\n', 21);
        MS_save(fid, Cs(:,:,i), rhos(i), 'eunit', 'pa', 'dunit', 'kgm3', 'Aij');
    end
    fclose(fid);
end

function [AVp, AVS_max] = anisotropic_data(C, rh)
      % This is basically copied from MS_plot. We should
      % pull this out into MS_something.

      grid_spacing = 6;
      
      % Set up inc-az grids...
      [INC,AZ] = meshgrid( 90:-grid_spacing:0 , 0:grid_spacing:360 ) ;
      
      % Invoke MS_phasevels to get wave velocities etc.
      [~,~,vs1,vs2,vp, S1P] = MS_phasevels(C,rh,...
        reshape(INC, ((360/grid_spacing+1)*(90/grid_spacing+1)),1), ...
        reshape(AZ,  ((360/grid_spacing+1)*(90/grid_spacing+1)),1));
    
      % reverse so sph2cart() works properly
      % AZ = -AZ;
      
      % Reshape results back to grids
      VS1 = reshape(vs1,61,16);
      VS2 = reshape(vs2,61,16);
      VP =  reshape(vp,61,16);
      % VS1_x = reshape(S1P(:,1),61,16);
      % VS1_y = reshape(S1P(:,2),61,16);
      % VS1_z = reshape(S1P(:,3),61,16);
        
      % aVS data
      dVS =  (VS1-VS2) ;
      VSmean = (VS1+VS2)./2.0 ;
      AVS = 100.0*(dVS./VSmean) ;
      
      % avp data
      AVS_max = max(max(AVS));
      AVp = (max(VP)- min(VP))/(max(VP)+min(VP)).*200.0;
      
end

