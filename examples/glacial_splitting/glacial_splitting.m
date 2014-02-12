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

% (C) Alan Baird and Andrew Walker, 2014
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
        
        % Single crystal elasticity of this layer
        [C, rho] = ice_cij(all_data(i).Tav);
        
        % Calculate poly xtal elasticity of this layer 
        [Cs(:,:,j), rhos(j)] = Cs_from_EBSD_file(C, rho, ...
                                             all_data(i).tex_file);

        report_layer(i, all_data(i).tex_file, Cs(:,:,j), rhos(j), ...
            all_data(i).dtb, dtt, all_data(i).Tav, plot_pole)
    end

    % The effective splitting calculation
    fprintf('\nEffective splitting calculation \n');
    fprintf('--------------------------------\n');
    % Set up inclination and azimuth and effective splitting
    % params. These two could be functions of spol, for eg.
    inc = 90.0;
    azi = 0.0;
    % Loop over source polarizations
    % and frequencies
    spol = 0:15:180; % Deg
    spol = min_azi:del_azi:max_azi; % Deg
    freq = [0.3, 3.0, 30.0]; % Hz
    fast_eff = zeros(length(freq),length(spol));
    tlag_eff = zeros(length(freq),length(spol));
    for f = 1:length(freq)
        for s = 1:length(spol)
    
            fprintf('\n');
            fprintf('For source polarization: %f (deg)\n', spol(s));
            fprintf('and frequency: %f (Hz)\n', freq(f));

            [fast_eff(f,s), tlag_eff(f,s)] = do_effective_splitting(Cs, ...
                     rhos, thickness, inc, azi, freq(f), spol(s), ...
                     plot_waves, eff_split_mode); 
    
            % FIXME: do we need to correct fast_eff here? We are working in
            % ray frame at the momenet.
            % fast(j) = MS_unwind_pm_90((azi+pol')) ; % geog. ref frame
            fprintf('Effective fast direction: %f (deg)\n', fast_eff(f,s));
            fprintf('Effective delay time:     %f (s)\n', tlag_eff(f,s));
        end
    end 

    % Make a figure to plot the results
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;

    % Plot lag times
    subplot(2,1,1)
    for f = 1:length(freq)
        plot(spol,tlag_eff(f,:))
        hold on
    end
    xlabel('Initial polarisation (deg)')
    ylabel('Lag time (s)')

    % Plot fast directions
    subplot(2,1,2)
    for f = 1:length(freq)
        plot(spol,fast_eff(f,:))
        hold on
    end
    xlabel('Initial polarisation (deg)')
    ylabel('Fast directions (deg)')
end


function [fast_eff, tlag_eff] = do_effective_splitting(Cs, rhos, ...
             thickness, inc, azi, freq, spol, plot_waves, eff_split_mode)

        % Header line for table...
        fprintf('time lag (s)     fast direction (deg)\n');
        % Loop over layers and calculate splitting parameters 
        fast = zeros(1,length(rhos));
        tlag = zeros(1,length(rhos));
        for j = 1:length(rhos)

            % Calculate the splitting parameters for this layer
            [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( Cs(:,:,j), ...
                rhos(j), inc, azi );
            fast(j) = pol; % FIXME: do we need to convert to geog ref? 
            tlag(j) = thickness(j)/vs2 - thickness(j)/vs1 ;

            fprintf('  %7f        %7f \n', tlag(j), fast(j));
        end
    
        % Calculate effective splitting
        if plot_waves
            [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,...
                tlag,'mode',eff_split_mode,'PlotWavelet');
        else 
            [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,...
                tlag,'mode',eff_split_mode);
        end
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

    % FIXME: AFB - A comment here as to where these layers 
    % come from is critical
    data(1) = struct('tex_file', fullfile(dat_dir, 'V97-248R.txt'), ...
                     'dtb', 450, ... % depth to base of layer m
                     'Tav', -25, ... % Average temperature (deg. C)
                     'azi', 0); % Horizonal rotation for texture
    data(2) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ...
                     'dtb', 2850, ... 
                     'Tav', -20, ... 
                     'azi', 0);    
    data(3) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ...
                     'dtb', 2900, ... 
                     'Tav', -15, ... 
                     'azi', 0); 
    data(4) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ...
                     'dtb', 3150, ... 
                     'Tav', -10, ... 
                     'azi', 0);    
    data(5) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ...
                     'dtb', 3225, ... 
                     'Tav', -9, ... 
                     'azi', 0); 
    data(6) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ...
                     'dtb', 3300, ... 
                     'Tav', -7, ... 
                     'azi', 0);    
    data(7) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ...
                     'dtb', 3350, ... 
                     'Tav', -5, ... 
                     'azi', 0); 
    data(8) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ...
                     'dtb', 3450, ... 
                     'Tav', -1, ... 
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
    fgetl(fid); % Header line - ignore
    data = fscanf(fid, '%f', [12 inf]);
    nxtl = length(data(1,:));
    eulers = zeros(3,nxtl);
    eulers(1,:) = data(3,:); % phi1
    eulers(2,:) = data(4,:); % phi1
    eulers(3,:) = data(5,:); % phi1
    fclose(fid);
end


function report_layer(layernum, filename, Cvrh, rho, dtb, dtt, tav, ...
    plot_pole)

    fprintf('\nLayer: %i\n', layernum);
    fprintf(['data file: %s, \ntop: %f m, base: %f m, \n'...
        'thickness: %f m, teperature: %f C\n'], ...
        filename, dtt, dtb, dtb-dtt, tav);
   
    if plot_pole
        MS_plot(Cvrh, rho, 'wtitle', filename, 'fontsize', 11, ...
            'avscontours', 0:0.2:10.16, 'pcontours', 3.60:0.01:3.80, ...
            'polsize', 0.18, 0.16, 2.0, 1.0, 'limitsonpol', 'silent');
    end
end

