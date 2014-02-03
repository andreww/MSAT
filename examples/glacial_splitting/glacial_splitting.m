% Comments

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

function [fast_eff,tlag_eff] = glacial_splitting()

    all_data = get_data(); % This should be optional - argument needed
    
    % Build Cij and thickness arrays - bottom to top
    thickness = zeros(1,length(all_data));
    Cs = zeros(6,6,length(all_data));
    rhos = zeros(1,length(all_data));
    fast = zeros(1,length(all_data));
    tlag = zeros(1,length(all_data));
    j = 0;
    for i=length(all_data):-1:1
        j = j + 1; % switch array order
        
        % Thickness of this layer
        if (i==1) 
            thickness(j) = (all_data(i).dtb)/1000.0; % In km!
            dtt = 0;
        else
            thickness(j) = (all_data(i).dtb - all_data(i-1).dtb)/1000.0;%km
            dtt = all_data(i-1).dtb;
        end
        
        % Single crystal elasticity of this layer
        [C, rho] = ice_cij(all_data(i).Tav);
        
        % Calculate poly xtal elasticity of this layer 
        [Cs(:,:,j), rhos(j)] = Cs_from_EBSD_file(C,rho,all_data(i).tex_file);
        
        report_layer(i, all_data(i).tex_file, Cs(:,:,j), rhos(j), ...
            all_data(i).dtb, dtt, all_data(i).Tav)
        
        % Set up inclination and azimuth
        inc=90.0;
        azi=0.0;
        
        [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( Cs(:,:,j), ...
            rhos(j), inc, azi );
        fast(j) = pol; % FIXME: do we need to convert to geog ref? 
        tlag(j) = thickness(j)/vs2 - thickness(j)/vs1 ;
        
    end
    
    % Calculate effective splitting
    freq = 30; % 30 Hz
    spol = 45;
    
    %[fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,tlag);
    [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,...
        tlag,'mode','GaussianWavelet','PlotWavelet');
    
    % FIXME: do we need to correct fast_eff here? We are working in
    % ray frame at the momenet.
    % fast(j) = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
    
    fprintf('\n');
    fprintf('Effective fast direction: %f (deg)\n', fast_eff);
    fprintf('Effective delay time:     %f (s)\n', tlag_eff);
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
    % FIXME: Do we need to use P too?
    % results are in GPa and kg/m^3
    
    % For now, just use one elasticty 
    % FIXME: add a sensible parameterisation
    % AMW to fix this.
    
    [C, r] = MS_elasticDB('ice');
    
end

function [data] = get_data()
    % Return default data structure
    % elements run from surface to base
    % (i.e. they need reversing for 
    % MS_effective_splitting_N)
    
    % Deal with file seperators in a cross platform way.
    dat_dir = fullfile('~andreww', 'Science', 'Ice_EBSD_data');
    
    data(1) = struct('tex_file', fullfile(dat_dir, 'V97-248R.txt'), ...
                     'dtb', 2850, ... % depth to base of layer m
                     'Tav', -10, ... % Average temperature (deg. C)
                     'azi', 0); % Horizonal rotation for texture
    data(2) = struct('tex_file', fullfile(dat_dir, 'V3311g.txt'), ...
                     'dtb', 3100, ... 
                     'Tav', -20, ... 
                     'azi', 10);             
    data(3) = struct('tex_file', fullfile(dat_dir, 'V3321c.txt'), ...
                     'dtb', 3600, ... 
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

function report_layer(layernum, filename, Cvrh, rho, dtb, dtt, tav)

    fprintf('Layer: %i\n', layernum);
    fprintf(['data file: %s, \ntop: %f m, base: %f m, \n'...
        'thickness: %f m, teperature: %f C\n'], ...
        filename, dtt, dtb, dtb-dtt, tav); 
    
    MS_plot(Cvrh, rho, 'wtitle', filename, 'fontsize', 11, ...
        'avscontours', 0:0.2:10.16, 'pcontours', 3.8:0.01:3.99, ...
        'polsize', 0.18, 0.16, 2.0, 1.0, 'limitsonpol');

end
