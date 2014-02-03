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
    j = 0;
    for i=length(all_data):-1:1
        j = j + 1; % switch array order
        
        % Thickness of this layer
        if (i==1) 
            thickness(j) = all_data(i).dtb;
        else
            thickness(j) = all_data(i).dtb - all_data(i-1).dtb;
        end
        
        % Single crystal elasticity of this layer
        [C, rho] = ice_cij(all_data(i).Tav);
        
        % Calculate poly xtal elasticity of this layer
        % TODO  
        [Cs(:,:,j), rhos(j)] = Cs_from_EBSD_file(C,rho,all_data(i).tex_file);
        
        % Set up inclination and azimuth
        inc=90.0;
        azi=0.0;
        
        [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( Cs(:,:,j), rhos(j), inc, azi );
        fast(j) = MS_unwind_pm_90((azi+pol')) ; % geog. reference frame
        tlag(j) = thickness(j)/vs2 - thickness(j)/vs1 ;
        
    end
    
    % Calculate effective splitting
    freq = 30; % 30 Hz
    spol = 0;
    
    [fast_eff,tlag_eff]=MS_effective_splitting_N(freq,spol,fast,tlag,'mode','GaussianWavelet','PlotWavelet');
    
end


function [C_poly, rho_poly] = Cs_from_EBSD_file(C_single, ...
                                                rho_single, filename)
                                            
    % TODO Read Eulers from file...
    % For now
    nxtls = 2;
    eulers = ones(3,2);
                                            
    Cs = zeros(6,6,nxtls);
    for i = 1:nxtls
        Cs(:,:,i) = MS_rotEuler(C_single, eulers(1,i), eulers(2,i), eulers(3,i));
    end
    rhos = ones(nxtls,1)*rho_single; % All crystals have the same density.
    vfs = ones(nxtls,1);             % Same volume fraction for each point 
                                     % - normalised by MS_VRH.
    [C_poly, rho_poly] = MS_VRH(vfs, Cs, rhos);
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
    data(1) = struct('tex_file', 'some_file_name', ...
                     'dtb', 10, ... % depth to base of layer m
                     'Tav', -10, ... % Average temperature (deg. C)
                     'azi', 0); % Horizonal rotation for texture
    data(2) = struct('tex_file', 'some_file_name', ...
                     'dtb', 100, ... 
                     'Tav', -20, ... 
                     'azi', 10);             
    data(3) = struct('tex_file', 'some_file_name', ...
                     'dtb', 400, ... 
                     'Tav', -1, ... 
                     'azi', 0); 
end