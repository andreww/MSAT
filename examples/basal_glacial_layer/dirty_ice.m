% DIRTY_ICE.M example to generate a ATRAK model
%
% This script create a ATRAK input file (to be used with "modsimple"
% that is supposed to represent a layer of "dirty ice" at the base of
% a glacier. This includes rock particles described as isotropic 
% elipoidal inclusions in an isotopic ice matrix. These inclusions are
% alligned, so the basal layer is elastically anisotropic. The rest of the
% glacier is considered to be isotropic, and we do not consider depth 
% dependent wave speeds.
% 
% See also: MS_effective_medium

% (C) Andrew Walker, University of Leeds 2014
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

function dirty_ice 

% Input for model - these are the things that will need to change 
% between different ATRAK runs.

% First we need the seismic properties of the ice. Note
% that the ice is isotropic, so we only need two velocities
% and a density.
Vp_ice = 3.9; % km/s
Vs_ice = 1.8; % km/s
density_ice = 900; % kg/m^3 

% Also the seismic properties of the rock fragments, also isotropic.
Vp_rock = 5.0; % km/s
Vs_rock = 3.0; % km/s
density_rock = 2600 ; % kg/m^3

% Define the number of rock particles (as a volume fraction) and
% their shape (as an aspect ratio). The volume fraction is a number
% between 0 (no rock) and 1 (100% rock) and should be fairly small
% (less than 0.5) or the method will fail. The aspect ratio is the 
% ratio of the length of unique axis to the length of the other two 
% axes (which are the same, the particles have rotational symmetry)
% so values less than 1 are like smarties and values more than 1 are
% like cigars. An aspect ratio of 1 gives a sphere, and an isotropic
% elastic model.
volume_fraction_rock = 0.1;
rock_aspect_ratio = 1.1 ;

% Geometrical properties - these are just used to set up the 
% modsimple input file.
depth_to_interface = 1000 ; % m

% Now calculate the (isotropic) elastic constants that represent the 
% ice and the rock (MSAT used GPa for all elasticity).
[C_ice] = MS_iso(Vp_ice, Vs_ice, density_ice);
[C_rock] = MS_iso(Vp_rock, Vs_rock, density_rock);

% Now we can claculate the elasticity and density of the dirty ice layer.
[C_dirty_ice,rho_dirty_ice] = MS_effective_medium('tandon', ...
    C_ice, density_ice, C_rock, density_rock,  ...
    rock_aspect_ratio, volume_fraction_rock);

% Check we have a sensible value - this will throw an error if C_dirty_ice
% is not physical 
MS_checkC(C_dirty_ice);

% Write out the result for the dirty ice layer (for the record)
report_elasticity(C_dirty_ice, rho_dirty_ice, 'dirty ice model')


end






%                       +++++++++++++++++++++
%                         SUPPORT FUNCTIONS
%                       +++++++++++++++++++++


function create_modsimple_rsp(filename, dtts, dtbs, Cs, rhos)
    % Create a modsimple input file with two layers each with a 
    % thickness, elasticity and density. Write to file "filename"

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
    fprintf(fid, '%d %f %f\n', 10, 0.0, 10000.0);
    fprintf(fid, '%d %f %f\n', 10, -5000.0, -5000.0);
    fprintf(fid, '%d %f %f\n', 10, 0.0, -dtbs(num_layers)); % Base of hole
    
    % Now add each layer in turn (remembering we need to go backwards)
    j = 0;
    for i = num_layers:-1:1
        j = j + 1;
        fprintf(fid, '%d %s\n', j, 'flat');
        fprintf(fid, '%f %f\n', -dtts(i), -dtts(i));
        fprintf(fid, '%s\n', 'ani');
        fprintf(fid, '%d\n', 0); 
        fprintf(fid, '%d\n', 21);
        MS_save(fid, Cs(:,:,i), rhos(i), 'eunit', 'mbar', 'dunit', ...
            'gcc', 'Aij');
    end
    fclose(fid);
end

function report_elasticity(C, rh, string)
    % Generate text and pole-figure of phase velocities given
    % an elasticity matrix and descriptive string.
    % This also shows how to make "nice" pole figures.

    MS_checkC(C);
    
    fprintf('\nInformation for %s:\n\n', string);
    fprintf('Elasticity matrix (Voigt notation, GPa):\n'); 
    fprintf([' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n', ...
             ' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n', ...
             ' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n', ...
             ' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n', ...
             ' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n',...
             ' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n'], ...
        C);
    fprintf('Density: %7.2f kg/m^3\n\n', rh);
    MS_plot(C, rh, 'wtitle', string, 'fontsize', 11, ...
        'polsize', 0.18, 0.16, 2.0, 1.0, 'limitsonpol');
end 