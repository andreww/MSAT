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

function ice_lpo 

%                       +++++++++++++++
%                         MODEL SETUP
%                       +++++++++++++++

% Input for model - these are the things that will need to change 
% between different ATRAK runs.

% First we need the seismic properties of the ice. Note
% that the ice is isotropic, so we only need two velocities
% and a density.
Vp_ice = 3.9; % km/s
Vs_ice = 1.8; % km/s
density_ice = 900; % kg/m^3 

% What is the file name for the texture file
filename = 'TEX_PH1_compress.OUT';

% Geometrical properties - these are just used to set up the 
% modsimple input file and get the ice texture in the right 
% orientation.

% First we need to rotate the elasticity into the 
% correct frame of reference. ATRAK works with X3 vertical, 
% we have the compression direction X1. Let's 
% assume that flow is in the X1 direction. If we have cigars, these
% will be alligned with flow, so we don't need to rotate
% (set rot_X2 = 0 and rot_X3 = 0). If we have
% smarties, the short direction could be vertical (in which case set
% rot_X2 = 90 and rot_X3 = 0) or perpendicular to flow (in which case set
% rot_X2 = 0 and rot_X3 = 90). It may be worth drawing a figure of this!
rot_X2 = 90.0; % Degrees
rot_X3 = 0.0; % Degrees. 

% This is just for modsimple
depth_to_interface = 1000 ; % m
depth_to_base = 1100;
modsimple_filename = 'modsimple.rsp';


%                       +++++++++++++++
%                         DO THE WORK
%                       +++++++++++++++

% Now calculate the (isotropic) elastic constants that represent the 
% ice and the rock (MSAT used GPa for all elasticity).
[C_ice] = MS_iso(Vp_ice, Vs_ice, density_ice);


% Calculation of the elasticity of textured ice
% First we need the elasticity of single crystal ice

% which we have in a database...
[C_ice_xtl, ~] = MS_elasticDB('ice');

% Now read the orientations of all ice crystals from the inpout 
% file, and build an array of elasticity tensors corresponding to each
% grain at the relevent orentation
[eulers, nxtls] = read_VPSC_file(filename);
Cs = zeros(6,6,nxtls);
for i = 1:nxtls
  Cs(:,:,i) = MS_rotEuler(C_ice_xtl, eulers(1,i), eulers(2,i), eulers(3,i));
end
rhos = ones(nxtls,1)*density_ice; % All crystals have the same density.
vfs = ones(nxtls,1);         % Same volume fraction for each point 
                             % - normalised by MS_VRH.
[C_ice_av, ~] = MS_VRH(vfs, Cs, rhos);

% Rotate ice into the correct orientation *note that I'm assuming 
% only 1 non zero rotation as rotations do not commute*.
C_ice_av = MS_rot3(C_ice_av, 0.0, rot_X2, rot_X3);

% Check we have a sensible value - this will throw an error if C_dirty_ice
% is not physical 
MS_checkC(C_ice_av);

% Write out the result for the dirty ice layer (for the record)
report_elasticity(C_ice_av, density_ice, 'textured ice model')

% Now generate a modsimple file. We are reusing old code here - which
% handles many layers - hence the arrays.
layer_Cs = zeros(6,6,2);
layer_Cs(:,:,1) = C_ice;
layer_Cs(:,:,2) = C_ice_av;
create_modsimple_rsp(modsimple_filename, [0 -depth_to_interface], ...
    [-depth_to_interface -depth_to_base], layer_Cs, ...
    [density_ice density_ice]);

end






%                       +++++++++++++++++++++
%                         SUPPORT FUNCTIONS
%                       +++++++++++++++++++++


function create_modsimple_rsp(filename, dtts, dtbs, Cs, rhos)
    % Create a modsimple input file with two layers each with a 
    % thickness, elasticity and density. Write to file "filename"
    % taking care of odd unit conversion.

    model_name = 'icesheet.mod';
    length_Cs = size(Cs);
    length_Cs = length_Cs(3);
    assert(length(dtts)==length(dtbs), 'Tops and bottoms must agree')
    assert(length(dtts)==length(rhos), 'Tops and density must agree')
    assert(length(dtts)==length_Cs, 'Tops and Cijs must agree')

    num_layers = length(dtts);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', model_name);
    fprintf(fid, '%d\n', num_layers);

    % Hard code X, Y and Z nodes here. We could make this
    % an argument, but I don't know what the rules are.
    fprintf(fid, '%d %f %f\n', 10, -10000.0, 10000.0);
    fprintf(fid, '%d %f %f\n', 10, -5000.0, 5000.0);
    fprintf(fid, '%d %f %f\n', 10, -1100, -dtbs(num_layers)); % Base of hole
    
    % Now add each layer in turn
    j = 0;
    for i = 1:num_layers
        fprintf(fid, '%d %s\n', i, 'flat');
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

function [eulers, nxtl] = read_VPSC_file(filename)
    % Reads data from a VPSC output or input texture file.
    % Must be Bunge convention and in degrees. Returns eulers, a (3,nxtl) 
    % array of Euler angles (Bunge convention, (1,:) = phi1, (2,:) = Phi
    % and (3,:) = phi2) and a scalar nxtl, specifying the number of 
    % crystals / Euler angle triples.

    % Read data from the file
    fid = fopen(filename); % Read - the default
    fgetl(fid); % Header line - ignore
    fgetl(fid); % Lengths of phase ellipsoid axes - ignore
    fgetl(fid); % Euler angles for phase ellipsoid - ignore
    L = sscanf(fgetl(fid), '%s %d'); % Convention and number of crystals
    E = fscanf(fid, '%f');
    fclose(fid);
    
    % Get hold of header info
    assert((char(L(1))=='B'), ... % Check Euler angle convention
        'Could not read VPSC file - not Bunge format\n');
    nxtl = L(2); % Number of crystals
    
    % Build Euler angles array.
    eulers = zeros(3,nxtl);
    eulers(1,:) = E(1:4:4*nxtl);
    eulers(2,:) = E(2:4:4*nxtl);
    eulers(3,:) = E(3:4:4*nxtl);
   
end
