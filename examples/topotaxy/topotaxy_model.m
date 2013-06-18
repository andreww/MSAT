% TOPOTAXY_MODEL.M - an example script for modelling pv/ppv topotaxy
% 
% This script simulates the effect of a texture preserving phase 
% transition (topotaxy) between MgSiO3 post-perovskite and MgSiO3
% perovskite on increasing temperature or decreasing pressure in 
% the lowermost mantle. Post-perovskite is assumed to deform by 
% simple shear in its stability field and this is modelled using 
% the VPSC code. The output of this modelling (a collection of 
% Euler angles) is the input to this script. The Euler angles
% are first used to generate a model of the elasticity of textured
% post-perovskite just before the phase transition. A rotation matrix
% is built to describe the topotaxy and, together with the Euler angles,
% this is used to generate a model of the elasticity of the perovskite
% following the phase transition. Finally, a set of Euler angles 
% representing the transformed perovskite is generated and these 
% can be used to begin a further VPSC calculation of the deformation 
% of pre-textured perovskite. This is the program used to generate
% the results in Figure 2 of Dobson et al. (2013). Further information 
% can be found in the supplementary information of that paper. 
% 
% The script is packaged as a function taking a number of optional
% arguments to give the source of the VPSC output and restart files
% but these are not needed. The script can most simply be run as 
% topotaxy_model() on the command line. The MTEX toolkit (Bachmann 
% et al. 2010) can be used to produce pole figures of the textures 
% but this is not required. The function can take one, two, or four 
% optional arguments:
%   * The first optional argument in the file name of a non-default 
%     input file to read the post-perovskite texture from. This 
%     defaults  to 'TEX_PH1.OUT' in the local directory. 
%   * The second is a file to be used to write the transformed 
%     perovskite texture in a format that can be read by the VPSC code.
%     If no argument is provided no file is produced. 
%   * The third and fourth arguments allow the post-perovskite
%     and perovskite elasticity to be written in a file for later reuse.
% 
% References:
%     D. P. Dobson, N. Miyajima, F. Nestola, M. Alvaro, N. Casati, 
%        C. Liebske, I. G. Wood and A. M. Walker (2013) "Strong texture 
%        inheritance between perovskite and post-perovskite in the D''
%        layer". Nature Geoscience. 10.1038/NGEO1844
%     F. Bachmann, R. Hielscher and H. Schaeben (2010) "Texture analysis 
%        with MTEX -- free and open source software toolbox". Solid State
%        Phenomena 160:63-68. 10.4028/www.scientific.net/SSP.160.63
%
% See also: MS_rotEuler, MS_rotR, MS_VRH

% (C) Andrew Walker, 2012, 2013
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

function topotaxy_model(varargin)

    % Set Defaults 
    % ============
    vpsc_in_fname = './TEX_PH1.OUT'; % read local copy of TEX_PH1 and 
    write_vpsc_out = 0; % don't generate a restart file.
    vpsc_out_fname = '';
    write_elastic = 0 ; % don't write elasticity to files.
    ppv_elastic_name = '';
    pv_elastic_name = '';
    
    % Process any input arguments
    % ===========================
    if ((length(varargin)) == 1)
        vpsc_in_fname = varargin{1};
    elseif ((length(varargin)) == 2)
        vpsc_in_fname = varargin{1};
        write_vpsc_out = 1;
        vpsc_out_fname = varargin{2};    
    elseif ((length(varargin)) == 4)
        vpsc_in_fname = varargin{1};
        write_vpsc_out = 1;
        vpsc_out_fname = varargin{2};  
        write_elastic = 1;
        ppv_elastic_name = varargin{3};
        pv_elastic_name = varargin{4};
    elseif ((length(varargin)) ~= 0)
        error('Only one, two, or four optional arguments.') ;
    end
        
    % Is MTEX installed 
    % =================
    try
        % check by generating a symmetry op.
        [~] = symmetry('-1');
        withMTEX = 1;
    catch e
        if strcmp(e.identifier, ...
            'MATLAB:UndefinedFunction')
            withMTEX = 0;
        end
    end
    
    % Write out the input arguments for reference
    % ===========================================
    fprintf('\n\n');
    fprintf('        Simulation of topotaxy\n');
    fprintf('        ======================\n\n');
    fprintf('Reading input (ppv) texture from: %s\n', vpsc_in_fname);
    if (write_vpsc_out)
        fprintf('Writing output (pv) texture to: %s\n', vpsc_out_fname);
    else
        fprintf('Not writing output (pv) texture.\n');
    end
    if (withMTEX)
        fprintf('Texture pole figures will be generated using MTEX.\n');
    else
        fprintf('MTEX not found - pole figures not generated.\n');
    end
    fprintf('\n');
    
    % Crystal parameters and elasticity 
    % =================================
    pv_a = 4.318; % Perovskite cell parameters
    pv_b = 4.595; 
    pv_c = 6.305;
    ppv_a = 2.455; % Post-erovskite cell parameters
    ppv_b = 8.061; 
    ppv_c = 6.099;
    % Perovskite and post-perovskite elasticity and density at 2800 K 
    % and 126/7 GPa (Wookey et al. Nature, 2005, 438:1004-1007; 
    % doi:10.1038/nature04345 )';
    [C_pv, rh_pv] = MS_elasticDB('llm_mgsio3pv');
    [C_ppv, rh_ppv] = MS_elasticDB('llm_mgsio3ppv');
    if (withMTEX)
        % We need the crystal symmetry to do MTEX plots
        CS_pv = symmetry('Pbnm', [pv_a, pv_b, pv_c]);
        CS_ppv = symmetry('Cmcm', [ppv_a, ppv_b, ppv_c]);
    end
    
    % Read ppv texture from VPSC file
    % ===============================
    
    [eulers, nxtls] = read_VPSC_file(vpsc_in_fname);
    fprintf('%i sets of Euler angles read from %s\n\n', ...
        [nxtls, vpsc_in_fname]);
    
    % Rotation matricies for transformation
    % =====================================
    [transA, transB] = getTopotaxyRotations(ppv_a, ppv_b, ppv_c);
    fprintf('Rotation matrix for [001]pv||[001]ppv & [100]pv||[-110]ppv:\n');
    fprintf(' %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n\n', ...
        transB);
    fprintf('Rotation matrix for [001]pv||[001]ppv & [100]pv||[110]ppv:\n');
    fprintf(' %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n', ...
        transA);
    
    % Elasticity calculation
    % ======================

    % Elasticity of textured ppv pre-transformeation
    Cs = zeros(6,6,nxtls);
    for i = 1:nxtls
        Cs(:,:,i) = MS_rotEuler(C_ppv, eulers(1,i), eulers(2,i), eulers(3,i));
    end
    rhos = ones(nxtls,1)*rh_ppv; % All crystals have the same density.
    vfs = ones(nxtls,1);         % Same volume fraction for each point 
                                 % - normalised by MS_VRH.
    [C_ppv_av, rh_ppv_av] = MS_VRH(vfs, Cs, rhos);
    report_elasticity(C_ppv_av, rh_ppv_av, ...
        'post-perovskite before transformation');
    if write_elastic
        % Save Cij in format for reflectivity calc.
        MS_save(ppv_elastic_name, C_ppv_av, rh_ppv_av, 'Aij', 'eunit', ...
            'Pa');    
    end
    
    % Elasticity of pv post-transformeation
    chooser = rand([nxtls,1]); % Randomly choose between two rotations
    for i = 1:nxtls
        Cs(:,:,i) = MS_rotEuler(C_pv, eulers(1,i), eulers(2,i), eulers(3,i));
        if (chooser(i) < 0.5)
            Cs(:,:,i) = MS_rotR(Cs(:,:,i),transA);
        else
            Cs(:,:,i) = MS_rotR(Cs(:,:,i),transB);
        end
    end
    rhos = ones(nxtls,1)*rh_pv; % All crystals have the same density.
    vfs = ones(nxtls,1);        % Same volume fraction for each point 
                                % - normalised by MS_VRH.
    [C_pv_av, rh_pv_av] = MS_VRH(vfs, Cs, rhos);
    report_elasticity(C_pv_av, rh_pv_av, ...
        'perovskite after transformation');
    if write_elastic
        % Save Cij in format that can be used for reflectivity calculation.
        MS_save(pv_elastic_name, C_pv_av, rh_pv_av, 'Aij', 'eunit', ...
            'Pa');    
    end
        
    % Texture of transformed material
    % ===============================
    if (write_vpsc_out||withMTEX)
        transformed_eulers = zeros(size(eulers));
        % Only do this if we need to for output or plotting
        for x = 1:nxtls
            if (chooser(x) < 0.5)
                g = rot_from_Euler(eulers(1,x), eulers(2,x), eulers(3,x));
                g = g*transA';
                [phi1, Phi, phi2] = rot_to_Euler(g);
            else
                g = rot_from_Euler(eulers(1,x), eulers(2,x), eulers(3,x));
                g = g*transB';
                [phi1, Phi, phi2] = rot_to_Euler(g);
            end
            transformed_eulers(1,x) = phi1;
            transformed_eulers(2,x) = Phi;
            transformed_eulers(3,x) = phi2;
        end
        
        % A test that this is correct is to see if these Euler angles give
        % C_pv_av...
        Cs = zeros(6,6,nxtls);
        for x = 1:nxtls
            Cs(:,:,x) = MS_rotEuler(C_pv, transformed_eulers(1,x),...
                transformed_eulers(2,x), transformed_eulers(3,x));
        end
        rhos = ones(nxtls,1)*rh_pv; % All crystals have the same density.
        vfs = ones(nxtls,1);         % Same volume fraction for each point 
                                     % - normalised by MS_VRH.
        [C_pv_av_chk, rh_pv_av_chk] = MS_VRH(vfs, Cs, rhos);
        try
            assertElementsAlmostEqual(C_pv_av_chk, C_pv_av, ...
                ['Elasticity from rotated tensor does not match ', ...
                'elasticity from rotated Euler angles!']);
            assertElementsAlmostEqual(rh_pv_av_chk, rh_pv_av, ...
                'Densities do not match.');
        catch e
            if strcmp(e.identifier, ...
                'MATLAB:UndefinedFunction')
                %XUnit not installed...
                fprintf(['\n Warning: XUnit is not installed. \n ', ...
                     'Proceeding with out checking for validity of ', ...
                     'rotated Euler angles\n']);
            else
                % Anything else is a problem.
                rethrow(e);
            end
        end
        fprintf('\n%i sets of transformed Euler angles calculated\n\n',...
            nxtls);
    end
    

    % Write out the transformed texture
    % =================================
    if (write_vpsc_out)
        write_VPSC_file(vpsc_out_fname, transformed_eulers, ...
            'Euler angles for perovskite transformed from post-perovskite');
        fprintf('\n%i sets of transformed Euler angles written to %s\n',...
        [nxtls, vpsc_out_fname]);
    end
    
    % If MTEX is installed draw pole figures
    % ======================================
    if (withMTEX)
        plot_pole_figure(eulers, CS_ppv, symmetry('-1'), ...
            'Texture of post-perovskite before transformation', [0, 12])
        plot_pole_figure(transformed_eulers, CS_pv, symmetry('-1'), ...
            'Texture of perovskite after transformation', [0, 12])
    end
    
end





%                       +++++++++++++++++++++
%                         SUPPORT FUNCTIONS
%                       +++++++++++++++++++++

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
        'avscontours', 0:0.1075:8.6, 'pcontours', 13.0:0.015:14.2, ...
        'polsize', 0.18, 0.16, 2.0, 1.0, 'limitsonpol');
    
    [ xi, phi ] = enforce_VTI_rot(C);
    fprintf('Parameters for imposed VTI symmetry: \n')
    fprintf('    xi  (Vsh^2/Vsv^2) = %7.4f     ln(xi)  = %7.2f %% \n', ...
        xi, log(xi)*100.0);
    fprintf('    phi (Vpv^2/Vph^2) = %7.4f     ln(phi) = %7.2f %% \n', ...
        phi, log(phi)*100.0);


end 


function [ xi, phi ] = enforce_VTI_rot(C)
    % Create an elasticity tensor with VTI symmetry by rotating a 
    % general 6*6 matrix around the X3 axis in 360 one degree steps
    % and calculating the VRH tensor from each rotation. Outputs 
    % the Xi (Vsh^2/Vsv^2) and Phi (Vpv^2/Vph^2) used in global 
    % anisotropic tomography.

    Cs = zeros(6,6,360);
    Vf = ones(360,1);
    rh = Vf*1000;

    for i = 1:360
        Cs(:,:,i) = MS_rot3(C, 0, 0, i);
    end
    C_VTI = MS_VRH(Vf, Cs, rh);
    [ pol, ~, vs1, vs2, vp_h] = MS_phasevels( C_VTI, 1000, 0.0, 20.0 );
    [ ~, ~, ~, ~, vp_v] = MS_phasevels( C_VTI, 1000, 90.0, 0.0 );
    if (abs(pol) < sqrt(eps)) 
        vs_v = vs1;
        vs_h = vs2;
    elseif ((abs(pol)-90.0) < sqrt(eps))
        vs_h = vs1;
        vs_v = vs2;
    else
       error('not VTI')
    end
    xi = (vs_h^2.0) / (vs_v^2.0);
    phi = (vp_v^2.0) / (vp_h^2.0);
    
end



function [transA, transB] = getTopotaxyRotations(ppv_a, ppv_b, ppv_c)
    % Return two rotation matrices to give the two possible orientations
    % of perovskite given a post-perovskite orientation.

    % We need the post-perovskite lattice vectors 
    ppv_lattice_a = [ppv_a 0.0   0.0  ];
    ppv_lattice_b = [0.0   ppv_b 0.0  ];
    ppv_lattice_c = [0.0   0.0   ppv_c];                                         
                                             
    % TOPOTAXY RULES
    % ==============

    % Easer to express for ppv -> pv. To do pv -> ppv just need to invert
    % (or transpose) trans matrix. Rules are [001]_pv||[001]_ppv, 
    % [100]_pv||[110]_ppv or [100]_pv||[-110]_ppv (twin). 

    % Get unit vectors pointing along the PPV axis frame
    ppv_X1 = ppv_lattice_a/norm(ppv_lattice_a);
    ppv_X2 = ppv_lattice_b/norm(ppv_lattice_b);
    ppv_X3 = ppv_lattice_c/norm(ppv_lattice_c);

    % Topotaxy rules...
    % c_pv = c_ppv
    pv_X3 = ppv_X3;
    % a_pv = -110_ppv
    absolute_pv_a = (ppv_lattice_b - ppv_lattice_a);
    pv_X1a = absolute_pv_a/norm(absolute_pv_a);
    % OR a_pv = 110_ppv
    absolute_pv_a = (ppv_lattice_b + ppv_lattice_a);
    pv_X1b = absolute_pv_a/norm(absolute_pv_a);
    % Third direction fixed by orthoganlity and RH frame...
    pv_X2a = cross(pv_X3,pv_X1a)/norm(cross(pv_X3,pv_X1a));
    pv_X2b = cross(pv_X3,pv_X1b)/norm(cross(pv_X3,pv_X1b));
    % Build the rotation matrix...
    transA = [dot(ppv_X1,pv_X1a) dot(ppv_X2,pv_X1a) dot(ppv_X3,pv_X1a); ...
              dot(ppv_X1,pv_X2a) dot(ppv_X2,pv_X2a) dot(ppv_X3,pv_X2a); ...
              dot(ppv_X1,pv_X3)  dot(ppv_X2,pv_X3)  dot(ppv_X3,pv_X3) ];
    transB = [dot(ppv_X1,pv_X1b) dot(ppv_X2,pv_X1b) dot(ppv_X3,pv_X1b); ...
              dot(ppv_X1,pv_X2b) dot(ppv_X2,pv_X2b) dot(ppv_X3,pv_X2b); ...
              dot(ppv_X1,pv_X3)  dot(ppv_X2,pv_X3)  dot(ppv_X3,pv_X3) ];

    % Check results
    % =============
    try
        % Check we have a proper rotation matrix (and not an impropper rotation
        % matrix) - need XUnit installed
        assertElementsAlmostEqual(transA', inv(transA), ...
            'TransA is not a rotation matrix: R^T != R^-1');
        assertElementsAlmostEqual(det(transA), 1.0, ...
            'TransA is not a rotation matrix: det(R) != 1.0');
        assertElementsAlmostEqual(transB', inv(transB), ...
            'TransB is not a rotation matrix: R^T != R^-1');
        assertElementsAlmostEqual(det(transB), 1.0, ...
            'TransB is not a rotation matrix: det(R) != 1.0');
    catch e
        if strcmp(e.identifier, ...
            'MATLAB:UndefinedFunction')
            %XUnit not installed...
            fprintf(['\n Warning: it seems XUnit is not installed. \n ', ...
                     'Proceeding with out checking for valididty of ', ...
                     'rotation matrices\n']);
        else
            % Anything else is a problem.
            rethrow(e);
        end
    end
    
end


% ============================
% Functions for VPSC in Matlab
% ============================
%
% Note that I have copies of these
% in a separate repo to ease reuse.

function [ g ] = rot_from_Euler(phi1, theta, phi2)
    % Given three Euler angles phi1, Phi, phi2 (Bunge notation, 
    % in degrees) return a rotation matrix, g, representing the 
    % rotation. This is the 'passive' rotation.
    
    % Pre-compute trig functions.
    cp1 = cos(phi1*pi/180.0);
    sp1 = sin(phi1*pi/180.0);
    cp2 = cos(phi2*pi/180.0);
    sp2 = sin(phi2*pi/180.0);
    cth = cos(theta*pi/180.0);
    sth = sin(theta*pi/180.0);

    % Form rotation matrix for Bunge convention of 
    % Euler angles, see eq 24 (pg81) of "Preferred
    % orientation in deformed meals and rocks... H-K Wenk (ed)"
    g = [cp1*cp2 - sp1*sp2*cth, sp1*cp2 + cp1*sp2*cth, sp2*sth ; ...
        -1*cp1*sp2 - sp1*cp2*cth, -1*sp1*sp2 + cp1*cp2*cth, cp2*sth; ...
        sp1*sth, -1*cp1*sth, cth];
    
   
end


function [phi1, Phi, phi2] = rot_to_Euler(R)
    % Given a rotation matrix, R, return the three 
    % Euler angles phi1, Phi, phi2 (Bunge notation, in degrees) 
    % representing the rotation. 

    if (R(3,3)>1.0-eps)
        Phi = 0.0;
        phi1 = atan2(R(1,2),R(1,1));
        phi2 = 0.0;
    else
        Phi = acos(R(3,3));
        phi1 = atan2( (R(3,1)/sin(Phi)) , (-1.0*(R(3,2)/sin(Phi))) );
        phi2 = atan2( (R(1,3)/sin(Phi)) , (R(2,3)/sin(Phi)) );
    end
    
    if (phi1 < 0.0) 
        phi1 = 2*pi + phi1;
    end
    if (phi2 < 0.0)
        phi2 = 2*pi + phi2;
    end
    
    phi1 = (180.0/pi) * phi1;
    Phi =  (180.0/pi) * Phi;
    phi2 = (180.0/pi) * phi2;

end


function plot_pole_figure(eulers, CS, SS, title, scale)
    % Use MTEX to plot a pole figure of a set of Euler 
    % angles. Euler angles are in Bunge notation and 
    % degrees. User must pass in cystal symmetry (CS)
    % and sample symmetry (SS) generated using MTEX symmetry
    % constructor along with a title (a string for the plot
    % header) and an array of two values - the mimimum and 
    % maximum value for the contour scheme.
   
    % Work in radians.
    eulers_r = eulers*degree;
    
    % Generate a list of orientation objects and an ebsd wrapper
    o = orientation('Euler', eulers_r(1,:), eulers_r(2,:), ...
        eulers_r(3,:), CS,SS, 'Bunge');
    ebsd = EBSD(o);
    
    % Fit an ODF
    odf = calcODF(ebsd,'HALFWIDTH', 10*degree, 'silent');
    
    % Make a new figure and plot three pole figures in it...
    h = figure;
    plotpdf(odf,[Miller(1,0,0),Miller(0,1,0),Miller(0,0,1)], ...
        'antipodal','silent','position',[100 100 1200 400], 'colorrange', ...
        scale);
    colorbar;
    set(h, 'Name', title);
    
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


function write_VPSC_file(filename, eulers, title)
    % Create a VPSC input or output texture file called <filename>.
    % Texture is described by a size(3,n) array of Euler angles
    % in <eulers> and these should be in degrees and in Bunge notation.
    % The string <title> is inserted into the header. At present, I assume
    % all crystals have the same volume fraction.

    % How many crystals
    n_xtals = size(eulers,2);
    % Open file for write (text mode)
    fid = fopen(filename,'wt');
    
    % Write the header lines
    fprintf(fid,'%s\n', title);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'B  %5i\n', n_xtals);
    
    % Write the data - assume that volume fraction is the same for
    % each crystal
    for i = 1:n_xtals
        fprintf(fid,'%6.2f  %6.2f  %6.2f  %12.4f\n',...
            eulers(1,i), eulers(2,i), eulers(3,i), (1/n_xtals));
    end
    
    % Close file
    fclose(fid);  
    
end
