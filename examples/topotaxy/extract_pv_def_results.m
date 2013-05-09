% EXTRACT_PV_DEF_RESULTS.M - cut down version for PV def after transform
% 
% This script is for the analysis of deformed perovskite after
% topotaxy_model has been run on deformed post_perovskite and the 
% VPSC code used to deform the resulting textured perovskite. It is
% included here for completeness.
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

function extract_pv_def_results(varargin)

    % Set Defaults 
    % ============
    vpsc_in_fname = './TEX_PH1.OUT'; % read local copy of TEX_PH1 and 
    write_elastic = 0 ; % don't write elasticity to files.
    pv_elastic_name = '';
    
    % Process any input arguments
    % ===========================
    if ((length(varargin)) == 1)
        vpsc_in_fname = varargin{1};
    elseif ((length(varargin)) == 2)
        vpsc_in_fname = varargin{1};
        write_elastic = 1;
        pv_elastic_name = varargin{2};
    elseif ((length(varargin)) ~= 0)
        error('Only one or two optional arguments.') ;
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
    fprintf('        Results for perovskite\n');
    fprintf('        ======================\n\n');
    fprintf('Reading input (pv) texture from: %s\n', vpsc_in_fname);
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
    % Perovskite and elasticity and density at 2800 K 
    % and 126/7 GPa (Wookey et al. Nature, 2005, 438:1004-1007; 
    % doi:10.1038/nature04345 )';
    [C_pv, rh_pv] = MS_elasticDB('llm_mgsio3pv');
    if (withMTEX)
        % We need the crystal symmetry to do MTEX plots
        CS_pv = symmetry('Pbnm', [pv_a, pv_b, pv_c]);
    end
    
    % Read ppv texture from VPSC file
    % ===============================
    
    [eulers, nxtls] = read_VPSC_file(vpsc_in_fname);
    fprintf('%i sets of Euler angles read from %s\n\n', ...
        [nxtls, vpsc_in_fname]);
    
    
    % Elasticity calculation
    % ======================

    % Elasticity of textured pv
    Cs = zeros(6,6,nxtls);
    for i = 1:nxtls
        Cs(:,:,i) = MS_rotEuler(C_pv, eulers(1,i), eulers(2,i), eulers(3,i));
    end
    rhos = ones(nxtls,1)*rh_pv; % All crystals have the same density.
    vfs = ones(nxtls,1);         % Same volume fraction for each point 
                                 % - normalised by MS_VRH.
    [C_pv_av, rh_pv_av] = MS_VRH(vfs, Cs, rhos);
    report_elasticity(C_pv_av, rh_pv_av, ...
        'perovskite at end of run');
    if write_elastic
        % Save Cij in format for reflectivity calc.
        MS_save(pv_elastic_name, C_pv_av, rh_pv_av, 'Aij', 'eunit', ...
            'Pa');    
    end
    
    % If MTEX is installed draw pole figures
    % ======================================
    if (withMTEX)
        plot_pole_figure(eulers, CS_pv, symmetry('-1'), ...
            'Texture of post-perovskite before transformation', [0, 12])
    end
    
end





%                       +++++++++++++++++++++
%                         SUPPORT FUNCTIONS
%                       +++++++++++++++++++++

function report_elasticity(C, rh, string)
    % Generate text and pole-figure of phase velocities given
    % an elasticity matrix and descriptive string.

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
        'scontours', 0:0.1075:8.6, 'pcontours', 13.0:0.015:14.2, ...
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




% ============================
% Functions for VPSC in Matlab
% ============================
%
% Note that I have copies of these
% in a separate repo to ease reuse.

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

