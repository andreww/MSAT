% INTERPOLATE_2D_EXAMPLE.m - example script showing two ways elasticity 
%                         matricies can be interpolated in MSAT.
%
% Interpolation of elastic constants can be a difficult problem. In this
% example we show how two approaches, both implemented in MSAT, give quite 
% different results. We imagine a region where the elasticity changes 
% smoothly from one described by single crystal olivine rotated 60 degrees 
% clockwise around its c-axis to one described by single crystal enstatie 
% rotated 120 degrees clockwise around its c-axis and seek to calculate 
% the elasticity of a points between these limits. One (commonly 
% used) approach is to take the average of each element of the two 
% elasticity tensors or the average of each element of the two compliance 
% tensors. This amounts to the Voigt, Reuss or Voigt-Reuss-Hill average 
% available via MS_VRH. If the interpolation is between identical rotated 
% tensors this approach is clearly sub-optimal and an alternative is 
% available in via MS_interpolate. This function by finding a common 
% orientation for the two tensors, using MS_axes, averaging the rotated 
% tensors, then moving the rotated tensor into the correct orientation. 
% This example script compares the two approaches.
%
% FIXME - sort out optional arguments below.
% This example function takes a handfull of optional arguments intended to
% allow robust timing of the two approaches. Turning off the generation of
% phase velocity plots and summary information ('no_velocities' and 
% 'no_summary') means that almost all time in this function is spent in the
% MS_VRH or MS_interpolate. Setting 'calc_times' gathers timing information
% for the avarage call (and implies 'no_velocities' and 'no_summary').
% MS_VRH based interpolation is typically ~200 times faster than
% MS_interpolate (0.0001 s per call versus 0.0023 s per call, on one
% machine).
% 
% See also: INTERPOLATE_1D_EXAMPLE MS_VRH, MS_INTERPOLATE, MS_AXES

% (C) James Wookey and Andrew Walker, 2015
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

function interpolate_2D_example(varargin)

    % Default options
    phase_vels = 1; % Plot phase velocities
    gather_results = 1; % Calculate and plot summary data
    report_time = 0; % Just report per-interpolation timeings

    % Process the optional arguments
    iarg = 1 ;
    while iarg<=(length(varargin))
       switch varargin{iarg}
          case 'no_velocities'
             phase_vels = 0 ;
             iarg = iarg + 1 ;
          case 'no_summary'
             gather_results = 0 ;
             iarg = iarg + 1 ;
          case 'calc_times'
             gather_results = 0 ;
             phase_vels = 0 ;
             report_time = 1 ;
             iarg = iarg + 1 ;
          otherwise 
             error('Unknown flag') ;   
       end
    end
    
    
    % Read full data set from textures and convert to elasticity.
    [C_single, rh] = MS_elasticDB('llm_mgsio3ppv');
    [Cs_full, X_full, Y_full] = setup_full_elasticity_grid(C_single, ...
        'vpsc_textures.tex');
    
    
    % Setup low resolution model
    [X_subset, Y_subset] = meshgrid(0:5:20, 0:5:20);
    [ni, nj] = size(X_subset);
    C_subset = zeros(6, 6, ni, nj);
    for i = 1:ni
        for j = 1:nj
            C_subset(:,:,i,j) = Cs_full(:,:,i*4,j*4);
        end
    end
    
    
    % Do the interpolation
    [ni_full, nj_full] = size(X_full);
    Cs_interp = zeros(6,6,ni_full,nj_full);
    Cs_voigt = zeros(6,6,ni_full,nj_full);
    r_interp = zeros(ni_full, nj_full);
    r_voigt = zeros(ni_full, nj_full);
    
    Vp_full = zeros(ni_full, nj_full);
    Vp_interp = zeros(ni_full, nj_full);
    Vp_voigt = zeros(ni_full, nj_full);
    P_ort_full = zeros(ni_full, nj_full);
    P_ort_interp = zeros(ni_full, nj_full);
    P_ort_voigt = zeros(ni_full, nj_full);
    
    inc = 0.0;
    azi = 0.0;
    
    for i = 1:ni_full
        for j = 1:nj_full
            
            x = X_full(i,j);
            y = Y_full(i,j);
            [dx, dy, illcnr, jllcnr] = point_finder(x, y, ...
                X_subset, Y_subset);
            
            [Cs_interp(:,:,i,j), r_interp(i,j)] = bilinear_interp(...
                                     C_subset(:,:,illcnr,jllcnr), rh, ...
                                     C_subset(:,:,illcnr,jllcnr+1), rh, ...
                                     C_subset(:,:,illcnr+1,jllcnr), rh, ...
                                     C_subset(:,:,illcnr+1,jllcnr+1), rh, ...
                                     dx/4, dy/4, 0);
                                 
            [Cs_voigt(:,:,i,j), r_voigt(i,j)] = bilinear_VRH(...
                                     C_subset(:,:,illcnr,jllcnr), rh, ...
                                     C_subset(:,:,illcnr,jllcnr+1), rh, ...
                                     C_subset(:,:,illcnr+1,jllcnr), rh, ...
                                     C_subset(:,:,illcnr+1,jllcnr+1), rh, ...
                                     dx/4, dy/4);
           
            [~, ~, ~, ~, Vp_full(i,j)] = ...
                MS_phasevels( Cs_full(:,:,i,j), rh, inc, azi );                   
            [~, ~, ~, ~, Vp_interp(i,j)] = ...
                MS_phasevels( Cs_interp(:,:,i,j), r_interp(i,j), inc, azi );
            [~, ~, ~, ~, Vp_voigt(i,j)] = ...
                MS_phasevels( Cs_voigt(:,:,i,j), r_voigt(i,j), inc, azi );
            
            [P_ort_full(i,j), ~,~] = summary_anisotropy(Cs_full(:,:,i,j));
            [P_ort_interp(i,j), ~,~] = summary_anisotropy(Cs_interp(:,:,i,j));
            [P_ort_voigt(i,j), ~,~] = summary_anisotropy(Cs_voigt(:,:,i,j));
        end
    end
    
    figure
    subplot(3, 1, 1);
    contourf(X_full, Y_full, Vp_full);
    colorbar;
    title('Full');
    subplot(3, 1, 2);
    contourf(X_full, Y_full, Vp_interp);
    colorbar;
    title('Interp');
    subplot(3, 1, 3);
    contourf(X_full, Y_full, Vp_voigt);
    colorbar;
    title('Voigt');
    
    figure
    subplot(3, 1, 1);
    contourf(X_full, Y_full, P_ort_full);
    colorbar;
    title('Full');
    subplot(3, 1, 2);
    contourf(X_full, Y_full, P_ort_interp);
    colorbar;
    title('Interpt')
    subplot(3, 1, 3);
    contourf(X_full, Y_full, P_ort_voigt);
    colorbar;
    title('Voigt');
    
        

%     
%     if report_time
%         fprintf('Time per Voigt interpolation = %f s\n', time_per_voigt)
%         fprintf('Time per common orientation interpolation = %f s\n', time_per_interp)
%     end
%     
    
end

function [Cxy, rxy] = bilinear_interp(C00, r00, C10, r10, C01, r01, C11,...
                                      r11, x, y, order)
     % Bilinear interpolation of elastic constants and density at point x,
     % y (between 0 and 1) given values at (x, y) of (0, 0), (1, 0), (0, 1)
     % and (1, 1) the symmetry aware interpolator.
     
     if order == 0
         % Interpolate along x at y=0
         [Cy0, ry0] = MS_interpolate(C00, r00, C10, r10, x);
         % Interpolate along x at y=1
         [Cy1, ry1] = MS_interpolate(C01, r01, C11, r11, x);
         % interpolate along y at x=x
         [Cxy, rxy] = MS_interpolate(Cy0, ry0, Cy1, ry1, y);
     else
         % Interpolate along y at x=0
         [Cx0, rx0] = MS_interpolate(C00, r00, C01, r01, y);
         % Interpolate along y at x=1
         [Cx1, rx1] = MS_interpolate(C10, r10, C11, r11, y);
         % interpolate along y at x=x
         [Cxy, rxy] = MS_interpolate(Cx0, rx0, Cx1, rx1, x);
     end
     
end


function [Cxy, rxy] = bilinear_VRH(C00, r00, C10, r10, C01, r01, C11,... 
                                      r11, x, y)
     % Bilinear interpolation of elastic constants and density at point x,
     % y (between 0 and 1) given values at (x, y) of (0, 0), (1, 0), (0, 1)
     % and (1, 1) the symmetry aware interpolator.
     
     % Interpolate along x at y=0
     [~, ry0, Cy0] = MS_VRH([x 1-x], C00, r00, C10, r10);
     % Interpolate along x at y=1
     [~, ry1, Cy1] = MS_VRH([x 1-x], C01, r01, C11, r11);
     % interpolate along y at x=x
     [~, rxy, Cxy] = MS_VRH([y 1-y], Cy0, ry0, Cy1, ry1); 
     
end


function [dx, dy, illcnr, jllcnr] = point_finder(x, y, Xgrid, Ygrid)
    % Given the point (x,y) and a grid of points Xgrid and Ygrid
    % which could have been produced by meshgrid, find the indicies of the
    % grids representing the 'lower left' corner of the grid - i.e. the 
    % point x,y is in the box defined by (Xgrid(i,j), Ygrid(i,j)), 
    % (Xgrid(i,j+1), Ygrid(i,j+1)), (Xgrid(i+1,j), Ygrid(i+1,j)),
    % and (Xgrid(i+1,j+1), Ygrid(i+1,j+1)). The distance from the corner
    % to the point is given by dx and dy.
    
    [ni, nj] = size(Xgrid);
    % Check most of the gird (brute force)
    for i = 1:ni-1
        for j = 1:nj-1
            xllcnr = Xgrid(i,j);
            yllcnr = Ygrid(i,j);
            xurcnr = Xgrid(i+1,j+1);
            yurcnr = Ygrid(i+1,j+1);
            if (x >= xllcnr) && (y >= yllcnr) ...
                    && (x < xurcnr) && (y < yurcnr)
                dx = x - xllcnr;
                dy = y - yllcnr;
                illcnr = i;
                jllcnr = j;
                return
            end        
        end
    end

    % Check the upper and right hand edges
    i = ni;
    for j = 1:nj-1
        xllcnr = Xgrid(i,j);
        yllcnr = Ygrid(i,j);
        yurcnr = Ygrid(i,j+1);
        if (x == xllcnr) && (y >= yllcnr) && (y < yurcnr)
            dx = x - xllcnr;
            dy = y - yllcnr;
            illcnr = i;
            jllcnr = j;
            return
        end        
    end
    j = nj;
    for i = 1:ni-1
        xllcnr = Xgrid(i,j);
        yllcnr = Ygrid(i,j);
        xurcnr = Xgrid(i+1,j);
        if (x >= xllcnr) && (y == yllcnr) && (x < xurcnr)
            dx = x - xllcnr;
            dy = y - yllcnr;
            illcnr = i;
            jllcnr = j;
            return
        end        
    end
    xllcnr = Xgrid(ni,nj);
    yllcnr = Ygrid(ni,nj);
    if (x == xllcnr) && (y == yllcnr)
        dx = 0;
        dy = 0;
        illcnr = ni;
        jllcnr = nj; 
        return
    end
    
    x
    y
    Xgrid
    Ygrid
    error('failed to find point in grid')

end


function [P_ort, P_low, Ua] = summary_anisotropy(C)

         % Calculate summary information about the anisotropy 
         % of an elasticity matrix for comparison. This returns
         % the fraction of the anisotropic part of the tensor that
         % is consistent with orthorhombic symmetry (P_ort) and the
         % anisotropic part that is of lower symmetry (P_low). We
         % also return the general measure of anisotropy propsed by
         % Ledbetter and Miglion (2006) - which is based on the maximum
         % and minimum S-wave velocities.

         CR = MS_axes(C) ;
         [C_iso,C_hex,C_tet,C_ort,C_mon,C_tri] = MS_decomp(CR);
         P = MS_norms(CR,C_iso,C_hex,C_tet,C_ort,C_mon,C_tri) ;
         P_ort = sum(P(1:4))./sum(P(1:6));
         P_low = sum(P(5:6))./sum(P(1:6));
         [Ua] = MS_anisotropy(C);

end


function [Cs, X, Y] = setup_full_elasticity_grid(C_single, filename)

    % Read Euler angles for each strain increment - these are stored
    % in cell arrays
    [eulers, nxtl] = read_VPSC_file(filename);
    num_blocks = length(nxtl) - 1; % We don't want the last block - which
                                   % is the last block repeted
    Cs_in = zeros([6 6 num_blocks]);
    for i = 1:num_blocks
        Cs_in(:, :, i) = elasticity_from_eulers(C_single, nxtl(i), ...
            eulers{i});
    end
    
    assert(num_blocks == 400, 'Grid assumes 400 strains');
    
    point_spacing = 1.0;
    Cs = zeros(6, 6, 20, 20);
    X = zeros(20, 20);
    Y = zeros(20, 20);
    % i is index along x (horizontal) and j is index aling y (vertical, 
    % increses downwards, towards CMB).
    for j = 1:20 
        for i = 1:20
          Cs(:,:,i,j) = Cs_in(:,:, i*j);
          X(i,j) = (i*point_spacing)-point_spacing;
          Y(i,j) = (j*point_spacing)-point_spacing;
        end
    end
    
end

function [C_poly] = elasticity_from_eulers(C_single, nxtls, eulers)
    % Elasticity of textured pv
    Cs = zeros(6,6,nxtls);
    for i = 1:nxtls
        Cs(:,:,i) = MS_rotEuler(C_single, eulers(1,i), eulers(2,i), ...
            eulers(3,i));
    end
    rhos = ones(nxtls,1); % All crystals have the same density.
                          % we end up ignoring this.
    vfs = ones(nxtls,1);  % Same volume fraction for each point 
                          % - normalised by MS_VRH.
    [C_poly, ~] = MS_VRH(vfs, Cs, rhos);

end


function [eulers, nxtl] = read_VPSC_file(filename)
    % Reads data from a VPSC output or input texture file.
    % Must be Bunge convention and in degrees. Returns eulers, a (3,nxtl) 
    % array of Euler angles (Bunge convention, (1,:) = phi1, (2,:) = Phi
    % and (3,:) = phi2) and a scalar nxtl, specifying the number of 
    % crystals / Euler angle triples.

    % Read data from the file
    fid = fopen(filename); % Read - the default
    
    blocks = 0; % Which number texture block are we on?

    while (true) % infinite loop (break with conditional)
        
        tmp_l = fgetl(fid);  % First header line - ignore...
        if (tmp_l == -1)     % if no strain info then we are at the end 
            break            % of the file, break out of read loop
        end
        
        fgetl(fid); % Lengths of phase ellipsoid axes - ignore
        fgetl(fid); % Euler angles for phase ellipsoid - ignore
        L = sscanf(fgetl(fid), '%s %d'); % Convention and number of crystals

        % Get hold of header info
        assert((char(L(1))=='B'), ... % Check Euler angle convention
            'Could not read VPSC file - not Bunge format\n');
        tmp_nxtl = L(2); % Number of crystals

        % loop over each grain and extract Eulers
        tmp_eulers = zeros(3,tmp_nxtl);
        for i = 1:tmp_nxtl
            E = sscanf(fgetl(fid),'%g %g %g %g');
        
            % Build Euler angles array.
            tmp_eulers(1,i) = E(1,:);
            tmp_eulers(2,i) = E(2,:);
            tmp_eulers(3,i) = E(3,:);
        end
                
        % Build output...
        blocks = blocks + 1;
        if blocks == 1
            eulers = tmp_eulers;
            nxtl = tmp_nxtl;
        elseif blocks == 2
            % Multiple textures, bundle into a cell array...
            eulers = {eulers};
            eulers(2) = {tmp_eulers};
            nxtl(blocks) = tmp_nxtl;
        else
            eulers(blocks) = {tmp_eulers};
            nxtl(blocks) = tmp_nxtl;
        end
        
    end
    fclose(fid);
   
end
