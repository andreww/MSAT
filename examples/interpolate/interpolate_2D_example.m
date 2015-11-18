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
    
    % Set up and plot values we 'know' - on one end we have rotated 
    % Fo90, on the other we have rotated fayalite
    
    % Forsterite: X = 0, Y=0
    [C_x0y0, rh_x0y0] = MS_elasticDB('olivine');
    C_x0y0 = MS_rot3(C_x0y0, 0, 0, 30);
    if phase_vels
       MS_plot(C_x0y0, rh_x0y0, 'wtitle', '100% forsterite', 'quiet');
    end 
    
    % Fayalite: X = 1, Y=0
    [C_x1y0, rh_x1y0] = MS_elasticDB('fayalite');
    C_x1y0 = MS_rot3(C_x1y0, 0, 0, 30);
    if phase_vels
        MS_plot(C_x1y0, rh_x1y0, 'wtitle', '100% fayalite', 'quiet');
    end
    
    % Diopside: X = 0, Y=1
    [C_x0y1, rh_x0y1] = MS_elasticDB('olivine');
    C_x0y1 = MS_rot3(C_x0y1, 0, 0, 50);
    if phase_vels
       MS_plot(C_x0y1, rh_x0y1, 'wtitle', '100% forsterite', 'quiet');
    end 
    
    % Enstatite: X = 1, Y=1
    [C_x1y1, rh_x1y1] = MS_elasticDB('fayalite');
    C_x1y1 = MS_rot3(C_x1y1, 0, 0, 50);
    if phase_vels
        MS_plot(C_x1y1, rh_x1y1, 'wtitle', '100% fayalite', 'quiet');
    end

    % Points for the interpolation
    [X, Y] = meshgrid(0:0.1:1, 0:0.1:1);
    C_interp = zeros([6 6 size(X)]);
    C_interp_yx = zeros([6 6 size(X)]);
    C_voigt = zeros([6 6 size(X)]);
    r_interp = zeros(size(X));
    r_interp_yx = zeros(size(X));
    r_voigt = zeros(size(X));
    Vp_interp = zeros(size(X));
    Vp_interp_yx = zeros(size(X));
    Vp_voigt = zeros(size(X));
    P_ort_interp = zeros(size(X));
    P_ort_interp_yx = zeros(size(X));
    P_ort_voigt = zeros(size(X));
    [ni, nj] = size(X);
    
    inc = 0.0;
    azi = 0.0;
    
    for i = 1:ni
        for j = 1:nj
            [C_interp(:,:,i,j), r_interp(i,j)] = bilinear_interp(C_x0y0,... 
                rh_x0y0, C_x1y0, rh_x1y0, C_x0y1, rh_x0y1, C_x1y1, ...
                rh_x1y1, X(i,j), Y(i,j), 0);
            [C_interp_yx(:,:,i,j), r_interp_yx(i,j)] = bilinear_interp(C_x0y0,... 
                rh_x0y0, C_x1y0, rh_x1y0, C_x0y1, rh_x0y1, C_x1y1, ...
                rh_x1y1, X(i,j), Y(i,j), 1);
            [C_voigt(:,:,i,j), r_voigt(i,j)] = bilinear_VRH(C_x0y0,... 
                rh_x0y0, C_x1y0, rh_x1y0, C_x0y1, rh_x0y1, C_x1y1, ...
                rh_x1y1, X(i,j), Y(i,j));
            
            [~, ~, ~, ~, Vp_interp(i,j)] = ...
                MS_phasevels( C_interp(:,:,i,j), r_interp(i,j), inc, azi );
            [~, ~, ~, ~, Vp_interp_yx(i,j)] = ...
                MS_phasevels( C_interp_yx(:,:,i,j), r_interp_yx(i,j), inc, azi );
            [~, ~, ~, ~, Vp_voigt(i,j)] = ...
                MS_phasevels( C_voigt(:,:,i,j), r_voigt(i,j), inc, azi );
            
            [P_ort_interp(i,j), ~,~] = summary_anisotropy(C_interp(:,:,i,j));
            [P_ort_interp_yx(i,j), ~,~] = summary_anisotropy(C_interp_yx(:,:,i,j));
            [P_ort_voigt(i,j), ~,~] = summary_anisotropy(C_voigt(:,:,i,j));
        end
    end
    
    figure
    subplot(3, 1, 1);
    contourf(X, Y, Vp_interp);
    colorbar;
    title('Interp')
    subplot(3, 1, 2);
    contourf(X, Y, Vp_interp_yx);
    colorbar;
    title('Voigt')
    subplot(3, 1, 3);
    contourf(X, Y, Vp_voigt);
    colorbar;
    title('Voigt')
    
        figure
    subplot(3, 1, 1);
    contourf(X, Y, P_ort_interp);
    colorbar;
    title('Interp')
    subplot(3, 1, 2);
    contourf(X, Y, P_ort_interp_yx);
    colorbar;
    title('Voigt')
    subplot(3, 1, 3);
    contourf(X, Y, P_ort_voigt);
    colorbar;
    title('Voigt')
    
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
    
%     
%     % Do the Voigt (element-wise) average
%     i = 0;
%     tic;
%     for x = xvals
%         i = i + 1;
%         [~, rh_interp, C_voigt] = MS_VRH([x 1-x], C_x0, rh_x0, C_x1, rh_x1);
%         if gather_results
%             [P_ort_voigt(i), P_low_voigt(i), lmA_voigt(i)] = ...
%                 summary_anisotropy(C_voigt);
%         end
%         
%         % Plot five values
%         if (mod(x,0.25) == 0) && phase_vels
%             plottitle = sprintf( ...
%                 '%d%% forsterite %d%% fayalite, elementwise (Voigt mean)', ...
%                 (x*100), (1-x)*100);          
%             MS_plot(C_voigt, rh_interp, 'wtitle', plottitle, 'quiet');
%         end
%         
%     end
%     time_per_voigt = toc/i;
%     
% 
%     % Do the MS_interpolate interpolation
%     i = 0;
%     tic;
%     for x = xvals
%         i = i + 1;            
%         C_interp = MS_interpolate(C_x0, rh_x0, C_x1, rh_x1, x);
%         if gather_results
%             [P_ort_interp(i), P_low_interp(i), lmA_interp(i)] = ...
%                 summary_anisotropy(C_interp);
%         end
%         
%         % Plot five values
%         if (mod(x,0.25) == 0) && phase_vels
%             plottitle = sprintf( ...
%                 '%d%% forsterite %d%% fayalite, (common orientation)', ...
%                 (x*100), (1-x)*100);
% 
%             MS_plot(C_interp, rh_interp, 'wtitle', plottitle, 'quiet');
%         end
%         
%     end
%     time_per_interp = toc/i;
%     
%     if gather_results
%         % Plot up the summary data
%         figure;
%         plot(xvals, P_ort_voigt*100, 'b-');
%         hold on;
%         plot(xvals, P_low_voigt*100, 'b--');
%         plot(xvals, P_ort_interp*100, 'g-');
%         plot(xvals, P_low_interp*100, 'g--');
%         title('Anisotropic decomposition for interpolated points');
%         legend('Consistent with orthorhombic symmetry, Voigt', ...
%             'Inconsistent with orthorhombic symmetry, Voigt', ...
%             'Consistent with orthorhombic symmetry, common orientation', ...
%             'Inconsistent with orthorhombic symmetry, common orientation', ...
%             'location', 'southoutside');
%         xlabel('Volume fraction forsterite')
%         ylabel('Norm of decomposed elasticity, \%')
%         hold off;
%         
%         figure
%         plot(xvals, lmA_voigt, 'b-');
%         hold on;
%         plot(xvals, lmA_interp, 'g-');
%         title('Universal anisotropy index for interpolated points')
%         legend('Voigt interpolation', 'Inperpolation based on common orientations')
%         xlabel('Volume fraction forsterite')
%         ylabel('Universal anisotropy index')
%         hold off;
%     end
%     
%     if report_time
%         fprintf('Time per Voigt interpolation = %f s\n', time_per_voigt)
%         fprintf('Time per common orientation interpolation = %f s\n', time_per_interp)
%     end
%     
% end

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
