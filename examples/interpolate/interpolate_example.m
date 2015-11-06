% INTERPOLATE_EXAMPLE.m - example script showing two ways elasticity 
%                         matricies can be interpolated in MSAT.
%
% Interpolation of elastic constants can be a difficult problem. In this
% example we show how two approaches, both implemented in MSAT, give quite 
% different results. We imagine a region where the elasticity changes 
% smoothly from one described by single crystal olivine rotated 30 degrees 
% clockwise around its c-axis to one described by single crystal enstatie 
% rotated 30 degrees anticlockwise around its c-axis and seek to calculate 
% the elasticity of a point mid-way between these limits. One (commonly 
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
% See also: MS_VRH, MS_INTERPOLATE

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

function interpolate_example(varargin)

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
    
    % Forsterite: X = 0
    [C_x0, rh_x0] = MS_elasticDB('olivine');
    C_x0 = MS_rot3(C_x0, 0, 0, 60);
    if phase_vels
       MS_plot(C_x0, rh_x0, 'wtitle', '100% forsterite', 'quiet');
    end 
    
    % Fayalite: X = 1
    [C_x1, rh_x1] = MS_elasticDB('fayalite');
    C_x1 = MS_rot3(C_x1, 0, 0, 120);
    if phase_vels
        MS_plot(C_x1, rh_x1, 'wtitle', '100% fayalite', 'quiet');
    end

    
    % We want to record some summary information
    % arrays for holding the info, and setting the 
    % x-spacing.
    xvals = 0.0:0.0025:1.0;
    P_ort_voigt = xvals * 0.0;
    P_low_voigt = P_ort_voigt;
    lmA_voigt = P_ort_voigt;
    P_ort_interp = P_ort_voigt;
    P_low_interp = P_ort_voigt;
    lmA_interp = P_ort_voigt;
    
    
    % Do the Voigt (element-wise) average
    i = 0;
    tic;
    for x = xvals
        i = i + 1;
        [~, rh_interp, C_voigt] = MS_VRH([x 1-x], C_x0, rh_x0, C_x1, rh_x1);
        if gather_results
            [P_ort_voigt(i), P_low_voigt(i), lmA_voigt(i)] = ...
                summary_anisotropy(C_voigt);
        end
        
        % Plot five values
        if (mod(x,0.25) == 0) && phase_vels
            plottitle = sprintf( ...
                '%d%% forsterite %d%% fayalite, elementwise (Voigt mean)', ...
                (x*100), (1-x)*100);          
            MS_plot(C_voigt, rh_interp, 'wtitle', plottitle, 'quiet');
        end
        
    end
    time_per_voigt = toc/i;
    

    % Do the MS_interpolate interpolation
    i = 0;
    tic;
    for x = xvals
        i = i + 1;            
        rh_interp = rh_x0*x + rh_x1*(1-x);
        C_interp = MS_interpolate(C_x0, C_x1, x);
        if gather_results
            [P_ort_interp(i), P_low_interp(i), lmA_interp(i)] = ...
                summary_anisotropy(C_interp);
        end
        
        % Plot five values
        if (mod(x,0.25) == 0) && phase_vels
            plottitle = sprintf( ...
                '%d%% forsterite %d%% fayalite, (common orientation)', ...
                (x*100), (1-x)*100);

            MS_plot(C_interp, rh_interp, 'wtitle', plottitle, 'quiet');
        end
        
    end
    time_per_interp = toc/i;
    
    if gather_results
        % Plot up the summary data
        figure;
        plot(xvals, P_ort_voigt*100, 'b-');
        hold on;
        plot(xvals, P_low_voigt*100, 'b--');
        plot(xvals, P_ort_interp*100, 'g-');
        plot(xvals, P_low_interp*100, 'g--');
        title('Anisotropic decomposition for interpolated points');
        legend('Consistent with orthorhombic symmetry, Voigt', ...
            'Inconsistent with orthorhombic symmetry, Voigt', ...
            'Consistent with orthorhombic symmetry, common orientation', ...
            'Inconsistent with orthorhombic symmetry, common orientation', ...
            'location', 'southoutside');
        xlabel('Volume fraction forsterite')
        ylabel('Norm of decomposed elasticity, \%')
        hold off;
        
        figure
        plot(xvals, lmA_voigt, 'b-');
        hold on;
        plot(xvals, lmA_interp, 'g-');
        title('Universal anisotropy index for interpolated points')
        legend('Voigt interpolation', 'Inperpolation based on common orientations')
        xlabel('Volume fraction forsterite')
        ylabel('Universal anisotropy index')
        hold off;
    end
    
    if report_time
        fprintf('Time per Voigt interpolation = %f s\n', time_per_voigt)
        fprintf('Time per common orientation interpolation = %f s\n', time_per_interp)
    end
    
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
         P_ort = sum(P(2:4))./sum(P(2:6));
         P_low = sum(P(5:6))./sum(P(2:6));
         [Ua] = MS_anisotropy(C);

end
