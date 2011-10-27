% SINGLE_CRYSTAL_EXAMPLE.M - example script for the analysis of
%                            experimental / computational data
%                            giving single crystal elastic constnats
%
% This example shows how MSAT can be used to perform analysis of single
% crystal elasticity data. In this case the data is taken from calculations
% using density functional theory stored in the file diopside_P.txt - see
% http://dx.doi.org/10.1016/j.pepi.2011.10.002. From this data the
% the properties of an aggregate, along with some measures of anisotropy 
% are calculated and the anisotropy plotted as a function of wave
% propogation direction. 

% (C) James Wookey and Andrew Walker, 2011

function diopside_example()

    % First load in the data.
    [P, C, r] = MS_load_list('diopside_P.txt');

    % Write table header
    % P - pressure
    % rho - density 
    % K - isotropic VRH estimate of bulk mod.
    % G - isotropic VRH estimate of shear mod.
    % VPav - isotropic P-wave velo.
    % VSav - isotropic S-wave velo.
    % uA - Universal anisotropy index
    % lmA - Ledbetter and Miglion anisotropy index
    % VS1 - Fast S-wave velocity in [010] direction
    % VS2 - Slow S-wave velocity in [010] direction
    % aVS - S-wave anisotropy in [010] direction
    % dt - SWS delay time for [010] direction 
    fprintf(' P   rho   K    G    VPav VSav  uA   lmA  VS1  VS2  aVS  dt \n');
    fprintf(' GPa kg/m3 GPa  GPa  km/s km/s  -    -    km/s km/s %%    s\n');
    fprintf('============================================================\n');
    
    % Loop over pressures
    for i = 1:length(P)
        
        % Calculate anisotropy index
        [ uA, lmA ] = MS_anisotropy( C(:,:,i) );
        
        % Calculate isotropic elasticity
        [ K_vrh, G_vrh ] = MS_polyaverage( C(:,:,i) );
        
        % Get isotropic phase velocities
        Ciso = MS_build_isotropic('lam', G_vrh, 'K', K_vrh);
        [ ~, ~, VSiso, ~, VPiso ] = MS_phasevels( Ciso, r(i), 0.0, 0.0 );
        
        % Calculate splitting for [010] direction (5km layer)
        % This is along the y axis.
        [ ~, aVS, VS1, VS2, ~ ] = MS_phasevels( C(:,:,i), r(i), 0.0, 90.0 );
        dt = (10/VS2 - 10/VS1);
        
        % Write data for this pressure
        fprintf('%4.1f %4.0f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.2f\n', ...
            P(i), r(i), K_vrh, G_vrh, VPiso, VSiso, uA, lmA, VS1, VS2, aVS, dt);
        
        % Plot polefigs
        title = sprintf('Diopside seismic anisotropy at %2.1f GPa', P(i));
        MS_plot(C(:,:,i),r(i), 'wtitle', title, 'fontsize', 10, 'quiet');
    
    end
    fprintf('============================================================\n');
end
