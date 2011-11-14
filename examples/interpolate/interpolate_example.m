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

% (C) James Wookey and Andrew Walker, 2011

function interpolate_example()

    % Set up and plot values we 'know' - on one end we have rotated 
    % olivine, on the other we have rotated enstatite
    [C0_ol, rh_ol] = MS_elasticDB('olivine');
    C30_ol = MS_rot3(C0_ol, 0, 0, 30);
    MS_plot(C30_ol, rh_ol, 'wtitle', '100% olivine', 'quiet');
    fprintf(['\n100%%' ' olivine\n'])
    MS_info(C30_ol, rh_ol);
    
    [C0_en, rh_en] = MS_elasticDB('enstatite');
    C30_en = MS_rot3(C0_en, 0, 0, -30);
    MS_plot(C30_en, rh_en, 'wtitle', '100% enstatite', 'quiet');
    fprintf(['\n100%%' ' enstatite\n'])
    MS_info(C30_en, rh_en);
    
    % Do the VRH (element-wise) mean.
    [C0_vrh, rh_vrh] = MS_VRH([0.5 0.5], C30_ol, rh_ol, C30_en, rh_en);
    MS_plot(C0_vrh, rh_vrh, 'wtitle', ...
        '50% olivine 50% enstatite (VRH)', 'quiet');
    fprintf(['\n50%%' ' VRH mean\n'])
    MS_info(C0_vrh, rh_vrh);
    
    % Do the MS_interpolate interpolation
    C0_int = MS_interpolate(C30_ol, C30_en, 0.5);
    MS_plot(C0_int, rh_vrh, 'wtitle', ...
        '50% olivine 50% enstatite (common orientation)', 'quiet');
    fprintf(['\n50%%' ' common orientation mean\n'])
    MS_info(C0_int, rh_vrh);
    
end

