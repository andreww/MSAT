% MS_POLYAVERAGE - Isotropic elasticity for polycrystal 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the bulk and shear moduli of an isotropic polycrystal. 
%
% [ K_vrh, G_vrh, ... ] = MS_polyavarage( C )
%
% Usage: 
%     [ K_vrh, G_vrh ] = MS_polyavarage( C )                    
%         Return Voigt-Reuss-Hill avarage bulk (K) and shear (G) moduli of 
%         an isotropic aggregate composed of randomly orentated ansiotropic
%         crystals described by the elasticity matrix C.
%
%     [ K_vrh, G_vrh, K_v, G_v, K_r, G_r ] = MS_polyavarage( C )
%         Also return Voigt and Reuss bounds on the bulk and shear moduli.
%
% Notes:
%     For the case of an isotropic solid formed from a random distribution
%     of elements, the bulk and shear moduli can be found directly without
%     the need to rotate the crystal. See Anderson, Theory of the Earth,
%     pp. 122.
%
%
% See also: MS_VRH, MS_ANISOTROPY

% (C) James Wookey and Andrew Walker, 2011

function [ K_vrh, G_vrh, K_v, G_v, K_r, G_r ] = MS_polyaverage( C )
    
    K_v = (1.0/9.0) * (C(1,1) + C(2,2) + C(3,3)) + ...
          (2.0/9.0) * (C(1,2) + C(1,3) + C(2,3));
    G_v = (1.0/15.0) * (C(1,1) + C(2,2) + C(3,3) - ...
                        C(1,2) - C(1,3) - C(2,3)) + ...
          (1.0/5.0) * (C(4,4) + C(5,5) + C(6,6));

    S = inv(C);
    K_r = 1.0 / (         (S(1,1) + S(2,2) + S(3,3)) + ...
                  ( 2.0 * (S(1,2) + S(1,3) + S(2,3))) );
    G_r = 15.0 / ( 4.0 * (S(1,1) + S(2,2) + S(3,3)) - ...
                   4.0 * (S(1,2) + S(1,3) + S(2,3)) + ...
                   3.0 * (S(4,4) + S(5,5) + S(6,6)) );
               
    K_vrh = (K_v+K_r)/2.0;           
    G_vrh = (G_v+G_r)/2.0;

end