% MS_CALC_ANISOTROPY - Simple measures of anisotropy 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the degree of ansiotropy of an elasticity matrix. 
%
% [ uA, ... ] = MS_calc_anisotropy( C )
%
% Usage: 
%     [ uA ] = MS_calc_anisotropy( C )                    
%         Return the Universal Elastic Anisotropy Index of Ranganathan
%         and Ostoja-Starzewski (2008). Valid for any elasticity matrix,
%         uA is zero for an isotropic case and increses for incresing 
%         anisotropy.
%
%     [ uA, zA ] = MS_calc_anisotropy( C )
%         Also return the Zenner (1948) measure of anisotropy. This is
%         only valid for cubic crystals (NaN is returned if C does not 
%         represent a cubic crystal). zA is 1 for an isotropic case and
%         increses or decreses with incresing anisotropy.
%
%     [ uA, zA, cbA ] = MS_calc_anisotropy( C )
%         Also return the Chung-Buessem (1967) anisotropy index. This
%         is a single valued measure of anisotropy derived from zA. Like
%         uA, this is zero for an isotropic case and increses for incresing 
%         anisotropy. Only valid for matricies representing cubic crystals.
%
% Notes:
%     These measures of anisotropy are independent of orientation. However,
%     the test for cubic symmetry assumes the matrix is in an ideal
%     orention. Use MS_AXES to reorentate the imput matrix for the general 
%     case. MS_NORMS can be used to provide an alternate measure of
%     anisotropy.
%
% References:
%
%     Zenner, C. (1948) Elasticity and Anelasticiy of Metals. University 
%     of Chicago. 
%
%     Chung, D. H. and W. R. Buessem (1967) Journal of Applied Physics 
%     vol.38 p.5 
%
%     Shivakumar, I. and M. Ostoja-Starzewski (2008) "Universal Elastic
%     Anisotropy Index" Physical Review Letters vol.101 art.num.055504.
%     http://dx.doi.org/10.1103/PhysRevLett.101.055504
%
% See also: MS_POLYAVERAGE, MS_NORMS, MS_AXES, MS_PHASEVELS

% (C) James Wookey and Andrew Walker, 2011

function [ uA, zA, cbA ] = MS_calc_anisotropy( C )

    uA = UniversalAnisotropyIndex(C);
    if (isCubic( C ))
        zA = ZenerAnisotropyIndex(C);
        cbA = (3*(zA-1)^2) / (3*(zA-1)^2 + 25*zA);
    else
        zA = NaN;
        cbA = NaN;
    end
end

function i = isCubic ( C )

    i = (C(1,1) == C(2,2)) ...
      & (C(1,1) == C(3,3)) ...
      & (C(1,2) == C(1,3)) ...
      & (C(1,2) == C(2,3)) ...
      & (C(4,4) == C(5,5)) ...
      & (C(4,4) == C(6,6)) ...
      & (C(1,4) == 0.0) ...
      & (C(1,5) == 0.0) ...
      & (C(1,6) == 0.0) ...
      & (C(2,4) == 0.0) ...
      & (C(2,5) == 0.0) ...
      & (C(2,6) == 0.0) ...
      & (C(3,4) == 0.0) ...
      & (C(3,5) == 0.0) ...
      & (C(3,6) == 0.0) ...
      & (C(4,5) == 0.0) ...
      & (C(4,6) == 0.0) ...
      & (C(5,6) == 0.0);
end

function zA = ZenerAnisotropyIndex( C )

    zA = C(4,4)^2 / (C(1,1)^2 - C(1,2)^2); 

end 

function uA = UniversalAnisotropyIndex( C )

    [~, ~, B_v, G_v, B_r, G_r] = MS_polyaverage( C );
    uA = (5*(G_v/G_r) + (B_v/B_r)) - 6;
    
end