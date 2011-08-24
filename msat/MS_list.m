% MS_LIST - Print elasticity matrix (in 'list' format)
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% CIJ_list(C,rho)
%
% Usage: 
%     CIJ_list(C,rho)
%         For a given 6x6 elasticity matrix, C, and density, rho, print 
%         the constants.
%
%  Notes: 
%    The 'list' format consists of three numbers per line. The first two 
%    are integers representing the location in the elasticity matrix in
%    Voigt notation, the third is the value of the elastic constant. After
%    the 36 lines of elastic constants a line containing the density is
%    printed with dummy indices '7 7'. This function exists to allow 
%    MSAT to be interfaced with legacy code.
%

% (C) James Wookey and Andrew Walker, 2011

function MS_list(C,rho)

    for i=1:6
       for j=i:6
          fprintf('%1i %1i %e\n',i,j,C(i,j))
       end
    end   
    fprintf('%1i %1i %f\n',7,7,rho)
end
