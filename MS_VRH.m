%  [Cave,rhave]=CIJ_VRH(VF,C1,rho1,C2,rho2,...)
%
%  n phase Voigt-Reuss-Hill average of elastic tensors (C1, C2 etc) and
%  densities (rho1, rho2 etc). VF is a vector of volume fractions. There should
%  be a C,rho argument pair for each entry in this vector. The sum of the 
%  elements of VF is normalised to 1 before averaging.  
%
%  (c) James Wookey, 2006
%
function [Cave,rhave]=CIJ_VRH(VF,varargin)

%--check the input
   if length(varargin)~=length(VF)*2
      error('There must be as many C,rho argument pairs as elements in VF') ;
   end   

%--normalise VF
   VF = VF ./ sum(VF) ;

   voigt_ave = zeros(6,6) ;
   reuss_ave = zeros(6,6) ;
   rhave = 0 ;

   n = length(VF) ;

   for i=1:n
      CR = varargin{(i-1).*2+1} ;
      rh = varargin{(i-1).*2+2} ;
      voigt_ave = voigt_ave + CR .* VF(i) ;      % stiffness average
      reuss_ave = reuss_ave + inv(CR) .* VF(i) ; % compliance average
      rhave = rhave + rh .* VF(i) ;
   end

   Cave = (voigt_ave + inv(reuss_ave)) ./ 2  ;

return
