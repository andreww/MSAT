%-------------------------------------------------------------------------------
%                  MSAT - Matlab Seismic Anisotropy Toolkit 
%-------------------------------------------------------------------------------
% MS_expand - expand a minimal set of elastic constants based on a specifed
%             symmetry
%-------------------------------------------------------------------------------
% 
% Usage: (some parts not yet implemented!)
%     [ Cf ] = MS_expand( C, nec, mode )
%         Fill out elastic tensor based on symmetry, defined by mode. This can
%         take the following values:
%            'auto' - assume symmetry based on number of Cijs specified 
%            'iso' - isotropic (nec=2) ; C33 and C66 must be specified.
%            'hex' - hexagonal (nec=5) ; C33, C44, C11, C66 and C13 must be
%                       specified, x3 is symmetry axis
%            'vti' - synonym for hexagonal
%            'cubic' - cubic (nec=3) ; C33, C66 and C12 must be specified

function [ C ] = MS_expand( Cin , nec, mode )

C = Cin ;

switch lower(mode)
case 'auto'
   switch nec
   case 2
      mode = 'iso' ;
   case 3
      mode = 'cubic' ;
   case 5
      mode = 'hex' ;
   otherwise
      error(['No automatic symmetry set for ' ...
         sprintf('%i',nec) ' elastic constants.']) ;
   end
otherwise
end

switch lower(mode)
case 'iso'
   C(1,1) = C(3,3) ; C(2,2) = C(3,3) ;
   C(5,5) = C(6,6) ; C(4,4) = C(6,6) ;
   C(1,2) = (C(3,3)-2.*C(4,4)) ;
   C(1,3) = C(1,2) ; C(2,3) = C(1,2) ;   
otherwise
   error('Unsupported symmetry.') ;
end

for i=1:6
   for j=(i+1):6
      C(j,i) = C(i,j) ;
   end
end

return 