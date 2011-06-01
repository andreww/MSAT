%-------------------------------------------------------------------------------
%                  MSAT - Matlab Seismic Anisotropy Toolkit 
%-------------------------------------------------------------------------------
% MS_expand - expand a minimal set of elastic constants based on a specifed
%             symmetry to a full Cij tensor. 
%-------------------------------------------------------------------------------
% 
% Usage: (some parts not yet implemented!)
%     [ Cf ] = MS_expand( C, mode )
%         Fill out elastic tensor C based on symmetry, defined by mode. This can
%         take the following values:
%            'auto' - assume symmetry based on number of Cijs specified 
%            'iso' - isotropic (nec=2) ; C33 and C66 must be specified.
%            'hex' - hexagonal (nec=5) ; C33, C44, C11, C66 and C13 must be
%                       specified, x3 is symmetry axis
%            'vti' - synonym for hexagonal
%            'cubic' - cubic (nec=3) ; C33, C66 and C12 must be specified
%
%     Cijs *not* specified in the appropriate symmetry should be zero in the 
%     input matrix. 
%

function [ C ] = MS_expand( Cin , mode )

try 
   MS_checkC(C,'fast') ; 
catch
   error('Bad input elasticity matrix.')
end   

C = zeros( 6 , 6 ) ;


nec = length( find( Cin ~= 0 ) ) ;

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
      error('MS:EXPANDnoautosym',['No automatic symmetry set for ' ...
         sprintf('%i',nec) ' elastic constants.']) ;
   end
otherwise
%  mode specifies symmetry
end

switch lower(mode)
case 'iso'
%  check nec
   if nec~=2, error('MS:EXPANDbadiso',...
      'Isotropic expansion requires C33 and C66 to be set') ;, end
%  check that C(1,1) and C(6,6) are set
   if Cin(3,3)==0 | Cin(6,6)==0
      error('MS:EXPANDbadiso',...
         'Isotropic expansion requires C33 and C66 to be set.')
   end
   C(3,3) = Cin(3,3) ;
   C(6,6) = Cin(6,6) ;
   
   C(1,1) = C(3,3) ; C(2,2) = C(3,3) ;
   C(5,5) = C(6,6) ; C(4,4) = C(6,6) ;
   C(1,2) = (C(3,3)-2.*C(4,4)) ;
   C(1,3) = C(1,2) ; C(2,3) = C(1,2) ;   
otherwise
   error('MS:EXPANDunsupportsymmetry','Unsupported symmetry.') ;
end

for i=1:6
   for j=(i+1):6
      C(j,i) = C(i,j) ;
   end
end

% check the resulting matrix
try 
   MS_checkC(C) ;
catch ME
   error('MS:EXPANDcheckfailed', ...
      ['MS_expand: Resulting Cij matrix failed checks with error: ' ME.message])
end

return 