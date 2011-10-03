% MS_EXPAND - fill out elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Expand a minimal set of elastic constants based on a specifed
%             symmetry to a full Cij tensor. 
%
% [ Cf ] = MS_expand( C, mode )
%
% Usage: 
%     Fill out elastic tensor C based on symmetry, defined by mode. This 
%     can take the following values:
%        'auto' - assume symmetry based on number of Cijs specified 
%        'iso' - isotropic (nec=2) ; C33 and C66 must be specified.
%        'hex' - hexagonal (nec=5) ; C33, C44, C11, C66 and C13 must be
%                specified, x3 is symmetry axis
%        'vti' - synonym for hexagonal
%        'cubic' - cubic (nec=3) ; C33, C66 and C12 must be specified
%
%     Cijs *not* specified in the appropriate symmetry should be zero in 
%     the input matrix. 
%
%
% See also: MS_ELASTICDB MS_LOAD MS_LOAD_LIST

% Copyright (c) 2011, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
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

function [ C ] = MS_expand( Cin, mode )

try 
   MS_checkC(Cin,'fast') ; 
catch ME
   error(ME.identifier, ['Bad input elasticity matrix:' ME.message])
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
         'Isotropic expansion requires C33 and C66 to be set.') ;
   end
   
   C(3,3) = Cin(3,3) ;
   C(6,6) = Cin(6,6) ;
   
   C(1,1) = C(3,3) ; C(2,2) = C(3,3) ;
   C(5,5) = C(6,6) ; C(4,4) = C(6,6) ;
   C(1,2) = (C(3,3)-2.*C(4,4)) ;
   C(1,3) = C(1,2) ; C(2,3) = C(1,2) ;

case 'cubic'
   if nec~=3, error('MS:EXPANDbadcubic',...
      'Cubic expansion requires C33, C66 and C12 to be set') ;, end
%  check that C(1,1) and C(6,6) are set
   if Cin(3,3)==0 | Cin(6,6)==0 | Cin(1,2) == 0
      error('MS:EXPANDbadcubic',...
         'Isotropic expansion requires C33, C66 and C12 to be set.') ;
   end
   
   C = Cin ;
   C(1,1) = C(3,3) ; C(2,2) = C(3,3) ;
   C(1,3) = C(1,2) ; C(2,3) = C(1,2) ;
   C(4,4) = C(6,6) ; C(5,5) = C(6,6) ;

case {'hex', 'vti'}
   if nec~=5, error('MS:EXPANDbadhex',...
      'Hexagonal expansion requires C11, C33, C44, C66 and C13 to be set.') ;, end
%  check that C(1,1) and C(6,6) are set
   if Cin(1,1)==0 | Cin(3,3)==0 | Cin(1,3) == 0 | Cin(4,4)==0 | Cin(6,6) == 0 
      error('MS:EXPANDbadhex',...
         'Isotropic expansion requires C11, C33, C44, C66 and C13 to be set.') ;
   end
   
   C = Cin ;
   
   C(1,2) = (C(1,1)-2.*C(6,6)) ;
   C(2,2) = C(1,1) ; C(2,3) = C(1,3) ;
   C(5,5) = C(4,4) ;
         
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