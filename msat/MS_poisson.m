% MS_POISSON - Poisson's ratio in anisotropic media.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate Poisson's ration from an elsticity matrix. 
%
%   nu = MS_poisson( C, m, n )
%
% Usage: 
%     nu = MS_poisson( C, m, n )                    
%         Calculate Poisson's ration for anisotropic media given an
%         elasticity matrix, C, longitudinal direction described by
%         a unit vector, n, and transverse direction given by the 
%         unit vector m.
%
% Notes:
%    Poisson's ratio, nu, is normally described as the negative ratio of 
%    a (typically) compressive strain normal a direction of extensional 
%    strain. In an isotropic material nu takes a single value between -1.0
%    an 0.5 (positive values around 0.2 or 0.3 are common). For the general
%    anisotropic case Poisson's ratio depends on the longitudinal extension
%    direction and the direction perpendicular to this (the transverse
%    direction) where the compression (positive nu) or extension (negative 
%    nu) is measured. This function allows nu to be calculated as a function
%    of these two directions, given as unit vectors. The approach follows 
%    equation 3 of Wojciechowski (2005)
%
% 
% Reference: Wojciechowski, K. W. (2005) "Poisson's ratio of anisotropic
%     systems" COMPUTATIONAL METHODS IN SCIENCE AND TECHNOLOGY 11, 73-79.
%
% See also: MS_PHASEVELS

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

function [ nu ] = MS_poisson( C, m, n )

    % check the inputs: C
    assert(MS_checkC(C)==1, 'MS:POISSON:badC',  ...
        'MS_checkC error MS_poisson') ;
    assert(MS_checkUnit(m)==1, 'MS:POISSON:badM',  ...
        'MS_checkUnit error for m in MS_poisson') ;
    assert(MS_checkUnit(n)==1, 'MS:POISSON:badN',  ...
        'MS_checkUnit error for n in MS_poisson') ;
    assert((abs(dot(m,n))<sqrt(eps)), 'MS:POISSON:NotOrthogonal',  ...
        'Unit vectors M and N are not orthogonal') ;

    S4 = MS_cij2cijkl(inv(C));
    
    Smmnn = 0.0;
    Snnnn = 0.0;
    for i = 1:3
      for j = 1:3
        for k = 1:3
          for l = 1:3
              Smmnn = Smmnn + m(i)*m(j)*n(k)*n(l)*S4(i,j,k,l);
              Snnnn = Snnnn + n(i)*n(j)*n(k)*n(l)*S4(i,j,k,l);
          end
        end
      end
    end
    nu = -1.0*(Smmnn/Snnnn);

end

function [ isgood ] = MS_checkUnit( V )

    errorval = sqrt(eps);
    % Is V a unit vector?
    if ~isnumeric(V)
         error('MS:CHECKUNITNotNumeric', ...
	      'Unit vector error: Appears not to be numeric.');
    end
    s = size(V);
    if (length(s) ~= 2) || ... 
        ~(((s(1) == 1) && (s(2) == 3)) || ((s(1) == 3) && (s(2) == 1)))
         error('MS:CHECKUNITNotVec', ...
	      'Unit vector error: Not a 3-d row or column vector.');
    end
    if (abs((sqrt(V(1)^2 + V(2)^2 + V(3)^2)) - 1.0) > errorval)
        error('MS:CHECKUNITNotUnit', ...
	      'Unit vector error: Length not equal to 1.');
    end
    isgood = 1;
    

end