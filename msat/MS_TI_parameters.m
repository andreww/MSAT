% MS_TI_parameters - TI parameters from C. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Returns various parameters describing TI anisotropy given an elasticity
% matrix (C) and density (rho) in GPa and kg/m^3. 
%
%  [loveA, loveC, loveL, loveN, loveF, vpv, vsv, eps, gam, del ...
%             vpa, vsa, xi, phi, eta] = MS_TI_parameters(C, rho)
%
%  Output:
%       vpa                               : isotropic average P-wave 
%                                           velocity, km/s 
%       vsa                               : isotropic average S-wave 
%                                           velocity, km/s
%       vpv                               : vertical P-wave velocity, km/s
%       vsv                               : vertical S-wave velocity, km/s
%       eps, gam, del                     : dimensionless Thomsen (1986) 
%                                           parameters
%       xi, phi, eta                      : dimensionless anisotropy 
%                                           parameters (e.g. Panning 
%                                           and Romanowicz, 2006)
%       loveA, loveC, loveL, loveN, loveF : Love's (1927) anisotropy 
%                                           parameters, GPa
%Notes
%~~~~~
%   Input units should be GPa for the elasticity matrix and kg/m^3 for the
%   density - this gives output parameters in units of km/s for velocity and 
%   GPa for elasticity. Choosing a different unit for the elasticity will 
%   change the units for the Love parameters and velocities but not the 
%   dimensionless parameters. Choosing an arbitrary value for the density
%   will not alter the dimensionless parameters or Love parameters, but the 
%   velocities will be meaningless. 
%
%References
%~~~~~~~~~~
%   Thomsen, L. (1986) "Weak elastic anisotropy" Geophysics 
%      vol.51 pp.1954-1966
%
%   Mark Panning and Barbara Romanowicz (2006) A three-dimensional radially 
%      anisotropic model of shear velocity in the whole mantle. Geophysical  
%      Journal International v167, 361–379. 
%      doi: 10.1111/j.1365-246X.2006.03100.x
%
%   Love, A.E.H., (1927). A Treatise on the Theory of Elasticity, 
%      Cambridge Univ. Press, Cambridge.
%
% See also: MS_TI, MS_polyaverage
 
% Copyright (c) 2011-2012, James Wookey, Andrew Walker and Alan Baird
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

function [loveA, loveC, loveL, loveN, loveF, vpv, vsv, eps, gam, del ...
             vpa, vsa, xi, phi, eta] = MS_TI_parameters(C, rho)

    % Check we have valid input...
    MS_checkC(C);
    if ~is_Ti(C)
       error('MS:BadTIelasticity',...
          'MS_TI_parameters requires VTI input') ;
    end
    
    % Pass off the calculations to subroutines.
    [eps, gam, del] = MS_thomsen_params(C);
    [loveA, loveC, loveL, loveN, loveF] = MS_love_params(C);
    [vpa, vsa, xi, phi, eta, vpv, vsv] = MS_panning_params(C, rho);

end

function [vpa, vsa, xi, phi, eta, vpv, vsv]=MS_panning_params(C, rho)

   % Convert to Pa
   C = C.*1e9 ;

   % Love params in Pa and density in Kg/m^3 give 
   % velocites in m/s
   [A, C, L, N, F] = MS_love_params(C);
   vph = sqrt(A/rho);
   vpv = sqrt(C/rho); % Thom VP
   vsv = sqrt(L/rho); % Thom VS
   vsh = sqrt(N/rho);
   
   %  note: xi and phi are defined oppositely; i.e.: 
   %  xi = vsh^2/vsv^2 and phi = vpv^2/vph^2
   xi = vsh^2 / vsv^2;
   phi = vpv^2 / vph^2;
   
   vpa = sqrt((vpv^2 + 4.0*vph.^2)/5.0) ;
   vsa = sqrt((2.0*vsv.^2 + vsh.^2)/3.0) ;
   
   eta = F/(A-2.0*L);
   
   %  convert to km/s
   vpa = vpa/1e3;
   vsa = vsa/1e3;
   vpv = vpv/1e3;
   vsv = vsv/1e3;

end

function [eps, gam, del] = MS_thomsen_params(C)
    % From a VTI (unique X3) elasticity matrix, return 
    % the three Thomsen parameters eps, gam and del

    % FIXME - do we need to do the 'vs too high' test here?
    
    eps=(C(1,1)-C(3,3))/(2.0*C(3,3));
    gam=(C(6,6)-C(4,4))/(2.0*C(4,4));
    del=((C(1,3)+C(4,4))^2.0-(C(3,3)-C(4,4))^2.0) / ...
        (2.0*C(3,3)*(C(3,3)-C(4,4)));

end


function [A, C, L, N, F] = MS_love_params(CC)
    % From a VTI (uneq X3) elasticity matrix, return 
    % the five Love parameters A, C, L, N and F 
    A = CC(1,1); C = CC(3,3) ;
    L = CC(4,4); N = CC(6,6) ;
    F = CC(1,3);
   
end

function [check] = is_Ti(C)
    % Check if input elasticity matrix is VTI with unique X3 axis
    % return 1 if it is and 0 if not.
    
    thresh = 1e-6 ;

    if ( (abs(C(1,4)) > thresh) || (abs(C(1,5)) > thresh) || ... 
         (abs(C(1,6)) > thresh) || (abs(C(2,4)) > thresh) || ...
         (abs(C(2,5)) > thresh) || (abs(C(2,6)) > thresh) || ...
         (abs(C(3,4)) > thresh) || (abs(C(3,5)) > thresh) || ...
         (abs(C(3,6)) > thresh) || (abs(C(4,5)) > thresh) || ...         
         (abs(C(4,6)) > thresh) || (abs(C(5,6)) > thresh) || ...
         (abs(C(1,1)-C(2,2)) > thresh) || ...
         (abs(C(4,4)-C(5,5)) > thresh) || ...
         (abs(C(2,3)-C(1,3)) > thresh) || ...
         (abs(C(1,2)-(C(1,1)-2.0*C(6,6))) > thresh) ) 
        check = 0;
    else
        check = 1;  
    end

end

