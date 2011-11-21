% MS_ROTM - Create a cartesian rotation matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [ M ] = MS_rotM( alp, bet, gam )
%
%  Inputs:
%     alp - clockwise about 1-axis (looking at origin, == yaw)
%     bet - clockwise about 2-axis (looking at origin, == -dip)
%     gam - clockwise about 3-axis (looking at origin, == azimuth)
%           angles are in ** degrees **, they should be scalars. 
%
%  Output: is M, a 3x3 rotation matrix
%
%     [M] = MS_rotM(...,'order',V)                    
%         order: V is a 3 element vector containing the order in which to apply the 3 rotations.
%                Default is [1 2 3].
%
%  Notes: 
%    The rotations are applied in order, ie: alpha, then beta then gamma (by default).
%

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

function [ M ] = MS_rotM( alp, bet, gam, varargin )


orderV = [1 2 3] ;

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'order' 
         orderV = varargin{iarg+1} ;
         if (length(find(orderV==1))+length(find(orderV==2))+length(find(orderV==3)))~=3
            error('MS:ROTM:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;   
         end
         if (length(orderV)~=3)
            error('MS:ROTM:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;  
         end                  
         iarg = iarg + 2 ;
      otherwise 
         error('MS:ROTM:UnknownOption',...
            ['Unknown option: ' varargin{iarg}]) ;   
   end   
end

a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

%  Make rotation matrices

R(1,:,:) = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R(2,:,:) = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R(3,:,:) = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

M =  squeeze(R(orderV(3),:,:)) * ...
         squeeze(R(orderV(2),:,:)) * ...
            squeeze(R(orderV(1),:,:));

return

