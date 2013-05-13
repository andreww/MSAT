% MS_SPLIT_TRACE - Apply splitting operator to a pair of orthogonal traces
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast,tlag)
%
%  Inputs:
%     time (vector) : Series of time points where amplitudes are known (s)
%     T00 (vector)  : Series of amplitudes of transverse wave #1 before
%                     splitting applied. 
%     T90 (vector)  : Series of amplitudes of transverse wave #2 before
%                     splitting applied.  
%     fast (scalar) : polarisation direction of fast wave (degrees)
%     tlag (scalar) : lag time splitting operator to be applied (s)
%
%  Outputs:
%     T00sp (vector)  : Series of amplitudes of transverse wave #1 after
%                       splitting applied.
%     T90sp (vector)  : Series of amplitudes of transverse wave #2 after
%                       splitting applied.
%
%  MS_split_trace applies a shear wave splitting operator (fast and tlag) 
%  to a time series of the amplitude of displacement in two orthogonal 
%  directions (e.g. of waves or wavelets) stored in time, T00 and T90. The
%  displacements are first rotated into the fast direction before being time
%  shifted by forwards (for the fast wave) or backwards (for the slow wave)
%  by half of the lag time and interpolated onto the same time series. The
%  resulting displacements are then rotated back into the original frame. 
%  The fast direction is measured clockwise from the direction of T00
%  with a 90 degree fast orientation being parallel with T90. 

% Copyright (c) 2013, James Wookey and Andrew Walker
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


function [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast,tlag)
    % Apply a splitting operator to a pair of orthogonal traces.
   
    % rotate to fast reference frame
    [ F, S ] = RF_rotate(T00,T90,fast) ;
   
    % shift both components tlag/2 in opposite directions
    F = tshift(time,F,+tlag/2) ;
    S = tshift(time,S,-tlag/2) ;
   
    % rotate back to original reference frame
    [ T00sp, T90sp ] = RF_rotate(F,S,-fast) ;
end

function [As] = tshift(time,A,tshift)
    % Apply an interpolated time shift to a trace
    As = pchip(time,A,time+tshift) ;
end

function [T00R,T90R] = RF_rotate(T00,T90,theta)
    % Apply a rotation to a pair of orthogonal traces (this is a reference 
    % frame rotation)

    % form 2 row matrix
    D = [T90 ; T00] ;
   
    % rotation matrix
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)] ;
   
    DR = R*D ;
   
    T00R = DR(2,:) ;
    T90R = DR(1,:) ; 
end
