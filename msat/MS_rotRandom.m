% MS_ROTRANDOM - Randomly rotate an elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Generate one or more randomly oriented elasticity matrices given a single
% input matrix.
%
% Usage:
%
%    [ CC, ... ] = MS_rotRandom( C, ... )
%
%     Rotate the elasticity matrix C to a random new orientation. CC is a
%     6x6 matrix
%
%    [ CC, ... ] = MS_rotRandom( C, 'number', n )
%
%     CC is a 6x6xn matrix where each orientation is drawn from a uniform
%     random distribution of orientations such that, for large n, the 
%     avarage of the matricies in CC is isotropic.
%
%    [ CC, ... ] = MS_rotRandom( C, 'stream', PRNgenerator )
%
%     As above, but allow the user to specify the random number generator 
%     to be used.
%
%    [ CC, eulers ] = MS_rotRandom( C, ... )
%
%     As above, but also output a 3xn array of Euler angles (Bunge notation
%     and in degrees) used to generate the rotated elasticity matrices.
%
% See also: MS_ROT3 MS_ROTR MS_ROTEULER

% Copyright (c) 2014 James Wookey and Andrew Walker
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

function [CC, eulers] = MS_rotRandom( C, varargin )

    % Default PRN stream and number of output tensors
    % Note that for older versions of Matlab, this should be
    % stream = RandStream.getDefaultStream;
    stream = RandStream.getGlobalStream;
    number = 1;
    
    % Process the optional arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
       switch lower(varargin{iarg})
          case 'number' 
             number = varargin{iarg+1} ;
             iarg = iarg + 2 ;
          case 'stream'
             stream = varargin{iarg+1} ;
             iarg = iarg + 2 ;
           otherwise 
             error('MS:ROTRAND:UnknownOption',...
                ['Unknown option: ' varargin{iarg}]) ;   
       end   
    end
      
    % Uniform random numbers between 0 and 1. We convert these
    % into the Euler angles.
    eulers = rand(stream,3,number);

    % Convert phi1 and phi2 to degrees between 0 and 360. Note
    % that if Phi = 0 the orientation is determined by the sum of
    % these two rotations, but the periodicity means the orentation
    % is still uniformly distributed.
    eulers(1,:) = eulers(1,:).*360.0;
    eulers(3,:) = eulers(3,:).*360.0;
     
    % Phi can be between -90 and 90 degrees and cos(Phi) must
    % be unformly distributed (the metric space is phi1, cos(Phi)
    % phi2. So first make the number between -1 and 1 then take the
    % arccos in degrees.
    eulers(2,:) = (eulers(2,:).*2.0)-1.0;
    eulers(2,:) = acosd(eulers(2,:));

    % Annoying transpose - see GitHub issue #7
    CC = MS_rotEuler(C, eulers(1,:)', eulers(2,:)', eulers(3,:)');
 
end
