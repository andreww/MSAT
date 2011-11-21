% MS_ROTEULER - Rotate an elasticity matrix using Bunge's Euler angles.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotate one or more elasticity matrices usign Bunge's Euler angles 
%
%  % [ CC ] = MS_rotEuler( C, phi1, theta, phi2 )
%
% Usage: 
%     For three angles in degrees (phi1, theta, phi2) representing a
%     rotation in Bunge convention rotate the elasticity matrix C. All
%     variables can be arrays or scalars but (1) the angles must all be the
%     same length and (2) either C or all the angles must be scalars unless
%     they are the same length.
%
% Notes:
%    In the context of this function the term 'Euler angle' is used in 
%    the sense common in texture analysis not computer graphics etc. where
%    the three angles are rotations around a fixes axis system (see 
%    MS_rot3 for this case). In the (Bunge) convention used here the object
%    (crystal) reference frame is first assumed to be alligned with the 
%    sample reference frame. The object is then rotated around the common
%    z axis by an angle phi1. The second rotation (theta) is about 
%    the object's x axis (which is no longer parallel to the external 
%    x axis). The third rotation, phi2, is about the samples z axis 
%    (which, by this point, is not paralell with the external z axis). See
%    Figure 5 of Bunge 1985 for a graphical example.
%
% References:
%    Bunge, H. J. (1985) "Representaton of Preferred Orientations" in H.-R.
%    Wenk (ed.) "Preferred Orientation in Deformed Metals and Rocks: 
%     An Introduction to Modern Texture Analysis" Academic Press inc.
%     Orlando.
%
% See also: MS_ROT3 MS_ROTR

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

function [ CC ] = MS_rotEuler( C, phi1, theta, phi2 )

    % How many Euler angles - are lists correct length. 
    numEs = size(phi1, 1);
    assert(numEs == size(theta, 1), 'MS:ListsMustMatch',...
        'The number of phi1 angles must be the same as the number of theta angles')
    assert(numEs == size(phi2, 1), 'MS:ListsMustMatch',...
        'The number of phi1 angles must be the same as the number of phi2 angles')
    
    % Convert all angles to radians
    phi1r = phi1*(pi/180.0);
    thetar = theta*(pi/180.0);
    phi2r = phi2*(pi/180.0);
    
    numdimCs = ndims(C);
    if ((numEs == 1) && (numdimCs == 2))
        % Two matrix scalar case
        CC = rot_Euler_scalar(C, phi1r, thetar, phi2r);
    elseif ((numEs > 1) && (numdimCs == 2))
        % Many rotations and one matrix case
        CC = zeros(6,6,numEs);
        for i = 1:numEs
            CC(:,:,i) = rot_Euler_scalar(C, phi1r(i), thetar(i), phi2r(i));
        end
    elseif ((numEs == 1) && (numdimCs == 3))
        % Many matrix and one rotation case
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C, 3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r, thetar, phi2r);
        end
    elseif ((numEs > 1) && (numdimCs == 3))
        % List of rotatiosn and matrices case
        assert((numEs == size(C, 3)), 'MS:ListsMustMatch', ...
            'The length of the lists of rotations and matrices must match');
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C,3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r(i), thetar(i), phi2r(i));
        end
    end
    
    
end


function [ CC ] = rot_Euler_scalar(C, phi1, theta, phi2)
    
    %FIXME - delegate this text to MS_rotR?
    assert(MS_checkC(C)==1, 'MS:Cinvalid', ...
        'Argument "C" must be a valid 6x6 elasticity matrix');
    
    % Pre-compute trig functions.
    cp1 = cos(phi1);
    sp1 = sin(phi1);
    cp2 = cos(phi2);
    sp2 = sin(phi2);
    cth = cos(theta);
    sth = sin(theta);

    % Form rotation matrix for Bunge convention of 
    % Euler angles, see eq 24 (pg81) of "Preferred
    % orientation in deformed meals and rocks... H-K Wenk (ed)"
    g = [cp1*cp2 - sp1*sp2*cth, sp1*cp2 + cp1*sp2*cth, sp2*sth ; ...
        -1*cp1*sp2 - sp1*cp2*cth, -1*sp1*sp2 + cp1*cp2*cth, cp2*sth; ...
        sp1*sth, -1*cp1*sth, cth];
    
    
    CC = MS_rotR(C, g);
   
end
