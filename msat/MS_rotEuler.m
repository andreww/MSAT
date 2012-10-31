% MS_ROTEULER - Rotate an elasticity matrix using Bunge's Euler angles.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotate one or more elasticity matrices using Bunge's Euler angles 
%
% Usage:
%
%    [ CC ] = MS_rotEuler( C, phi1, Phi, phi2 )
%
%     For three angles in degrees (phi1, theta, phi2) representing a
%     rotation in Bunge convention rotate the elasticity matrix C. All
%     variables can be arrays or scalars but (1) the angles must all be the
%     same length and (2) either C or all the angles must be scalars unless
%     they are the same length.
%
%    [ CC ] = MS_rotEuler( C, phi1, Phi, phi2, 'sense', 'active' )
%    [ CC ] = MS_rotEuler( C, phi1, Phi, phi2, 'sense', 'passive' )
%
%     As above, but allow user choice between 'active' (the default) or 
%     'passive' rotations (the usual description for texture analysis).
%     See notes, below, for details.
%
% Notes:
%    In the context of this function the term 'Euler angle' is used in 
%    the sense common in texture analysis not computer graphics etc. where
%    the three angles are rotations around a fixes axis system (see 
%    MS_rot3 for this case). In the (Bunge) convention used here the
%    relationship between the orientation of a set of "crystal" axes is 
%    considered with respect to an external "sample" axis system. This 
%    is described by considering a moving reference frame which starts out 
%    parallel to the sample axes. This moving reference frame is first 
%    rotated anticlockwise about its z axis by phi1 degrees. A second 
%    anticlockwise rotation of Phi degrees is made about the frame's x axis 
%    (which, at this point, is no longer parallel to the sample x axis).
%    Finally, a second anticlockwise rotation of phi2 degrees about the 
%    frame's z axis is made to bring the moving reference into line with
%    the crystal axes. See Figure 5 of Bunge 1985 for a graphical example.
% 
%    As described, the Euler angles are passive rotations showing how the
%    reference frame must be rotated. However, for most applications one
%    knows the elasticity matrix of the crystal on the crystal axis system
%    and wants to calculate the elasticity matrix represented on the global
%    "sample" reference frame. In this case the corresponding active
%    rotation must be applied to the crystal's elasticity matrix. By 
%    default, the MS_rotEuler function applies the active rotation.  
%    The optional arguments 'sense' and 'passive' can be used to reverse
%    the sense of rotation. This can be used if, for example, the
%    elasticity has been measured on an experimental frame which is offset
%    from the crystal frame more usually used.
%
% References:
%    Bunge, H. J. (1985) "Representation of Preferred Orientations" in H.-R.
%    Wenk (ed.) "Preferred Orientation in Deformed Metals and Rocks: 
%     An Introduction to Modern Texture Analysis" Academic Press inc.
%     Orlando.
%
% See also: MS_ROT3 MS_ROTR

% Copyright (c) 2011, 2012 James Wookey and Andrew Walker
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

function [ CC ] = MS_rotEuler( C, phi1, theta, phi2, varargin )

    transpose = 1;
    % Process the optional arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
       switch lower(varargin{iarg})
          case 'sense' 
             senseval = varargin{iarg+1} ;
             switch lower(senseval)
                 case 'active'
                     transpose = 1;
                 case 'passive'
                     transpose = 0;
                 otherwise
                error('MS:ROTE:UnknownOption',...
                ['Unknown option for sense: ' senseval]) ;
             end
             iarg = iarg + 2 ;
           otherwise 
             error('MS:ROTE:UnknownOption',...
                ['Unknown option: ' varargin{iarg}]) ;   
       end   
    end


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
        CC = rot_Euler_scalar(C, phi1r, thetar, phi2r, transpose);
    elseif ((numEs > 1) && (numdimCs == 2))
        % Many rotations and one matrix case
        CC = zeros(6,6,numEs);
        for i = 1:numEs
            CC(:,:,i) = rot_Euler_scalar(C, phi1r(i), thetar(i), ...
                phi2r(i), transpose);
        end
    elseif ((numEs == 1) && (numdimCs == 3))
        % Many matrix and one rotation case
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C, 3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r, thetar, ...
                phi2r, transpose);
        end
    elseif ((numEs > 1) && (numdimCs == 3))
        % List of rotatiosn and matrices case
        assert((numEs == size(C, 3)), 'MS:ListsMustMatch', ...
            'The length of the lists of rotations and matrices must match');
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C,3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r(i), thetar(i), ...
                phi2r(i), transpose);
        end
    end
    
    
end


function [ CC ] = rot_Euler_scalar(C, phi1, theta, phi2, transpose)
    
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
    % Note that this is for the passive roation.
    g = [cp1*cp2 - sp1*sp2*cth, sp1*cp2 + cp1*sp2*cth, sp2*sth ; ...
        -1*cp1*sp2 - sp1*cp2*cth, -1*sp1*sp2 + cp1*cp2*cth, cp2*sth; ...
        sp1*sth, -1*cp1*sth, cth];
    
    if transpose
        % Active rotation
        CC = MS_rotR(C, g');
    else
        % Passive rotation
        CC = MS_rotR(C, g);
    end
end
