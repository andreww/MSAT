% MS_CIJ2CIJKL - Convert from Voigt elasticity matrix to tensor
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Converts between a a 6x6 Voigt representation and a 3x3x3x3 tensor 
%     representation of anisotropic elasticity. 
%
%  % [CC] = MS_cij2cijkl(C)
%
% Usage: 
%     C must be a rank 2 array with size (6,6). The returned array C will 
%     be rank 4 with size (3,3,3,3). 
%
% Notes:
%     Do not use this function for the elastic compliance as additional
%     terms are needed in this case.
%
% See also: MS_cijkl2cij

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

% 2005/07/04 - fixed Vera Schulte-Pelkum's bug

function [CC] = MS_cij2cijkl(C)
 CC = zeros(3,3,3,3);
 CC(1,1,1,1) = C(1,1)         ;
 CC(2,2,2,2) = C(2,2)         ;
 CC(3,3,3,3) = C(3,3)         ;
 CC(2,3,2,3) = C(4,4)         ;
 CC(3,2,3,2) = CC(2,3,2,3)     ;
 CC(2,3,3,2) = CC(2,3,2,3)     ;
 CC(3,2,2,3) = CC(2,3,2,3)     ;
 CC(1,3,1,3) = C(5,5)         ;
 CC(3,1,1,3) = CC(1,3,1,3)     ;
 CC(1,3,3,1) =CC(1,3,1,3)     ;
 CC(3,1,3,1) =CC(1,3,1,3)     ;
 CC(1,1,2,2) = C(1,2)         ;
 CC(2,2,1,1) =CC(1,1,2,2)     ;
 CC(1,1,3,3) = C(1,3)         ;
 CC(3,3,1,1) =CC(1,1,3,3)     ;
 CC(1,1,2,3) = C(1,4)         ;
 CC(1,1,3,2) =CC(1,1,2,3)     ;
 CC(2,3,1,1) =CC(1,1,2,3)     ;
 CC(3,2,1,1) =CC(1,1,2,3)     ;
 CC(1,1,1,3) = C(1,5)         ;
 CC(1,1,3,1) =CC(1,1,1,3)     ;
 CC(1,3,1,1) =CC(1,1,1,3)     ;
 CC(3,1,1,1) =CC(1,1,1,3)     ;
 CC(1,1,1,2) = C(1,6)         ;
 CC(1,1,2,1) =CC(1,1,1,2)     ;
 CC(1,2,1,1) =CC(1,1,1,2)     ;
 CC(2,1,1,1) =CC(1,1,1,2)     ;
 CC(2,2,3,3) = C(2,3)         ;
 CC(3,3,2,2) =CC(2,2,3,3)     ;
 CC(2,2,2,3) = C(2,4)         ;
 CC(2,2,3,2) =CC(2,2,2,3)     ;
 CC(2,3,2,2) =CC(2,2,2,3)     ;
 CC(3,2,2,2) =CC(2,2,2,3)     ;
 CC(2,2,1,3) = C(2,5)         ;
 CC(2,2,3,1) =CC(2,2,1,3)     ;
 CC(1,3,2,2) =CC(2,2,1,3)     ;
 CC(3,1,2,2) =CC(2,2,1,3)     ;
 CC(2,2,1,2) = C(2,6)         ;
 CC(2,2,2,1) =CC(2,2,1,2)     ;
 CC(1,2,2,2) =CC(2,2,1,2)     ;
 CC(2,1,2,2) =CC(2,2,1,2)     ;
 CC(3,3,2,3) = C(3,4)         ;
 CC(3,3,3,2) = CC(3,3,2,3)    ;
 CC(2,3,3,3) = CC(3,3,2,3)    ;
 CC(3,2,3,3) = CC(3,3,2,3)    ;
 CC(3,3,1,3) = C(3,5)         ;
 CC(3,3,3,1) = CC(3,3,1,3)    ;
 CC(1,3,3,3) = CC(3,3,1,3)    ;
 CC(3,1,3,3) = CC(3,3,1,3)    ;
 CC(3,3,1,2) = C(3,6)         ;
 CC(3,3,2,1) = CC(3,3,1,2)    ;
 CC(1,2,3,3) = CC(3,3,1,2)    ;
 CC(2,1,3,3) = CC(3,3,1,2)    ;
 CC(2,3,1,3) = C(4,5)         ;
 CC(3,2,1,3) =CC(2,3,1,3)     ;
 CC(1,3,3,2) =CC(2,3,1,3)     ;
 CC(1,3,2,3) =CC(2,3,1,3)     ;
 CC(2,3,3,1) =CC(2,3,1,3)     ;
 CC(3,2,3,1) =CC(2,3,1,3)     ;
 CC(3,1,2,3) =CC(2,3,1,3)     ;
 CC(3,1,3,2) =CC(2,3,1,3)     ;
 CC(2,3,1,2) = C(4,6)         ;
 CC(3,2,1,2) =CC(2,3,1,2)     ;
 CC(1,2,2,3) =CC(2,3,1,2)     ;
 CC(1,2,3,2) =CC(2,3,1,2)     ;
 CC(2,3,2,1) =CC(2,3,1,2)     ;
 CC(3,2,2,1) =CC(2,3,1,2)     ;
 CC(2,1,2,3) =CC(2,3,1,2)     ;
 CC(2,1,3,2) =CC(2,3,1,2)     ;
 CC(1,3,1,2) = C(5,6)         ;
 CC(3,1,1,2) =CC(1,3,1,2)     ;
 CC(1,2,1,3) =CC(1,3,1,2)     ;
 CC(1,2,3,1) =CC(1,3,1,2)     ;
 CC(1,3,2,1) =CC(1,3,1,2)     ;
 CC(3,1,2,1) =CC(1,3,1,2)     ;
 CC(2,1,1,3) =CC(1,3,1,2)     ;
 CC(2,1,3,1) =CC(1,3,1,2)     ;
 CC(1,2,1,2) = C(6,6)         ;
 CC(2,1,1,2) =CC(1,2,1,2)     ;
 CC(1,2,2,1) =CC(1,2,1,2)     ;
 CC(2,1,2,1) =CC(1,2,1,2)     ;
return


