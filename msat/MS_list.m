% MS_LIST - Print elasticity matrix (in 'list' format)
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% CIJ_list(C,rho)
%
% Usage: 
%     CIJ_list(C,rho)
%         For a given 6x6 elasticity matrix, C, and density, rho, print 
%         the constants.
%
%  Notes: 
%    The 'list' format consists of three numbers per line. The first two 
%    are integers representing the location in the elasticity matrix in
%    Voigt notation, the third is the value of the elastic constant. After
%    the 36 lines of elastic constants a line containing the density is
%    printed with dummy indices '7 7'. This function exists to allow 
%    MSAT to be interfaced with legacy code.
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

function MS_list(C,rho)

    for i=1:6
       for j=i:6
          fprintf('%1i %1i %e\n',i,j,C(i,j))
       end
    end   
    fprintf('%1i %1i %f\n',7,7,rho)
end
