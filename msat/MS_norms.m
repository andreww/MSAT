% MS_NORMS - Browaeys and Chevrot analysis of the elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  Calculate the percentages of the original matrix norm represented by the
%  individual decompositions of an elasticity tensor C, after: 
%     Browaeys and Chevrot (GJI, v159, 667-678, 2004)
%  
%  [P] = MS_norms(C,Ciso)
%     Isotropic part (scalar)
%
%  [P] = MS_norms(C,Ciso,Chex)
%     Isotropic, and hexagonal parts of the elastic tensor (vector)
%
%  [P] = MS_norms(C,Ciso,Chex,Ctet,Cort,Cmon,Ctri)
%     All parts of the elastic tensor (vector)
%     
%
% Notes:
%     Input components of of C must first be calculated using MS_decomp.
%     The output vector P will have as may elements as the number of
%     decomposed parts of the input matrix that are provided.
%
% References:
%     Browaeys, J. T. and S. Chevrot (2004) Decomposition of the elastic
%         tensor and geophysical applications. Geophysical Journal 
%         international v159, 667-678.
%
% See also: MS_AXES, MS_DECOMP

% Copyright (c) 2011-2015 James Wookey and Andrew Walker
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
   
function [ P ] = MS_norms( varargin )

    Cref = varargin{1} ;
    [Xref]=C2X(Cref) ;

    N = sqrt(Xref'*Xref) ;

    Ctot=zeros(6,6);
    P(1:nargin-1) = 0.0;

    for i=2:nargin
        Ctot = Ctot + varargin{i} ;
        [Xtot]=C2X(Ctot) ;
        XD=Xref-Xtot;
        P(i-1)=1-(sqrt(XD'*XD)/N) ;
    end

    % transform P
    PP = [0 P] ;
    P(1:nargin-1) = PP(2:nargin)-PP(1:nargin-1) ;

end

function [X]=C2X(C)
%  after Browaeys and Chevrot (GJI, 2004)
	X = zeros(21,1) ;
	
	X(1)  = C(1,1) ;
	X(2)  = C(2,2) ;
	X(3)  = C(3,3) ;
	X(4)  = sqrt(2).*C(2,3) ;
	X(5)  = sqrt(2).*C(1,3) ;
	X(6)  = sqrt(2).*C(1,2) ;
	X(7)  = 2.*C(4,4) ;
	X(8)  = 2.*C(5,5) ;
	X(9)  = 2.*C(6,6) ;
	X(10) = 2.*C(1,4) ;
	X(11) = 2.*C(2,5) ;
	X(12) = 2.*C(3,6) ;
	X(13) = 2.*C(3,4) ;
	X(14) = 2.*C(1,5) ;
	X(15) = 2.*C(2,6) ;
	X(16) = 2.*C(2,4) ;
	X(17) = 2.*C(3,5) ;
	X(18) = 2.*C(1,6) ;
	X(19) = 2.*sqrt(2).*C(5,6) ;
	X(20) = 2.*sqrt(2).*C(4,6) ;
	X(21) = 2.*sqrt(2).*C(4,5) ;
	
end
