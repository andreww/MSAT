% MS_SPLITTING_MISFIT - Calculate the misfit between a pair of splitting
%                       operators.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Usage: 
%     [misfit] = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,spol,dfreq)         
%         Calculate the misfit between splitting operators [fast1,tlag1]
%         and [fast2,tlag2], using the default mode (see below). SPOL is
%         the initial source polarisation, dfreq is the dominant frequency.
%
%
% Modes of operation:
%     [misfit] = MS_splitting_misfit(...,'mode','lam2')
%        Calculate misfit based on the residual second eigenvalue of applying
%        [f1,t1] to a test wavelet, then removing [f2,t2]. [DEFAULT].
%
%
% See also: MS_EFFECTIVE_SPLITTING_N

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

%===============================================================================
function [misfit] = MS_splitting_misfit(fast1,tlag1,...
                                        fast2,tlag2,spol,dfreq,varargin) ;
%===============================================================================


end