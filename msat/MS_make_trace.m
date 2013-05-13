% MS_MAKE_TRACE - Create a pair of orthogonal traces
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [time, T00, T90] = MS_make_trace(spol, dfreq, max_tlag)
%
%  Inputs:
%     dfreq (scalar)    : Dominant frequency of trace (Hz)  
%     spol (scalar)     : source polarisation direction (degrees)
%     max_tlag (scalar) : maximum lag time to be accommodated on trace (s)
%
%  Outputs:
%     time (vector) : Series of time points where amplitudes are known (s)
%     T00 (vector)  : Series of amplitudes of transverse wave #1. 
%     T90 (vector)  : Series of amplitudes of transverse wave #2.
%
%  MS_make_trace creates a time series of the amplitude of displacement in 
%  two orthogonal directions (e.g. of waves or wavelets) stored in time, 
%  T00 and T90. In the current version only a (first-derivative) Gaussian
%  wavelet can be generated. Future additional optional arguments will 
%  allow other choices.

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

function [time, T00, T90] = MS_make_trace(spol, dfreq, max_tlag)

    % Further options should be added here for other wavelets.

    [time, T00, T90] = FDGaussian_wavelet(spol,dfreq,max_tlag);

end


function [time,amp0,amp90] = FDGaussian_wavelet(spol,dfreq,max_tlag)
    % generate a (first-derivative) Gaussian wavelet centred on 0, time base is
    % set so the maximum tlag can be accomodated. This is defined as
    % -(max_tlag+2*T):T/100:(max_tlag+2*T) ;

    % calculate time base
    T = 1/dfreq ;
    dt = T/100 ;
    time = -(max_tlag+2*T):dt:(max_tlag+2*T) ;

    % calculate wavelet
    sig = T/6 ;
    amp = -(time./(sig.^3*sqrt(2.*pi))).*exp(-time.^2./(2.*sig.^2)) ;
   
    % normalise and project amplitude
    amp = amp./max(amp) ;
    amp0 = amp.*cosd(spol) ;
    amp90 = amp.*sind(spol) ;
end