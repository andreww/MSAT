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

   % default mode
   modeStr = 'lam2' ;
   
   % process the optional arguments
   iarg = 1 ;
   while iarg <= (length(varargin))
      switch lower(varargin{iarg})
      case 'mode'  % parameter definition (value required)
         modeStr = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      otherwise 
         error('MS:SPLITTING_MISFIT:UnknownOption',...
            ['Unknown option: ' varargin{iarg}]) ;   
      end
   end

   % call the appropriate mode function
   switch lower(modeStr)
   case 'lam2'
      misfit=MS_splitting_misfit_lam2(fast1,tlag1,fast2,tlag2,spol,dfreq) ;
   otherwise
      error('MS:SPLITTING_MISFIT:UnknownMode',...
         ['Unknown mode: ' varargin{iarg}]) ;      
   end
   
end
%===============================================================================

%===============================================================================
function [misfit] = MS_splitting_misfit_lam2(fast1,tlag1,fast2,tlag2,spol,dfreq)
%===============================================================================

   % generate a test wavelet.
   
   
   [time,T00,T90] = FDGaussian_wavelet(spol,dfreq,max([tlag1 tlag2])) ;

   % apply splitting operator 1 
   [T00sp,T90sp] = split(time,T00,T90,fast1,tlag1) ;
   
   %plot_wavelet(time,T00sp,T90sp)
   
   % apply inverse of splitting operator 2
   [T00sp2,T90sp2] = split(time,T00sp,T90sp,fast2,-tlag2) ;   

   %plot_wavelet(time,T00sp2,T90sp2)
   
   
   %% measure second eignvalue of the resulting wavelet
   % calculate the covariance matrix
   COVM = cov(T90sp2,T00sp2) ;

   % take the eigenvalues   
   [V,D] = eig(COVM) ;
      
   % calculate normalised misfit.   
   misfit = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;   
   
end
%===============================================================================

%===============================================================================
function [T00sp,T90sp] = split(time,T00,T90,fast,tlag)
%===============================================================================
% Apply a splitting operator to a pair of orthogonal traces. 
   
   % rotate to fast reference frame
   [ F, S ] = RF_rotate(T00,T90,fast) ;
   
   % shift both components tlag/2 in opposite directions
   F = tshift(time,F,+tlag/2) ;
   S = tshift(time,S,-tlag/2) ;
   
   % rotate back to original reference frame
   [ T00sp, T90sp ] = RF_rotate(F,S,-fast) ;
   
end
%===============================================================================

%===============================================================================
function [As] = tshift(time,A,tshift)
%===============================================================================
% Apply an interpolated time shift to a trace
   
   As = pchip(time,A,time+tshift) ;
   
end
%===============================================================================

%===============================================================================
function [T00R,T90R] = RF_rotate(T00,T90,theta)
%===============================================================================
% Apply a rotation to a pair of orthogonal traces (this is a reference frame
% rotation)

   % form 2 row matrix 
   D = [T90 ; T00] ;
   
   % rotation matrix 
   R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)] ;
   
   DR = R*D ;
   
   T00R = DR(2,:) ;
   T90R = DR(1,:) ;
   
end
%===============================================================================

%===============================================================================
function [time,amp0,amp90] = FDGaussian_wavelet(spol,dfreq,max_tlag)
%===============================================================================
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
   
   %  normalise and project amplitude
   amp = amp./max(amp) ;  
   amp0 = amp.*cosd(spol) ;
   amp90 = amp.*sind(spol) ;
   
   % done
end
%===============================================================================

%===============================================================================
function [] = plot_wavelet(time,T00,T90)
%===============================================================================
% generate a (first-derivative) Gaussian wavelet centred on 0, time base is 
% set so the maximum tlag can be accomodated. This is defined as
% -(max_tlag+2*T):T/100:(max_tlag+2*T) ;

figure

plot(time,T00,'b-')
hold on
plot(time,T90,'r-')


end
%===============================================================================
