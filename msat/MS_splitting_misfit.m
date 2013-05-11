% MS_SPLITTING_MISFIT - Calculate the misfit between a pair of splitting
%                       operators.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Usage: 
%     [misfit] = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,spol,dfreq)         
%         Calculate the misfit between splitting operators [fast1,tlag1]
%         and [fast2,tlag2], using the default mode (see below). spol is
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
   case 'lam2s'
      misfit=MS_splitting_misfit_lam2S(fast1,tlag1,fast2,tlag2,spol,dfreq) ;
   otherwise
      error('MS:SPLITTING_MISFIT:UnknownMode',...
         ['Unknown mode: ' modeStr]) ;      
   end
   
end
%===============================================================================

%===============================================================================
function [misfit] = MS_splitting_misfit_lam2(fast1,tlag1,fast2,tlag2,spol,dfreq)
%===============================================================================

   % generate a test wavelet.
   
   
   [time,T00,T90] = MS_make_trace(spol,dfreq,max([tlag1 tlag2])) ;

   % apply splitting operator 1 
   [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast1,tlag1) ;
   
   %plot_wavelet(time,T00sp,T90sp)
   
   % apply inverse of splitting operator 2
   [T00sp2,T90sp2] = MS_split_trace(time,T00sp,T90sp,fast2,-tlag2) ;   

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
function [misfit] = MS_splitting_misfit_lam2S(fast1,tlag1,fast2,tlag2,spol,dfreq)
%===============================================================================

% foward 

   % generate a test wavelet.
   [time,T00,T90] = MS_make_trace(spol,dfreq,max([tlag1 tlag2])) ;

   % apply splitting operator 1 
   [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast1,tlag1) ;
      
   % apply inverse of splitting operator 2
   [T00sp2,T90sp2] = MS_split_trace(time,T00sp,T90sp,fast2,-tlag2) ;   
   
   %% measure second eignvalue of the resulting wavelet
   % calculate the covariance matrix
   COVM = cov(T90sp2,T00sp2) ;

   % take the eigenvalues   
   [V,D] = eig(COVM) ;
      
   % calculate normalised misfit.   
   misfit1 = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;   

   %% reverse
   
   [time,T00,T90] = MS_make_trace(spol,dfreq,max([tlag1 tlag2])) ;

   % apply splitting operator 2 
   [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast2,tlag2) ;
      
   % apply inverse of splitting operator 1
   [T00sp2,T90sp2] = MS_split_trace(time,T00sp,T90sp,fast1,-tlag1) ;   
   
   %% measure second eignvalue of the resulting wavelet
   % calculate the covariance matrix
   COVM = cov(T90sp2,T00sp2) ;

   % take the eigenvalues   
   [V,D] = eig(COVM) ;
      
   % calculate normalised misfit.   
   misfit2 = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;     
   
%  result is the average of the two
   misfit = (misfit1+misfit2)/2 ;

   
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
