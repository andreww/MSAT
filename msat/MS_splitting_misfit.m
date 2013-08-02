% MS_SPLITTING_MISFIT - Calculate the misfit between a pair of splitting
%                       operators.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Usage: 
%     [misfit] = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,spol,dfreq)         
%        Calculate the misfit between splitting operators [fast1,tlag1]
%        and [fast2,tlag2], using the default mode (see below). SPOL is
%        the initial source polarisation, dfreq is the dominant frequency.
%
%
% Modes of operation:
%     [misfit] = MS_splitting_misfit(...,'mode','lam2')
%        Calculate misfit based on the residual second eigenvalue of applying
%        [f1,t1] to a test wavelet, then removing [f2,t2]. [DEFAULT].
%
%     [misfit] = MS_splitting_misfit(...,'mode','lam2S')
%        Calculate misfit based on the residual second eigenvalue of applying
%        [f1,t1] to a test wavelet, then removing [f2,t2]. This is then 
%        reversed, and the final misfit is the average of the two. 
%
%     [misfit] = MS_splitting_misfit(...,'mode','intensity')
%        Calculates the misfit as the difference between the splitting intensity
%        (see, e.g., Long and Silver, 2009) for the two operators. 
%
%     [misfit] = MS_splitting_misfit(...,'mode','simple')
%        A simple (normalised) arithmetic misfit. 
%
%  References: 
%      Long, M. D., & Silver, P. G. (2009). "Shear Wave Splitting and Mantle 
%        Anisotropy: Measurements, Interpretations, and New Directions". Surveys 
%        in Geophysics, 30(4-5), 407–461. doi:10.1007/s10712-009-9075-1
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
function [varargout] = MS_splitting_misfit(fast1,tlag1,...
                                        fast2,tlag2,spol,dfreq,varargin)
%===============================================================================

   % default mode
   modeStr = 'lam2' ;
   idebug = 0 ;
   max_tlag = NaN ; % autoselect max tlag
   
   % process the optional arguments
   iarg = 1 ;
   while iarg <= (length(varargin))
      switch lower(varargin{iarg})
      case 'debug'  % parameter definition (value required)
         idebug = 1 ;
         iarg = iarg + 1 ;
      case 'mode'  % parameter definition (value required)
         modeStr = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'max_tlag'  % parameter definition (value required)
         max_tlag = varargin{iarg+1} ;
         iarg = iarg + 2 ;         
      otherwise 
         error('MS:SPLITTING_MISFIT:UnknownOption',...
            ['Unknown option: ' varargin{iarg}]) ;   
      end
   end
   
   % check outputs
   if nargout>1
      error('MS:SPLITTING_MISFIT:TooManyOutputs','Too many outputs specified')
   end

   % call the appropriate mode function
   switch lower(modeStr)
   case 'simple'
       if isscalar(fast1)
           misfit = MS_splitting_misfit_simple(...
               fast1,tlag1,fast2,tlag2,spol,dfreq) ;
       else
           misfit = 0.0;
           for i = 1:length(fast1)
               misfit = misfit + MS_splitting_misfit_simple(...
                   fast1(i),tlag1(i),fast2(i),tlag2(i),spol,dfreq) ; 
           end
       end
   case 'intensity'
      [misfit]=MS_splitting_misfit_intensity(fast1,tlag1,fast2,tlag2,spol,dfreq) ;      
   case 'lam2'
      [misfit]=MS_splitting_misfit_lam2( ...
         fast1,tlag1,fast2,tlag2,spol,dfreq,max_tlag,idebug) ;
   case 'lam2s'
      [misfit]=MS_splitting_misfit_lam2S( ...
         fast1,tlag1,fast2,tlag2,spol,dfreq,max_tlag,idebug) ;
   otherwise
      error('MS:SPLITTING_MISFIT:UnknownMode',...
         ['Unknown mode: ' modeStr]) ;      
   end
   
   % assign outputs
   switch nargout
      case 1
         varargout = {misfit} ;
   end
      
   
end
%===============================================================================

%===============================================================================
function [misfit] = ...
   MS_splitting_misfit_simple(fast1,tlag1,fast2,tlag2,spol,dfreq)
%===============================================================================
% Simple numerical misfit. Average of normalised fast angular difference and 
% time lag. TLAG misfit is normalised by 1/4 of the dominant period. Also,
% apply a weighting factor for the polarisation.

   fast_misfit = abs(MS_unwind_pm_90(fast1-fast2))./90 ;
   tlag_misfit = abs(tlag1-tlag2)./(0.25/dfreq) ;
   
   misfit = mean([fast_misfit tlag_misfit]) ;
      
end
%===============================================================================

%===============================================================================
function [misfit] = ...
   MS_splitting_misfit_intensity(fast1,tlag1,fast2,tlag2,spol,dfreq)
%===============================================================================
% Calculate the misfit between the splitting intensities predicted by the 
% two [fast,tlag] operators combined with the initial source polarisations. 

   beta1 = abs(MS_unwind_pm_90(fast1-spol)) ;
   S1 = tlag1 * sind(2*beta1) ;
   
   beta2 = abs(MS_unwind_pm_90(fast2-spol)) ;
   S2 = tlag2 * sind(2*beta2) ;
 
   misfit = abs(S1-S2) ; 
      
end
%===============================================================================

%===============================================================================
function [misfit] = ...
   MS_splitting_misfit_lam2(fast1,tlag1,fast2,tlag2,spol,dfreq,max_tlag,idebug)
%===============================================================================

   % generate a test wavelet.
   if isnan(max_tlag), max_tlag = max([tlag1 tlag2]); , end
   
   [time,T00,T90] = MS_make_trace(spol,dfreq,max_tlag) ;

   if (idebug), plot_wavelet(time,T00,T90), end

   % apply splitting operator 1 
   [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast1,tlag1) ;

   
   if (idebug), plot_wavelet(time,T00sp,T90sp), end
   
   % apply inverse of splitting operator 2
   [T00sp2,T90sp2] = MS_split_trace(time,T00sp,T90sp,fast2,-tlag2) ;   

   if (idebug), plot_wavelet(time,T00sp2,T90sp2), end
   
   
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
function [misfit] = ...
   MS_splitting_misfit_lam2S(fast1,tlag1,fast2,tlag2,spol,dfreq,max_tlag,idebug)
%===============================================================================

% foward 

   % generate a test wavelet.
   if isnan(max_tlag), max_tlag = max([tlag1 tlag2]); , end

   [time,T00,T90] = MS_make_trace(spol,dfreq,max_tlag) ;

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
   
   [time,T00,T90] = MS_make_trace(spol,dfreq,max_tlag) ;

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

