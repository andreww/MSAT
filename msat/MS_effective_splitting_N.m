% MS_EFFECTIVE_SPLITTING_N - N-layer effective splitting operator calculation.
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [fast_eff,tlag_eff]=MS_effective_splitting_N(f,spol,fast,tlag)
%  
%  Inputs:
%     f (scalar) : Dominant frequency of wave 
%     spol (scalar) : initial source polarisation
%     fast (scalar/vector) : fast direction(s) of layer(s) to be included
%     tlag (scalar/vector) : lag time(s) of layer(s) to be included
%
%  Fast and tlag can be scalars or vectors but must all be the same length.  
%
%  CAUTION! The method gives unstable results when
%  the source polarisation is near the effective fast direction. This would,
%  however probably be seen as a null result anyway.
%
%  Reference: 
%      Silver, P. G. and Savage, M. K. 1994 "The interpretation of shear-wave
%      splitting parameters in the presence of two anisotropic layers" 
%      Geophysical Journal International, v119 pp 949-963.  

% Copyright (c) 2011, James Wookey and Andrew Walker
% Copyright (c) 2005-2011, James Wookey
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
function [fast_eff,tlag_eff]=MS_effective_splitting_N(f,spol,fast,tlag)
%===============================================================================
%  check inputs are the same size   
   if ~isequal(length(fast),length(tlag))
      error('MS:ListsMustMatch', ...
          'Input fast and tlag vectors must be of the same length')            
   end

%  check for just one layer
   if length(fast)==1
      fast_eff = fast ;
      tlag_eff = tlag ;
      return
   end

%  aggregate layers which differ by 0 or 90 (+/- 1 degree)
   [ fast , tlag ] = aggregate( fast , tlag, 'threshold' , 1.0 ) ;

%  check (again) for just one layer
   if length(fast)==1
      fast_eff = fast ;
      tlag_eff = tlag ;
      return
   end

%  remove any layers with 0 tlag time; these break the calculation
   ind = find(tlag~=0.0) ;
   fast2 = fast(ind) ;
   tlag2 = tlag(ind) ;
   fast = fast2 ;
   tlag = tlag2 ;

%  unwind the fast directions
   fast = MS_unwind_pm_90(fast) ;
   
%**process   
   w = 2.0 .* pi .* f ;
   th = w .* tlag ./ 2. ;
   al = 2.*(fast-spol) .* pi/180.0 ;
         
   S = prod(cos(th)) ;
   Cc = S.*sum((tan(th).*cos(al))) ;
   Cs = S.*sum((tan(th).*sin(al))) ;
   
   n=length(fast) ;
   
   ap = 0;
   app = 0;
   
   for i=1:n-1
      for j=i+1:n
         ap = ap + (tan(th(i)).*tan(th(j)).*cos(al(i)-al(j))) ;
         app = app + (tan(th(i)).*tan(th(j)).*sin(al(i)-al(j))) ;
      end
   end
   
   ap = S.*(1-ap) ;
   app = S.*app ;
   ala = atan ( (app.^2.+Cs.^2.) ./ (app.*ap + Cs.*Cc) ) ;
   tha = atan ( (app) ./ (Cs.*cos(ala)-Cc.*sin(ala)) ) ;
   
   fast_eff = spol + (ala.*180./pi) ./ 2. ;
   tlag_eff = 2.*tha./w ;                                     
   
   fast_eff = MS_unwind_pm_90(fast_eff) ;

%**handle zero tlags      
   if (tlag_eff < 0)
      fast_eff = MS_unwind_pm_90(fast_eff+90) ;
      tlag_eff = abs(tlag_eff) ;
   end   

end
%===============================================================================

%===============================================================================
function [ fastA , tlagA ] = aggregate( fast , tlag, varargin )
%  Agglomerate splitting operators which are (near) parallel or
%  perpendicular. 
%
   
   threshold = 1 ;
   
%  process the optional arguments
   iarg = 1 ;
   while iarg <= (length(varargin))
      switch lower(varargin{iarg})
                case 'threshold'  % parameter definition (value required)
                  threshold = varargin{iarg+1} ;
                  iarg = iarg + 2 ;
                otherwise 
                  error('MS:EFFECTIVE_SPLITTING_N:UnknownOption',...
                     ['Unknown option: ' varargin{iarg}]) ;   
      end
   end   

   % new arrays
   fast_agg = zeros(size(fast))+NaN ;
   tlag_agg = zeros(size(fast))+NaN ;
   
   N = length(fast) ;
   
   iRef = 1 ;
   iAgg = 1 ;
   iComp = 2 ;
   
   fast_agg(1) = fast(1) ;
   tlag_agg(1) = tlag(1) ;
   
   % process the fast directions. If the reference direction and the 
   % comparison  direction are (nearly) identical, sum the tlags. If they 
   % differ by (nearly) 90 degrees, take the difference. 
   while iComp<=N
      if abs(fast(iRef)-fast(iComp))<threshold
         % add the tlag
         tlag_agg(iAgg) = tlag_agg(iAgg) + tlag(iComp) ;
         iComp = iComp + 1 ;
      elseif abs((abs(fast(iRef)-fast(iComp))-90))<threshold
         % subtract the tlag
         tlag_agg(iAgg) = tlag_agg(iAgg) - tlag(iComp) ;
         iComp = iComp + 1 ;
      else
         % move onto the next element
         iRef = iComp ;
         iComp = iComp + 1 ;
         iAgg = iAgg + 1;
         fast_agg(iAgg) = fast(iRef) ;
         tlag_agg(iAgg) = tlag(iRef) ;
      end
   end
   
   % remove the NaNs.
   ind = find(~isnan(fast_agg)) ;
   fastA = fast_agg(ind) ;
   tlagA = tlag_agg(ind) ;
   
   % finally, see if we have ended up with any negative numbers. 
   % if so, flip the sign and add 90. 
   ind = find(tlagA<0) ;
   tlagA(ind) = abs(tlagA(ind)) ;
   fastA(ind) = fastA(ind)+90 ; 
   
end
