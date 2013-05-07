% MS_EFFECTIVE_SPLITTING_N - N-layer effective splitting operator calculation.
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [fast_eff,tlag_eff]=MS_effective_splitting_N(f, spol, fast, tlag, ...)
%
%  Common Inputs:
%     f (scalar) : Dominant frequency of wave (Hz)
%     spol (scalar) : initial source polarisation (degrees)
%     fast (scalar/vector) : fast direction(s) of layer(s) to be 
%                            included (s)
%     tlag (scalar/vector) : lag time(s) of layer(s) to be included
%                            (degrees)
%
%  Fast and tlag can be scalars or vectors but must all be the same length.
% 
%  Usage:
%   [fast_eff,tlag_eff]=MS_effective_splitting_N(f, spol, fast, tlag)
%     or 
%   [fast_eff,tlag_eff]=MS_effective_splitting_N(f, spol, fast, tlag, ...
%   'mode', 'S&S')
%
%     Calculate effective splitting parameters using the method of Silver
%     and Savage (1994). CAUTION! This is quick but is known to fail if the delay
%     time becomes large compared to the dominant frequence of the wave. It
%     also gives unstable results when the source polarisation is near the 
%     effective fast direction. This would, however probably be seen as a 
%     null result anyway.
% 
%   [fast_eff,tlag_eff]=MS_effective_splitting_N(f, spol, fast, tlag, ...
%   'mode', 'GaussianWavelet')
%
%     Calculate effective splitting parameters by first applying the 
%     splitting operators in sequence to an unsplit Gaussian wavelet then
%     searching for the effective splitting operator that minimises the
%     second eigenvalue of the covariance matrix of the result. This is the 
%     ray theory approach described by Bonnin et al. (2012). 
%
%   [fast_eff,tlag_eff]=MS_effective_splitting_N(f, spol, fast, tlag, ...
%   'mode', 'GaussianWavelet', 'PlotWavelet')
%
%     Plot the resulting wavelet when in GaussianWavelet mode.
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
function [fast_eff,tlag_eff]=MS_effective_splitting_N(f,spol,fast,tlag, varargin)
%===============================================================================
%  check inputs are the same size   
   if ~isequal(length(fast),length(tlag))
      error('MS:ListsMustMatch', ...
          'Input fast and tlag vectors must be of the same length')            
   end
   
   mode = 's&s';
   plotwave = 0;
   
   iarg = 1 ;
   while iarg <= (length(varargin))
       switch lower(varargin{iarg})
           case 'mode'
              mode = varargin{iarg+1} ; 
              iarg = iarg + 2;
           case 'plotwavelet'
               plotwave = 1;
               iarg = iarg + 1;
       end
   end
   
   if strcmp(mode, 's&s') 
       [fast_eff,tlag_eff]=MS_effective_splitting_N_SS(f,spol,fast,tlag);
   elseif strcmp(mode, 'gaussianwavelet')
       [fast_eff,tlag_eff]=MS_effective_splitting_N_GW(f,spol,...
           fast,tlag,plotwave);
   else
       error('MS:effectiveUnkownMode', ...
           ['unknown mode in MS_effective_splitting_N: ', mode])
   end
end

function [fast_eff,tlag_eff]=MS_effective_splitting_N_SS(f,spol,fast,tlag)
   
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

function [fast_eff,tlag_eff]=MS_effective_splitting_N_GW(f, spol, fast, ...
    tlag, plotwave)

    [time,T00,T90] = FDGaussian_wavelet(spol,f, sum(tlag)) ;

    if length(fast)==1
        % So we can test - should get [fast_eff,tlag_eff] == [fast, tlag]
        [T00,T90] = split(time,T00,T90,fast,tlag);
    else
        for i =1:length(fast)
            % Split by each op in turn...
            [T00,T90] = split(time,T00,T90,fast(i),tlag(i));
        end
    end
    
    [fast_eff, tlag_eff] = measure_splitting_function(time, T00, T90, ...
        15, 0.5.*min(tlag), sum(tlag));
 
    fast_eff = MS_unwind_pm_90(fast_eff) ;
    if (tlag_eff < 0)
       fast_eff = MS_unwind_pm_90(fast_eff+90) ;
       tlag_eff = abs(tlag_eff) ;
    end   
    
    if plotwave
        plot_splitting(time, T00, T90, fast_eff, tlag_eff);
    end
end

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


function [fast, tlag] = measure_splitting_function(time, T00, T90, ...
    d_fast, d_tlag, max_tlag)
    % Given a two traces (T00 and T90) on a common time axis (time) find
    % the splitting operator (fast and tlag) that minimises the second
    % eigenvalue of the covariance matrix. This is done in two steps. First
    % using an approximate grid search (spacing and time limits set by 
    % d_fast, d_tlag and max_tlag) then by using the Matlab built
    % in simplex optimiser from the best grid point.

    % Do the grid search
    best_misfit = 100;
    for fast = 0.0:d_fast:180.0
        for tlag = 0.0:d_tlag:max_tlag
            split_op = [fast, tlag];
            misfit = splitting_function(split_op, time, T00, T90);
            if (misfit < best_misfit)
                best_split_op = split_op;
                best_misfit = misfit;
            end
        end
    end
    
    % Simplex minimisation
    [new_split_op, misfit] = fminsearch(...
                @(split_op) splitting_function(split_op, time, ...
                T00, T90), best_split_op, optimset('Display','notify'));
    if (misfit < best_misfit)
        best_split_op = new_split_op;
    end
            
    % Unpack results
    fast = best_split_op(1);
    tlag = best_split_op(2);        
end

function misfit = splitting_function(split_op, time, T00, T90)
    % Apply a spliting operator to two traces and calculate the misfit in a
    % way that can be used by measure_splitting_function
    
    % Unpack args
    fast = split_op(1);
    tlag = split_op(2);
    % Calculate the traces with the revered split
    [T00sp,T90sp] = split(time,T00,T90,fast,-tlag) ;
    % measure second eignvalue of the resulting wavelet
    % calculate the covariance matrix
    COVM = cov(T90sp,T00sp) ;
    % take the eigenvalues
    [~,D] = eig(COVM) ;
    % calculate normalised misfit.
    misfit = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;     
end

function [T00sp,T90sp] = split(time,T00,T90,fast,tlag)
    % Apply a splitting operator to a pair of orthogonal traces.
   
    % rotate to fast reference frame
    [ F, S ] = RF_rotate(T00,T90,fast) ;
   
    % shift both components tlag/2 in opposite directions
    F = tshift(time,F,+tlag/2) ;
    S = tshift(time,S,-tlag/2) ;
   
    % rotate back to original reference frame
    [ T00sp, T90sp ] = RF_rotate(F,S,-fast) ;
end

function [As] = tshift(time,A,tshift)
    % Apply an interpolated time shift to a trace
    As = pchip(time,A,time+tshift) ;
end

function [T00R,T90R] = RF_rotate(T00,T90,theta)
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


function [] = plot_splitting(time, T00, T90, fast, tlag)
    % Plot some possibly useful diagnostics;

    [T00sp,T90sp] = split(time,T00,T90,fast,-tlag) ;

    figure
    subplot(221);
    plot(time,T00,'b-')
    hold on
    plot(time,T90,'r-')
    title('Split wavelet')
    subplot(222);
    plot(time,T00sp,'b-')
    hold on
    plot(time,T90sp,'r-')
    title('After effective splitting removed')
    subplot(223);
    plot(T00, T90,'b-');
    subplot(224);
    plot(T00sp, T90sp,'b-');
  
end

 
