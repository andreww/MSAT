function test_suite = test_MS_splitting_misfit
   initTestSuite;
end

function test_MS_splitting_misfit_Zero_cases
   % check misfit between two identical splitting operators is zero for different
   % polarisations, frequencies and operators. 

   fast = [0 30 45 60 90]; 
   tlag = [0.003 0.01 0.4 3 20 ];
   spol = [0 60 135 150 180 ] ;
   freq = [100 10 1.0 0.1 0.01] ;
   
   for i=1:length(fast)
      misfit = MS_splitting_misfit(fast(i),tlag(i),fast(i),tlag(i),spol(i),freq(i)) ;
      assertElementsAlmostEqual(misfit, 0, 'absolute',0.001) ;
   end

end 

function test_MS_splitting_misfit_lam2S
   % Check that the symmetrical mode works properly. 
   
   fast1 = 30  ;
   tlag1 = 1.5 ;
   
   fast2 = 75  ;
   tlag2 = 0.8 ;
   
   misfit1 = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,0.0,0.1,'mode','lam2S') ;
   misfit2 = MS_splitting_misfit(fast2,tlag2,fast1,tlag1,0.0,0.1,'mode','lam2S') ;

   assertElementsAlmostEqual(misfit1, misfit2, 'absolute',0.001) ;

end 

function test_MS_splitting_misfit_intensity
   % Check that the intensity mode works properly. This has a 
   % different periodicity to the other misfit measures.
   
   fast1 = 90  ;
   tlag1 = 1.5 ;
   
   fast2 = 0  ;
   tlag2 = 1.5 ;
   
   misfit = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,0.0,0.1,'mode','intensity') ;

   assertElementsAlmostEqual(misfit,0,'absolute',0.001) ;

end 

function test_MS_splitting_misfit_simple
   % Check that the simple mode works properly. 
   
   fast1 = 90  ;
   tlag1 = 1.5 ;
   
   fast2 = 90  ;
   tlag2 = 1.5 ;
   
   misfit = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,0.0,0.1,'mode','simple') ;

   assertElementsAlmostEqual(misfit,0,'absolute',0.001) ;

end 

function test_MS_splitting_misfit_Worst_case
   % check misfit between two identical splitting operators is unity for
   % opposite parameters

   fast1 = 90 ; tlag1 = 5;
   fast2 = 00 ; tlag2 = 5;
   spol = 45 ;
   freq = 0.2 ;
   
   misfit = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,spol,freq) ;
   assertElementsAlmostEqual(misfit, 1, 'absolute',0.001) ;
end 

function test_MS_splitting_misfit_analytical
   % Check that the misfit for a given parameter pair and a null result is
   % the same as the 2nd eigenvalue from an analytically constructed split
   % wave. 

   % construct an effectively split shear wave
   fast = 30 ; padf=1; T=5 ; dt= 0.05 ;
   
   N = [ FDGaussian(dt,T,padf).*cosd(fast) FDGaussian(dt,T,padf).*cosd(fast+90) ] ;
   E = [ FDGaussian(dt,T,padf).*sind(fast) FDGaussian(dt,T,padf).*sind(fast+90) ] ;
   
   M = cov(E,N) ;

   % take the eigenvalues   
   [V,D] = eig(M) ;
   
   % calculate normalised lam2
   lam2 = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;

   misfit = MS_splitting_misfit(fast,10,90,0.0,75,0.2) ;
   assertElementsAlmostEqual(misfit, lam2, 'absolute',0.001) ;

end 

function [amp] = FDGaussian(dt,T,padf)
%function [A,t] = GaussWin(dt, T)
%
%
%	dt - sample interval
%	T - dominant period
%  padf - padding factor
%
%	Returns
%	amp - First-derivative gaussian waveform
%	time - time series for waveform which will be -T*pf:dt:T*pf
%
   
   time = -T*padf:dt:T*padf ;

%  peak period 
   sig = T/6 ;
   amp = -(time./(sig.^3*sqrt(2.*pi))).*exp(-time.^2./(2.*sig.^2)) ;

%  normalise amplitude
   amp = amp./max(amp) ;
   
end


%function test_MS_splitting_misfit_reversed
%   % check misfit between two identical splitting operators is zero for different
%   % polarisations, frequencies and operators. 
%
%   % NOTE, this test does not give identical results. I think that this is the
%     non-commutative nature of splitting operators rearing its head. 
%
%   fast1 = 30 ; tlag1 = 1.5 ;
%   fast2 = 60 ; tlag2 = 2.0 ;
%   spol = 45 ;
%   freq = 0.1 ;
%   
%   misfit1 = MS_splitting_misfit(fast1,tlag1,fast2,tlag2,spol,freq) ;
%   misfit2 = MS_splitting_misfit(fast2,tlag2,fast1,tlag1,spol,freq) ;
%
%   assertElementsAlmostEqual(misfit1, misfit2, 'absolute',0.001) ;
%
%end 


