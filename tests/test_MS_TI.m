function test_suite = test_MS_TI
initTestSuite;
end

function test_MS_TI_crosscheck

%% Check we get the same results from the two VTI routines. 
    
%% generate elasticities from thomsen parameters.
   rh = 4000 ; 
   vpv = 8 ;
   vsv = 5 ;
   
   C = MS_TI(vpv,vsv,4000,0.05,0.05,0.00,'thomsen') ;
   
%% calculate parameters for VTI2
   vph = sqrt(C(1,1)*1e9./rh)./1e3 ;
   vsh = sqrt(C(6,6)*1e9./rh)./1e3 ;
   
%% average velocities
   vpa = sqrt((vpv.^2+4.*vph.^2)./5) ;
   vsa = sqrt((2.*vsv.^2+vsh.^2)./3) ;
   
   xi = vsh.^2 / vsv.^2 ;
   phi = vpv.^2 / vph.^2 ; 
   eta = C(1,3) ./ (C(1,1) - 2.*C(4,4)) ;
   
   Cl = MS_TI(vpa,vsa,rh,xi,phi,eta,'panning') ;
   
   assertElementsAlmostEqual(C, Cl) ;

end 
