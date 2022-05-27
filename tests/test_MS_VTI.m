%% test_MS_VTI_crosscheck

% Check we get the same results from the two VTI routines. 
    
% generate elasticities from thomsen parameters.
   rh = 4000 ; 
   vpv = 8 ;
   vsv = 5 ;
   
   C = MS_VTI(vpv,vsv,4000,0.05,0.05,0.00) ;
   
% calculate parameters for VTI2
   vph = sqrt(C(1,1)*1e9./rh)./1e3 ;
   vsh = sqrt(C(6,6)*1e9./rh)./1e3 ;
   
% average velocities
   vpa = sqrt((vpv.^2+4.*vph.^2)./5) ;
   vsa = sqrt((2.*vsv.^2+vsh.^2)./3) ;
   
   xi = vsh.^2 / vsv.^2 ;
   phi = vpv.^2 / vph.^2 ; 
   eta = C(1,3) ./ (C(1,1) - 2.*C(4,4)) ;
   
   Cl = MS_VTI2(vpa,vsa,rh,xi,phi,eta) ;
   
   assertElementsAlmostEqual(C, Cl) ;
