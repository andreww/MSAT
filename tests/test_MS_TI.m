function test_suite = test_MS_TI
initTestSuite;
end

function test_MS_TI_crosscheck

% Check we get the same results from the two VTI routines. 
    
% generate elasticities from thomsen parameters.
   rh = 4000 ; 
   vpv = 8 ;
   vsv = 5 ;
   
   C = MS_TI(vpv,vsv,4000,0.05,0.05,0.00,'thomsen') ;
   
% calculate parameters for VTI2
   vph = sqrt(C(1,1)*1e9./rh)./1e3 ;
   vsh = sqrt(C(6,6)*1e9./rh)./1e3 ;
   
% average velocities
   vpa = sqrt((vpv.^2+4.*vph.^2)./5) ;
   vsa = sqrt((2.*vsv.^2+vsh.^2)./3) ;
   
   xi = vsh.^2 / vsv.^2 ;
   phi = vpv.^2 / vph.^2 ; 
   eta = C(1,3) ./ (C(1,1) - 2.*C(4,4)) ;
   
   Cl = MS_TI(vpa,vsa,rh,xi,phi,eta,'panning') ;
   
   assertElementsAlmostEqual(C, Cl) ;

end 

function test_MS_TI_parameters_1

% Check we get the same results from the two VTI routines. 
    
% generate elasticities from thomsen parameters.
   rh_in = 4000 ; 
   vpv_in = 8 ;
   vsv_in = 5 ;
   eps_in = 0.05;
   gam_in = 0.05;
   del_in = 0.00;
   
   C = MS_TI(vpv_in, vsv_in, rh_in, eps_in, gam_in, del_in,'thomsen'); 
 
   [loveA, loveC, loveL, loveN, loveF, vpv, vsv, eps, gam, del ...
             vpa, vsa, xi, phi, eta] = MS_TI_parameters(C, 4000);
          
   assertElementsAlmostEqual(vpv_in, vpv);
   assertElementsAlmostEqual(vsv_in, vsv);
   assertElementsAlmostEqual(eps_in, eps);
   assertElementsAlmostEqual(gam_in, gam);
   assertElementsAlmostEqual(del_in, del);
   
   C2 = MS_TI(loveA, loveC, loveL, loveN, loveF,'love') ;
   assertElementsAlmostEqual(C, C2);
   
   C3 = MS_TI(vpa, vsa, rh_in, xi, phi, eta,'panning') ;
   assertElementsAlmostEqual(C, C3);

end

function test_MS_TI_parameters_errors

    [C, rh] = MS_elasticDB('stishovite');
    
    f = @()MS_TI_parameters(C, rh);
        assertExceptionThrown(f, 'MS:BadTIelasticity');
    
    
    [C, rh] = MS_elasticDB('apatite');
    
    f = @()MS_TI_parameters(MS_rot3(C, 0.0, 90.0, 0.0), rh);
        assertExceptionThrown(f, 'MS:BadTIelasticity');
end


function test_MS_TI_apatite

    [C, rh] = MS_elasticDB('apatite');
    
    [loveA, loveC, loveL, loveN, loveF, vpv, vsv, eps, gam, del ...
        vpa, vsa, xi, phi, eta] = MS_TI_parameters(C, rh);
    
    C1 = MS_TI(vpv, vsv, rh, eps, gam, del,'thomsen');
    C2 = MS_TI(loveA, loveC, loveL, loveN, loveF,'love') ;
    C3 = MS_TI(vpa, vsa, rh, xi, phi, eta,'panning') ;
    
    assertElementsAlmostEqual(C, C1);
    assertElementsAlmostEqual(C, C2);
    assertElementsAlmostEqual(C, C3);
    
end

function test_MS_TI_apatite_units

    [C, rh] = MS_elasticDB('apatite');
    
    [~, ~, ~, ~, ~, ~, ~, eps, gam, del ...
        ~, ~, xi, phi, eta] = MS_TI_parameters(C, rh);
    
    [~, ~, ~, ~, ~, ~, ~, eps1, gam1, del1 ...
        ~, ~, xi1, phi1, eta1] = MS_TI_parameters(C/10, rh);
    
    assertElementsAlmostEqual(xi,  xi1);
    assertElementsAlmostEqual(phi, phi1);
    assertElementsAlmostEqual(eta, eta1);
    assertElementsAlmostEqual(eps, eps1);
    assertElementsAlmostEqual(gam, gam1);
    assertElementsAlmostEqual(del, del1);
    
    [~, ~, ~, ~, ~, ~, ~, eps1, gam1, del1 ...
        ~, ~, xi1, phi1, eta1] = MS_TI_parameters(C*10, rh);
    
    assertElementsAlmostEqual(xi,  xi1);
    assertElementsAlmostEqual(phi, phi1);
    assertElementsAlmostEqual(eta, eta1);
    assertElementsAlmostEqual(eps, eps1);
    assertElementsAlmostEqual(gam, gam1);
    assertElementsAlmostEqual(del, del1);
    
    [~, ~, ~, ~, ~, ~, ~, eps1, gam1, del1 ...
        ~, ~, xi1, phi1, eta1] = MS_TI_parameters(C, 3000);
    
    assertElementsAlmostEqual(xi,  xi1);
    assertElementsAlmostEqual(phi, phi1);
    assertElementsAlmostEqual(eta, eta1);
    assertElementsAlmostEqual(eps, eps1);
    assertElementsAlmostEqual(gam, gam1);
    assertElementsAlmostEqual(del, del1);
end