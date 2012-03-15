function test_suite = test_MS_effective_medium
initTestSuite;
end

function test_MS_effective_medium_TW

   Ceff = [2.7412    1.4338    1.4338         0         0         0 ; ...
           1.4338    6.0575    2.5015         0         0         0 ; ...
           1.4338    2.5015    6.0575         0         0         0 ; ...
                0         0         0    1.7780         0         0 ; ...
                0         0         0         0    1.1864         0 ; ...
                0         0         0         0         0    1.1864 ] ;
   
   rheff = 1900 ;
   
   vpm = 2.0 ; vsm = 1.0 ; rhm = 2000.0 ;
   vpi = 1.5 ; vsi = 0.001 ; rhi = 1000 ;
   
   del = 0.1 ;
   vf = 0.1 ;
   
   Cm = MS_iso(vpm,vsm,rhm) ;
   Ci = MS_iso(vpi,vsi,rhi) ;
   
   [Ceff1, rheff1] = MS_effective_medium('tandon',Cm,rhm,Ci,rhi,del,vf) ;
   [Ceff2, rheff2] = MS_effective_medium('tandon',vpm,vsm,rhm,vpi,vsi,rhi,del,vf) ;
   
   assertElementsAlmostEqual(Ceff, Ceff1, 'absolute',0.001) ;
   assertElementsAlmostEqual(Ceff, Ceff2, 'absolute',0.001) ;
   assertElementsAlmostEqual(rheff, rheff1, 'absolute',0.001) ;
   assertElementsAlmostEqual(rheff, rheff2, 'absolute',0.001) ;

end

function test_MS_effective_medium_Hudson

%  water saturated cracks (HCS1) case from Crampin, 1984 (GJRAS, 76, pp 135-145)
%  see table 1, page 139
   
   C11 = 87.464 - 0.142 ;
   C22 = 87.464 - 0.016 ;
   C33 = C22 ;
   C23 = 29.142 - 0.016 ;
   C13 = 29.142 - 0.047 ;
   C12 = C13 ;
   C44 = 29.161 ;
   C55 = 29.161 - 6.666 + 0.745 ;
   C66 = C55 ;
   
   Ceff = [C11 C12 C13  0   0   0   ; ...
           C12 C22 C23  0   0   0   ; ...
           C13 C23 C33  0   0   0   ; ...
            0   0   0  C44  0   0   ; ...
            0   0   0   0  C55  0   ; ...
            0   0   0   0   0  C66] ;
   
   rheff = 2599.9497 ;
   
   rhm = 2600.0 ;
   vpm = 5.8000 ;
   vsm = 3.3490 ;
   aspr = 0.0001 ;
   cden = 0.1 ;
   
   vpc = 1.5000 ;
   vsc = 0.0001 ;
   rhc = 1000.0  ;

   [Ceff1, rheff1] = MS_effective_medium('hudson',vpm,vsm,rhm,vpc,vsc,rhc,aspr,cden) ;
      
   assertElementsAlmostEqual(Ceff, Ceff1, 'absolute',0.001) ;
   assertElementsAlmostEqual(rheff, rheff1, 'absolute',0.01) ;
   
%  note, changing vsc so it doesn't break MS_iso. This requires the tolerances
%  to be lowered a bit, since the Crampin numbers assume vs=0, whereas this is
%  unsupported in MSAT.
   
   vsc = 0.001 ;
   
   Cm = MS_iso(vpm,vsm,rhm) ;
   Cc = MS_iso(vpc,vsc,rhc) ;
   
   [Ceff2, rheff2] = MS_effective_medium('hudson',Cm,rhm,Cc,rhc,aspr,cden) ;
   
   assertElementsAlmostEqual(Ceff, Ceff2, 'absolute',0.01) ;
   assertElementsAlmostEqual(rheff, rheff2, 'absolute',0.01) ;

end
