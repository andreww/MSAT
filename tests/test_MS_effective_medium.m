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


