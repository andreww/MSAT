function test_suite = test_MS_effective_medium
initTestSuite;
end



function test_MS_effective_medium_TW

   Ceff = [6.4667    3.5875    3.5875         0         0         0 ; ...
           3.5875    7.3026    3.7466         0         0         0 ; ...
           3.5875    3.7466    7.3026         0         0         0 ; ...
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

function test_MS_effective_medium_Backus
   
   [h,vp,vs,rh] = example_layering()  ;
   
   n = length(h) ;
   
   C = zeros(6,6,n) ;
   
   for i=1:n
      C(:,:,i) = MS_iso(vp(i),vs(i),rh(i)) ;
   end
   
   [Ceff1,rheff1] = MS_effective_medium('backus',h,vp,vs,rh) ;
   
   [Ceff2,rheff2] = MS_effective_medium('backus',h,C,rh) ;
   
   
   Vp0 =   3360.728 ;
   Vs0 =   1754.414 ;
   Rho =   2439.401 ;

   Eps =   0.122594 ;
   Gam =   0.165415 ;
   Del =  -0.007099 ;
   
   Ceffr = MS_TI(Vp0./1e3,Vs0./1e3,Rho,Eps,Gam,Del,'thomsen') ;
   
   assertElementsAlmostEqual(Ceffr, Ceff1, 'absolute',0.001) ;
   assertElementsAlmostEqual(Ceffr, Ceff2, 'absolute',0.001) ;
   assertElementsAlmostEqual(Rho, rheff1, 'absolute',0.001) ;
   assertElementsAlmostEqual(Rho, rheff2, 'absolute',0.001) ;
   
end

function [h,vp,vs,rh] = example_layering()

   data = [ 5.0  2659.99    1319.64    2265.49 ; ...
            5.0  2666.03    1551.17    2120.92 ; ...
            5.0  4666.90    2059.59    2551.36 ; ...
            5.0  3718.53    2130.02    2483.37 ; ...
            5.0  5169.32    2793.55    2732.48 ; ...
            5.0  3328.20    1937.93    2435.49 ; ...
            5.0  3183.62    1650.54    2282.05 ; ...
            5.0  5193.23    2718.97    2788.21 ; ...
            5.0  5447.86    2881.82    2823.21 ; ...
            5.0  2931.93    1490.97    2282.22 ; ...
            5.0  3017.37    1681.12    2217.13 ; ...
            5.0  3115.85    1517.10    2393.99 ; ...
            5.0  3451.60    1794.47    2390.58 ; ...
            5.0  3083.05    1561.44    2216.66 ; ...
            5.0  3317.40    1733.62    2283.42 ; ...
            5.0  4734.10    2578.15    2588.18 ; ...
            5.0  4238.52    1963.57    2662.89 ; ...
            5.0  2484.67    1302.08    2247.81 ; ...
            5.0  3642.96    2261.31    2387.66 ; ...
            5.0  4385.17    2107.17    2693.91 ; ...
            5.0  3234.26    1545.98    2380.38 ] ; 

   h = data(:,1) ;
   vp = data(:,2)./1e3 ;
   vs = data(:,3)./1e3 ;
   rh = data(:,4) ;

end

