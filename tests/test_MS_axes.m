function test_suite = test_MS_axes
initTestSuite;
end

function test_MS_axes_ortho
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C, r] = MS_elasticDB('olivine');
   for i=1:25
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end
end

function test_MS_axes_triclin_weak
   % Test the interpolator works for constant
   % matrices rotated about axes
   close all
   [C, r] = MS_elasticDB('albite');
   [Ciso] = MS_decomp(C) ;
   
   [Cweak] = MS_VRH([0.05 0.95],C,r,Ciso,r) ;
   
   C=Cweak ;
    
   for i=1:25
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end
end

function test_MS_axes_triclin_strong
    % Test the interpolator works for constant
    % matrices rotated about axes
 
   [C, r] = MS_elasticDB('albite');

   for i=1:25
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end
end

