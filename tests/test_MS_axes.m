%% test_MS_axes_reference
    % Test the olivine example from page 671 of B&C
    % This checks that (for this case) the orentation
    % is correct (no permutation of axes).
    C_ref = [ 192.0  66.0  60.0   0.0   0.0   0.0 ; ...
               66.0 160.0  56.0   0.0   0.0   0.0 ; ...
               60.0  56.0 272.0   0.0   0.0   0.0 ; ...
                0.0   0.0   0.0  60.0   0.0   0.0 ; ...
                0.0   0.0   0.0   0.0  62.0   0.0 ; ...
                0.0   0.0   0.0   0.0   0.0  49.0 ];
    % NB: The reference matrix (above) has the second
    % stiffest direction in the X1 direction. Because our
    % d eigen-vecors are sorted we have stiffnesses X3>X2>X1
    % Rotate C_ref into our frame...
    C_ref = MS_rot3(C_ref,0,0,90);
    % Shouldn't matter for future decoposition (T, H and I maticies
    % are not changed by this, only O.
            
    for i=1:25
      a1 = rand.*360 ;
      b1 = rand.*360 ; 
      g1 = rand.*360 ;
      
      if 0
          fprintf('Rotations: [%6.1f %6.1f %6.1f]\n',a1,b1,g1)
      end
      
      C1 = MS_rot3(C_ref,a1,b1,g1) ;

      assertElementsAlmostEqual(MS_axes(C1),C_ref);
    end


%% test_MS_axes_ortho
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C, r] = MS_elasticDB('olivine');
   for i=1:10
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      if 0
          fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      end
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end


%% test_MS_axes_triclin_weak
   % Test the interpolator works for constant
   % matrices rotated about axes
   close all
   [C, r] = MS_elasticDB('albite');
   [Ciso] = MS_decomp(C) ;
   
   [Cweak] = MS_VRH([0.05 0.95],C,r,Ciso,r) ;
   
   C=Cweak ;
    
   for i=1:10
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      if 0
          fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      end
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end


%% test_MS_axes_triclin_strong
% Test the interpolator works for constant
% matrices rotated about axes
 
   [C, r] = MS_elasticDB('albite');

   for i=1:10
      a1 = rand.*360 ; a2 = rand.*360 ;
      b1 = rand.*360 ; b2 = rand.*360 ;
      g1 = rand.*360 ; g2 = rand.*360 ;
      
      if 0
          fprintf('Rotations: [%6.1f %6.1f %6.1f] [%6.1f %6.1f %6.1f]\n',a1,b1,g1,a2,b2,g2)
      end
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      C2 = MS_rot3(C,a2,b2,g2) ;

      assertElementsAlmostEqual(MS_axes(C1),MS_axes(C2));
    
   end


%% test_MS_axes_triclin_strong_rots
    % Test the interpolator works for constant
    % matrices rotated about axes
 
   [C, r] = MS_elasticDB('albite');

   for i=1:25
      a1 = rand.*360 ;
      b1 = rand.*360 ; 
      g1 = rand.*360 ; 
      
      if 0
          fprintf('Rotations: [%6.1f %6.1f %6.1f] \n',a1,b1,g1)
      end
      
      C1 = MS_rot3(C,a1,b1,g1) ;
      [CC, RR] = MS_axes(C1);
      C2 = MS_rotR(C1,RR) ;
      C3 = MS_rotR(CC,RR');
      assertElementsAlmostEqual(CC,C2);
      assertElementsAlmostEqual(C1,C3);
    
   end


%% test_MS_axes_errors

    [C, ~] = MS_elasticDB('olivine');
    f = @()MS_axes(C, 'NotAnArgument');
    assertExceptionThrown(f, 'MS:AXES:UnknownOption')

    % This is open coded so that we can check the output
    % ... cannot work out how to do this by passing a function
    % handle into assertExceptionThrown. Maybe a patch to xunit
    % is in order
    val = 0;  
    try
        [out1, out2, out3] = MS_axes(C);
    catch e
        assert(strcmp(e.identifier, ...
      'MS:AXES:BadOutputArgs'), 'Wrong error thrown')
       val = 1;
    end
    assert(val==1, 'No error thrown')

