%% test_MS_poisson_isotropic
    
    Ciso = MS_build_isotropic('e', 300,'nu',0.2);
    % Do we always get the answer we expect for an isotropic case?
    assertElementsAlmostEqual(MS_poisson( Ciso, [1 0 0], [0 0 1] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [1 0 0], [0 1 0] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [1; 0; 0], [0; 0; 1] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [1; 0; 0], [0; 1; 0] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [0 1 0], [0 0 1] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [0 1 0], [1 0 0] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [0 0 1], [0 1 0] ), 0.2);
    assertElementsAlmostEqual(MS_poisson( Ciso, [0 0 1], [1 0 0] ), 0.2);

%% test_MS_poisson_errors

    Ciso = MS_build_isotropic('e', 300,'nu',0.2);
    % N not a unit vector
    f = @()MS_poisson( Ciso, [1 0 0], [0 0 10]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [1 0 0], [0 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
    % M not a unit vector
    f = @()MS_poisson( Ciso, [10 0 0], [0 0 1]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [0 0], [0 0 1]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
    % N and M not a normal
    f = @()MS_poisson( Ciso, [1 0 0], [1 0 0]);
    assertExceptionThrown(f, 'MS:POISSON:NotOrthogonal');
    % N not a unit vector
    f = @()MS_poisson( Ciso, [1 ;0 ;0], [0 ;0 ;10]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [1; 0; 0], [0; 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
    % M not a unit vector
    f = @()MS_poisson( Ciso, [10; 0; 0], [0; 0; 1]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [0; 0], [0; 0; 1]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
    % N and M not a normal
    f = @()MS_poisson( Ciso, [1; 0; 0], [1; 0; 0]);
    assertExceptionThrown(f, 'MS:POISSON:NotOrthogonal');
    % Odd cases...
    f = @()MS_poisson( Ciso, [0; 0; 0], [1; 0; 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [1; 0; 0], [0; 0; 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [0; 0; 0], [0; 0; 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotUnit');
    f = @()MS_poisson( Ciso, [1 0 0], [0 0 1 0]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
    f = @()MS_poisson( Ciso, [1 0 0 0], [0 0 1]);
    assertExceptionThrown(f, 'MS:CHECKUNITNotVec');
