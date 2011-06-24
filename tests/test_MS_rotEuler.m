function test_suite = test_MS_rotEuler
initTestSuite;
end

function test_MS_rot_Euler_scalar_vector

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs_in = zeros(6,6,3);
    Cs_in(:,:,1) = C_in;
    Cs_in(:,:,2) = C_in;
    Cs_in(:,:,3) = C_in;
    phi1 = 90.0;
    theta = 0.0;
    phi2 = 0.0;
    phi1s = [phi1 ; phi1; phi1];
    thetas = [theta ; theta ; theta];
    phi2s = [phi2 ; phi2 ; phi2];
    
    % Scalar test
    assertElementsAlmostEqual(C_in, MS_rotEuler(C_in, phi1, theta, phi2))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1, theta, phi2))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(C_in, phi1s, thetas, phi2s))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1s, thetas, phi2s))
end

function test_MS_rot_Euler_errors

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs_in = zeros(6,6,3);
    Cs_in(:,:,1) = C_in;
    Cs_in(:,:,2) = C_in;
    Cs_in(:,:,3) = C_in;
    Cs_in(:,:,4) = C_in;
    phi1 = 90.0;
    theta = 0.0;
    phi2 = 0.0;
    phi1s = [phi1 ; phi1; phi1];
    thetas = [theta ; theta ; theta];
    phi2s = [phi2 ; phi2 ; phi2];
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rotEuler(Cs_in, phi1s, thetas, phi2s);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rotEuler(Cs_in, phi1s, theta, phi2s);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rotEuler(Cs_in, phi1s, thetas, phi2);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
end