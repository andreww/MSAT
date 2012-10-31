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

function test_MS_rot_Euler_scalar_vector_opts

    % Note that for these rotations passive and active rotatiosn give the 
    % same result.
   
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
    assertElementsAlmostEqual(C_in, MS_rotEuler(C_in, phi1, theta, phi2, ...
        'sense', 'active'))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1, theta, phi2, ...
        'sense', 'active'))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(C_in, phi1s, thetas, phi2s, ...
        'sense', 'active'))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1s, thetas, phi2s, ...
        'sense', 'active'))
    
    % Scalar test
    assertElementsAlmostEqual(C_in, MS_rotEuler(C_in, phi1, theta, phi2, ...
        'sense', 'passive'))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1, theta, phi2, ...
        'sense', 'passive'))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(C_in, phi1s, thetas, phi2s, ...
        'sense', 'passive'))
    assertElementsAlmostEqual(Cs_in, MS_rotEuler(Cs_in, phi1s, thetas, phi2s, ...
        'sense', 'passive'))
    
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
    
    % Opt args
    f = @()MS_rotEuler(Cs_in, phi1s, thetas, phi2, 'bob');
    assertExceptionThrown(f, 'MS:ROTE:UnknownOption')
    f = @()MS_rotEuler(Cs_in, phi1s, thetas, phi2, 'sense', 'bob');
    assertExceptionThrown(f, 'MS:ROTE:UnknownOption')
    
end

function test_MS_rot_Euler_VPSC

    % The VPSC code (version 6c) allows the calculation of texture (given
    % by sets of Euler angles) following deformation. It also calculates
    % the Voigt (and SC) elasticity of the aggregate at the end of
    % defomation. Here we use the results of a VPSC calculation on 10
    % crystals in simple shear to check that MS_rot_Euler works.

    % Single crystal elastic constants (GPa)
    C_single_xtl =  [ 1225.0   409.0   328.0   000.0   000.0   000.0 ; ...                                  
                       409.0   928.0   484.0   000.0   000.0   000.0 ; ...                                 
                       328.0   484.0  1211.0   000.0   000.0   000.0 ; ...                                 
                       000.0   000.0   000.0   281.0   000.0   000.0 ; ...                                 
                       000.0   000.0   000.0   000.0   260.0   000.0 ; ...                                 
                       000.0   000.0   000.0   000.0   000.0   404.0 ];   

    % Euler angles (Bunge notation) at end of VPSC run
    euler_angles = [  129.78  145.70   -9.06 ; ...
                      -76.71   23.15   34.51 ; ...
                       77.12   16.58  125.27 ; ...
                      145.53  145.62   -3.27 ; ...
                      142.69  142.05  -28.88 ; ...
                      -12.35   17.86 -108.23 ; ...
                       -6.76   26.27   -8.20 ; ...
                      178.81  161.46  -23.48 ; ...
                      -37.69   29.70  169.14 ; ...
                     -113.62  153.84 -140.89 ];

    % Voigt avarage elasticity from end of VPSC run - should reproduce this
    C_vpsc_voigt =  [1139.88  382.70  409.74  -18.03   16.17  -12.94 ; ...
                      382.70 1063.68  440.82   -1.80  -18.70   -3.63 ; ...
                      409.74  440.82 1135.92  -36.72  -14.86   16.33 ; ...
                      -18.03   -1.80  -36.72  289.68    0.82   -9.29 ; ...
                       16.17  -18.70  -14.86    0.82  308.29   -3.17 ; ...
                      -12.94   -3.63   16.33   -9.29   -3.17  359.30];
    
    
    % Calculate avatage elasticity
    % First do some setup
    n_xtls = size(euler_angles,1); % Number of crystals
    Cs = zeros(6,6,n_xtls);        % Array for each crystal stiffness
    rh = ones(n_xtls,1);           % density (we do not use this)
    VF = ones(n_xtls,1)/n_xtls;    % Volume fraction
    % Then rotate the elasticites.
    for i = 1:n_xtls
        Cs(:,:,i) = MS_rotEuler( C_single_xtl, euler_angles(i,1), ...
                                 euler_angles(i,2), euler_angles(i,3) );
    end
    % Calculate avarage.
    [~,~,C_voigt,~] = MS_VRH(VF, Cs, rh);
                  

    % Check result - only know VPSC result to 0.01 GPa
    assertElementsAlmostEqual(C_vpsc_voigt, C_voigt, 'absolute', 0.01)

    
    % Calculate avatage elasticity again with options!
    % First do some setup
    n_xtls = size(euler_angles,1); % Number of crystals
    Cs = zeros(6,6,n_xtls);        % Array for each crystal stiffness
    rh = ones(n_xtls,1);           % density (we do not use this)
    VF = ones(n_xtls,1)/n_xtls;    % Volume fraction
    % Then rotate the elasticites.
    for i = 1:n_xtls
        Cs(:,:,i) = MS_rotEuler( C_single_xtl, euler_angles(i,1), ...
                                 euler_angles(i,2), euler_angles(i,3), ...
                                 'sense', 'active');
    end
    % Calculate avarage.
    [~,~,C_voigt,~] = MS_VRH(VF, Cs, rh);
                  

    % Check result - only know VPSC result to 0.01 GPa
    assertElementsAlmostEqual(C_vpsc_voigt, C_voigt, 'absolute', 0.01)

end

