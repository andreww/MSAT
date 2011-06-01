function test_suite = test_MS_rot_mtex_rotation
initTestSuite;
end

function test_MS_rot_mtex_rotation_scalar_vector

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
    o = rotation('Euler',90*degree,0*degree,0*degree,'Roe');
    os = [o ; o ; o];
    
    % Scalar test
    assertElementsAlmostEqual(C_in, MS_rot_mtex_rotation(C_in, o))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rot_mtex_rotation(Cs_in, o))
    assertElementsAlmostEqual(Cs_in, MS_rot_mtex_rotation(C_in, os))
    assertElementsAlmostEqual(Cs_in, MS_rot_mtex_rotation(Cs_in, os))
end

function test_MS_rot_mtex_rotation_errors

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
    o = rotation('Euler',90*degree,0*degree,0*degree,'Roe');
    os = [o ; o ; o];
    noto = 'not a rotation object';
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot_mtex_rotation(Cs_in, os);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % Not a rotation object
    f = @()MS_rot_mtex_rotation(Cs_in, noto);
    assertExceptionThrown(f, 'MS:notarot')
    
end