function test_suite = test_MS_rot3
initTestSuite;
end

function test_MS_rot_3_triclinic

    [C,r] = MS_elasticDB('alb') ;
    C1 = MS_rot3(C,30,45,60) ;
    C2 = MS_rot3(C1,-30,-45,-60,'order',[3 2 1]) ;
    
    M = MS_rotM(30,45,60) ;
    
    C3=MS_rotR(C,M) ;
    
    assertElementsAlmostEqual(C, C2) ;
    assertElementsAlmostEqual(C1, C3) ;

end

function test_MS_rot_3_scalar_vector

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
    a = 90.0;
    b = 0.0;
    g = 0.0;
    as = [a ; a; a];
    bs = [b ; b ; b];
    gs = [g ; g ; g];
    
    % Scalar test
    assertElementsAlmostEqual(C_in, MS_rot3(C_in, a, b, g))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rot3(Cs_in, a, b, g))
    assertElementsAlmostEqual(Cs_in, MS_rot3(C_in, as, bs, gs))
    assertElementsAlmostEqual(Cs_in, MS_rot3(Cs_in, as, bs, gs))
end

function test_MS_rot_3_errors

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    badC_in = [  49.67  13.18   13.18   000.0   000.0   000.0
                 31.18  49.67   13.18   000.0   000.0   000.0
                 13.18  13.18   49.67   000.0   000.0   000.0
                 000.0  000.0   000.0   12.78   000.0   000.0
                 000.0  000.0   000.0   000.0   12.78   000.0
                 000.0  000.0   000.0   000.0   000.0   12.78];
    a = 90.0;
    b = 0.0;
    g = 0.0;
    
    f = @()MS_rot3(badC_in, a, b, g);
    assertExceptionThrown(f, 'MS:ROT3BadInputMatrix')
        
    f = @()MS_rot3(C_in, a, b, g, 'order',[1 1 1 1]);
    assertExceptionThrown(f, 'MS:ROT3:BadOrder')
    
    f = @()MS_rot3(C_in, a, b, g, 'order',[1 2 4]);
    assertExceptionThrown(f, 'MS:ROT3:BadOrder')
    
    f = @()MS_rot3(C_in, a, b, g, 'bob',[1 2 3]);
    assertExceptionThrown(f, 'MS:ROT3:UnknownOption')
    
   
end

function test_MS_rot_3_vecerrors

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
    a = 90.0;
    b = 0.0;
    g = 0.0;
    as = [a ; a; a];
    bs = [b ; b ; b];
    gs = [g ; g ; g];
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, bs, gs);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, b, gs);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, bs, g);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
end

