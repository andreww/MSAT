%% test_MS_info_string

    [C, rh] = MS_elasticDB('ol');
    res = sprintf('\n---ROTATIONAL INFORMATION \n\n   Elasticity matrix appears unrotated.\n');  
    S = MS_info(C);
    assert(strcmp(res,S));
    
    res = sprintf([res, '\n---VELOCITY INFORMATION \n\n',...
        '   Isotropic velocities: Vp =  8.415 ,  Vs =  4.869\n\n',...
        '   Direction      VP    VS1    VS2  VSAni \n',...
        '          X1   9.774  4.843  4.791  1.092 \n',...
        '          X2   7.653  4.843  4.368 10.329 \n',...
        '          X3   8.343  4.791  4.368  9.240 \n']);
    S = MS_info(C, rh);
    assert(strcmp(res,S));

%% test_MS_info_args

    [C, rh] = MS_elasticDB('ol');
    res = sprintf('\n---ROTATIONAL INFORMATION \n\n   Elasticity matrix appears unrotated.\n');  
    S = MS_info(C, 'mode', 'rot');
    assert(strcmp(res,S));
    
    res = sprintf([res, '\n---VELOCITY INFORMATION \n\n',...
        '   Isotropic velocities: Vp =  8.415 ,  Vs =  4.869\n\n',...
        '   Direction      VP    VS1    VS2  VSAni \n',...
        '          X1   9.774  4.843  4.791  1.092 \n',...
        '          X2   7.653  4.843  4.368 10.329 \n',...
        '          X3   8.343  4.791  4.368  9.240 \n']);
    S = MS_info(C, rh, 'mode', 'all');
    assert(strcmp(res,S));
    
    res = sprintf(['\n---VELOCITY INFORMATION \n\n',...
        '   Isotropic velocities: Vp =  8.415 ,  Vs =  4.869\n\n',...
        '   Direction      VP    VS1    VS2  VSAni \n',...
        '          X1   9.774  4.843  4.791  1.092 \n',...
        '          X2   7.653  4.843  4.368 10.329 \n',...
        '          X3   8.343  4.791  4.368  9.240 \n']);
    S = MS_info(C, rh, 'mode', 'vel');
    assert(strcmp(res,S));

%% test_MS_info_errors

    [C, rh] = MS_elasticDB('ol');
    
    f = @()MS_info(C, rh, 'mode', 'bob');
    assertExceptionThrown(f, 'MS:INFO:UnsupportedMode');
    
    f = @()MS_info(C, rh, 'notanopt');
    assertExceptionThrown(f, 'MS:INFO:UnknownOption');
    
    f = @()MS_info(C, 'mode', 'bob');
    assertExceptionThrown(f, 'MS:INFO:UnsupportedMode');
    
    f = @()MS_info(C, 'notanopt');
    assertExceptionThrown(f, 'MS:INFO:UnknownOption');
    
    f = @()MS_info(C, 'mode', 'vel');
    assertExceptionThrown(f, 'MS:INFO:NoVelocityInfo');

