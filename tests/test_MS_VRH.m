%% test_MS_VRH_onephase
    % Test the Voigt-Reuss-Hill avarage does not
    % change the result for a single phase 
    % calculation
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    assertElementsAlmostEqual(MS_VRH([0.5,0.5],C_in,10,C_in,10), C_in);
    assertElementsAlmostEqual(MS_VRH([0.75,0.25],C_in,10,C_in,10), C_in);
    assertElementsAlmostEqual(MS_VRH([50,50],C_in,10,C_in,10), C_in);
    assertElementsAlmostEqual(MS_VRH([9987,932],C_in,10,C_in,10), C_in);
    assertElementsAlmostEqual(MS_VRH([9987 932 5],C_in,10,C_in,10,C_in,10), C_in);

%% test_MS_VRH_badinput
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    f = @()MS_VRH([0.333,0.333,0.333],C_in,10,C_in,10);
    assertExceptionThrown(f, 'MS:VRH:args');

%% test_MS_VRH_onephase_vectors
    % Test the Voigt-Reuss-Hill avarage does not
    % change the result for a single phase 
    % calculation - vector form
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs = zeros(6,6,2);
    Cs(:,:,1) = C_in;
    Cs(:,:,2) = C_in;
    assertElementsAlmostEqual(MS_VRH([9987 932],Cs, [10 10]), C_in);
    assertElementsAlmostEqual(MS_VRH([0.75 0.25],Cs, [10 10]), C_in);
    assertElementsAlmostEqual(MS_VRH([50 50], Cs, [10 10]), C_in);
    assertElementsAlmostEqual(MS_VRH([9987 932], Cs, [10 10]), C_in);
    Cs = zeros(6,6,3);
    Cs(:,:,1) = C_in;
    Cs(:,:,2) = C_in;
    Cs(:,:,3) = C_in;
    assertElementsAlmostEqual(MS_VRH([9987 932 5],Cs,[10 10 10]), C_in);

%% test_MS_VRH_badinput_vectors
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs = zeros(6,6,3);
    Cs(:,:,1) = C_in;
    Cs(:,:,2) = C_in;
    Cs(:,:,3) = C_in;
    f = @()MS_VRH([0.333,0.333,0.333],Cs,[10 10]);
    assertExceptionThrown(f, 'MS:VRH:args');
    f = @()MS_VRH([0.333,0.333],Cs,[10 10 10]);
    assertExceptionThrown(f, 'MS:VRH:args');
    f = @()MS_VRH([0.333,0.333,0.333],C_in,[10 10 10]);
    assertExceptionThrown(f, 'MS:VRH:args');
    f = @()MS_VRH([0.333,0.333,0.333],ones(6,5,3),[10 10 10]);
    assertExceptionThrown(f, 'MS:VRH:args');
    f = @()MS_VRH([0.333,0.333,0.333],ones(5,6,3),[10 10 10]);
    assertExceptionThrown(f, 'MS:VRH:args');

