function test_suite = test_MS_VRH
initTestSuite;
end

function test_MS_VRH_onephase
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
end

function test_MS_VRH_badinput
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    f = @()MS_VRH([0.333,0.333,0.333],C_in,10,C_in,10);
    assertExceptionThrown(f, 'MS:VRHargs');
end
