function test_suite = test_MS_interpolate
initTestSuite;
end

function [Cint] = mock_MS_interpolate(C1, C2, frac)
    % Mockup of MS_interpolate that does not care about density.
    [ Cint, ~ ] = MS_interpolate(C1, 1.0, C2, 1.0, frac);
end

function test_MS_interpolate_angs_ortho
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C, r] = MS_elasticDB('olivine');
    C90x = MS_rot3(C, 90, 0, 0);
    C45x = MS_rot3(C, 45, 0, 0);
    assertElementsAlmostEqual(mock_MS_interpolate(C, C90x, 0.5), C45x);
    assertElementsAlmostEqual(mock_MS_interpolate(C90x, C, 0.5), C45x);
    C40y = MS_rot3(C, 0, 40, 0);
    C30y = MS_rot3(C, 0, 30, 0);
    C31y = MS_rot3(C, 0, 31, 0);
    C35y = MS_rot3(C, 0, 35, 0);
    C39y = MS_rot3(C, 0, 39, 0);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.5), C35y);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.9), C31y);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.1), C39y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.9), C39y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.5), C35y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.1), C31y);
    C10z = MS_rot3(C, 0, 0, 10);
    C110z = MS_rot3(C, 0, 0, 110);
    C60z = MS_rot3(C, 0, 0, 60);
    assertElementsAlmostEqual(mock_MS_interpolate(C10z, C110z, 0.5), C60z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 0.5), C60z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 1.0), C110z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 0.0), C10z);
end

function test_MS_interpolate_angs_triclin
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C, ~] = MS_elasticDB('albite');
    C90x = MS_rot3(C, 90, 0, 0);
    C45x = MS_rot3(C, 45, 0, 0);
    assertElementsAlmostEqual(mock_MS_interpolate(C, C90x, 0.5), C45x);
    assertElementsAlmostEqual(mock_MS_interpolate(C90x, C, 0.5), C45x);
    C40y = MS_rot3(C, 0, 40, 0);
    C30y = MS_rot3(C, 0, 30, 0);
    C31y = MS_rot3(C, 0, 31, 0);
    C35y = MS_rot3(C, 0, 35, 0);
    C39y = MS_rot3(C, 0, 39, 0);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.5), C35y);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.9), C31y);
    assertElementsAlmostEqual(mock_MS_interpolate(C30y, C40y, 0.1), C39y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.9), C39y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.5), C35y);
    assertElementsAlmostEqual(mock_MS_interpolate(C40y, C30y, 0.1), C31y);
    C10z = MS_rot3(C, 0, 0, 10);
    C110z = MS_rot3(C, 0, 0, 110);
    C60z = MS_rot3(C, 0, 0, 60);
 
    C20z = MS_rot3(C, 0, 0, 20);
    C30z = MS_rot3(C, 0, 0, 30);
    C45z = MS_rot3(C, 0, 0, 45);
    C90z = MS_rot3(C, 0, 0, 90);
    assertElementsAlmostEqual(mock_MS_interpolate(C, C90z, 0.5), C45z);
    assertElementsAlmostEqual(mock_MS_interpolate(C10z, C30z, 0.5), C20z);

    assertElementsAlmostEqual(mock_MS_interpolate(C10z, C110z, 0.5), C60z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 0.5), C60z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 1.0), C110z);
    assertElementsAlmostEqual(mock_MS_interpolate(C110z, C10z, 0.0), C10z);
end

function test_MS_interpolate_forst_fay
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C0, r0] = MS_elasticDB('olivine');
    [C1, r1] = MS_elasticDB('fayalite');
    
    C1_90x = MS_rot3(C1, 90, 0, 0);
    
    for x = 0.0:0.1:1.0
        [~, r_int, C_voigt] = MS_VRH([x 1-x], C0, r0, C1, r1);
        rot = 90*(1-x);
        C_voigt_rot = MS_rot3(C_voigt, rot, 0, 0);
        
        [C_common, rho_common] = MS_interpolate(C0, r0, C1_90x, r1, x);
        assertElementsAlmostEqual(C_common, C_voigt_rot);
        assertElementsAlmostEqual(rho_common, r_int);
        
    end
end
