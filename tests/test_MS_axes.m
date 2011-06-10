function test_suite = test_MS_axes
initTestSuite;
end


function test_MS_axes_ortho
    % Test the interpolator works for constant
    % matrices rotated about axies
 
    [C, r] = MS_elasticDB('olivine');
    C90x = MS_rot3(C, 90, 0, 0);
    C45x = MS_rot3(C, 45, 0, 0);
    C40y = MS_rot3(C, 0, 40, 0);
    C30y = MS_rot3(C, 0, 30, 0);
    C31y = MS_rot3(C, 0, 31, 0);
    C35y = MS_rot3(C, 0, 35, 0);
    C39y = MS_rot3(C, 0, 39, 0);
    C10z = MS_rot3(C, 0, 0, 10);
    C110z = MS_rot3(C, 0, 0, 110);
    C60z = MS_rot3(C, 0, 0, 60);
 
    C20z = MS_rot3(C, 0, 0, 20);
    C30z = MS_rot3(C, 0, 0, 30);
    C45z = MS_rot3(C, 0, 0, 45);
    C90z = MS_rot3(C, 0, 0, 90);
    
    Cmulti = MS_rot3(C, 15, 67, 23);

    assertElementsAlmostEqual(MS_axes(C),MS_axes(C90z));
    assertElementsAlmostEqual(MS_axes(C90z),MS_axes(Cmulti));

end

function test_MS_axes_triclin
    % Test the interpolator works for constant
    % matrices rotated about axes
 
    [C, r] = MS_elasticDB('albite');
    C90x = MS_rot3(C, 90, 0, 0);
    C45x = MS_rot3(C, 45, 0, 0);
    C40y = MS_rot3(C, 0, 40, 0);
    C30y = MS_rot3(C, 0, 30, 0);
    C31y = MS_rot3(C, 0, 31, 0);
    C35y = MS_rot3(C, 0, 35, 0);
    C39y = MS_rot3(C, 0, 39, 0);
    C10z = MS_rot3(C, 0, 0, 10);
    C110z = MS_rot3(C, 0, 0, 110);
    C60z = MS_rot3(C, 0, 0, 60);
 
    C20z = MS_rot3(C, 0, 0, 20);
    C30z = MS_rot3(C, 0, 0, 30);
    C45z = MS_rot3(C, 0, 0, 45);
    C90z = MS_rot3(C, 0, 0, 90);

    assertElementsAlmostEqual(MS_axes(C),MS_axes(C90z));

end

