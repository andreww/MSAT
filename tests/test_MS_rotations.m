function test_suite = test_MS_rotations
initTestSuite;
end

function test_MS_rotations_triclinic

    [C,r] = MS_elasticDB('alb') ;
    C1 = MS_rot3(C,30,45,60) ;
    C2 = MS_rot3(C1,-30,-45,-60,'order',[3 2 1]) ;
    
    M = MS_rotM(30,45,60) ;
    
    C3=MS_rotR(C,M) ;
    
    assertElementsAlmostEqual(C, C2) ;
    assertElementsAlmostEqual(C1, C3) ;

end

