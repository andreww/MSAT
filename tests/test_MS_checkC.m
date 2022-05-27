%% test_MS_checkC_isNum
    f = @()MS_checkC('My matrix');
    assertExceptionThrown(f, 'MS:CHECKCNotNumeric');
    f = @()MS_checkC({23,22});
    assertExceptionThrown(f, 'MS:CHECKCNotNumeric');
    



%% test_MS_checkC_size
    f = @()MS_checkC([1 1 1 1 1 1]);
    assertExceptionThrown(f, 'MS:CHECKCnot6x6');
    f = @()MS_checkC(ones(6,6,6));
    assertExceptionThrown(f, 'MS:CHECKCnot6x6');
    f = @()MS_checkC(ones(5,5));
    assertExceptionThrown(f, 'MS:CHECKCnot6x6');
    f = @()MS_checkC(ones(7,7));
    assertExceptionThrown(f, 'MS:CHECKCnot6x6');


%% test_MS_checkC_psodef

    f = @()MS_checkC(-1.*ones(6,6));
    assertExceptionThrown(f, 'MS:CHECKCnotposdef');


%% test_MS_checkC_zeros

    f = @()MS_checkC(zeros(6,6));
    assertExceptionThrown(f, 'MS:CHECKCbadzeros');
    


%% test_MS_checkC_supression

    assert(MS_checkC(zeros(6,6), 'nozerochk', 'nopdefchk') == 1)
    assert(MS_checkC(-1*ones(6,6), 'nopdefchk') == 1)
    assert(MS_checkC(zeros(6,6), 'fast') == 1)
    assert(MS_checkC(-1*ones(6,6), 'fast') == 1)

