%% test_MS_measure_no_splitting_dfreq
  
   % If we measure the splitting of an unsplit wave it should always
   % be 0 for tlag... this checks the dominant freq
   [time, T00, T90] = MS_make_trace(0.0, 0.2, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;

   [time, T00, T90] = MS_make_trace(0.0, 0.02, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;
   
   [time, T00, T90] = MS_make_trace(0.0, 2.0, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;

%% test_MS_measure_no_splitting_pol
  
   % If we measure the splitting of an unsplit wave it should always
   % be 0 for tlag... this checks the init_pol freq
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;
   
   [time, T00, T90] = MS_make_trace(30.0, 0.2, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;
   
   [time, T00, T90] = MS_make_trace(45.0, 0.2, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;
   
   [time, T00, T90] = MS_make_trace(90.0, 0.2, 4.0);
   [~, tlag] = MS_measure_trace_splitting(time, T00, T90, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([0.0 tlag], [0.0 0.0]) ;

%% test_MS_measure_small_splitting_pol_1
  
   % We should recover the splitting we impose
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 25.0, 0.2);
   [fast, tlag] = MS_measure_trace_splitting(time, T00r, T90r, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([fast tlag], [25.0 0.2], 'absolute',0.0001) ;
   
%% test_MS_measure_small_splitting_pol_2
  
   % We should recover the splitting we impose
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 25.0, 0.02);
   [fast, tlag] = MS_measure_trace_splitting(time, T00r, T90r, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([fast tlag], [25.0 0.02], 'absolute',0.1) ;
   
%% test_MS_measure_small_splitting_pol_3
  
   % We should recover the splitting we impose
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 75.0, 0.2);
   [fast, tlag] = MS_measure_trace_splitting(time, T00r, T90r, 5.0, ...
      0.5, 4.0);
   assertElementsAlmostEqual([fast tlag], [75.0 0.2], 'absolute',0.0001) ;
   
%% test_MS_measure_small_splitting_pol_4
  
   % We should recover the splitting we impose
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 16.0, 0.1);
   [fast, tlag] = MS_measure_trace_splitting(time, T00r, T90r, 5.0, ...
      0.5, 4.0);
   % Errors go up close to null.
   assertElementsAlmostEqual([fast tlag], [16.0 0.1], 'absolute',0.1) ;
   
%% test_MS_remove_splitting

   % If we apply and remove splitting from a trace we should
   % get back what we started
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 75.0, 0.2);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 75.0, -0.2);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 75.0, 2.0);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 75.0, -2.0);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 5.0, 0.2);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 5.0, -0.2);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 15.0, 0.2);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 15.0, -0.2);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 0.0, 0.2);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 0.0, -0.2);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 75.0, 0.1);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 75.0, -0.1);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
   
   [time, T00, T90] = MS_make_trace(15.0, 0.2, 4.0);
   [T00r T90r] = MS_split_trace(time, T00, T90, 82.0, 0.2);
   [T00r T90r] = MS_split_trace(time, T00r, T90r, 82.0, -0.2);
   assertElementsAlmostEqual([T00 T90], [T00r, T90r]);
