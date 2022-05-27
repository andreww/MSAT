%% test_MS_load_default

   [Cp,rp] = C_Ol() ;
   
   % create the output file to load.
   fid = fopen('/tmp/test_MS_load_default.txt','wt') ;
   
   fprintf(fid,'%% Default format           \n') ;
   fprintf(fid,'1 1 320.5000                \n') ;
   fprintf(fid,'1 2 68.1000                 \n') ;
   fprintf(fid,'1 3 71.6000                 \n') ;
   fprintf(fid,'2 2 196.5000                \n') ;
   fprintf(fid,'2 3 76.8000                 \n') ;
   fprintf(fid,'3 3 233.5000                \n') ;
   fprintf(fid,'4 4 64.0000                 \n') ;
   fprintf(fid,'5 5 77.0000                 \n') ;
   fprintf(fid,'6 6 78.7000                 \n') ;
   fprintf(fid,'7 7 3355.0 %% density, kg/m3\n') ;
   
   fclose(fid) ;
   
   [Cl,rl] = MS_load('/tmp/test_MS_load_default.txt') ;

   assertElementsAlmostEqual(Cp, Cl, 'absolute',0.001) ;
   assertElementsAlmostEqual(rp, rl, 'absolute',0.001) ;

%% test_MS_load_Aij

   [Cp,rp] = C_Ol() ;
   
   % create the output file to load.
   fid = fopen('/tmp/test_MS_load_Aij.txt','wt') ;
   
   fprintf(fid,'% Aij               \n') ;
   fprintf(fid,'1 1 9.552906e+07    \n') ;
   fprintf(fid,'1 2 2.029806e+07    \n') ;
   fprintf(fid,'1 3 2.134128e+07    \n') ;
   fprintf(fid,'2 2 5.856930e+07    \n') ;
   fprintf(fid,'2 3 2.289121e+07    \n') ;
   fprintf(fid,'3 3 6.959762e+07    \n') ;
   fprintf(fid,'4 4 1.907601e+07    \n') ;
   fprintf(fid,'5 5 2.295082e+07    \n') ;
   fprintf(fid,'6 6 2.345753e+07    \n') ;
   fprintf(fid,'7 7 3355.0 % density\n') ;
   
   fclose(fid) ;
   
   [Cl,rl] = MS_load('/tmp/test_MS_load_Aij.txt','Aij','eunit','Pa') ;

   assertElementsAlmostEqual(Cp, Cl, 'absolute',0.001) ;
   assertElementsAlmostEqual(rp, rl, 'absolute',0.001) ;

%% test_MS_load_ematrix

   [Cp,rp] = C_Ol() ;
   
   % create the output file to load.
   fid = fopen('/tmp/test_MS_load_ematrix.txt','wt') ;
   
   fprintf(fid,'Olivine elastic constants (Abramson et al, 1996) in Mbar\n') ;
   fprintf(fid,'This is the ematrix file format.                        \n') ;
   fprintf(fid,'  0.0000  0.0000  0.0000 90.0000 90.0000 90.0000        \n') ;
   fprintf(fid,'  3.2050  0.6810  0.7160  0.0000  0.0000  0.0000        \n') ;
   fprintf(fid,'  0.6810  1.9650  0.7680  0.0000  0.0000  0.0000        \n') ;
   fprintf(fid,'  0.7160  0.7680  2.3350  0.0000  0.0000  0.0000        \n') ;
   fprintf(fid,'  0.0000  0.0000  0.0000  0.6400  0.0000  0.0000        \n') ;
   fprintf(fid,'  0.0000  0.0000  0.0000  0.0000  0.7700  0.0000        \n') ;
   fprintf(fid,'  0.0000  0.0000  0.0000  0.0000  0.0000  0.7870        \n') ;
   
   fclose(fid) ;
   
   [Cl] = MS_load('/tmp/test_MS_load_ematrix.txt','format','ematrix') ;

   assertElementsAlmostEqual(Cp, Cl, 'absolute',0.001) ;

%% test_MS_load_iso

   C33 = 237.5533 ;
   C66 = 79.5400 ;
   rho = 3355 ;
   
   % create the output file to load.
   fid = fopen('/tmp/test_MS_load_iso.txt','wt') ;
   
   fprintf(fid,'%% Isotropic format\n') ;
   fprintf(fid,'3 3 237.5533\n') ;
   fprintf(fid,'6 6 79.5400\n') ;
   fprintf(fid,'7 7 3355.0 % density, kg/m3\n') ;
   
   fclose(fid) ;
   
   [Cl,rl] = MS_load('/tmp/test_MS_load_iso.txt','symmetry','iso') ;

   assertElementsAlmostEqual(C33, Cl(3,3), 'absolute',0.001) ;
   assertElementsAlmostEqual(C66, Cl(6,6), 'absolute',0.001) ;
   assertElementsAlmostEqual(rho, rl, 'absolute',0.001) ;

function [C,r] = C_Ol()

C = [ 320.5000   68.1000   71.6000         0         0         0 ;
       68.1000  196.5000   76.8000         0         0         0 ;
       71.6000   76.8000  233.5000         0         0         0 ;
             0         0         0   64.0000         0         0 ;
             0         0         0         0   77.0000         0 ;
             0         0         0         0         0   78.7000 ] ;
             
r = 3355 ;              

end