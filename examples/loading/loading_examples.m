% example script to demonstrate loading elastic constants in MSAT. 

%  ** Load a file in default format :
      C = MS_load('Olivine.txt') ;

%  ** ^^ with density
      [C,rh] = MS_load('Olivine.txt') ;

%  ** density normalised elasticities, where the pre-normalised elasticities are
%     in Pa.
      [C,rh] = MS_load('Olivine.Aij','Aij','eunit','Pa') ;
      
%  ** in Ematrix format      
      [C] = MS_load('Olivine.Ematrix','format','ematrix') ;

%  ** Isotropic file
      [C,rh] = MS_load('Olivine.iso','symmetry','iso') 

