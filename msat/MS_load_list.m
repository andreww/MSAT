% MS_LOAD_LIST - Load a list of CIJs, in a specific format
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [x,C,rho] = MS_load_list(fname)
%
% Usage:
%    As above, x is the independent variable (the first column in the file, 
%    of length n) C is an array 6*6*n where n is the number of tensors 
%    loaded rho is a vector of length n, the last column in the file
%    currently, only 21 constant elastic files can be loaded. Lines in 
%    the file should be of the form:
%        x, c11,c12,...,c16,c22,...,c26,c33,...,c66,rho 
%
% See also: MS_LOAD, MS_ELASTICDB

% (C) James Wookey and Andrew Walker, 2011


function [x,C,rho] = MS_load_list(fname) 


ecdata=load(fname) ;

[n nec] = size(ecdata) ;

nec = nec - 2; % two columns for 

if nec~=21
   error('Only 21 elastic constants allowed') ;
end

C = zeros(6,6,n) ;
x = zeros(n,1) ;
rho = zeros(n,1) ;

for irow = 1:n
   ipos = 1 ;
   x(irow) = ecdata(irow,1) ;
   for i=1:6 ;
      for j=i:6 ;
         ipos = ipos + 1 ;
         C(i,j,irow) = ecdata(irow,ipos) ;
         if (i ~= j) ; % Populate the lower half of the matrix.
             C(j, i, irow) = C(i, j, irow);
         end
      end
   end
   rho(irow) = ecdata(irow,23) ;
   
end

return
