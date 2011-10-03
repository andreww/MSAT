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

% Copyright (c) 2011, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
