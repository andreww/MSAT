% MS_NORMS - Browaeys and Chevrot analysis of the elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  Calculate the percentages of the original matrix norm represented by the
%  individual decompositions of an elasticity tensor C, after: 
%     Browaeys and Chevrot (GJI, v159, 667-678, 2004)
%  
%  [P] = MS_norms(C,Ciso)
%     Isotropic part (scalar)
%
%  [P] = MS_norms(C,Ciso,Chex)
%     Isotropic, and hexagonal parts of the elastic tensor (vector)
%
%  [P] = MS_norms(C,Ciso,Chex,Ctet,Cort,Cmon,Ctri)
%     All parts of the elastic tensor (vector)
%     
%
% Notes:
%     Input components of of C must first be calculated using MS_decomp.
%     The output vector P will have as may elements as the number of
%     decomposed parts of the input matrix that are provided.
%
% References:
%     Browaeys, J. T. and S. Chevrot (2004) Decomposition of the elastic
%         tensor and geophysical applications. Geophysical Journal 
%         international v159, 667-678.
%
% See also: MS_AXES, MS_DECOMP

% (C) James Wookey and Andrew Walker, 2011   
   
function [ P ] = MS_norms( varargin )
      
%if (nargin~=(nargout+1)), ...
%   error('Must be one more input than output, see help.'), end

Cref = varargin{1} ;
[Xref]=C2X(Cref) ;

N = sqrt(Xref'*Xref) ;

Ctot=zeros(6,6);
   
for i=2:nargin
   C = varargin{i} ;
   Ctot = Ctot + varargin{i} ;
   [Xtot]=C2X(Ctot) ;
   XD=Xref-Xtot;
   P(i-1)=1-(sqrt(XD'*XD)/N) ;
end

% transform P
PP = [0 P] ;    
P(1:nargin-1) = PP(2:nargin)-PP(1:nargin-1) ;

return

function Xn=mynorm(X)
%   Xn=sqrt(sum(X.^2)) ; % Euclidean norm?
    Xn=sqrt(dot(X,X)) ;
    Xn=norm(X) ;
return

function [X]=C2X(C)
%  after Browaeys and Chevrot (GJI, 2004)
	X = zeros(21,1) ;
	
	X(1)  = C(1,1) ;
	X(2)  = C(2,2) ;
	X(3)  = C(3,3) ;
	X(4)  = sqrt(2).*C(2,3) ;
	X(5)  = sqrt(2).*C(1,3) ;
	X(6)  = sqrt(2).*C(1,2) ;
	X(7)  = 2.*C(4,4) ;
	X(8)  = 2.*C(5,5) ;
	X(9)  = 2.*C(6,6) ;
	X(10) = 2.*C(1,4) ;
	X(11) = 2.*C(2,5) ;
	X(12) = 2.*C(3,6) ;
	X(13) = 2.*C(3,4) ;
	X(14) = 2.*C(1,5) ;
	X(15) = 2.*C(2,6) ;
	X(16) = 2.*C(2,4) ;
	X(17) = 2.*C(3,5) ;
	X(18) = 2.*C(1,6) ;
	X(19) = 2.*sqrt(2).*C(5,6) ;
	X(20) = 2.*sqrt(2).*C(4,6) ;
	X(21) = 2.*sqrt(2).*C(4,5) ;
	
return
