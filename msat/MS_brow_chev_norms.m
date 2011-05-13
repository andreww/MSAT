function [P]=CIJ_brow_chev_norms(varargin)
%
%  Calculate the percentages of the original matrix norm represented by the
%  individual decompositions of an elasticity tensor. 
%
%  See Browaeys and Chevrot (GJI, v159, 667-678, 2004), and function
%     CIJ_brow_chev_decomp()
%
%  [P] = CIJ_brow_chev_norms(C,Ciso)
%     Isotropic part (scalar)
%
%  [P] = CIJ_brow_chev_norms(C,Ciso,Chex)
%     Isotropic, and hexagonal parts of the elastic tensor (vector)
%
%  [P] = CIJ_brow_chev_norms(C,Ciso,Chex,Ctet,Cort,Cmon,Ctri)
%     All parts of the elastic tensor (vector)
%   
   
   
   
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
