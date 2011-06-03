% MS_AXES - Reorient elasticity matrix for optimal decomposition.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the principle axes of elasticity matrix C, after:
%     Browaeys and Chevrot (2004) (GJI, v159, 667-678, 2004)
%
%  % [ CR, ... ] = MS_axes(C)
%
% Usage: 
%     [CR] = MS_axes(C)                    
%         Return a rotated elasticity matrix minimising the number of 
%         distinct elements. This is the orientation which, on further 
%         decomposition using MS_NORMS, will maximise the high symmetry
%         components of the matrix. 
%
%     [CR, RR] = MS_axes(C)
%         In addition, return the rotation matrix, RR, used to perform
%         the rotation to generate C from CR.
%
% Notes:
%     If the input matrix has isotropic, hexagonal or tetragonal 
%     symmetry there are multiple orentations of the principle axes.
%     in the isotropic case CR not rotated with respect to C (and RR
%     is the identity matrix). In the hexagonal and tetragonal cases, 
%     X3 is defined by the distinct eigenvalue (see Browaeys and Chevrot).
%     For the monoclinic or triclinic cases we have to make a 'best-guess'
%     and following Browraeys and Chevrot we use the bisectrix of each of 
%     the eigenvectors of d and its closest match in v.
%
% References:
%     Browaeys, J. T. and S. Chevrot (2004) Decomposition of the elastic
%         tensor and geophysical applications. Geophysical Journal 
%         international v159, 667-678.
%     Cowin, S. C. and M. M. Mehrabadi (1987) On the identification of 
%         material symmetry for anisotropic elastic materials. Quartely
%         Journal of Mechanics and Applied Mathematics v40, 451-476.

% (C) James Wookey and Andrew Walker, 2011.
% 
%

function [ CR RR ] = MS_axes( C )

det_thresh = 0.01 ; % threshold on flagging an error on the orthogonality 
                    % of the best guess axes 

% calculate D and V matrices   
d= [... 
   (C(1,1)+C(1,2)+C(1,3)) (C(1,6)+C(2,6)+C(3,6)) (C(1,5)+C(2,5)+C(3,5)) ; ...
   (C(1,6)+C(2,6)+C(3,6)) (C(1,2)+C(2,2)+C(3,2)) (C(1,4)+C(2,4)+C(3,4)) ; ...
   (C(1,5)+C(2,5)+C(3,5)) (C(1,4)+C(2,4)+C(3,4)) (C(1,3)+C(2,3)+C(3,3)) ; ...
] ;

v=[...
   (C(1,1)+C(6,6)+C(5,5)) (C(1,6)+C(2,6)+C(4,5)) (C(1,5)+C(3,5)+C(4,6)) ; ...
   (C(1,6)+C(2,6)+C(4,5)) (C(6,6)+C(2,2)+C(4,4)) (C(2,4)+C(3,4)+C(5,6)) ; ...
   (C(1,5)+C(3,5)+C(4,6)) (C(2,4)+C(3,4)+C(5,6)) (C(5,5)+C(4,4)+C(3,3)) ; ...
];
   
% calculate eigenvectors and eigenvalues of D and V 
[vecd,val]=eig(d) ; vald = [val(1,1) val(2,2) val(3,3)] ;    
[vecv,val]=eig(v) ; valv = [val(1,1) val(2,2) val(3,3)] ;

% count number of distinct eigenvalues. Maximum allowable difference is set
% 1/1000th of the norm of the matrix. 
[nud,id] = ndistinct(valv,norm(valv)./1000) ;
[nuv,iv] = ndistinct(vald,norm(vald)./1000) ;

% set rotation flag
irot = 0 ;

% use the number of distinct eigenvalues to decide the symmetry of the medium
switch nud
case 1
%  isotropic tensor, choice is arbitrary, so leave it the way it is!
   X1=[1 0 0] ;
   X2=[0 1 0] ;
   X3=[0 0 1] ;
   irot = 0 ;
case 2
% hexagonal or tetragonal tensor, only one is defined, other ones are 
% any two othogonal vectors. X3 is defined by the distinct eigenvalue
% (see Browaeys and Chevrot).
  X3=vecd(:,id)' ;
% check that X3 has changed
   if (X3(1)==0 & X3(2)==0)   
      irot=0;
   else      
%     now set the other two vectors.    
%     we want X2 to be horizontal ...
      X2=cross(X3,[X3(1) X3(2) 0]) ; X2=X2./norm(X2);
%     now have a definition for X1
      X1=cross(X3,X2) ;  X1=X1./norm(X1) ;     
      irot=1;
   end   
case 3
%  orthorhombic or lower symmetry
%  first figure out how many common vectors there are   
   neq = 0;
   for i=1:3
      for j=1:3
         neq = neq + veceq(vecd(:,i),vecv(:,j),0.01) ;
      end
   end
% for if neq=3, then symmetry is orthorhombic   
   if (neq==3)
%  significant axes are the three eigenvectors
      X1=vecd(:,1)' ;
      X2=vecd(:,2)' ;
      X3=vecd(:,3)' ;
      irot = 1 ;
      if (X3(1)==0 & X3(2)==0), irot=0;, end
   else
% monoclinic or triclinic. Here we have to make a 'best-guess'. Following
% Browraeys and Chevrot we use the bisectrix of each of the eigenvectors of
% d and its closest match in v.
      D1=vecd(:,1) ; D2=vecd(:,2) ; D3=vecd(:,3) ;
      V1=vecv(:,1) ; V2=vecv(:,2) ; V3=vecv(:,3) ;

%      plotvec(D1,'g-'); plotvec(D2,'g-'); plotvec(D3,'g-'); 
%      plotvec(V1,'r-'); plotvec(V2,'r-'); plotvec(V3,'r-'); 
      
      [dum,ind] = max([dot(D1,V1) dot(D1,V2) dot(D1,V3)]) ;
      X1 = bisectrix(D1,vecv(:,ind))' ;
      [dum,ind] = max([dot(D2,V1) dot(D2,V2) dot(D2,V3)]) ;
      X2 = bisectrix(D2,vecv(:,ind))' ;
      [dum,ind] = max([dot(D3,V1) dot(D3,V2) dot(D3,V3)]) ;
      X3 = bisectrix(D3,vecv(:,ind))' ;

%      plotvec(X1,'b-'); plotvec(X2,'b-'); plotvec(X3,'b-'); 
%      axis([-1 1 -1 1 -1 1]); daspect([1 1 1]) ;
%      fprintf('Perp. check: %f %f\n',dot(X1,X2),dot(X1,X3));

      irot = 1 ;
      if (X3(1)==0 & X3(2)==0), irot=0;, end      
   end
otherwise
end

% Now apply the necessary rotation. The three new axes define the rotation
% matrix which turns the 3x3 unit matrix (original axes) into the best projection.
% So the rotation matrix we need to apply to the elastic tensor is the inverse
% of this. 

if irot
%  check axes
   dps = abs([dot(X1,X2) dot(X1,X3) dot(X2,X3)]) 
   if (length(find(dps>det_thresh))>0)
      warning('MS_axes: Improper rotation matrix resulted, not rotating.') ;
      CR=C ;
      return
   end

%  fix up the axes, for safety, by redefining X2 and X3   
   X3 = cross(X1,X2) ;
   X2 = cross(X1,X3) ;
      
%  construct forward rotation matrix
   R1 = [X1' X2' X3'] ;

%  calculate reverse rotation
   RR = inv(R1) ;

   % check rotation matrix
   if (abs(det(RR))-1)>det_thresh ;
       warning('MS_axes: Improper rotation matrix resulted, not rotating.') ;
       fprintf('Determinant = %20.18f\n',det(RR))
       RR = [1 0 0 ; 0 1 0 ; 0 0 1];
       CR = C ;
       return 
   end
   
   % apply to the input elasticity matrix
   CR = MS_rotR(C,RR) ;
else
   warning('No rotation was deemed necessary.')
   RR = [1 0 0 ; 0 1 0 ; 0 0 1];
   CR = C;
end
return

%%%
%%%   SUBFUNCTIONS
%%%

function [C]=bisectrix(A,B)
% return the unit length bisectrix of 3-vectors A and B
     C=(A+B);
     C=C./norm(C) ;
return

function [i]=veceq(x,y,thresh)
% return the number and indices of distinct entries in a 3x3 element vector, ignoring a 
% difference of thresh
   i=1 ;
   if length(find(abs(x-y)>thresh))>0, i=0 ;, end   
return

function [nd,i]=ndistinct(x,thresh)
% return the number and indices of distinct entries in a 3 element vector, ignoring a 
% difference of thresh
   if (abs(x(1)-x(2))<thresh) & (abs(x(1)-x(3))<thresh) % all the same
      nd = 1 ; i = 0 ;
   elseif (abs(x(1)-x(2))>thresh) & (abs(x(1)-x(3))>thresh) % all different
      nd = 3 ; i = 0 ;
   else % one is different
      nd = 2 ; 
      if (abs(x(1)-x(2))<thresh)
         i = 3 ;
      elseif (abs(x(1)-x(3))<thresh)
         i = 2 ;
      else 
         i = 1 ;
      end  
   end   
return
 
function plotvec(V,spec)
   plot3([0 V(1)],[0 V(2)],[0 V(3)],spec)
   hold on
return
