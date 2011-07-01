% MS_AXES - Reorient elasticity matrix for optimal decomposition.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  Calculate the principle axes of elasticity tensor C, after: 
%     Browaeys and Chevrot (GJI, v159, 667-678, 2004)
%  
% [ CR, ... ] = MS_axes(C)
%
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
%  [ ... ] = MS_axes( C, 'nowarn' )
%     Suppress all warnings.
%
%  [ ... ] = MS_axes( C, 'X3_stiff' )
%     For hexagonal or tetragonal tensors, make X3 the stiffest direction,
%     not the distinct direction.
%
%  [ ... ] = MS_axes( C, 'debug' )
%     Enable debugging plots and messages. These are quite messy.  
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
%
% See also: MS_NORMS, MS_INTERPOLATE

% (C) James Wookey and Andrew Walker, 2011

function [ varargout ] = MS_axes( C, varargin )

warn = 1 ;
debug = 0 ;
X3_distinct = 1;

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'nowarn' % flag (i.e., no value required)
               warn = 0 ;
               iarg = iarg + 1 ;
            case 'x3_stiff'
               X3_distinct = 0;
               iarg = iarg + 1 ;
            case 'debug'
               debug = 1;
               iarg = iarg + 1 ;
            otherwise 
               error('MS:AXES:UnknownOption',...
                  ['Unknown option: ' varargin{iarg}]) ;   
         end   
      end

det_thresh = 100 * sqrt(eps) ; % threshold on flagging an error on the orthogonality 
                               % of the best guess axes 

% construct D and V matrices   
d= [... 
   (C(1,1)+C(1,2)+C(1,3)) (C(1,6)+C(2,6)+C(3,6)) (C(1,5)+C(2,5)+C(3,5)) ; ...
   (C(1,6)+C(2,6)+C(3,6)) (C(1,2)+C(2,2)+C(3,2)) (C(1,4)+C(2,4)+C(3,4)) ; ...
   (C(1,5)+C(2,5)+C(3,5)) (C(1,4)+C(2,4)+C(3,4)) (C(1,3)+C(2,3)+C(3,3)) ; ...
] ;

v=[...
   (C(1,1)+C(6,6)+C(5,5)) (C(1,6)+C(2,6)+C(4,5)) (C(1,5)+C(3,5)+C(4,6)) ; ...
   (C(1,6)+C(2,6)+C(4,5)) (C(6,6)+C(2,2)+C(4,4)) (C(2,4)+C(3,4)+C(5,6)) ; ...
   (C(1,5)+C(3,5)+C(4,6)) (C(2,4)+C(3,4)+C(5,6)) (C(5,5)+C(4,4)+C(3,3)) ; ...
] ;
   
% calculate eigenvectors and eigenvalues of D and V 
[vecd,val]=eig(d) ; vald = [val(1,1) val(2,2) val(3,3)] ;    
[vecv,val]=eig(v) ; valv = [val(1,1) val(2,2) val(3,3)] ;

D1=vecd(:,1) ; D2=vecd(:,2) ; D3=vecd(:,3) ;
V1=vecv(:,1) ; V2=vecv(:,2) ; V3=vecv(:,3) ;

% check orthogonality of these
if ~isortho(D1,D2,D3) ;
   warning('D basis is not othogonal') ;
end

if ~isortho(V1,V2,V3) ;
   warning('D basis is not othogonal') ;
end

% count number of distinct eigenvalues. Maximum allowable difference is set
% to be 1/1000th of the norm of the matrix. 
[nud,id] = ndistinct(valv,norm(valv)./1000) ;
[nuv,iv] = ndistinct(vald,norm(vald)./1000) ;

% set rotation flag
irot = 0 ;
isMonoOrTri = 0 ;
% use the number of distinct eigenvalues to decide the symmetry of the medium
switch nud
case 1
%  isotropic tensor, choice is arbitrary, so leave it the way it is!
   X1=[1 0 0] ;
   X2=[0 1 0] ;
   X3=[0 0 1] ;
   irot = 0 ;
case 2
   if X3_distinct
%  ** hexagonal or tetragonal tensor, only one is defined, other ones are 
%     any two othogonal vectors. X3 is defined by the distinct eigenvalue
%     (see Browaeys and Chevrot).
      X3=vecd(:,id)' ;
%  ** check that X3 has changed
      if (X3(1)==0 & X3(2)==0)   
         irot=0;
      else      
%     ** now set the other two vectors.    
%        we want X2 to be horizontal ...
         X2=cross(X3,[X3(1) X3(2) 0]) ; X2=X2./norm(X2);

%        now have a definition for X1
         X1=cross(X3,X2) ;  X1=X1./norm(X1) ;     
         irot=1;
      end
   else
      % Eigen vectors are sorted, so treat like orthorhombic
      X1=vecd(:,1)' ;
      X2=vecd(:,2)' ;
      X3=vecd(:,3)' ;
      if (all(all([X1' X2' X3'] == eye(3))))
        irot = 0;
      end
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
      
      [X1,X2,X3] = makeRH(X1,X2,X3) ; % make sure they form a RH set.

      irot = 1 ;
      if (all(all([X1' X2' X3'] == eye(3))))
          irot = 0;
      end
   else
% monoclinic or triclinic. Here we have to make a 'best-guess'. Following
% Browraeys and Chevrot we use the bisectrix of each of the eigenvectors of
% d and its closest match in v.
       
      isMonoOrTri = 1 ;
      
      [X1,X2,X3]=estimate_basis(D1,D2,D3,V1,V2,V3) ;

      irot = 1 ;
      if (X3(1)==0 & X3(2)==0), irot=0; end      
   end
otherwise
% not possible.
end

% Now apply the necessary rotation. The three new axes define the rotation
% matrix which turns the 3x3 unit matrix (original axes) into the best projection.
% So the rotation matrix we need to apply to the elastic tensor is the inverse
% of this. 

if irot
%  check axes
   dps = abs([dot(X1,X2) dot(X1,X3) dot(X2,X3)]) ;
   if (~isempty(find(dps>det_thresh)))
      if warn
         warning('MS_axes: Determined axes not orthogonal: DPS=') ;
         dps
      end   
   end
   
%  construct forward rotation matrix
   R1 = [X1' X2' X3'] ;
   
%  calculate reverse rotation
   RR = inv(R1) ;
   
   if ((det(RR)+1) < det_thresh)
%   ** Fix rotoinversion
       if 0
          warning('MS_axes: Fixing up rotoinversion.') ;
          fprintf('Determinant = %20.18f\n',det(RR))
       end
       R1 = [-X1' X2' X3'] ;
       RR = inv(R1) ;
   end

   % check rotation matrix
   if (abs(det(RR))-1)>det_thresh ;
       if warn
          warning('MS_axes: Improper rotation matrix resulted, not rotating.') ;
          fprintf('Determinant = %20.18f\n',det(RR))
       end
       RR = eye(3);
       CR=C ;
   else
%  ** apply to the input elasticity matrix
      CR = MS_rotR(C,RR) ;

%  ** if crystal is monoclinic or triclinic, we have some more work to do. We need to make
%     sure that differently rotated crystals return the same result. So they are consistent
%     we try 180 degree rotations around all the principle axes. For each flip, we choose
%     the result that gives the fastest P-wave velocity in the [1 1 1] direction
      if isMonoOrTri
         [CR,RR] = tryrotations(CR,RR) ;
      end   
   end
else
   if warn, warning('No rotation was deemed necessary.');, end
   CR=C;
   RR = eye(3) ;
end

   switch nargout
   case 0
      varargout{1} = CR ;
   case 1
      varargout{1} = CR ;
   case 2
      varargout{1} = CR ;
      varargout{2} = RR ;
   otherwise   
      error('MS:INFO:BadOutputArgs','Requires 1 or 2 output arguments.')
   end
      

return

%%%
%%%   SUBFUNCTIONS
%%%

function [ CRR, RRR ] = tryrotations( CR, RR ) ;
%  ** try all combinations of 180 rotations around principle axes
%     select the one that gives the highest Vp in the [1 1 1] direction
      x1r = [000 180 000 000 180 000] ;
      x2r = [000 000 180 000 180 180] ;
      x3r = [000 000 000 180 000 180] ;
      vp = zeros(1,6) ;
      for i=1:6
         [ ~, ~, ~, ~, vp(i), ~, ~ ] = ...
            MS_phasevels( MS_rot3(CR,x1r(i),x2r(i),x3r(i)), 2000, 45, 45 ) ;
      end
      [~,ind] = sort(vp) ; % select fastest
      
      CRR = MS_rot3(CR,x1r(ind(end)),x2r(ind(end)),x3r(ind(end))) ;

      a = x1r(ind(end)) * pi/180. ; R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
      b = x2r(ind(end)) * pi/180. ; R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
      g = x3r(ind(end)) * pi/180. ; R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

      RRR =  R3 * R2 * R1 * RR ;

return

function [C1,C2,C3]=estimate_basis(A1,A2,A3,B1,B2,B3)
%
%     Estimate the best basis vectors for decomposition. This is bisectrices
%     of vectors A1 and B1. 
%           
      [~, ind] = max(abs([dot(A1,B1) dot(A1,B2) dot(A1,B3)])) ;
      eval(sprintf('X=B%1.1i;',ind)) ;
      if dot(A1,X)<0
         C1 = bisectrix(A1,-X)' ;
      else   
         C1 = bisectrix(A1, X)' ;
      end   
      
      [dum, ind] = max(abs([dot(A2,B1) dot(A2,B2) dot(A2,B3)])) ;
      eval(sprintf('X=B%1.1i;',ind)) ;
      if dot(A2,X)<0
         C2 = bisectrix(A2,-X)' ;
      else   
         C2 = bisectrix(A2, X)' ;
      end   

      [dum, ind] = max(abs([dot(A3,B1) dot(A3,B2) dot(A3,B3)])) ;
      eval(sprintf('X=B%1.1i;',ind)) ;
      if dot(A3,X)<0
         C3 = bisectrix(A3,-X)' ;
      else   
         C3 = bisectrix(A3, X)' ;
      end   

%
%     Enforce orthoganality, and renormalise
%
      C3 = cross(C1,C2) ;
      C2 = cross(C1,C3) ;
      C1 = C1 ./ norm(C1) ;
      C2 = C2 ./ norm(C2) ;
      C3 = C3 ./ norm(C3) ;
return

function [C]=bisectrix(A,B)
% return the unit length bisectrix of 3-vectors A and B
     C=(A+B);
     C=C./norm(C) ;
return

function [i]=veceq(x,y,thresh)
% return the number and indices of distinct entries in a 3 element vector, 
% ignoring a difference of thresh
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

function [Y1,Y2,Y3] = makeRH(X1,X2,X3) ;
   if (dot(cross(X1,X2),X3))<0
      Y1 = X1; Y2 = X2; Y3 = -X3 ;
   else
      Y1 = X1; Y2 = X2; Y3 = X3 ;
   end   
return

%  vector plotting functions for debug purposes
function plotvec(V,spec,str)
   plot3([0 V(1)],[0 V(2)],[0 V(3)],spec)
   text([V(1)],[V(2)],[V(3)],str)
   hold on
return

function plot2vec(V,spec,str)
   if V(3)<0, spec = [spec '-'];, end
   plot([0 V(1)],[0 V(2)],spec)
   text([V(1)],[V(2)],str)
   hold on
return

function [iortho] = isortho(X1,X2,X3)
   thresh = 10 * sqrt(eps) ;
   dps = abs([dot(X1,X2) dot(X1,X3) dot(X2,X3)]) ;
   iortho = isempty(find(dps>thresh)) ;
return
   
   
