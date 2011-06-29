% MS_VROT3 - Rotate a (set of) 3-vector(s) by 3 angles.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Short description of functionality. 
%
% [VR] = MS_rot3(V,alp,bet,gam)
%     alp = clockwise about 1-axis (looking at origin, == yaw)
%     bet = clockwise about 2-axis (looking at origin, == -dip)
%     gam = clockwise about 3-axis (looking at origin, == azimuth)
%
%  [VR] = MS_rot3(V,alp,bet,gam)
%
%
%  Notes:
%     The rotations are applied in order, ie: alpha, then beta then gamma
%  One dimension of the input V matrix should be 3. The code tries to infer the 
%  orientation of the matix by testing which this is. If the matrix is 3x3,
%  it is assumed that it contains row vectors (i.e., 
%  [x1,y1,z1; x2,y2,y2; x3 y3 y3])

% Change log / copyright statements.
% 
%

%===============================================================================
   function [VR] = MS_Vrot3(Vin,alp,bet,gam,varargin)
%===============================================================================

      reverse = 0 ;

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'reverse' % flag (i.e., no value required)
               reverse = 1 ;
               iarg = iarg + 1 ;
            otherwise 
               error('MS:VROT3:UnknownOption',...
                  ['Unknown option: ' varargin{iarg}]) ;   
         end   
      end

      [nr nc] = size(Vin) ; % store this to produce output
      itran = 0 ;
      if nr==3 & nc==3 
%     ** row vector
         itran = 1 ;
         V = transpose(Vin) ;
      elseif nr~=3 & nc==3 
%     ** row vector
         itran = 1 ;
         V = transpose(Vin) ;
      elseif nr==3 & nc~=3
%     ** column vector
         V = Vin ;
      else
         error('MS:VROT3BadInputMatrix', ...
            'Input vector matrix must be 3xN or Mx3') ;
      end      

%  ** Make rotation matrix
      a = alp * pi/180. ;
      b = bet * pi/180. ;
      g = gam * pi/180. ;
      
      R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
      R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
      R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;
      
      if reverse
          RR = R1 * R2 * R3 ;
      else    
          RR =  R3 * R2 * R1;
      end
%  ** apply it
      VR = RR * V;
      
%  ** put it back in the original shape      
      if itran, VR=transpose(VR);, end ;
 
%===============================================================================
return
%===============================================================================


