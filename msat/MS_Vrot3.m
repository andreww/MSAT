% MS_VROT3 - Rotate a (set of) 3-vector(s) by 3 angles.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotates 3-vector(s). 
%
% [VR] = MS_Vrot3(V,alp,bet,gam)
%     alp = clockwise about 1-axis (looking at origin, == yaw)
%     bet = clockwise about 2-axis (looking at origin, == -dip)
%     gam = clockwise about 3-axis (looking at origin, == azimuth)
%
%     [VR] = MS_Vrot3(...,'order',V)                    
%         order: V is a 3 element vector containing the order in which to apply the 3 rotations.
%                Default is [1 2 3].
%
%
%  Notes:
%     The rotations are applied in order, ie: alpha, then beta then gamma
%  One dimension of the input V matrix should be 3. The code tries to infer the 
%  orientation of the matix by testing which this is. If the matrix is 3x3,
%  it is assumed that it contains row vectors (i.e., 
%  [x1,y1,z1; x2,y2,y2; x3 y3 y3])
%
%  See also: MS_ROT3

% (C) James Wookey and Andrew Walker, 2011 
%

%===============================================================================
   function [VR] = MS_Vrot3(Vin,alp,bet,gam,varargin)
%===============================================================================

orderV = [1 2 3] ;

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'order' 
         orderV = varargin{iarg+1} ;
         if (length(find(orderV==1))+length(find(orderV==2))+length(find(orderV==3)))~=3
            error('MS:VROT3:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;   
         end
         if (length(orderV)~=3)
            error('MS:VROT3:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;  
         end                  
         iarg = iarg + 2 ;
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
      

R(1,:,:) = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R(2,:,:) = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R(3,:,:) = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  squeeze(R(orderV(3),:,:)) * ...
         squeeze(R(orderV(2),:,:)) * ...
            squeeze(R(orderV(1),:,:));
            
%  ** apply it
      VR = RR * V;
      
%  ** put it back in the original shape      
      if itran, VR=transpose(VR);, end ;
 
%===============================================================================
return
%===============================================================================


