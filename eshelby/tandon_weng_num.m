%  [CC,rh]=MS_tandon_and_weng(vp,vs,rho,del,c,vpi,vsi,rhoi)
%
%  Calculate anisotropic elastic constants for an isotropic material containing
%  spheroidal inclusions with a symmetry axis aligned with the x1-direction.
%
%  From the theory of Tandon and Weng (1984), based on original FORTRAN code by Mike
%  Kendall.  
%
%  Input parameters:
%       vp,vs,rho : isotropic parameters of the matrix (km/s, kg/m3)
%             del : aspect ratio of spheroids:
%                      del < 1 = smarties
%                      del = 1 = spheres
%                      del > 1 = cigars
%               c : volume fractions of inclusions       
%  vpi, vsi, rhoi : isotropic parameters of the inclusions
%
%  Output parameters:
%     CC : Elastic constants (GPa)
%     rh : aggregate density (kg/m3)
%

 
function [esh]=tandon_weng_num(vp,vs,rho,del,c,vpi,vsi,rhoi) ;

%  weighted average density
   rh = (1.0-c)*rho + c*rhoi ;

   vp = vp * 1e3 ; % convert to m/s 
   vs = vs * 1e3 ; % convert to m/s

   amu = vs*vs*rho ;
   amui = vsi*vsi*rhoi ;
   alam = vp*vp*rho - 2.0*amu ;
   alami = vpi*vpi*rhoi - 2.0*amui ;
   bmi = alami + amui*2.0/3.0 ;
   bmps = alam + amu ;
%  Young's modulus for matrix
   E0 = amu*(3.0*alam + 2.0*amu)/(alam + amu) ;
%  Poisson's ratio of the matrix.
   anu = alam/(2.0*(alam + amu)) ;

%  Some time saving terms
   t1 = del*del - 1.0 ;
   t2 = 1.0 - anu ;
   t3 = 1.0 - 2.0*anu ;
   t4 = 3.0*del*del ;
   t5 = 1.0 - del*del ;
%
% D1, D2 and D3 from Tandon and Weng (1984) (just before equation (18)).
   D1 = 1.0 + 2.0*(amui - amu)/(alami - alam) ;
   D2 = (alam + 2.0*amu)/(alami - alam) ;
   D3 = alam/(alami-alam) ;
%
% g and g' terms (appendix of Tandon and Weng 1984). g is for spheroidal
% inclusions (del>1), whilst g' is for disc-like inclusions (del<1).
%
      if (del >= 1)
       acshdel = log(del + sqrt(t1)) ;
       g =(del*sqrt(t1) - acshdel)*del/(sqrt(t1)^3) ;
      else
%      g' below
       g =(acos(del) - del*sqrt(t5))*del/(sqrt(t5)^3) ;
      end
%
% Eshelby's Sijkl tensor (appendix of Tandon and Weng 1984).
%
       s11 = (t3 + (t4-1.0)/t1 - (t3 + t4/t1)*g)/(2.0*t2)                ;
       s22 = (t4/(t1*2.0) + (t3 - 9.0/(4.0*t1))*g)/(4.0*t2)              ;
       s33 = s22                                                         ;
       s23 = (del*del/(2.0*t1) - (t3 + 3.0/(4.0*t1))*g)/(4.0*t2)         ;
       s32 = s23                                                         ;
       s21 = (-2.0*del*del/t1 + (t4/t1 - t3)*g)/(4.0*t2)                 ;
       s31 = s21                                                         ;
       s12 = (-1.0*(t3 + 1.0/t1) + (t3 + 3.0/(2.0*t1))*g)/(2.0*t2)       ;
       s13 = s12                                                         ;
       s44 = (del*del/(2.0*t1) + (t3 - 3.0/(4.0*t1))*g)/(4.0*t2)         ;
       s66 = (t3 - (t1+2.0)/t1 - (t3 - 3.0*(t1+2.0)/t1)*g/2.0)/(4.0*t2)  ;
       s55 = s66  ;

       esh = [ s11 s12 s13 0.0 0.0 0.0 ; ...
               s21 s22 s23 0.0 0.0 0.0 ; ...
               s31 s32 s33 0.0 0.0 0.0 ; ...
               0.0 0.0 0.0 s44 0.0 0.0 ; ...
               0.0 0.0 0.0 0.0 s55 0.0 ;
               0.0 0.0 0.0 0.0 0.0 s66];
