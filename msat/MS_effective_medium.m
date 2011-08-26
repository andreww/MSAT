% [Ceff,rh]=MS_effective_medium(theory, ...)
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%  
%  Generate elastic constants for the effective medium parameters based on 
%  various theories, identified by the string 'theory'. Subsequent required 
%  arguments depend on the theory invoked. 
%
%  Currently available theories:
%
%  ** Tandon and Weng, 84 ('tandon' or 't&w')
%     (Isotropic host matrix with unidirectionally aligned isotropic spheroid 
%      inclusions) 
%
%     [Ceff,rh]=MS_effective_medium('t&w',vpm,vsm,rhm,vpi,vsi,rhi,del,f) or
%        Input parameters:
%           vpm,vsm,rhm : isotropic parameters of the matrix (km/s, kg/m3)
%                   del : aspect ratio of spheroids:
%                            del < 1 = smarties
%                            del = 1 = spheres
%                            del > 1 = cigars
%                    f : volume fractions of inclusions       
%        vpi, vsi, rhi : isotropic parameters of the inclusions
%
%     [Ceff,rh]=MS_effective_medium('t&w',Cm,rhm,Ci,rhi,del,f) or
%        Input parameters:
%             Cm,rh : elasticity and density of the matrix (GPa, kg/m3)
%                   del : aspect ratio of spheroids:
%                            del < 1 = smarties
%                            del = 1 = spheres
%                            del > 1 = cigars
%                    f : volume fractions of inclusions       
%         Ci, rhi : isotropic parameters of the inclusions
%
%       Output parameters:
%          Ceff : Elastic constants (GPa) (symmetry in X1 direction)
%          rh : aggregate density (kg/m3)
%
%
%  References:
% - Tandon, GP and Weng, GJ. The Effect of Aspect Ratio of Inclusions on the 
%   Elastic Properties of Unidirectionally Aligned Composites. Polymer 
%   Composites, 5, pp 327-333, 1984. 
%
%
% See also: MS_elasticDB

function [Ceff,rh]=MS_effective_medium(theory, varargin) ;

if ~ischar(theory)
   error('MS:EFFECTIVE_MEDIUM:BadTString', ...
   'A string specifying the theory to use is required.') ;
end

switch lower(theory)
%-------------------------------------------------------------------------------
   case {'tandon', 't&w'}
      if length(varargin)~=8 & length(varargin)~=6
         error('MS:EFFECTIVE_MEDIUM:TWWrongArgs', ...
         'Tandon and Weng (1984) requires 6 or 8 input parameters.') ;
      end
      if length(varargin)==6 % elasticity matrix form
         Cm = varargin{1} ; rhm = varargin{2} ; 
         Ci = varargin{3} ; rhi = varargin{4} ;
         del = varargin{5} ; f = varargin{6} ;
%     ** check the matrices
         MS_checkC(Cm) ;
         if MS_anisotropy( Cm) > 10*sqrt(eps) % not isotropic           
            error('MS:EFFECTIVE_MEDIUM:TWBadC', ...
               'Tandon and Weng (1984) requires isotropic inputs matrices.') ;
         end
         MS_checkC(Ci) ;
         if MS_anisotropy( Ci ) > 10*sqrt(eps) % not isotropic           
            error('MS:EFFECTIVE_MEDIUM:TWBadC', ...
               'Tandon and Weng (1984) requires isotropic input matrices.') ;
         end
%     ** unload the 
         vpm = sqrt(Cm(3,3)*1e3./rhm) ;  vsm = sqrt(Cm(6,6)*1e3./rhm) ;
         vpi = sqrt(Ci(3,3)*1e3./rhi) ;  vsi = sqrt(Ci(6,6)*1e3./rhi) ;
      else % velocity form
         vpm = varargin{1} ; vsm = varargin{2} ; rhm = varargin{3} ; 
         vpi = varargin{4} ; vsi = varargin{5} ; rhi = varargin{6} ;
         del = varargin{7} ; f = varargin{8} ;
      end
      [Ceff,rh]=tandon_and_weng(vpm,vsm,rhm,vpi,vsi,rhi,del,f) ;
%-------------------------------------------------------------------------------
   otherwise
      error('MS:EFFECTIVE_MEDIUM:UnknownTheory', ...
         'Specified theory is not supported.') ;
%-------------------------------------------------------------------------------
end % of switch

return



