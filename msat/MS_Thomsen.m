
function [eps,gam,del] = MS_Thomsen(C)
% 	MS_THOMSEN   
% 		[EPS,GAM,DEL] = MS_THOMSEN(C,RH)
% 
% 	Calculates Thomsen parameters from the elasticity matrix
% 	
% 	Created by Alan Baird on 2012-02-07.
% 	Copyright (c)  . All rights reserved.

eps=(C(1,1)-C(3,3))/(2.0*C(3,3));
gam=(C(6,6)-C(4,4))/(2.0*C(4,4));
del=((C(1,3)+C(4,4))^2.0-(C(3,3)-C(4,4))^2.0)/(2.0*C(3,3)*(C(3,3)-C(4,4)));

end %  function

