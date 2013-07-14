

function [ C ] = build_isotropic( K, G )

    % Get Lam\'e constants
    lam = K-(2/3)*G; 
    mu = G;

    % build C
    c1 = 2*mu+lam;
    c2 = lam;
    c3 = mu;

    C = [c1  c2  c2  0.0 0.0 0.0; ...
         c2  c1  c2  0.0 0.0 0.0; ...
         c2  c2  c1  0.0 0.0 0.0; ...
         0.0 0.0 0.0 c3  0.0 0.0; ...
         0.0 0.0 0.0 0.0 c3  0.0; ...
         0.0 0.0 0.0 0.0 0.0 c3 ];

     
end