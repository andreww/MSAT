function [ceff_hud, ceff_pade] = cheng_core(gamma, mu, phi, alpha, fluid_k)

    % Crack density, below cheng e.q. 10, eps in the paper
    cd = (3.0*phi)/(4.0*pi*alpha);

    % Cheng e.q. 7
    c_0_11 = gamma + 2.0 * mu;
    c_0_33 = c_0_11;
    c_0_13 = gamma;
    c_0_44 = mu;
    c_0_66 = c_0_44;

    % Assume weak fluid for now
    [U1, U3] = hudson_weak_fluid(gamma, mu, fluid_k, alpha);

    % Chrng e.q. 8
    c_1_11 = -1.0*((gamma^2)/mu)*cd*U3;
    c_1_13 = -1.0*((gamma*(gamma+2.0*mu))/mu)*cd*U3;
    c_1_33 = -1.0*(((gamma*2.0*mu)*2.0)/mu)*cd*U3; % Should that last 2.0 be a power?
    c_1_44 = -1.0*mu*cd*U1;
    c_1_66 = 0.0;

    % Cheng e.q. 10
    q = 15.0*(gamma^2.0/mu^2) + 15.0*(gamma/mu) + 28.0;

    % Cheng e.q. 9
    c_2_11 = (q/15.0)*(gamma^2/(gamma+2.0*mu))*(cd*U3)^2;
    c_2_13 = (q/15.0)*gamma*(cd*U3)^2;
    c_2_33 = (q/15.0)*(gamma+2.0*mu)*(cd*U3)^2;
    c_2_44 = (2/15.0)*((mu*(3.0*gamma+8.0*mu))/(gamma+2.0*mu))*(cd*U1)^2; % Is that really 2 not q at the start? It's 2 in the paper
    c_2_66 = 0.0;

    % Cheng equation 6
    ceff_hud = zeros(6,6);
    ceff_hud(1,1) = c_0_11 + c_1_11 + c_2_11;
    ceff_hud(1,3) = c_0_13 + c_1_13 + c_2_13;
    ceff_hud(3,3) = c_0_33 + c_1_33 + c_2_33;
    ceff_hud(4,4) = c_0_44 + c_1_44 + c_2_44;
    ceff_hud(6,6) = c_0_66 + c_1_66 + c_2_66;
    % Symmetry
    ceff_hud(2,2) = ceff_hud(1,1);
    ceff_hud(5,5) = ceff_hud(4,4);
    ceff_hud(2,3) = ceff_hud(1,3);
    ceff_hud(1,2) = ceff_hud(1,1) - 2.0*ceff_hud(6,6);
    ceff_hud(3,1) = ceff_hud(1,3);
    ceff_hud(3,2) = ceff_hud(2,3);
    ceff_hud(2,1) = ceff_hud(1,2);

    % eq 21
    ceff_pade = zeros(6,6);
    ceff_pade(1,1) = pade(c_0_11, c_1_11, c_2_11, cd);
    ceff_pade(1,3) = pade(c_0_13, c_1_13, c_2_13, cd);
    ceff_pade(3,3) = pade(c_0_33, c_1_33, c_2_33, cd);
    ceff_pade(4,4) = pade(c_0_44, c_1_44, c_2_44, cd);
    ceff_pade(6,6) = c_0_66;
    % Symmetry
    ceff_pade(2,2) = ceff_pade(1,1);
    ceff_pade(5,5) = ceff_pade(4,4);
    ceff_pade(2,3) = ceff_pade(1,3);
    ceff_pade(1,2) = ceff_pade(1,1) - 2.0*ceff_pade(6,6);
    ceff_pade(3,1) = ceff_pade(1,3);
    ceff_pade(3,2) = ceff_pade(2,3);
    ceff_pade(2,1) = ceff_pade(1,2);

end

function cij = pade(c0ij, c1ij, c2ij, cd)

    b = c2ij/(c1ij*cd);
    a = (c1ij/(c0ij*cd)) - b;

    cij = c0ij*((1.0-a*cd)/(1.0+b*cd));

end

function [U1, U3] = hudson_fluid_filled(gamma, mu)

    % Cheng e.q. 11
    [U1, ~] = hudson_dry(gamma, mu);
    U3 = 0.0;

end

function [U1, U3] = hudson_dry(gamma, mu)

    % Cheng e.q. 12
    U1 = (16.0*(gamma+2.0*mu)) / (3.0*(3.0*gamma + 4.0*mu));
    U3 = (4.0*(gamma+2.0*mu)) / (3.0*(3.0*gamma+4.0*mu));

end

function [U1, U3] = hudson_weak_fluid(gamma, mu, kf, alpha)

    % Cheng e.q. 14
    K = (kf*(gamma+2.0*mu))/(pi*alpha*mu*(gamma+mu));

    [U1, U3] = hudson_dry(gamma, mu);

    U3 = U3 * (1.0/(1.0+K));

end