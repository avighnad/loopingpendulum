function dYdt = odefun(t, Y, Ms, m, Mu_d, g, lambda, rr, L, omega,M)
    theta = Y(1);
    dtheta = Y(2);
    y = Y(3);
    dy = Y(4);
    %over here, the second order differential equations have been
    %represented as a vector where the first index is the angular position,
    %second as its angualar acceleration (with respect to time), and the
    %same for the y position (vertical position and acceleration) 
   
    %below all the functions have been listed out, with which the ode45 in
    %built matlab function solving these equations


    % Calculate l
    l = L + y - rr*(3.14159 + theta);
    
    % Calculate dl_dt
    dl_dt = dy - rr*dtheta;
    
    % First equation
    d2theta = (-3/(2*l)) * g * sin(theta) * ((Ms + 2*m)/(Ms + 3*m)) - (2/l) * dtheta * dl_dt;
    
    % Second equation
    d2y = (Mu_d*m*g*cos(theta)*exp(Mu_d*(pi+theta)) + Mu_d*m*l*omega^2*exp(Mu_d*(pi+theta)) - Mu_d*M*g) / ...
          (M*Mu_d + lambda*rr*exp(Mu_d*(pi+theta)) - lambda*rr + m*Mu_d*exp(Mu_d*(pi+theta)));
    
    dYdt = [dtheta; d2theta; dy; d2y];
    %avighna daruka, ST YAU 2024.
end

%AVIGHNA DARUKA ST YAU 2024 RESEARCH COMPETITION