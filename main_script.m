% Define parameters (these may need further tuning)
Ms = 0.0325; % Mass of slider
m = 0.0014; % Mass of moving part
Mu_d = 0.257; % Friction coefficient
g = 9.81; % Gravity
lambda = 0.00059; % Some constant
rr = 0.0005; % Radius
L = 0.880; % Length
omega = 3.13; % Angular velocity
M = 0.0325;

% Set initial conditions, the ode45 requires an initial condition to
% integrate with and this is stored in the XO vector.

theta0 = 1.42; % Initial angle
dtheta0 = 1.42;
y0 = 0;
dy0 = 0;
X0 = [theta0; dtheta0; y0; dy0];

% Set time span
tspan = [0 5]; % Simulation time

% Solve ODE using ode45
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
%RelTol (Relative Tolerance) Controls the relative error tolerance for the solver. It sets the accuracy with which the solution is computed. Smaller values make the solver work harder for higher accuracy.
%AbsTol (Absolute Tolerance): Controls the absolute error tolerance. It
%specifies the threshold below which the error is acceptable. Smaller
%values ensure higher accuracy in absolute terms. in this, the value has
%been standardized. 

[t, X] = ode45(@(t,X) odefun(t, X, Ms, m, Mu_d, g, lambda, rr, L, omega, M), tspan, X0, options);

% Extract solutions
theta = X(:,1);
y = X(:,3);

% Calculate the position of the lighter mass
l = L + y - rr*(pi + theta);
x = l .* sin(theta);
y_mass = -l .* cos(theta);

coordinates = [x, y];
writematrix(coordinates, 'coordinates.txt');

% Plot results
figure;
plot(x, y_mass);
title('Trajectory of Lighter Mass');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
hold on;
plot(x(1), y_mass(1), 'ro', 'MarkerSize', 10);
plot(x(end), y_mass(end), 'go', 'MarkerSize', 10);
plot([0, x(end)], [0, y_mass(end)], 'k--');
legend('Trajectory', 'Start', 'End', 'Final Rod Position', 'Location', 'best');
xlim([-1 1]);
ylim([-0.5 0.3]);

% Display final values
fprintf('Final theta: %.4f rad\n', theta(end));
fprintf('Final y: %.4f m\n', y(end));
fprintf('Final x position of mass: %.4f m\n', x(end));
fprintf('Final y position of mass: %.4f m\n', y_mass(end));

% Print out coordinates
fprintf('\nCoordinates (x, y):\n');
for i = 1:length(x)
    fprintf('(%.4f, %.4f)\n', x(i), y_mass(i));
end

