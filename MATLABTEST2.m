% Defining the parameters
Ms = 0.0004; % Mass of string
m = 0.0039; % Mass of moving part (lighter bob m)
Mu_d = 0.257; % Friction coefficient (dynamic friction)
g = 9.80665; % Gravity (constant)
lambda = 0.00059; % Linear density of the string 
rr = 0.0003; % Radius of the cylindrical rod
L = 0.500; % Length of the moving bob from the rod
M = 0.0212; % Mass of heavier pendulum bob 
Mm_ratio = M / m; % Mass ratio of heavier and lighter mass
omega=1.5;
% Setting initial conditions
theta0 = 1.57079632679; % Initial angle
dtheta0 = 0; % Initial angular velocity
y0 = 0; % Initial y position
dy0 = 0; % Initial y velocity
X0 = [theta0; dtheta0; y0; dy0]; % Initial conditions vector

% Set time span
tspan = [0 0.37]; % Simulation time (changing as mentioned in the paper)

% Solve ODE using ode45 (built in matlab function to solve the odes)
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, X] = ode45(@(t,X) odefun(t, X, Ms, m, Mu_d, g, lambda, rr, L, omega, M), tspan, X0, options);

% Extract solutions 
theta = X(:,1);
y = X(:,3);

% Calculate the position of the lighter mass using trignometric functions -
l = L + y - rr * (pi + theta);
x = l .* sin(theta);
y_mass = -l .* cos(theta);

% Calculate the vertical displacement of the heavier bob (for the research
% question)
y_heavy = y;

% Save coordinates to file (so that it can be imported into excel)
coordinates = [x, y_mass];
writematrix(coordinates, 'coord.xlsx');
y_final = y(end); % Final vertical position of the lighter mass (final valye that was used in the data table)

heavier_mass_coordinates = [t, y_heavy];
writematrix(heavier_mass_coordinates, 'heaviermass.xlsx');
% Plot the trajectory of the lighter mass
figure;
subplot(2, 1, 1); % Create a subplot for the lighter mass
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

% Plot the vertical displacement of the heavier mass
subplot(2, 1, 2); % Create a subplot for the heavier mass
plot(t, y_heavy);
title('Vertical Displacement of Heavier Mass');
xlabel('Time (s)');
ylabel('Vertical Displacement (m)');
grid on;

% Display final values
fprintf('Final theta: %.4f rad\n', theta(end));
fprintf('Final y: %.4f m\n', y(end));
fprintf('Final x position of mass: %.4f m\n', x(end));
fprintf('Final y position of mass: %.4f m\n', y_mass(end));
fprintf('Vertical distance traveled by heavier mass M: %.4f m\n', y_final);

% Print out coordinates
fprintf('\nCoordinates (x, y):\n');
for i = 1:length(x)
    fprintf('(%.4f, %.4f)\n', x(i), y_mass(i));

end

%AVIGHNA DARUKA ST YAU 2024 RESEARCH COMPETITION
