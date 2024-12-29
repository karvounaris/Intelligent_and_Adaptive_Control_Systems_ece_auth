% =========================================================================
% We use adaptive control in order to create an adaptive controller to make
% the system to behave like a desired system with desired characteristics
% in the closed loop situation. We have a multi-input multi-output
% linear system in the form x_dot = A*x + B*Lamda*(u + f(x));
%
% This script is the main control simulation for the system that uses the
% following control law: u = theta^T + theta_dot^T * phi
%
% After the system reach balance, a disturbance d(t) applies
% to the system for a total duration on 5 seconds
%
% VARIALBES:
%   - M:  (N*m*s^2)/rad
%   - G:  N*m
%   - C:  (N*m*s)/rad
%   - x1: rad
%   - x2: rad/s
%   - u:  N*m
%
% Formulas used:
%   - x_ref_dot = A_ref * x_ref + B_ref * r(t)
%   - y_ref = C_ref.' * x_ref
%   - u = theta' * omega + theta_dot' * phi
%   - x_dot = A*x + B*Lamda*(u + f(x))
%   - x_dot = A*x + B*Lamda*(u + f(x) + d(t))
%   - y = C.' * x
%   - omega_1_dot = F * omega_1 + g * u
%   - omega_2_dot = F * omega_2 + g * y
%   - omega = [omega_1 omega_2 y r(t)]
%   - phi_dot = - P_0 * phi + omega
%   - theta_dot = - Gamma * e * phi * sgn(k_p/k_m)
%
% =========================================================================

%% Reference Model

clear;
close all;

% Simulation parameters
timestep = 0.001;
simulation_duration = 30;           % Total duration of simulation
N = simulation_duration / timestep; % Number of steps
t = linspace(0, simulation_duration, N);

% Reference Model parameters
A_ref = [0 1; -1 -2];
B_ref = [0; 1];
C_ref = [1; 0];
x_ref(:,1) = [0 0];
y_ref = C_ref(:,1).' * x_ref(:,1);
k_m = 1; 

% System's input
r = @(t) 2 * sin(1 * t);

% Calculate the x_ref_dot and then use this value to calculate x_ref
for i = 1:(length(t) - 1)
    x_ref_dot = A_ref * x_ref(:,i) + B_ref * r(t(i+1));
    x_ref(:,i+1) = x_ref(:,i) + x_ref_dot * timestep;
    y_ref(i+1) = C_ref.' * x_ref(:,i+1);
end

%% Real System

% Real System parameters
M = 1;
G = 10;
C = 1;
A = [0 1;
     0 -C/M];
B = [0; 1];
C = [1; 0];
Lamda = 1/M;
theta_star = -G;
phi_sin = @(t) sin(t);
phi_function = @(t) theta_star*phi_sin(t);
x(:,1) = [0 0];
y = C.' * x(:,1);
k_p = 1/M;                           % Y(s) = [1/M]*[1/(s^2+Cs/M+G/M)]U(s)

% Initialize equations' parameters. Changable parameters that impact the
% simulation.
G_param = 1000;
Gamma = diag([G_param G_param G_param G_param]);
P_0 = 1;                            % 0 < P_0 < 2
lamda_0 = 5;
F = -lamda_0;                       % F = -lamda_0(s), lamda_0 = constant, lamda_0 > 0

% Initialize equations' parameters. Not changable.
phi(:,1) = [0 0 0 0];
theta(:,1) = [0 0 0 0];
theta_dot = [0; 0; 0; 0];
omega(:,1) = [0 0 0 0];
omega_1 = 0;
omega_2 = 0;
omega(:,1) = [omega_1 omega_2 y r(t(1))];
u = 0;
g = 1;
e = y - y_ref;

% % Disturbance
% w = 2;
% a = 50;
% d = @(t) a * sin(w*t);

% Ad = 50;
% fs = 1000; 
% td = 0:1/fs:5; 
% fd = 2;           % frequency
% Td = 1/fd;        % period
% d = Ad*sign(sin(2*pi*fd*td));

d = 10;

% Time management parameters
disturbance_timestep_counter = 0;
disturbance_timestep_duration = 5 / timestep;   % disturbance duration in timesteps
balance_timestep_counter = 0;
balance_timestep_duration = 2 / timestep;       % duration of balance needed before activate disturbances in timesteps
error_limit = 0.002;
start_counting_parameter = 120;                 % exclude the first buch of iterations from the balance counter, due to initialization
stop_counting_parameter = 0;                    % when it is equal to 0 keep counting , if it is equal to 1 stop counting

% Implement the formulas for the simulation
for i = 1:(length(t) - 1)
    if i>start_counting_parameter  && abs(e(i))<error_limit && abs(e(i-1))<error_limit && abs(e(i-2))<error_limit && stop_counting_parameter==0
        balance_timestep_counter = balance_timestep_counter + 1;
    end

    omega_1_dot = F * omega_1(i) + g * u;
    omega_1(i+1) = omega_1(i) + omega_1_dot * timestep;
    
    omega_2_dot = F * omega_2(i) + g * y(i);
    omega_2(i+1) = omega_2(i) + omega_2_dot * timestep;
    
    omega(:,i+1) = [omega_1(i+1) omega_2(i+1) y(i) r(t(i+1))];
    
    phi_dot = - P_0 * phi(:,i) + omega(:,i+1);
    phi(:,i+1) = phi(:,i) + phi_dot * timestep;
    
    theta_dot = - Gamma * e(i) * phi(:,i+1) * sign(k_p/k_m);
    theta(:,i+1) = theta(:,i) + theta_dot * timestep;
    
    u = theta(:,i+1)' * omega(:,i+1) + theta_dot' * phi(:,i+1);
    
    if balance_timestep_counter==balance_timestep_duration && disturbance_timestep_counter~=disturbance_timestep_duration
        stop_counting_parameter = 1;
        % x_dot = A*x(:,i) + B*Lamda*(u + phi_function(x(1,i)) + d(t(i)));
        disturbance_timestep_counter = disturbance_timestep_counter + 1;
        % x_dot = A*x(:,i) + B*Lamda*(u + phi_function(x(1,i)) + d(disturbance_timestep_counter));
        x_dot = A*x(:,i) + B*Lamda*(u + phi_function(x(1,i)) + d);
    else
        x_dot = A*x(:,i) + B*Lamda*(u + phi_function(x(1,i)));
    end

    x(:,i+1) = x(:,i) + x_dot * timestep;
    y(i+1) = C.' * x(:,i+1);
    e(i+1) = y(i+1) - y_ref(i+1);

end

%% Plot Section

% % Plot ouput of reference model and real system in the same plot
% figure()
% plot(t, y_ref(1,:), t, y(1,:));
% title('Rotation Angle, Real-Reference');
% ylabel('angle (rad)');
% xlabel('time (s)');
% legend('y\_ref', 'y');

% Plot output's errors 
figure()
plot(t, e);
title('Error y');
ylabel('angle (rad)'); 
xlabel('time (s)');
legend('y\_error');
