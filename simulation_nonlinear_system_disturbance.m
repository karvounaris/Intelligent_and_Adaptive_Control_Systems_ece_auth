% =========================================================================
% We use adaptive control in order to create an adaptive controller to make
% the system to behave like a desired system with desired characteristics
% in the closed loop situation. We have a multi-input multi-output
% nonlinear system in the form x_dot = A*x + B*Lamda*(u + f(x)), f(x) is a
% nonlinear vector function: f(x) = Theta^(T) * Phi(x)
%
% This script is the control simulation for the system that uses the
% following control law: u = - K_hat^(T)*x - L_hat^(T)*r - Theta_hat^(T)*Phi(x)
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
%   - x_ref_dot = A_ref*x_ref + B_ref*r(t)
%   - x_dot = A_ref*x + B_ref*r(t) + B*Lamda*(-k_bar.'*x - l_bar*r(t) - theta_bar*sin(x1))
%   - x_dot = A_ref*x + B_ref*r(t) + B*Lamda*(-k_bar.'*x - l_bar*r(t) - theta_bar*sin(x1) + d(t))
%   - u = - K_hat^(T)*x - L_hat^(T)*r - Theta_hat^(T)*Phi(x)
%   - x_dot = A*x + B*Lamda*(u + f(x))
%   - x_dot = A*x + B*Lamda*(u + f(x) + d(t))
%   - f(x) = theta_start^T * Phi(x)
%   - k_hat_dot = gamma_1*B.'*P*e*x
%   - l_hat_dot = gamma_2*B.'*P*e*r(t)
%   - theta_hat_dot = gamma_3*B.'*P*e(:,i)*sin(x1)
%   - k_bar = k_hat - k_star
%   - l_bar = l_hat - l_star
%   - theta_bar = theta_hat - theta_star
%   - e = x - x_ref
%
% =========================================================================

%% Reference Model

clear;
close all
% Simulation parameters
timestep = 0.001;
simulation_duration = 30;           % Total duration of simulation
N = simulation_duration / timestep; % Number of steps
t = linspace(0, simulation_duration, N);

% Reference Model parameters 
A_ref = [0 1; -1 -1.4];
B_ref = [0; 1];
x_ref(:,1) = [0 0];

% System's input
r = @(t) 4 * sin(t);

% Calculate the x_ref_dot and then use this value to calculate x_ref
for i = 1:(length(t) - 1)
    x_ref_dot = A_ref * x_ref(:,i) + B_ref * r(t(i+1));
    x_ref(:,i+1) = x_ref(:,i) + x_ref_dot * timestep;
end

%% Real System: 

% Real System parameters
M = 1;
G = 10;
C = 1;
x(:,1) = [0 0];
A = [0 1;
     0 -C/M];
B = [0; 1];
Lamda = 1/M;
theta_star = -G;
k_star = [M; 1.4*M-C];
l_star = -M;
phi_sin = @(t) sin(t);
phi_function = @(t) theta_star*phi_sin(t);

% Initialize estimations
k_hat = [0; 0];
l_hat = 0;
theta_hat = 0;

% Initialize equations' parameters. Changable parameters that impact the
% simulation.
gamma_1 = 1000;
gamma_2 = 1000;
gamma_3 = 1000;
Q = [150 0;
     0 150];
P = lyap(A_ref', Q);

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
balance_timestep_duration = 5 / timestep;       % duration of balance needed before activate disturbances in timesteps
error_limit = 0.001;
start_counting_parameter = 120;                 % exclude the first buch of iterations from the balance counter, due to initialization
stop_counting_parameter = 0;                    % when it is equal to 0 keep counting , if it is equal to 1 stop counting

% Initialize bar values
k_bar(:,1) = k_hat(:,1) - k_star;
l_bar = l_hat - l_star;
theta_bar = theta_hat - theta_star;

% Use the initialized parameters to find the second value of x and error
e(:,1) = x(:,1) - x_ref(:,1);
% x_dot = A_ref*x(:,1) + B_ref*r(t(1)) + B*Lamda*(-k_bar(:,1).'*x(:,1) - l_bar*r(t(1)) - theta_bar*phi_sin(x(1,1))); 
u = - k_hat(:,1).'*x(:,1) - l_hat(1)*r(t(1)) - theta_hat(1)*phi_sin(x(1,1));
x_dot = A*x(:,1) + B*Lamda*u;
x(:,2) = x(:,1) + x_dot * timestep;
e(:,2) = x(:,2) - x_ref(:,2);

% Implement the formulas for the simulation
for i = 2:(length(t) - 1)
    if i>start_counting_parameter  && abs(e(i))<error_limit && abs(e(i-1))<error_limit && abs(e(i-2))<error_limit && stop_counting_parameter==0
        balance_timestep_counter = balance_timestep_counter + 1;
    end

    k_hat_dot = gamma_1 * B.' * P * e(:,i) * x(:,i);
    k_hat(:,i) = k_hat(:,i-1) + k_hat_dot * timestep;
    
    l_hat_dot = gamma_2 * B.' * P * e(:,i) * r(t(i));
    l_hat(i) = l_hat(i-1) + l_hat_dot * timestep;
    
    theta_hat_dot = gamma_3 * B.' * P * e(:,i) * phi_sin(x(1,i));
    theta_hat(i) = theta_hat(i-1) + theta_hat_dot * timestep;
    
    % k_bar(:,i) = k_hat(:,i) - k_star;
    % l_bar(i) = l_hat(i) - l_star;
    % theta_bar(i) = theta_hat(i) - theta_star;
    
    u = - k_hat(:,i).'*x(:,i) - l_hat(i)*r(t(i)) - theta_hat(i)*phi_sin(x(1,i));

    % if balance_timestep_counter==balance_timestep_duration && disturbance_timestep_counter~=disturbance_timestep_duration
    %     stop_counting_parameter = 1;
    %     x_dot = A_ref*x(:,i) + B_ref*r(t(i)) + B*Lamda*(-k_bar(:,i).'*x(:,i) - l_bar(i)*r(t(i)) - theta_bar(i)*phi_sin(x(1,i)) + d(t(i))); 
    %     disturbance_timestep_counter = disturbance_timestep_counter + 1;
    % else
    %     x_dot = A_ref*x(:,i) + B_ref*r(t(i)) + B*Lamda*(-k_bar(:,i).'*x(:,i) - l_bar(i)*r(t(i)) - theta_bar(i)*phi_sin(x(1,i))); 
    % end

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
    e(:,i+1) = x(:,i+1) - x_ref(:,i+1);
end

%% Plotting section

% Plot angle's errors 
figure();
subplot(2, 1, 1); % This creates a subplot grid of 2 rows and 1 column, and selects the 1st subplot
plot(t, e(1,:));
title('Rotation Angle Error');
ylabel('angle (rad)'); 
xlabel('time (s)');
legend('e_1');

% Plot angular velocity's errors 
subplot(2, 1, 2); % This selects the 2nd subplot in the grid
plot(t, e(2,:));
title('Angular Velocity Error');
ylabel('angular velocity (rad/s)');
xlabel('time (s)');
legend('e_2');

% % Plot rotation angle of reference model and real system in the same plot
% figure()
% plot(t, x_ref(1,:), t, x(1,:));
% title('Rotation Angle, Real-Reference');
% ylabel('angle (rad)');
% xlabel('time (s)');
% legend('x1\_ref', 'x1');
% 
% % Plot rotation angle of reference model and real system in the same plot
% figure()
% plot(t, x_ref(2,:), t, x(2,:));
% title('Angular Velocity, Real-Reference');
% ylabel('angular velocity (rad/s)');
% xlabel('time (s)');
% legend('x2\_ref', 'x2');
