% Initialisation
clear 
close all

student_ID = 11540435;

%% Defining Parameters
M1 = 9000; % Mass 1
M2 = 10000; % Mass 2
ba1 = 0.5; % Aerodynamic Drag Coefficient 1
ba2 = 0.45;% Aerodynamic Drag Coefficient 2
kc = 5e4; % Coupler's Linear Spring Constant
dc = 1e5; % Coupler's Linear Damping Constant 
knl = 1e4; % Coupler's Nonlinear Spring Constant
Lc = 10; % Natural Length of Coupler

%% Operating Point
% Choosing 360 km/hr as conservative estimate based on current operating mag
% levs

vstar = 360 * 1000 / 3600;                    % Convert speed from km/hr to m/s
vlin = 70;                                    % Linearisation point
max_thrust = M1 * 1 + (ba1 + ba2) * vstar^2;  % 1 m/s² * M1 + drag

%%
% %% Defining States (original 4 states)
% % A matrix
% A = [0, 1, 0, 0;
%      -kc/M1, -(2*ba1*vstar + dc)/M1, kc/M1, dc/M1;
%      0, 0, 0, 1;
%      kc/M2, dc/M2, -kc/M2, -(2*ba2*vstar + dc)/M2];
% % B matrix
% B = [0; 1/M1; 0; 0];
% % C matrix - velocity measurements only
% C = [0,1,0,0;0,0,0,1];
% % D matrix
% D = [0;0];
% 
% %% Building State Space Model
% sys_unreduced = ss(A, B, C, D);
% % Checking controllability
% Co = ctrb(sys_unreduced);
% disp('Controllability rank:'); disp(rank(Co))
% % Rank = 4 Therefore controllable
% % Checking observability
% Ob = obsv(sys_unreduced);
% disp('Observability rank:'); disp(rank(Ob))
% % Rank = 3 Therefore not fully observable - we need to change state
% % definitions to reduce states to 3 (removing absolute positioins for
% % relative)

%% States were reduced by hand and are reimplemented below
A_r = [0,1,-1;
        (-kc/M1), -(2*ba1*vlin+dc)/M1,dc/M1;
        kc/M2,dc/M2,-(2*ba2*vlin+dc)/M2];

B_r = [0;1/M1;0];

C_r = [0,1,0;
       0,0,1];

D_r = [0;0];

sys = ss(A_r, B_r, C_r, D_r);

%% LQR Design

% Bryson's rule - weights based on maximum allowable deviations
max_delta = 0.5;     % max coupler deflection (m)
max_x2    = 0.3333;  % max velocity error carriage 1 (m/s)
max_x4    = 0.3333;  % max velocity error carriage 2 (m/s)

% Feedforward thrust to maintain cruise speed
ustar = (ba1 + ba2) * vstar^2;
max_u = 1.1 * ustar;  % max allowable thrust (N)

Q = diag([(1/max_delta^2)*1000, 1/max_x2^2, 1/max_x4^2]);

R = (1/max_u^2)*600;

%% Discretising the System

% Sampling time chosen to minimise error propagation without instability
Ts = 0.01;
disp('Sample time (s):'); disp(Ts)

% Discretising system with 
sys_rd = c2d(ss(A_r, B_r, C_r, D_r), Ts, 'zoh');
A_d = sys_rd.A;
B_d = sys_rd.B;
C_d = sys_rd.C;

%% Discrete LQR
K_d = dlqr(A_d, B_d, Q, R);
disp('Discrete LQR gain K_d:'); disp(K_d)

%% Discrete Kalman Filter
Q_kf = 0.5*diag([1e-4, 1e-4, 1e-4]);
R_kf = diag([0.01, 0.01]);

% Note - used to show observer poles (not for simulink implementation)
[L_d, ~, ~] = dlqe(A_d, eye(3), C_d, Q_kf, R_kf); 
disp('Discrete Kalman gain L_d:'); disp(L_d)

%% Pole Separation Check
% Observer poles should be faster (smaller magnitude) than controller poles
disp('Controller poles:'); disp(eig(A_d - B_d*K_d))
disp('Observer poles:');   disp(eig(A_d - L_d*C_d))

%% Integral Action

IntGain = 200; % Integral gain, iterated until it stabilised velocity within 2 mins 
SatLimit = 0.12 * ustar;  % Integral contributes max 12% of thrust to prevent overshoot from excessive windup

%% Simulation Parameters
zstar = [0; vstar; vstar];  % desired reduced state [delta*, v1*, v2*]
Tsim  = 120;                 % simulation duration (s)
perturb = 10;                % Initial Speed error
x0 = [0.1; vstar-perturb; vstar-perturb];   % [Lc, v1, v2] - both carriages at cruise, separated by Lc
z_hat0 = [0.1; vstar-perturb; vstar-perturb];   % initial observer guess [delta, v1, v2]

sim_data = sim('Train_Model_Real.slx', 'StopTime', num2str(Tsim));

%% Extract simulation data

% Defining number of frames to cut at so all data has same dimensions
sim_frame_len = size(sim_data.tout, 1) - 8; % Ensuring all data has same length
sim_time = sim_data.tout(2:sim_frame_len);
speed1 = squeeze(sim_data.y(1, 2:sim_frame_len))';
speed2 = squeeze(sim_data.y(2, 2:sim_frame_len))';

U = sim_data.u;
U_time = linspace(sim_time(1), sim_time(end), length(U));
Delta    = speed1 - speed2;           % relative velocity difference

%% Extract xhat
xhat = sim_data.xhat;
xhat_time = (linspace(sim_time(1), sim_time(end), size(xhat, 1)))';

delta_hat  = -interp1(xhat_time, xhat(:,1), sim_time);
speed1_hat = interp1(xhat_time, xhat(:,2), sim_time);
speed2_hat = interp1(xhat_time, xhat(:,3), sim_time);

%% Acceleration using F=MA with clean xhat states

win = 700; % Size of the window to average over

% Calculating accelerations based on filtered speed estimates
accel1 = gradient(movmean(speed1_hat, [win 0]), sim_time);
accel2 = gradient(movmean(speed2_hat, [win 0]), sim_time);

%% Figure 1 — Carriage Speeds
figure(1); clf;
subplot(2,1,1)
plot(sim_time, speed1, 'Color', [0.5 0.7 1], 'LineWidth', 0.8)      
hold on
plot(sim_time, speed1_hat, 'b', 'LineWidth', 1.5)                    
plot(sim_time, ones(size(sim_time))*vstar, 'color', [1, 0.5, 0], 'LineWidth', 1.5,'LineStyle',':')  
xlabel('Time (s)'); ylabel('Speed (m/s)')
title('Carriage 1 Speed vs Cruise')
legend('v1 measured', 'v1 estimated', sprintf('v* = %.0f m/s', vstar)); grid on

subplot(2,1,2)
plot(sim_time, speed2, 'Color', [1 0.6 0.6], 'LineWidth', 0.8)       
hold on
plot(sim_time, speed2_hat, 'r', 'LineWidth', 1.5)                    
plot(sim_time, ones(size(sim_time))*vstar, 'color', [1, 0.5, 0], 'LineWidth', 1.5,'LineStyle',':')   
xlabel('Time (s)'); ylabel('Speed (m/s)')
title('Carriage 2 Speed vs Cruise')
legend('v2 measured', 'v2 estimated', sprintf('v* = %.0f m/s', vstar)); grid on

sgtitle('Closed-Loop Speed Regulation')

%% Figure 2 — Control Input (Thrust)
figure(2); clf;
plot(U_time, U, 'r', 'LineWidth', 1.5)
yline(ustar, 'g--', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('Thrust u (N)')
title('Control Input (Thrust Force)')
legend('u(t) [N]', sprintf('u* = %.0f (Force required to maintain cruise speed [N])', ustar)); grid on

%% Figure 3 — Carriage Accelerations (Observer Estimate)
figure(3); clf;
subplot(2,1,1)
plot(sim_time, accel1, 'b', 'LineWidth', 1.5)
yline(0, 'k--')
xlabel('Time (s)'); ylabel('Acceleration (m/s²)')
title('Carriage 1 Acceleration')
legend('a1'); grid on

subplot(2,1,2)
plot(sim_time, accel2, 'r', 'LineWidth', 1.5)
yline(0, 'k--')
xlabel('Time (s)'); ylabel('Acceleration (m/s²)')
title('Carriage 2 Acceleration')
legend('a2'); grid on

sgtitle('Closed-Loop Acceleration')

%% Figure 4 - Coupler Deflection From Equilibrium (x1-x3-Lc)
figure(4); clf;
plot(sim_time, delta_hat, 'r', 'LineWidth', 1.5)
yline(0, 'k--')
xlabel('Time (s)'); ylabel('Coupler Deflection (m)')
title('Coupler Deflection \delta = x_1 - x_3 - L_c')
legend('\delta hat'); grid on

