% Initialisation
clear all
close all

student_ID = 11540435;

%% Defining Parameters
M1 = 9000;
M2 = 10000;
ba1 = 0.5;
ba2 = 0.45;
kc = 5e4;
dc = 1e5;
knl = 1e4;
Lc = 10;

%% Operating Point
% Choosing 400km/hr as conservative estimate based on current operating mag
% levs

vstar = 400 * 1000 / 3600; % Convert speed from km/hr to m/s

%% Defining States
% A matrix
A = [0, 1, 0, 0;
     -kc/M1, -(2*ba1*vstar + dc)/M1, kc/M1, dc/M1;
     0, 0, 0, 1;
     kc/M2, dc/M2, -kc/M2, -(2*ba2*vstar + dc)/M2];
% B matrix
B = [0; 1/M1; 0; 0];
% C matrix - velocity measurements only
C = [0,1,0,0;0,0,0,1];
% D matrix
D = [0;0];

%% Building State Space Model
sys = ss(A, B, C, D);
% Checking controllability
Co = ctrb(sys);
disp('Controllability rank:'); disp(rank(Co))
% Rank = 4 Therefore controllable
% Checking observability
Ob = obsv(sys);
disp('Observability rank:'); disp(rank(Ob))
% Rank = 3 Therefore not fully observable - we need to change state
% definitions to reduce states to 3 (removing absolute positioins for
% relative)

%% Reducing States to Ensure Full Observability
% Absolute position is unobservable from velocity sensors alone
% Redefine states as [delta, v1, v2] where delta = x1-x3 (coupler deflection)
% Velocity regulation only requires relative position - absolute position is irrelevant
T = [ 1,  0, -1,  0;   % delta = x1 - x3
      0,  1,  0,  0;   % v1
      0,  0,  0,  1];  % v2

% Project system matrices into reduced coordinates
A_r = T * A * pinv(T);   % 3x3
B_r = T * B;             % 3x1
C_r = C * pinv(T);       % becomes [0,1,0; 0,0,1]

% Verify full observability on reduced system
disp('Reduced observability rank:'); disp(rank(obsv(A_r, C_r)))
% Rank = 3 Therefore fully observable

%% LQR Design
% Bryson's rule - weights based on maximum allowable deviations
max_delta = 0.5;     % max coupler deflection (m)
max_x2    = 0.3333;  % max velocity error carriage 1 (m/s)
max_x4    = 0.3333;  % max velocity error carriage 2 (m/s)

% Feedforward thrust to maintain cruise speed
ustar = (ba1 + ba2) * vstar^2;
max_u = 1.1 * ustar;  % max allowable thrust (N)

Q = diag([(1/max_delta^2)*10, 1/max_x2^2, 1/max_x4^2]);
R = (1/max_u^2)*80;

%% Discretising the System
% Sample time chosen as 20x faster than the fastest continuous pole
fastest_pole = max(abs(real(eig(A_r))));
Ts = 1 / (20 * fastest_pole);
disp('Sample time (s):'); disp(Ts)

sys_rd = c2d(ss(A_r, B_r, C_r, zeros(2,1)), Ts, 'zoh');
A_d = sys_rd.A;
B_d = sys_rd.B;
C_d = sys_rd.C;

%% Discrete LQR
K_d = dlqr(A_d, B_d, Q, R);
disp('Discrete LQR gain K_d:'); disp(K_d)

%% Discrete Kalman Filter
% Bryson's rule - weights based on expected noise magnitudes
max_delta_noise = 0.01;  % coupler position uncertainty (m)
max_v1_noise    = 0.05;  % velocity process noise (m/s)
max_v2_noise    = 0.05;
Q_kf = diag([1/max_delta_noise^2, 1/max_v1_noise^2, 1/max_v2_noise^2]);

% Velocity sensor accuracy
max_v1_sensor = 0.08;    % (m/s)
max_v2_sensor = 0.08;
R_kf = diag([1/max_v1_sensor^2, 1/max_v2_sensor^2]);

[L_d, ~, ~] = dlqe(A_d, eye(3), C_d, Q_kf, R_kf);
disp('Discrete Kalman gain L_d:'); disp(L_d)

%% Pole Separation Check
% Observer poles should be faster (smaller magnitude) than controller poles
disp('Controller poles:'); disp(eig(A_d - B_d*K_d))
disp('Observer poles:');   disp(eig(A_d - L_d*C_d))

%% Simulation Parameters
zstar = [Lc; vstar; vstar];  % desired reduced state [delta*, v1*, v2*]
Tsim  = 45;                 % simulation duration (s)


%% Running the simulink file
% Initial Conditions
x0 = [Lc; 0; 0];   % [Lc, v1, v2] - both carriages at cruise, separated by Lc
z_hat0 = [Lc; vstar; vstar];   % initial observer guess [delta, v1, v2]
sim_data = sim('Train_Model_Real.slx', 'StopTime', num2str(Tsim));


%% Extract simulation data

sim_time = sim_data.tout;
Speed1   = squeeze(sim_data.y(1,:))';
Speed2   = squeeze(sim_data.y(2,:))';
U        = squeeze(sim_data.u(1,:))';
Xhat     = squeeze(sim_data.xhat(1,:))';   % [delta_hat, v1_hat, v2_hat]

Delta    = Speed1 - Speed2;           % relative velocity difference
Delta_hat = Xhat(:,1);               % estimated coupler deflection

%% Figure 1 — Carriage Speeds
figure(1); clf;
subplot(2,1,1)
plot(sim_time, Speed1, 'b', sim_time, ones(size(sim_time))*vstar, 'k--', 'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('Speed (m/s)')
title('Carriage 1 Speed vs Cruise')
legend('v1', 'v* = 75 m/s'); grid on

subplot(2,1,2)
plot(sim_time, Speed2, 'r', sim_time, ones(size(sim_time))*vstar, 'k--', 'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('Speed (m/s)')
title('Carriage 2 Speed vs Cruise')
legend('v2', 'v* = 75 m/s'); grid on

sgtitle('Closed-Loop Speed Regulation')

%% Figure 4 — Control Input (Thrust)
figure(4); clf;
plot(sim_time, U, 'k', 'LineWidth', 1.5)
yline(ustar, 'b--', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('Thrust u (N)')
title('Control Input (Thrust Force)')
legend('u(t)', sprintf('u* = %.0f N (feedforward)', ustar)); grid on

