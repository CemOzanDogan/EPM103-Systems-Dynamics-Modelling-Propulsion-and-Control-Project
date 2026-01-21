%% =========================================================================
%% QUANSER HELICOPTER CONTROL PROJECT - GROUP 3
%% Integrated Design and Simulation of All Parts
%% City, University of London - Dr. Abdelhafid Zenati
%% =========================================================================
%%   Part 1: Classical Roll and Pitch Controller Design
%%   Part 2: Travel Controller and DOFB Design
%%   Part 3: Flight Mission Simulation (Bonus)
%% =========================================================================

clear all; close all; clc;

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║     QUANSER HELICOPTER CONTROL PROJECT                     ║\n');
fprintf('║     Integrated Design and Simulation of All Parts          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% =========================================================================
%% SECTION 1: PHYSICAL PARAMETERS AND SYSTEM DEFINITION
%% =========================================================================

fprintf('========================================\n');
fprintf('INITIALIZING PHYSICAL PARAMETERS\n');
fprintf('========================================\n');

% Physical parameters from Table 1 (Page 3)
params = struct();
params.m = 1.15;
params.M = 3.57;
params.l_boom = 0.66;
params.l_phi = 0.004;
params.l_theta = 0.014;
params.l_h = 0.177;
params.I_xx = 0.036;
params.I_yy = 0.93;
params.K_tau = 4.25e-3;
params.theta_rest = -25;
params.g = 9.81;

% Uncertain damping coefficients (ranges given in Table 1)
params.L_p_range = [0.02, 0.2];   % roll damping coefficient range
params.M_q_range = [0.1, 0.9];    % pitch damping coefficient range
params.L_p_nominal = 0.11;        % nominal roll damping (mid-range)
params.M_q_nominal = 0.5;         % nominal pitch damping (mid-range)

% Additional parameters for Part 2 (Page 7)
params.I_zz = 0.93;
params.K_D = 0;

% Derived parameters (trim thrust and travel coefficient)
params.omega_coll_bar = (params.M * params.g * params.l_theta * ...
                        sind(params.theta_rest)) / ...
                        (params.K_tau * params.l_boom);  % trim collective thrust
params.K_v = 0.0125 * params.omega_coll_bar * ...
             params.K_tau * params.l_boom;               % travel rate thrust coefficient

fprintf('Physical Parameters Loaded:\n');
fprintf('  Roll Damping (L_p):   %.3f (range: [%.2f, %.2f])\n', ...
        params.L_p_nominal, params.L_p_range(1), params.L_p_range(2));
fprintf('  Pitch Damping (M_q):  %.3f (range: [%.1f, %.1f])\n', ...
        params.M_q_nominal, params.M_q_range(1), params.M_q_range(2));
fprintf('  Trim Collective:      %.3f N\n', params.omega_coll_bar);
fprintf('  Travel Coefficient:   %.6f\n', params.K_v);
fprintf('\n');

%% =========================================================================
%% PART 1: CLASSICAL CONTROLLER DESIGN
%% =========================================================================

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                      PART 1                                ║\n');
fprintf('║         CLASSICAL ROLL AND PITCH CONTROLLERS               ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Section 1.1: Transfer Function Derivation
fprintf('--- Step 1: Deriving Transfer Functions ---\n');

% Motor dynamics (linear model from equations 7 and 8, Page 3)
M_cyc  = tf(780, [1 6]);   % cyclic motor
M_coll = tf(540, [1 6]);   % collective motor
fprintf('  M_cyc(s)  = 780/(s+6)\n');
fprintf('  M_coll(s) = 540/(s+6)\n');

% Roll transfer function Q_R(s): omega_cyc -> phi
num_QR = params.K_tau * params.l_h;
den_QR = [params.I_xx, params.L_p_nominal, params.m * params.g * params.l_phi];
Q_R = tf(num_QR, den_QR);
fprintf('  Q_R(s) = %.6f / (%.3f*s^2 + %.3f*s + %.6f)\n', ...
        num_QR, den_QR(1), den_QR(2), den_QR(3));

% Pitch transfer function Q_P(s): omega_coll -> theta
K_theta = params.M * params.g * params.l_theta * cosd(params.theta_rest);
num_QP = params.K_tau * params.l_boom;
den_QP = [params.I_yy, params.M_q_nominal, K_theta];
Q_P = tf(num_QP, den_QP);
fprintf('  Q_P(s) = %.6f / (%.3f*s^2 + %.3f*s + %.6f)\n', ...
        num_QP, den_QP(1), den_QP(2), den_QP(3));

% Full plants from voltage to angle
G_roll  = M_cyc  * Q_R;  % V_cyc  -> phi
G_pitch = M_coll * Q_P;  % V_coll -> theta

fprintf('  G_roll(s)  = M_cyc(s)  * Q_R(s)\n');
fprintf('  G_pitch(s) = M_coll(s) * Q_P(s)\n\n');

%% Section 1.2: Roll Controller Design (PD)
fprintf('--- Step 2: Designing Roll Controller ---\n');

% PD controller for roll: We choose gains for fast response and <10%% error @ 1 rad/s
Kp_roll = 6.0;
Kd_roll = 1.0;
Gc_roll = tf([Kd_roll, Kp_roll], 1);

fprintf('  Controller: PD (Proportional-Derivative)\n');
fprintf('  Gc_roll(s) = %.2f + %.2f*s\n', Kp_roll, Kd_roll);

% Closed-loop behaviour
L_roll = Gc_roll * G_roll;     % open loop
T_roll = feedback(L_roll, 1);  % closed loop (unit feedback)

[Gm_roll, Pm_roll, ~, ~] = margin(L_roll);
info_roll = stepinfo(T_roll);

fprintf('  Phase Margin: %.2f° (Requirement: ≥20°)\n', Pm_roll);
fprintf('  Gain Margin: %.2f dB\n', 20*log10(Gm_roll));
fprintf('  Rise Time: %.3f s\n', info_roll.RiseTime);
fprintf('  Overshoot: %.2f%%\n', info_roll.Overshoot);

% Frequency-domain tracking at 1 rad/s
T_roll_1rad = evalfr(T_roll, 1i*1);
tracking_error_roll = abs(1 - abs(T_roll_1rad)) * 100;
fprintf('  Tracking Error @ 1 rad/s: %.2f%% (Requirement: <10%%)\n', tracking_error_roll);

if Pm_roll >= 20 && tracking_error_roll < 10
    fprintf('  Roll controller meets all specifications.\n\n');
else
    fprintf('  Roll controller does not meet all specifications; consider retuning.\n\n');
end

%% Section 1.3: Pitch Controller Design (PID)
fprintf('--- Step 3: Designing Pitch Controller ---\n');

% PID controller: integral for small steady-state error, derivative to improve phase
Kp_pitch = 4.5;
Ki_pitch = 0.4;
Kd_pitch = 6.0;
Gc_pitch = tf([Kd_pitch, Kp_pitch, Ki_pitch], [1, 0]);

fprintf('  Controller: PID (Proportional-Integral-Derivative)\n');
fprintf('  Gc_pitch(s) = %.2f + %.2f/s + %.2f*s\n', Kp_pitch, Ki_pitch, Kd_pitch);

% Closed-loop behaviour
L_pitch = Gc_pitch * G_pitch;
T_pitch = feedback(L_pitch, 1);

[Gm_pitch, Pm_pitch, ~, ~] = margin(L_pitch);
info_pitch = stepinfo(T_pitch);

fprintf('  Phase Margin: %.2f° (Requirement: ≥20°)\n', Pm_pitch);
fprintf('  Gain Margin: %.2f dB\n', 20*log10(Gm_pitch));
fprintf('  Rise Time: %.3f s\n', info_pitch.RiseTime);
fprintf('  Overshoot: %.2f%%\n', info_pitch.Overshoot);

% Steady-state and frequency-domain tracking performance. 
% We do it from 0 to 60 seconds
t_ss = 0:0.01:60;
[y_ss, ~] = step(T_pitch, t_ss);
ss_error_pitch = abs(1 - y_ss(end)) * 100;
T_pitch_1rad = evalfr(T_pitch, 1i*1);
tracking_error_pitch = abs(1 - abs(T_pitch_1rad)) * 100;

fprintf('  Steady-State Error: %.4f%% (Requirement: ~0%%)\n', ss_error_pitch);
fprintf('  Tracking Error @ 1 rad/s: %.2f%% (Requirement: <10%%)\n', tracking_error_pitch);

if Pm_pitch >= 20 && tracking_error_pitch < 10 && ss_error_pitch < 1
    fprintf('  Pitch controller meets all specifications.\n\n');
else
    fprintf('  Pitch controller does not meet all specifications; consider retuning.\n\n');
end

%% Section 1.4: Generating Plots for Part 1
fprintf('--- Generating Part 1 Plots ---\n');

% Combined figure for Part 1 (Bode, margins, step responses, robustness)
fig1 = figure('Name', 'Part 1: Classical Controller Design', ...
              'Position', [50 50 1400 900]);

% Roll analysis (top row plots)
subplot(2,4,1);
bode(G_roll);
grid on;
title('Roll Plant Bode');

subplot(2,4,2);
margin(L_roll);
grid on;
title('Roll Open-Loop');

subplot(2,4,3);
step(T_roll * 20 * pi/180, 10);
grid on;
title('Roll Step Response (20°)');
xlabel('Time (s)');
ylabel('Roll (rad)');

subplot(2,4,4);
% Robustness across L_p range
L_p_test = linspace(params.L_p_range(1), params.L_p_range(2), 5);
hold on;
for i = 1:length(L_p_test)
    den_QR_test = [params.I_xx, L_p_test(i), params.m * params.g * params.l_phi];
    Q_R_test = tf(num_QR, den_QR_test);
    G_roll_test = M_cyc * Q_R_test;
    T_roll_test = feedback(Gc_roll * G_roll_test, 1);
    step(T_roll_test * 20 * pi/180, 10);
end
grid on;
title('Roll Robustness (L_p range)');
xlabel('Time (s)');
ylabel('Roll (rad)');

% Pitch analysis (bottom row plots)
subplot(2,4,5);
bode(G_pitch);
grid on;
title('Pitch Plant Bode');

subplot(2,4,6);
margin(L_pitch);
grid on;
title('Pitch Open-Loop');

subplot(2,4,7);
step(T_pitch * 20 * pi/180, 20);
grid on;
title('Pitch Step Response (20°)');
xlabel('Time (s)');
ylabel('Pitch (rad)');

subplot(2,4,8);
% Robustness across M_q range
M_q_test = linspace(params.M_q_range(1), params.M_q_range(2), 5);
hold on;
for i = 1:length(M_q_test)
    den_QP_test = [params.I_yy, M_q_test(i), K_theta];
    Q_P_test = tf(num_QP, den_QP_test);
    G_pitch_test = M_coll * Q_P_test;
    T_pitch_test = feedback(Gc_pitch * G_pitch_test, 1);
    step(T_pitch_test * 20 * pi/180, 20);
end
grid on;
title('Pitch Robustness (M_q range)');
xlabel('Time (s)');
ylabel('Pitch (rad)');

fprintf('  Part 1 plots generated.\n\n');

fprintf('========================================\n');
fprintf('PART 1 COMPLETE\n');
fprintf('========================================\n');
fprintf('  Roll controller designed and validated.\n');
fprintf('  Pitch controller designed and validated.\n\n');

%% =========================================================================
%% PART 2: TRAVEL CONTROLLER AND DOFB DESIGN
%% =========================================================================

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                      PART 2                                ║\n');
fprintf('║         TRAVEL CONTROLLER AND DOFB DESIGN                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Section 2.1: Augmented Plant Dynamics for Travel
fprintf('--- Step 1: Deriving Augmented Travel Plant ---\n');

% Closed-loop roll system (inner loop)
T_roll_inner = feedback(Gc_roll * G_roll, 1);
fprintf('  Inner loop: T_roll(s) = phi(s)/phi_c(s)\n');

% Travel dynamics from roll: psi(s)/phi(s) = K_v / (I_zz * s^2)
G_travel_from_roll = tf(params.K_v, [params.I_zz, 0, 0]);
fprintf('  Travel from roll: psi(s)/phi(s) = K_v / (I_zz * s^2)\n');
fprintf('  K_v = %.6f (negative, sign handled in controller design)\n', params.K_v);

% Augmented travel plant from travel command to travel angle
G_travel_augmented = T_roll_inner * G_travel_from_roll;
fprintf('  Augmented plant: G_travel(s) = T_roll(s) * [psi/phi](s)\n\n');

%% Section 2.2: Classical Travel Controller Design
fprintf('--- Step 2: Designing Classical Travel Controller ---\n');

% PD controller for travel (double integrator + inner loop).
% Because K_v < 0, stabilizing control requires negative gains here.
Kp_travel = -5.0;
Kd_travel = -8.0;
Gc_travel = tf([Kd_travel, Kp_travel], 1);

fprintf('  Controller: PD (Proportional-Derivative)\n');
fprintf('  Gc_travel(s) = %.2f + %.2f*s\n', Kp_travel, Kd_travel);

L_travel = Gc_travel * G_travel_augmented;
T_travel = feedback(L_travel, 1);
[Gm_travel, Pm_travel, ~, ~] = margin(L_travel);

fprintf('  Phase Margin: %.2f°\n', Pm_travel);
fprintf('  Gain Margin: %.2f dB\n', 20*log10(Gm_travel));

info_travel = stepinfo(T_travel);
fprintf('  Rise Time: %.3f s\n', info_travel.RiseTime);
fprintf('  Settling Time: %.3f s\n', info_travel.SettlingTime);
fprintf('  Overshoot: %.2f%%\n\n', info_travel.Overshoot);

%% Section 2.3: DOFB Pitch Controller Design
fprintf('--- Step 3: DOFB Pitch Controller Design ---\n');

% Convert pitch plant to state-space (including motor dynamics)
[A_p, B_p, C_p, D_p] = ssdata(G_pitch);
n_states = size(A_p, 1);

fprintf('  State-space representation: %d states\n', n_states);

% Controllability and observability checks
Co = ctrb(A_p, B_p);
Ob = obsv(A_p, C_p);
fprintf('  Controllable: %s\n', ternary(rank(Co) == n_states, 'Yes', 'No'));
fprintf('  Observable:   %s\n', ternary(rank(Ob) == n_states, 'Yes', 'No'));

% LQR design for state-feedback K
Q_lqr = diag([100, 50, 1]);   % state weights (angle, rate, motor)
R_lqr = 1;                    % control effort weight
K_lqr = lqr(A_p, B_p, Q_lqr, R_lqr);
fprintf('  LQR Gain K = [%.4f, %.4f, %.4f]\n', K_lqr(1), K_lqr(2), K_lqr(3));

% LQE (Kalman filter) observer design for L
Q_n = 0.001 * eye(n_states);  % process noise covariance
R_n = 0.01;                   % measurement noise covariance
G_n = eye(n_states);          % noise input matrix
[L_lqe, ~, ~] = lqe(A_p, G_n, C_p, Q_n, R_n);
fprintf('  LQE Gain L = [%.4f; %.4f; %.4f]\n', L_lqe(1), L_lqe(2), L_lqe(3));

% DOFB controller realisation: state-feedback + observer + reference scaling
N_bar = -inv(C_p * inv(A_p - B_p*K_lqr) * B_p);  % reference scaling for unity gain
A_ctrl = A_p - B_p*K_lqr - L_lqe*C_p;
B_ctrl = L_lqe;
C_ctrl = -K_lqr;
D_ctrl = 0;

[num_dofb, den_dofb] = ss2tf(A_ctrl, B_ctrl, C_ctrl, D_ctrl);
Gc_dofb = tf(num_dofb, den_dofb) * N_bar;
fprintf('  Reference Scaling: N_bar = %.4f\n', N_bar);

% Closed-loop DOFB performance
L_dofb = Gc_dofb * G_pitch;
T_dofb = feedback(L_dofb, 1);
[Gm_dofb, Pm_dofb, ~, ~] = margin(L_dofb);

fprintf('  DOFB Phase Margin: %.2f°\n', Pm_dofb);
fprintf('  DOFB Gain Margin: %.2f dB\n', 20*log10(Gm_dofb));

info_dofb = stepinfo(T_dofb);
fprintf('  DOFB Rise Time: %.3f s\n', info_dofb.RiseTime);
fprintf('  DOFB Settling Time: %.3f s\n', info_dofb.SettlingTime);
fprintf('  DOFB Overshoot: %.2f%%\n\n', info_dofb.Overshoot);

%% Section 2.4: Generate Plots for Part 2
fprintf('--- Generating Part 2 Plots ---\n');

% Combined figure for Part 2
fig2 = figure('Name', 'Part 2: Advanced Control Design', ...
              'Position', [50 50 1400 900]);

% Travel controller plots
subplot(2,4,1);
bode(Gc_travel);
grid on;
title('Travel Controller Bode');

subplot(2,4,2);
margin(L_travel);
grid on;
title('Travel Open-Loop');

subplot(2,4,5);
step(T_travel * 20 * pi/180, 40);
grid on;
title('Travel Step Response (20°)');
xlabel('Time (s)');
ylabel('Travel (rad)');

subplot(2,4,6);
t_comp = 0:0.01:40;
u_comp = (20*pi/180) * (t_comp >= 0);
y_comp = lsim(T_travel, u_comp, t_comp);
plot(t_comp, u_comp*180/pi, 'b--', 'LineWidth', 1.5);
hold on;
plot(t_comp, y_comp*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
title('Travel Complete Response');
xlabel('Time (s)');
ylabel('Angle (°)');
legend('Command', 'Response', 'Location', 'best');

% DOFB vs classical pitch controller plots 
subplot(2,4,3);
bode(Gc_dofb);
grid on;
title('DOFB Controller Bode');

subplot(2,4,4);
margin(L_dofb);
grid on;
title('DOFB Open-Loop');

subplot(2,4,7);
t_comp = 0:0.01:15;
[y_classical, ~] = step(T_pitch * 15 * pi/180, t_comp);
[y_dofb, ~] = step(T_dofb * 15 * pi/180, t_comp);
plot(t_comp, y_classical*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_comp, y_dofb*180/pi, 'r--', 'LineWidth', 1.5);
plot(t_comp, 15*ones(size(t_comp)), 'k:', 'LineWidth', 1);
grid on;
title('DOFB vs Classical (15° step)');
xlabel('Time (s)');
ylabel('Pitch (°)');
legend('Classical PID', 'DOFB', 'Command', 'Location', 'best');

subplot(2,4,8);
% Simple normalised performance comparison for checking the system
categories = {'Rise Time', 'Settling', 'Overshoot', 'Phase Margin'};
classical_metrics = [info_pitch.RiseTime, info_pitch.SettlingTime/10, ...
                     info_pitch.Overshoot/10, Pm_pitch/10];
dofb_metrics = [info_dofb.RiseTime, info_dofb.SettlingTime/10, ...
                info_dofb.Overshoot/10, Pm_dofb/10];
x = 1:4;
bar(x, [classical_metrics; dofb_metrics]');
set(gca, 'XTickLabel', categories, 'XTick', 1:4);
xtickangle(45);
legend('Classical', 'DOFB', 'Location', 'best');
title('Performance Comparison (normalised)');
grid on;

fprintf('  Part 2 plots generated.\n\n');

fprintf('========================================\n');
fprintf('PART 2 COMPLETE\n');
fprintf('========================================\n');
fprintf('  Travel controller designed.\n');
fprintf('  DOFB pitch controller designed.\n\n');

%% =========================================================================
%% PART 3: FLIGHT MISSION (BONUS) - SIMULATION
%% =========================================================================

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                      PART 3                                ║\n');
fprintf('║       FLIGHT MISSION - TRAJECTORY SIMULATION               ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Section 3.1: Define Mission Waypoints
fprintf('--- Step 1: Defining Mission Waypoints ---\n');

% Mission logic:
%  - Two solid buildings (psi 60–90, 270–300, theta -25 to 15): fly over at ~20°
%  - Tunnel (psi 150–210) with gap theta -10 to 5: fly through at about -2.5°
%  - Surveillance point and payload hover around psi = 330°

% Columns: [psi (deg), theta (deg), hover_time (s)]
waypoints = [
    % --- LEG 1: START TO SURVEILLANCE ---
    0,    -25,   0;
    20,   -25,   0;
    40,    0,    0;
    50,    20,   0;

    % Over Building 1
    60,    20,   0;
    90,    20,   0;

    % Approach tunnel
    110,   20,   0;
    130,   5,    0;
    140,   -2.5, 0;

    % Through tunnel
    150,   -2.5, 0;
    180,   -2.5, 0;
    210,   -2.5, 0;

    % Approach Building 2
    230,   -2.5, 0;
    250,   10,   0;
    260,   20,   0;

    % Over Building 2
    270,   20,   0;
    300,   20,   0;

    % Surveillance point
    330,   20,   0;

    % --- LEG 2: PAYLOAD PICKUP ---
    330,   0,    0;
    330,   -20,  15;   % hover ~15 s at payload

    % --- LEG 3: RETURN HOME (reverse path) ---
    330,   20,   0;
    300,   20,   0;
    270,   20,   0;

    260,   20,   0;
    230,   -2.5, 0;
    210,   -2.5, 0;
    150,   -2.5, 0;

    140,   -2.5, 0;
    110,   20,   0;
    90,    20,   0;
    60,    20,   0;

    % Landing
    50,    20,   0;
    30,    0,    0;
    10,    -25,  0;
    0,     -25,  0;
];

fprintf('  Total waypoints: %d\n', size(waypoints, 1));
fprintf('  Hover at payload: waypoint 21 (12 s commanded)\n\n');

%% Section 3.2: Generate Smooth Trajectory
fprintf('--- Step 2: Generating Trajectory ---\n');

% Maximum command rates (deg/s) for trajectory generation
v_travel_max = 1.0;   % travel rate limit
v_pitch_max  = 3.0;   % pitch rate limit

dt_traj = 0.05;
psi_traj = [];
theta_traj = [];
t_traj = [];
t_current = 0;

for i = 1:size(waypoints,1)-1
    psi_start   = waypoints(i, 1);
    theta_start = waypoints(i, 2);
    psi_end     = waypoints(i+1, 1);
    theta_end   = waypoints(i+1, 2);
    hover_time  = waypoints(i+1, 3);

    d_psi   = psi_end   - psi_start;
    d_theta = theta_end - theta_start;

    t_psi   = abs(d_psi)   / v_travel_max;
    t_theta = abs(d_theta) / v_pitch_max;
    t_segment = max(t_psi, t_theta) * 1.2;  % 20% safety margin on duration

    if t_segment > 0
        n_points = max(3, ceil(t_segment / dt_traj) + 1);
        t_seg = linspace(0, t_segment, n_points);

        psi_seg   = linspace(psi_start,   psi_end,   n_points);
        theta_seg = linspace(theta_start, theta_end, n_points);

        if i == 1
            psi_traj   = [psi_traj,   psi_seg];
            theta_traj = [theta_traj, theta_seg];
            t_traj     = [t_traj,     t_current + t_seg];
        else
            psi_traj   = [psi_traj,   psi_seg(2:end)];
            theta_traj = [theta_traj, theta_seg(2:end)];
            t_traj     = [t_traj,     t_current + t_seg(2:end)];
        end
        t_current = t_current + t_segment;
    end

    % Add hover segment (constant psi/theta)
    if hover_time > 0
        n_hover = max(2, ceil(hover_time / dt_traj));
        psi_traj   = [psi_traj,   psi_end   * ones(1, n_hover)];
        theta_traj = [theta_traj, theta_end * ones(1, n_hover)];
        t_traj     = [t_traj,     linspace(t_current + dt_traj, ...
                                           t_current + hover_time, n_hover)];
        t_current = t_current + hover_time;
    end
end

fprintf('  Trajectory duration: %.1f seconds\n', t_traj(end));
fprintf('  Number of trajectory points: %d\n', length(t_traj));

% Convert trajectory and waypoints from degrees to radians (for dynamics)
fprintf('  Converting trajectory from degrees to radians...\n');
psi_traj_rad   = psi_traj   * pi/180;
theta_traj_rad = theta_traj * pi/180;
waypoints_rad  = waypoints  * pi/180;
fprintf('  Conversion complete.\n\n');

%% Section 3.3: Closed-Loop Simulation
fprintf('--- Step 3: Running Closed-Loop Simulation ---\n');

dt_sim = 0.01;
t_sim = 0:dt_sim:t_traj(end);
n_sim = length(t_sim);

% Interpolate commands (in radians)
psi_cmd   = interp1(t_traj, psi_traj_rad,   t_sim, 'linear', 'extrap');
theta_cmd = interp1(t_traj, theta_traj_rad, t_sim, 'linear', 'extrap');

% State histories (all in radians / rad/s)
psi_actual   = zeros(1, n_sim);
theta_actual = zeros(1, n_sim);
phi_actual   = zeros(1, n_sim);
dpsi   = zeros(1, n_sim);
dtheta = zeros(1, n_sim);
dphi   = zeros(1, n_sim);

V_cyc_history  = zeros(1, n_sim);
V_coll_history = zeros(1, n_sim);
phi_cmd_history = zeros(1, n_sim);

% Initial conditions
psi_actual(1)   = waypoints_rad(1, 1);
theta_actual(1) = waypoints_rad(1, 2);
theta_int_error = 0;

% Simple Euler integration of the coupled roll–pitch–travel model
for i = 2:n_sim
    psi_c   = psi_cmd(i);
    theta_c = theta_cmd(i);

    % Travel outer loop: psi_c -> roll command phi_c
    psi_error = psi_c - psi_actual(i-1);
    phi_c = Kp_travel * psi_error - Kd_travel * dpsi(i-1);
    phi_c = max(-35*pi/180, min(35*pi/180, phi_c));  % roll command limit
    phi_cmd_history(i) = phi_c;

    % Roll inner loop: phi_c -> V_cyc
    phi_error = phi_c - phi_actual(i-1);
    V_cyc = Kp_roll * phi_error - Kd_roll * dphi(i-1);
    V_cyc = max(-5, min(5, V_cyc));                  % voltage saturation
    V_cyc_history(i) = V_cyc;

    % Pitch loop: theta_c -> V_coll (with integral action)
    theta_error = theta_c - theta_actual(i-1);
    theta_int_error = theta_int_error + theta_error * dt_sim;
    theta_int_error = max(-10, min(10, theta_int_error));  % anti-windup

    V_coll = Kp_pitch * theta_error + ...
             Ki_pitch * theta_int_error - ...
             Kd_pitch * dtheta(i-1);
    V_coll = max(-5, min(5, V_coll));
    V_coll_history(i) = V_coll;

    % Approximate motor dynamics (first-order) to rotor speeds
    omega_cyc  = 780 * V_cyc_history(i)  / (1 + 6*dt_sim);
    omega_coll = 540 * V_coll_history(i) / (1 + 6*dt_sim);

    % Second-order dynamics for phi, theta, and psi
    ddphi = (params.K_tau * params.l_h * omega_cyc - ...
             params.m * params.g * params.l_phi * phi_actual(i-1) - ...
             params.L_p_nominal * dphi(i-1)) / params.I_xx;

    ddtheta = (params.K_tau * params.l_boom * omega_coll - ...
               K_theta * theta_actual(i-1) - ...
               params.M_q_nominal * dtheta(i-1)) / params.I_yy;

    ddpsi = params.K_v * phi_actual(i-1) / params.I_zz;

    % Integrate to obtain rates and angles
    dphi(i)   = dphi(i-1)   + ddphi   * dt_sim;
    dtheta(i) = dtheta(i-1) + ddtheta * dt_sim;
    dpsi(i)   = dpsi(i-1)   + ddpsi   * dt_sim;

    phi_actual(i)   = phi_actual(i-1)   + dphi(i)   * dt_sim;
    theta_actual(i) = theta_actual(i-1) + dtheta(i) * dt_sim;
    psi_actual(i)   = psi_actual(i-1)   + dpsi(i)   * dt_sim;
end

fprintf('  Simulation complete.\n\n');

%% Section 3.4: Obstacle Detection and Scoring
fprintf('--- Step 4: Obstacle Detection and Scoring ---\n');

% Convert actual trajectory and commands back to degrees for scoring
psi_actual_deg   = psi_actual   * 180/pi;
theta_actual_deg = theta_actual * 180/pi;
psi_cmd_deg      = psi_cmd      * 180/pi;
theta_cmd_deg    = theta_cmd    * 180/pi;
phi_actual_deg   = phi_actual   * 180/pi;
phi_cmd_history_deg = phi_cmd_history * 180/pi;

% Obstacle definitions (deg)
bldg1 = struct('psi', [60, 90],   'theta', [-25, 15]);
tunnel_walls = struct('psi', [150, 210], 'theta_low', -10, 'theta_high', 5);
bldg2 = struct('psi', [270, 300], 'theta', [-25, 15]);

% Building 1: hit if inside psi window and inside theta rectangle
hit_building_one_mask = (psi_actual_deg >= bldg1.psi(1) & psi_actual_deg <= bldg1.psi(2)) & ...
                        (theta_actual_deg >= bldg1.theta(1) & theta_actual_deg <= bldg1.theta(2));

% Tunnel: hit if within psi range but outside the vertical gap
hit_tunnel_mask = (psi_actual_deg >= tunnel_walls.psi(1) & psi_actual_deg <= tunnel_walls.psi(2)) & ...
                  (theta_actual_deg > tunnel_walls.theta_high | ...
                   theta_actual_deg < tunnel_walls.theta_low);

% Building 2: similar to Building 1
hit_building_two_mask = (psi_actual_deg >= bldg2.psi(1) & psi_actual_deg <= bldg2.psi(2)) & ...
                        (theta_actual_deg >= bldg2.theta(1) & theta_actual_deg <= bldg2.theta(2));

% Ceiling / ground violations
hit_ceiling = theta_actual_deg > 25;
% Slight tolerance on ground to avoid numerical artefacts at start
hit_ground  = theta_actual_deg < -26;

% Flags for any collisions
hit_building_one = any(hit_building_one_mask);
hit_tunnel       = any(hit_tunnel_mask);
hit_building_two = any(hit_building_two_mask);

% Time-weighted penalties (arbitrary scaling for scoring)
penalty_building_one = sum(hit_building_one_mask) * dt_sim * 1000;
penalty_tunnel       = sum(hit_tunnel_mask)       * dt_sim * 1000;
penalty_building_two = sum(hit_building_two_mask) * dt_sim * 1000;
penalty_ceiling      = sum(hit_ceiling)           * dt_sim * 1000;
penalty_ground       = sum(hit_ground)            * dt_sim * 1000;

total_obstacle_penalty = penalty_building_one + penalty_tunnel + ...
                         penalty_building_two + penalty_ceiling + ...
                         penalty_ground;

% Hover performance around payload location (psi ≈ 330, theta ≈ -20)
hover_mask = abs(psi_actual_deg - 330) < 3 & ...
             abs(theta_actual_deg - (-20)) < 3;
d_hover = diff([0, hover_mask, 0]);
hover_starts = find(d_hover == 1);
hover_ends   = find(d_hover == -1);

if isempty(hover_starts)
    hover_time_actual = 0;
    fprintf('  No hover segment detected at the payload.\n');
else
    hover_durations   = (hover_ends - hover_starts) * dt_sim;
    hover_time_actual = max(hover_durations);
    fprintf('  Hover detected at payload: %.1f seconds\n', hover_time_actual);
end

% Penalise missing hover time (target = 10 s)
hover_quality = min(hover_time_actual / 10, 1);
penalty_hover = max(0, (10 - hover_time_actual)) * 10;

% Integral of absolute tracking error in degrees
psi_error_integral   = sum(abs(psi_actual_deg   - psi_cmd_deg))   * dt_sim;
theta_error_integral = sum(abs(theta_actual_deg - theta_cmd_deg)) * dt_sim;

% Total score (lower is better)
total_score = total_obstacle_penalty + penalty_hover + ...
              psi_error_integral + theta_error_integral;

fprintf('\nPenalty Breakdown:\n');
fprintf('  Building One:     %8.0f\n', penalty_building_one);
fprintf('  Tunnel:           %8.0f\n', penalty_tunnel);
fprintf('  Building Two:     %8.0f\n', penalty_building_two);
fprintf('  Ceiling:          %8.0f\n', penalty_ceiling);
fprintf('  Ground:           %8.0f\n', penalty_ground);
fprintf('  -------------------------------\n');
fprintf('  Obstacle Total:   %8.0f\n\n', total_obstacle_penalty);

fprintf('  Hover penalty:    %8.2f (%.1f/10 s)\n', penalty_hover, hover_time_actual);
fprintf('  Travel error:     %8.2f\n', psi_error_integral);
fprintf('  Pitch error:      %8.2f\n', theta_error_integral);
fprintf('  -------------------------------\n');
fprintf('  TOTAL SCORE:      %8.2f\n\n', total_score);

if total_obstacle_penalty == 0
    fprintf('  No obstacles were hit.\n\n');
elseif total_obstacle_penalty < 2000
    fprintf('  Some obstacle contacts occurred (low penalty).\n\n');
else
    fprintf('  Multiple obstacle collisions (high penalty).\n\n');
end

%% Section 3.5: Mission Visualization
fprintf('--- Step 5: Generating Mission Visualization ---\n');

fig3 = figure('Name', 'Part 3: Flight Mission', ...
              'Position', [50 50 1600 900]);

% Mission map (travel vs pitch)
subplot(2,3,[1 4]);
hold on;

% Building 1
fill([60 90 90 60], [-25 -25 15 15], 'r', 'FaceAlpha', 0.3, ...
     'EdgeColor', 'r', 'LineWidth', 2);
text(75, -5, 'Bldg 1', 'Color', 'r', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Tunnel (upper and lower blocks)
fill([150 210 210 150], [5 5 25 25], 'r', 'FaceAlpha', 0.3, ...
     'EdgeColor', 'r', 'LineWidth', 2);     % upper
fill([150 210 210 150], [-25 -25 -10 -10], 'r', 'FaceAlpha', 0.3, ...
     'EdgeColor', 'r', 'LineWidth', 2);     % lower
text(180, 0, 'Tunnel Gap', 'Color', 'b', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Building 2
fill([270 300 300 270], [-25 -25 15 15], 'r', 'FaceAlpha', 0.3, ...
     'EdgeColor', 'r', 'LineWidth', 2);
text(285, -5, 'Bldg 2', 'Color', 'r', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Ceiling and ground lines
plot([0, 360], [25, 25], 'k--', 'LineWidth', 2);
text(180, 26.5, 'CEILING', 'FontSize', 10, ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
plot([0, 360], [-25, -25], 'k--', 'LineWidth', 2);
text(180, -26.5, 'GROUND', 'FontSize', 10, ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Trajectories and waypoints
plot(psi_actual_deg, theta_actual_deg, 'g-', 'LineWidth', 2.5);
plot(psi_cmd_deg,    theta_cmd_deg,    'b--', 'LineWidth', 1, ...
     'Color', [0.5 0.5 1]);
plot(waypoints(:,1), waypoints(:,2), 'mo', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'm');

% Key points
plot(0, -25, 'bs', 'MarkerSize', 18, 'MarkerFaceColor', 'b', 'LineWidth', 2);
text(5, -25, 'START/FINISH', 'FontSize', 10, 'FontWeight', 'bold');

plot(330, 20, 'r^', 'MarkerSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 2);
text(325, 22, 'Surveillance', 'FontSize', 10, 'FontWeight', 'bold');

plot(330, -20, 'gd', 'MarkerSize', 18, 'MarkerFaceColor', [0 0.7 0], ...
     'LineWidth', 2);
text(335, -20, 'Payload', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Travel Angle \psi (deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pitch Angle \theta (deg)', 'FontSize', 12, 'FontWeight', 'bold');
title('Mission Map - Complete Trajectory', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis([0 360 -30 30]);
legend('Bldg1', 'Tunnel Upper', 'Tunnel Lower', 'Bldg2', ...
       'Ceiling/Ground', '', 'Actual Path', 'Planned Path', ...
       'Waypoints', 'Location', 'northwest', 'FontSize', 8);

% Score annotation box
annotation('textbox', [0.02 0.45 0.3 0.1], 'String', ...
          sprintf('TOTAL SCORE: %.1f\nObstacle Hits: %d\nHover: %.1f/10 s', ...
          total_score, sum([hit_building_one hit_tunnel hit_building_two]), ...
          hover_time_actual), ...
          'FontSize', 12, 'FontWeight', 'bold', ...
          'BackgroundColor', 'yellow', ...
          'EdgeColor', 'black', 'LineWidth', 2);

% Travel vs time
subplot(2,3,2);
plot(t_sim, psi_actual_deg, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_sim, psi_cmd_deg, 'r--', 'LineWidth', 1);
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Travel \psi (deg)', 'FontWeight', 'bold');
title('Travel Angle vs Time');
legend('Actual', 'Command');
grid on;

% Pitch vs time
subplot(2,3,3);
plot(t_sim, theta_actual_deg, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_sim, theta_cmd_deg, 'r--', 'LineWidth', 1);
plot([0, t_sim(end)], [25, 25], 'k--', 'LineWidth', 1.5);
plot([0, t_sim(end)], [-25, -25], 'k--', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Pitch \theta (deg)', 'FontWeight', 'bold');
title('Pitch Angle vs Time');
legend('Actual', 'Command', 'Ceiling', 'Ground');
grid on;

% Roll vs time
subplot(2,3,5);
plot(t_sim, phi_actual_deg, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_sim, phi_cmd_history_deg, 'r--', 'LineWidth', 1);
plot([0, t_sim(end)], [35, 35], 'k--');
plot([0, t_sim(end)], [-35, -35], 'k--');
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Roll \phi (deg)', 'FontWeight', 'bold');
title('Roll Angle vs Time');
legend('Actual', 'Command', 'Limits');
grid on;

% Control voltages vs time
subplot(2,3,6);
plot(t_sim, V_cyc_history, 'b-', 'LineWidth', 1);
hold on;
plot(t_sim, V_coll_history, 'r-', 'LineWidth', 1);
plot([0, t_sim(end)], [5, 5], 'k--');
plot([0, t_sim(end)], [-5, -5], 'k--');
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontWeight', 'bold');
title('Control Voltages');
legend('V_{cyc} (Roll)', 'V_{coll} (Pitch)', 'Limits');
grid on;
ylim([-6 6]);

fprintf('  Visualization complete.\n\n');

fprintf('========================================\n');
fprintf('PART 3 COMPLETE\n');
fprintf('========================================\n');

%% =========================================================================
%% FINAL SUMMARY
%% =========================================================================

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              PROJECT COMPLETION SUMMARY                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('PART 1 - Classical Controllers:\n');
fprintf('  Roll (PD):   Kp = %.2f, Kd = %.2f\n', Kp_roll, Kd_roll);
fprintf('    Phase Margin: %.2f°\n', Pm_roll);
fprintf('    Tracking Error @ 1 rad/s: %.2f%%\n', tracking_error_roll);
fprintf('\n');
fprintf('  Pitch (PID): Kp = %.2f, Ki = %.2f, Kd = %.2f\n', ...
        Kp_pitch, Ki_pitch, Kd_pitch);
fprintf('    Phase Margin: %.2f°\n', Pm_pitch);
fprintf('    Steady-State Error (step): %.4f%%\n', ss_error_pitch);
fprintf('    Tracking Error @ 1 rad/s: %.2f%%\n\n', tracking_error_pitch);

fprintf('PART 2 - Advanced Controllers:\n');
fprintf('  Travel (PD): Kp = %.2f, Kd = %.2f\n', Kp_travel, Kd_travel);
fprintf('    Phase Margin: %.2f°\n', Pm_travel);
fprintf('    Note: sign accounts for K_v < 0.\n\n');

fprintf('  DOFB Pitch:  Q = diag([%.0f, %.0f, %.0f]), R = %.0f\n', ...
        Q_lqr(1,1), Q_lqr(2,2), Q_lqr(3,3), R_lqr);
fprintf('    Phase Margin: %.2f°\n', Pm_dofb);
fprintf('    K = [%.4f, %.4f, %.4f]\n\n', K_lqr(1), K_lqr(2), K_lqr(3));

fprintf('PART 3 - Mission Performance:\n');
fprintf('  Mission Duration: %.1f seconds\n', t_sim(end));
fprintf('  Final Score: %.2f points (lower is better)\n', total_score);
fprintf('  Obstacles Clear: %s\n', ternary(total_obstacle_penalty==0, 'YES', 'NO'));
fprintf('  Hover Performance: %.0f%% of target time\n', hover_quality*100);
fprintf('  Average Travel Error: %.2f deg\n', psi_error_integral/t_sim(end));
fprintf('  Average Pitch Error: %.2f deg\n\n', theta_error_integral/t_sim(end));

if total_score < 500
    fprintf('Overall performance rating: excellent.\n');
elseif total_score < 1500
    fprintf('Overall performance rating: good.\n');
elseif total_score < 3000
    fprintf('Overall performance rating: acceptable; further optimisation is possible.\n');
else
    fprintf('Overall performance rating: needs improvement; consider refining the trajectory or gains.\n');
end

fprintf('\n');
fprintf('Visualization Summary:\n');
fprintf('  Figure 1: Part 1 - Classical Controllers (8 plots)\n');
fprintf('  Figure 2: Part 2 - Advanced Control (8 plots)\n');
fprintf('  Figure 3: Part 3 - Mission Results (6 plots)\n\n');

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                  END OF CONTROL PROJECT                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

% Helper function for compact conditional printing
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

