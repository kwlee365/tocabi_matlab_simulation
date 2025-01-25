addpath('casadi-3.6.7-linux64-matlab2018b')
import casadi.*

PARA = PARA;

H = PARA.H;
state_length = PARA.state_length;
input_length = PARA.input_length;

m = SX.sym('m');
g = SX.sym('g');
h = SX.sym('h');
I = SX.sym('I', 3, 3);
I_inv = inv(I);

mu = SX.sym('mu');  % FRICTION COEFFICIENT
dT = SX.sym('dT');  % MPC TIMESTEP

H = PARA.H;
state_length = PARA.state_length;
input_length = PARA.input_length;

rL_ref_horizon = SX.sym('rL_ref_horizon', 3, H);
rR_ref_horizon = SX.sym('rR_ref_horizon', 3, H);

theta_ref_horizon = SX.sym('theta_ref_horizon', 3, H);

etaL_ref_horizon = SX.sym('etaL_ref_horizon', H);
etaR_ref_horizon = SX.sym('etaR_ref_horizon', H);

% STATE AND CONTROL INPUTS
X     = SX.sym('X',     state_length * H, 1);
X_ref = SX.sym('X_ref', state_length * H, 1);
U     = SX.sym('U',     input_length * H, 1);
U_ref = SX.sym('U_ref', input_length * H, 1);

v = [X; U];

W_Q = SX.sym('W_Q', H * state_length, 1);
W_R = SX.sym('W_R', H * input_length, 1);

J = (X - X_ref)' * diag(W_Q) * (X - X_ref) + (U - U_ref)' * diag(W_R) * (U - U_ref);

J_v = jacobian(J, v);
J_vv = hessian(J, v);

% SRBD DYNAMICS
ceq1 = [];

x0 = SX.sym('x0', state_length, 1);  % CURRENT ROBOT STATE
x_k = x0;
for i = 1:H
    
    u_k      = U((input_length * (i-1) + 1):(input_length * i));

    rL = rL_ref_horizon(:, i);
    rR = rR_ref_horizon(:, i);
    theta = theta_ref_horizon(:, i);

    etaL = etaL_ref_horizon(i);
    etaR = etaR_ref_horizon(i);

    T = SX.zeros(3, 3);
    roll = theta(1); pitch = theta(2); yaw = theta(3);

    T(1, 1) = cos(pitch) * cos(yaw);
    T(1, 2) = -sin(yaw);
    T(2, 1) = cos(pitch) * sin(yaw);
    T(2, 2) = cos(yaw);
    T(3, 1) = -sin(pitch);
    T(3, 3) = 1.0;

    T_inv = inv(T);

    % CONTINUOUS SYSTEM
    A = SX.zeros(state_length, state_length);
    B = SX.zeros(state_length, input_length);
    d = SX.zeros(state_length, 1);

    A(1:3, 7:9) = T_inv;
    % A(1:3, 7:9) = T;
    A(4:6, 10:12) = SX.eye(3);

    B(7:9, 1:3)   = etaL * I_inv;
    B(7:9, 4:6)   = etaL * I_inv * skew(rL);
    B(7:9, 7:9)   = etaR * I_inv;
    B(7:9, 10:12) = etaR * I_inv * skew(rR);

    B(10:12, 4:6)   = etaL * SX.eye(3) / m;
    B(10:12, 10:12) = etaR * SX.eye(3) / m;

    d(12) = -g;

    % DISCRETE SYSTEM
    % Explicit Euler
    % Ad = SX.eye(state_length) + A * dT;
    % Bd = B * dT;
    % dd = d * dT;
    % 
    % x_k_next = Ad * x_k + Bd * u_k + dd;
    % 
    % 4-th Runge-Kutta Method 
    k1 = A *  x_k + B * u_k + d;
    k2 = A * (x_k + 0.5 * dT * k1) + B * u_k + d;
    k3 = A * (x_k + 0.5 * dT * k2) + B * u_k + d;
    k4 = A * (x_k + dT * k3) + B * u_k + d;

    x_k_next = x_k + (dT / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    ceq1_sub = X((state_length * (i-1) + 1):(state_length * i)) - x_k_next;
    ceq1     = [ceq1; ceq1_sub];

    x_k = x_k_next;
end

% CAPTURABILITY CONSTRAINTS
ceq2 = []

b = sqrt(h/g);
dCOM_x = X(state_length * (H-1) + 10);
dCOM_y = X(state_length * (H-1) + 11);
dCOM_z = X(state_length * (H-1) + 12);

fL_x = U(input_length * (H-1) + 4);
fL_y = U(input_length * (H-1) + 5);
fL_z = U(input_length * (H-1) + 6);

fR_x = U(input_length * (H-1) + 10);
fR_y = U(input_length * (H-1) + 11);
fR_z = U(input_length * (H-1) + 12);

ceq2_x_sub = dCOM_x + b * (fL_x / m + fR_x / m);
ceq2_y_sub = dCOM_y + b * (fL_y / m + fR_y / m);
ceq2_z_sub = dCOM_z + b * (fL_z / m + fR_z / m- g);

ceq2 = [ceq2; ceq2_x_sub; ceq2_y_sub; ceq2_z_sub];

% FRICTION CONE
cineq1_max = []; cineq1_min = [];
cineq2_max = []; cineq2_min = [];
cineq3_max = []; cineq3_min = [];
cineq4_max = []; cineq4_min = [];
cineq5_max = []; cineq5_min = [];

f_z_max = SX.sym('f_z_max');
f_z_min = SX.sym('f_z_min');

footX = SX.sym('footX');    
footY = SX.sym('footY');    

for i = 1:H
    mL_x = U(input_length * (i-1) + 1);
    mL_y = U(input_length * (i-1) + 2);
    mL_z = U(input_length * (i-1) + 3);
    fL_x = U(input_length * (i-1) + 4);
    fL_y = U(input_length * (i-1) + 5);
    fL_z = U(input_length * (i-1) + 6);

    mR_x = U(input_length * (i-1) + 7);
    mR_y = U(input_length * (i-1) + 8);
    mR_z = U(input_length * (i-1) + 9);
    fR_x = U(input_length * (i-1) + 10);
    fR_y = U(input_length * (i-1) + 11);
    fR_z = U(input_length * (i-1) + 12);

    cineq1_max_l_sub = fL_z - f_z_max;
    cineq1_max_r_sub = fR_z - f_z_max;
    cineq1_min_l_sub =-fL_z + f_z_min;
    cineq1_min_r_sub =-fR_z + f_z_min;
    cineq1_max = [cineq1_max; cineq1_max_l_sub; cineq1_max_r_sub];
    cineq1_min = [cineq1_min; cineq1_min_l_sub; cineq1_min_r_sub];

    cineq2_max_l_sub = fL_x - mu * fL_z;
    cineq2_max_r_sub = fR_x - mu * fR_z;
    cineq2_min_l_sub =-fL_x - mu * fL_z;
    cineq2_min_r_sub =-fR_x - mu * fR_z;
    cineq2_max = [cineq2_max; cineq2_max_l_sub; cineq2_max_r_sub];
    cineq2_min = [cineq2_min; cineq2_min_l_sub; cineq2_min_r_sub];

    cineq3_max_l_sub = fL_y - mu * fL_z;
    cineq3_max_r_sub = fR_y - mu * fR_z;
    cineq3_min_l_sub =-fL_y - mu * fL_z;
    cineq3_min_r_sub =-fR_y - mu * fR_z;
    cineq3_max = [cineq3_max; cineq3_max_l_sub; cineq3_max_r_sub];
    cineq3_min = [cineq3_min; cineq3_min_l_sub; cineq3_min_r_sub];

    cineq4_max_l_sub = mL_x - footY * fL_z;
    cineq4_max_r_sub = mR_x - footY * fR_z;
    cineq4_min_l_sub =-mL_x - footY * fL_z;
    cineq4_min_r_sub =-mR_x - footY * fR_z;
    cineq4_max = [cineq4_max; cineq4_max_l_sub; cineq4_max_r_sub];
    cineq4_min = [cineq4_min; cineq4_min_l_sub; cineq4_min_r_sub];

    cineq5_max_l_sub = mL_y - footX * fL_z;
    cineq5_max_r_sub = mR_y - footX * fR_z;
    cineq5_min_l_sub =-mL_y - footX * fL_z;
    cineq5_min_r_sub =-mR_y - footX * fR_z;
    cineq5_max = [cineq5_max; cineq5_max_l_sub; cineq5_max_r_sub];
    cineq5_min = [cineq5_min; cineq5_min_l_sub; cineq5_min_r_sub];
end

ceq1_v = jacobian(ceq1, v);
ceq2_v = jacobian(ceq2, v);

cineq1_max_v = jacobian(cineq1_max, v);
cineq1_min_v = jacobian(cineq1_min, v);
cineq2_max_v = jacobian(cineq2_max, v);
cineq2_min_v = jacobian(cineq2_min, v);
cineq3_max_v = jacobian(cineq3_max, v);
cineq3_min_v = jacobian(cineq3_min, v);
cineq4_max_v = jacobian(cineq4_max, v);
cineq4_min_v = jacobian(cineq4_min, v);
cineq5_max_v = jacobian(cineq5_max, v);
cineq5_min_v = jacobian(cineq5_min, v);

%--- Function generation
disp('FUNCTION GENERATION START')

output_dir = strcat(folder, '/Function/');

J_v_func = Function('J_v_func',   {X, U, X_ref, U_ref, W_Q, W_R}, {J_v});
J_vv_func = Function('J_vv_func', {X, U, X_ref, U_ref, W_Q, W_R}, {J_vv});
 
ceq1_func = Function('ceq1_func', {x0, X, U, m, g, I, dT, rL_ref_horizon, rR_ref_horizon, theta_ref_horizon, etaL_ref_horizon, etaR_ref_horizon}, {ceq1});
ceq1_v_func = Function('ceq1_v_func', {x0, X, U, m, g, I, dT, rL_ref_horizon, rR_ref_horizon, theta_ref_horizon, etaL_ref_horizon, etaR_ref_horizon}, {ceq1_v});

ceq2_func = Function('ceq2_func', {X, U, m, g, h}, {ceq2});
ceq2_v_func = Function('ceq2_v_func', {X, U, m, g, h}, {ceq2_v});

cineq1_max_func = Function('cineq1_max_func', {U, f_z_max}, {cineq1_max});
cineq1_min_func = Function('cineq1_min_func', {U, f_z_min}, {cineq1_min});
cineq2_max_func = Function('cineq2_max_func', {U, mu}, {cineq2_max});
cineq2_min_func = Function('cineq2_min_func', {U, mu}, {cineq2_min});
cineq3_max_func = Function('cineq3_max_func', {U, mu}, {cineq3_max});
cineq3_min_func = Function('cineq3_min_func', {U, mu}, {cineq3_min});
cineq4_max_func = Function('cineq4_max_func', {U, footY}, {cineq4_max});
cineq4_min_func = Function('cineq4_min_func', {U, footY}, {cineq4_min});
cineq5_max_func = Function('cineq5_max_func', {U, footX}, {cineq5_max});
cineq5_min_func = Function('cineq5_min_func', {U, footX}, {cineq5_min});

cineq1_max_v_func = Function('cineq1_max_v_func', {U, f_z_max}, {cineq1_max_v});
cineq1_min_v_func = Function('cineq1_min_v_func', {U, f_z_min}, {cineq1_min_v});
cineq2_max_v_func = Function('cineq2_max_v_func', {U, mu}, {cineq2_max_v});
cineq2_min_v_func = Function('cineq2_min_v_func', {U, mu}, {cineq2_min_v});
cineq3_max_v_func = Function('cineq3_max_v_func', {U, mu}, {cineq3_max_v});
cineq3_min_v_func = Function('cineq3_min_v_func', {U, mu}, {cineq3_min_v});
cineq4_max_v_func = Function('cineq4_max_v_func', {U, footY}, {cineq4_max_v});
cineq4_min_v_func = Function('cineq4_min_v_func', {U, footY}, {cineq4_min_v});
cineq5_max_v_func = Function('cineq5_max_v_func', {U, footX}, {cineq5_max_v});
cineq5_min_v_func = Function('cineq5_min_v_func', {U, footX}, {cineq5_min_v});
disp('FUNCTION GENERATION END')

%--- Function generation
opts=struct('main',true,'mex',true,'with_header',true);

% Cost function -> Function('name', {Inputs,...}, {Output})
disp('HESSIAN CODE GENERATION START')
tic; 
J_vv_func.generate('J_vv_func.c', opts);
mex J_vv_func.c;
elapsedTime = toc;
disp(['HESSIAN CODE GENERATION COMPLETED IN ', num2str(elapsedTime), ' SECONDS'])

disp('GRADIENT CODE GENERATION START')
tic; 
J_v_func.generate('J_v_func.c', opts);
mex J_v_func.c
elapsedTime = toc;
disp(['GRADIENT CODE GENERATION COMPLETED IN ', num2str(elapsedTime), ' SECONDS'])

disp('CONSTRAINTS CODE GENERATION START')
tic; 
ceq1_func.generate('ceq1_func.c', opts);
mex ceq1_func.c
ceq1_v_func.generate('ceq1_v_func.c', opts);
mex ceq1_v_func.c
ceq2_func.generate('ceq2_func.c', opts);
mex ceq2_func.c
ceq2_v_func.generate('ceq2_v_func.c', opts);
mex ceq2_v_func.c

cineq1_max_func.generate('cineq1_max_func.c', opts);
mex cineq1_max_func.c
cineq1_min_func.generate('cineq1_min_func.c', opts);
mex cineq1_min_func.c
cineq2_max_func.generate('cineq2_max_func.c', opts);
mex cineq2_max_func.c
cineq2_min_func.generate('cineq2_min_func.c', opts);
mex cineq2_min_func.c
cineq3_max_func.generate('cineq3_max_func.c', opts);
mex cineq3_max_func.c
cineq3_min_func.generate('cineq3_min_func.c', opts);
mex cineq3_min_func.c
cineq4_max_func.generate('cineq4_max_func.c', opts);
mex cineq4_max_func.c
cineq4_min_func.generate('cineq4_min_func.c', opts);
mex cineq4_min_func.c
cineq5_max_func.generate('cineq5_max_func.c', opts);
mex cineq5_max_func.c
cineq5_min_func.generate('cineq5_min_func.c', opts);
mex cineq5_min_func.c

cineq1_max_v_func.generate('cineq1_max_v_func.c', opts);
mex cineq1_max_v_func.c
cineq1_min_v_func.generate('cineq1_min_v_func.c', opts);
mex cineq1_min_v_func.c
cineq2_max_v_func.generate('cineq2_max_v_func.c', opts);
mex cineq2_max_v_func.c
cineq2_min_v_func.generate('cineq2_min_v_func.c', opts);
mex cineq2_min_v_func.c
cineq3_max_v_func.generate('cineq3_max_v_func.c', opts);
mex cineq3_max_v_func.c
cineq3_min_v_func.generate('cineq3_min_v_func.c', opts);
mex cineq3_min_v_func.c
cineq4_max_v_func.generate('cineq4_max_v_func.c', opts);
mex cineq4_max_v_func.c
cineq4_min_v_func.generate('cineq4_min_v_func.c', opts);
mex cineq4_min_v_func.c
cineq5_max_v_func.generate('cineq5_max_v_func.c', opts);
mex cineq5_max_v_func.c
cineq5_min_v_func.generate('cineq5_min_v_func.c', opts);
mex cineq5_min_v_func.c
elapsedTime = toc;
disp(['CONSTRAINTS CODE GENERATION COMPLETED IN ', num2str(elapsedTime), ' SECONDS'])
% ---