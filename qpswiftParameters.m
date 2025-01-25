function [P, c, A, b, G, h] = qpswiftParameters(x0, theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, ...
                                                rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon)

H = PARA.H;
state_length = PARA.state_length;
input_length = PARA.input_length;

m = PARA.m_all;
g = PARA.g;
zc = PARA.zc;
I = PARA.I;
dT = PARA.dt_MPC;

% Gain vectors
gain_state_horizon = zeros(H*state_length, 1);
gain_input_horizon = zeros(H*input_length, 1);
for i = 1:H
    gain_state_horizon((i-1)*state_length + 1 : (i-1)*state_length + 3,  1) = PARA.Q_theta;
    gain_state_horizon((i-1)*state_length + 4,  1) = PARA.Q_COM_x;
    gain_state_horizon((i-1)*state_length + 5,  1) = PARA.Q_COM_y;
    gain_state_horizon((i-1)*state_length + 6,  1) = PARA.Q_COM_z;
    gain_state_horizon((i-1)*state_length + 7 : (i-1)*state_length + 9,  1) = PARA.Q_w;
    gain_state_horizon((i-1)*state_length + 10: (i-1)*state_length + 12, 1) = PARA.Q_dCOM;

    gain_input_horizon((i-1)*input_length + 1 : (i-1)*state_length + 3,  1) = PARA.R_mL;
    gain_input_horizon((i-1)*input_length + 4 : (i-1)*state_length + 6,  1) = PARA.R_fL;
    gain_input_horizon((i-1)*input_length + 7 : (i-1)*state_length + 9,  1) = PARA.R_mR;
    gain_input_horizon((i-1)*input_length + 10: (i-1)*state_length + 12, 1) = PARA.R_fR;   
end

% MPC Reference
X = zeros(H*state_length, 1);
X_ref = zeros(H*state_length, 1);
U = zeros(H*input_length, 1);
U_ref = zeros(H*input_length, 1);
for i = 1:H
    X_ref((i-1)*state_length + 1 : (i-1)*state_length + 3)  = theta_ref_horizon(:, i);
    X_ref((i-1)*state_length + 4 : (i-1)*state_length + 6)  = COM_ref_horizon(:,i);
    X_ref((i-1)*state_length + 7 : (i-1)*state_length + 9)  = w_ref_horizon(:,i);
    X_ref((i-1)*state_length + 10: (i-1)*state_length + 12) = dCOM_ref_horizon(:,i);
end

c = J_v_func(X, U, X_ref, U_ref, gain_state_horizon, gain_input_horizon);
P = J_vv_func(X, U, X_ref, U_ref, gain_state_horizon, gain_input_horizon);

ceq1   = ceq1_func(x0, X, U, m, g, I, dT, rL_ref_horizon, rR_ref_horizon, theta_ref_horizon, etaL_ref_horizon, etaR_ref_horizon);
ceq1_v = ceq1_v_func(x0, X, U, m, g, I, dT, rL_ref_horizon, rR_ref_horizon, theta_ref_horizon, etaL_ref_horizon, etaR_ref_horizon);

ceq2   = ceq2_func(X, U, m, g, zc);
ceq2_v = ceq2_v_func(X, U, m, g, zc);

cineq1_max = cineq1_max_func(U, PARA.f_z_max);
cineq1_min = cineq1_min_func(U, PARA.f_z_min);
cineq2_max = cineq2_max_func(U, PARA.mu);
cineq2_min = cineq2_min_func(U, PARA.mu);
cineq3_max = cineq3_max_func(U, PARA.mu);
cineq3_min = cineq3_min_func(U, PARA.mu);
cineq4_max = cineq4_max_func(U, PARA.Foot_width / 2.0);
cineq4_min = cineq4_min_func(U, PARA.Foot_width / 2.0);
cineq5_max = cineq5_max_func(U, PARA.Foot_length/ 2.0);
cineq5_min = cineq5_min_func(U, PARA.Foot_length/ 2.0);

cineq1_max_v = cineq1_max_v_func(U, PARA.f_z_max);
cineq1_min_v = cineq1_min_v_func(U, PARA.f_z_min);
cineq2_max_v = cineq2_max_v_func(U, PARA.mu);
cineq2_min_v = cineq2_min_v_func(U, PARA.mu);
cineq3_max_v = cineq3_max_v_func(U, PARA.mu);
cineq3_min_v = cineq3_min_v_func(U, PARA.mu);
cineq4_max_v = cineq4_max_v_func(U, PARA.Foot_width / 2.0);
cineq4_min_v = cineq4_min_v_func(U, PARA.Foot_width / 2.0);
cineq5_max_v = cineq5_max_v_func(U, PARA.Foot_length/ 2.0);
cineq5_min_v = cineq5_min_v_func(U, PARA.Foot_length/ 2.0);
    
% A = [ceq1_v];
% b = (-1).*[ceq1];
% 
A = [ceq1_v;
     ceq2_v];
b = (-1).*[ceq1;
           ceq2];

G = [cineq1_max_v;
     cineq1_min_v;
     cineq2_max_v;
     cineq2_min_v;
     cineq3_max_v;
     cineq3_min_v;
     cineq4_max_v;
     cineq4_min_v;
     cineq5_max_v;
     cineq5_min_v];
h = (-1).*[cineq1_max;
           cineq1_min;
           cineq2_max;
           cineq2_min;
           cineq3_max;
           cineq3_min;
           cineq4_max;
           cineq4_min;
           cineq5_max;
           cineq5_min];
end