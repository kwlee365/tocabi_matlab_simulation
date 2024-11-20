addpath('casadi-3.6.7-linux64-matlab2018b')
import casadi.*

% Orientation-Aware Model Predictive Control with Footstep Adaptation for Dynamic Humanoid Walking

PARA = PARA;

H = PARA.H;
state_length = PARA.state_length;
input_length = PARA.input_length;

Q = SX.sym('Q', H * state_length);
R = SX.sym('R', H * input_length);

M = SX.sym('M');
g = SX.sym('g');
I = SX.sym('I', 3, 3);
dT_mpc = SX.sym('dT_mpc'); 

% p_c    = SX.sym('p_c', 3);      % com position w.r.t. global frame
% dp_c   = SX.sym('dp_c', 3);     % com velocity w.r.t. global frame
% theta  = SX.sym('theta', 3);    % body orientation w.r.t. glboal frame
% dtheta = SX.sym('dtheta', 3);   % body angular velocity w.r.t. global frame
% c      = SX.sym('c', 3);        % foot contact location

dc = SX.sym('dc', 3); % footstep adaptation
eta = SX.sym('eta', H);

dc_x_max = SX.sym('dc_x_max');
dc_y_max = SX.sym('dc_y_max');
dc_x_min = SX.sym('dc_x_min');
dc_y_min = SX.sym('dc_y_min');

r_x_max = SX.sym('r_x_max');
r_x_min = SX.sym('r_x_max');
r_y_max = SX.sym('r_y_max');
r_y_min = SX.sym('r_y_max');

mu = SX.sym('mu');
Foot_length = SX.sym('Foot_length');
Foot_width = SX.sym('Foot_width');

x_init = SX.sym('x_init', state_length);

X  = SX.sym('X',     H*state_length);   % 1 ~ H
Xd = SX.sym('X_ref', H*state_length);   % 0 ~ H-1
U  = SX.sym('U',     H*input_length);
Ud = SX.sym('U_ref', H*input_length);

v = [X;U;dc];

% Cost function
J = 0;
for i = 1:H
    x  = X([(i-1)*state_length + 1 : i*state_length - 3]); % state (1 ~ 12)
    xd =Xd([(i-1)*state_length + 1 : i*state_length - 3]);
    u  = U([(i-1)*input_length + 1 : i*input_length]); %input (1 ~ 6)
    ud =Ud([(i-1)*input_length + 1 : i*input_length]); 

    c  = X([i*state_length - 2 : i*state_length]);  % state (13 ~ 15) 
    cd =Xd([i*state_length - 2 : i*state_length]);  

    J = J + (xd - x)'      * diag(Q((i-1)*state_length + 1:i*state_length - 3)) * (xd - x) ...
          + (ud - u)'      * diag(R((i-1)*input_length + 1:i*input_length))     * (ud - u) ...
          + (cd - c - dc)' * diag(Q(i*state_length - 2    :i*state_length))     * (cd - c - dc);
end

J_v = jacobian(J, v);
J_vv = hessian(J, v);

% Equality Constraint - discrete dynamics
ceq1 = [];
f = SX.sym('f', state_length); 
x = x_init;
for i = 1:H
    % state
    p_c    = x(1:3);     
    dp_c   = x(4:6);
    theta  = x(7:9);
    dtheta = x(10:12);
    c      = x(13:15);
    r = c - p_c;

    % control input
    u = U((i-1)*input_length + 1 :(i-1)*input_length + 6);
    F = u(1:3);
    m = u(4:6);
     
    % SRBD dynamics
    f([1:3])   = dp_c;
    f([4:6])   = (1/M) * F + [0;0;-g];
    f([7:9])   = dtheta;
    f([10:12]) = I\(cross(r,F) + m);
    f([13:15]) = zeros(3,1);

    % state space model
    df_dx = jacobian(f,x);
    df_du = jacobian(f,u);
    A = eye(state_length, state_length) + dT_mpc * df_dx;
    B = dT_mpc * df_du;
    d = dT_mpc * (f - A * x - B * u);
    Ac= A(:,[13:15]);

    % equality constraint
    ceq1_sub = X((i-1)*state_length + 1 :(i-1)*state_length + 15) - (A * x + B * u + d + eta(i) * Ac * dc);
    ceq1 = [ceq1 ; ceq1_sub];

    % next step
    x = X((i-1)*state_length+1:(i-1)*state_length+15); 
end

ceq1_v = jacobian(ceq1, v);

% Inequality constraint (1) - footstep location constraint
cineq1_max = [dc(1) - dc_x_max;
              dc(2) - dc_y_max;
              dc(3) - 0.0];
cineq1_min =[-dc(1) + dc_x_min;
             -dc(2) + dc_y_min;
             -dc(3) + 0.0];
cineq1_max_v = jacobian(cineq1_max, v);
cineq1_min_v = jacobian(cineq1_min, v);

% Inequality constraint (2) - kinematic reacheability constraint (length between CoM and foot projected on the ground.) 
c = [];
cineq2_max = [];
for i = 1:H
    p_c = X((i-1)*state_length+1:(i-1)*state_length+3);
    c   = X((i-1)*state_length+13:(i-1)*state_length+15); 
    
    r = (c + dc - p_c);
    rx = r(1);
    ry = r(2);
    
    cineq2_max_x_sub = rx - r_x_max;
    cineq2_max_y_sub = ry - r_y_max;

    cineq2_max = [cineq2_max; cineq2_max_x_sub;cineq2_max_y_sub];
end

cineq2_min = [];
for i = 1:H
    p_c = X((i-1)*state_length+1:(i-1)*state_length+3);
    c   = X((i-1)*state_length+13:(i-1)*state_length+15); 
    
    r = (c + dc - p_c);
    rx = r(1);
    ry = r(2);
    
    cineq2_min_x_sub =-rx + r_x_min;
    cineq2_min_y_sub =-ry + r_y_min;

    cineq2_min = [cineq2_min; cineq2_min_x_sub;cineq2_min_y_sub];
end
cineq2_max_v = jacobian(cineq2_max, v);
cineq2_min_v = jacobian(cineq2_min, v);

cineq3_max = [];
for i = 1:H
    u = U((i-1)*input_length + 1 :(i-1)*input_length + 6);
    fx = u(1);
    fy = u(2);
    fz = u(3);
    mx = u(4);
    my = u(5);
    mz = u(6);

    cineq3_max_sub = [fx - mu * fz;
                      fy - mu * fz;
                     -fz;
                      mx - Foot_width * fz;
                      my - Foot_length / 2.0 * fz];

    cineq3_max = [cineq3_max; cineq3_max_sub];
end

cineq3_min = [];
for i = 1:H
    u = U((i-1)*input_length + 1 :(i-1)*input_length + 6);
    fx = u(1);
    fy = u(2);
    fz = u(3);
    mx = u(4);
    my = u(5);
    mz = u(6);

    cineq3_min_sub = [-fx + mu * fz;
                      -fy + mu * fz;
                       fz;
                      -mx + Foot_width * fz;
                      -my + Foot_length / 2.0 * fz];

    cineq3_min = [cineq3_min; cineq3_min_sub];
end

cineq1_max_v = jacobian(cineq1_max, v);
cineq1_min_v = jacobian(cineq1_min, v);
cineq2_max_v = jacobian(cineq2_max, v);
cineq2_min_v = jacobian(cineq2_min, v);
cineq3_max_v = jacobian(cineq3_max, v);
cineq3_min_v = jacobian(cineq3_min, v);

%--- CasADi Function Generation and Compilation ---
opts = struct('main', true, 'mex', true, 'with_header', true);

% cost function generation
J_v_func = Function('J_v_func', {Q, R, X, Xd, U, Ud, dc}, {J_v});
J_vv_func = Function('J_vv_func', {Q, R, X, Xd, U, Ud, dc}, {J_vv});

J_v_func.generate('J_v_func.c', opts);
mex J_v_func.c
J_vv_func.generate('J_vv_func.c', opts);
mex J_vv_func.c

% equality constraint generation
ceq1_func = Function('ceq1_func', {x_init, X, U, M, g, I, dT_mpc, eta, dc}, {ceq1});
ceq1_v_func = Function('ceq1_v_func', {x_init, X, U, M, g, I, dT_mpc, eta, dc}, {ceq1_v});

ceq1_func.generate('ceq1_func.c', opts);
mex ceq1_func.c
ceq1_v_func.generate('ceq1_v_func.c', opts);
mex ceq1_v_func.c

cineq1_max_func = Function('cineq1_max_func', {dc, dc_x_max, dc_y_max}, {cineq1_max});
cineq1_min_func = Function('cineq1_min_func', {dc, dc_x_min, dc_y_min}, {cineq1_min});

cineq2_max_func = Function('cineq2_max_func', {X, dc, r_x_max, r_y_max}, {cineq2_max});
cineq2_min_func = Function('cineq2_min_func', {X, dc, r_x_min, r_y_min}, {cineq2_min});

cineq3_max_func = Function('cineq3_max_func', {U, mu, Foot_width, Foot_length}, {cineq3_max});
cineq3_min_func = Function('cineq3_min_func', {U, mu, Foot_width, Foot_length}, {cineq3_min});

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

cineq1_max_v_func = Function('cineq1_max_v_func', {dc, dc_x_max, dc_y_max}, {cineq1_max_v});
cineq1_min_v_func = Function('cineq1_min_v_func', {dc, dc_x_min, dc_y_min}, {cineq1_min_v});

cineq2_max_v_func = Function('cineq2_max_v_func', {X, dc, r_x_max, r_y_max}, {cineq2_max_v});
cineq2_min_v_func = Function('cineq2_min_v_func', {X, dc, r_x_min, r_y_min}, {cineq2_min_v});

cineq3_max_v_func = Function('cineq3_max_v_func', {U, mu, Foot_width, Foot_length}, {cineq3_max_v});
cineq3_min_v_func = Function('cineq3_min_v_func', {U, mu, Foot_width, Foot_length}, {cineq3_min_v});

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