%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL PREDICTIVE CONTROLLER (SRBD)
%%%
%%% Author: Kwanwoo Lee (kwlee365@snu.ac.kr)
%%% Date: 2025. 01. 21. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Path setting
clc; clear all; close all;
restoredefaultpath;
folder = fileparts(which('main.m'));
addpath(genpath(folder));
%---

%--- Robot setting
tocabi = importrobot('dyros_tocabi.urdf');
tocabi.DataFormat = 'column';  
tocabi.Gravity = [0; 0; -9.81];
q_dim = 40
q = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    
    0.0;        % L_HipYaw_Joint
    0.0;        % L_HipRoll_Joint
   -0.28;       % L_HipPitch_Joint
    0.6;        % L_Knee_Joint
   -0.32;       % L_AnklePitch_Joint
    0.0;        % L_AnkleRoll_Joint

    0.0;        % R_HipYaw_Joint
    0.0;        % R_HipRoll_Joint
   -0.28;       % R_HipPitch_Joint
    0.6;        % R_Knee_Joint
   -0.32;       % R_AnklePitch_Joint
    0.0;        % R_AnkleRoll_Joint

    0.0;        % Waist1_Joint
    0.0;        % Waist2_Joint
    0.0;        % Upperbody_Joint

    0.3;        % L_Shoulder1_Joint
    0.174533;   % L_Shoulder2_Joint
    1.22173;    % L_Shoulder3_Joint
   -1.27;       % L_Armlink_Joint
   -1.57;       % L_Elbow_Joint
    0.0;        % L_Forearm_Joint
   -1.0;        % L_Wrist1_Joint
    0.0;        % L_Wrist2_Joint

    0.0;        % Neck_Joint
    0.0;        % Head_Joint

   -0.3;        % R_Shoulder1_Joint
   -0.174533;   % R_Shoulder2_Joint
   -1.22173;    % R_Shoulder3_Joint
    1.27;       % R_Armlink_Joint
    1.57;       % R_Elbow_Joint
    0.0;        % R_Forearm_Joint
    1.0;        % R_Wrist1_Joint
    0.0         % R_Wrist2_Joint
];
%---

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step information
number_of_step = 10;             % Number of steps
step_length = 0.3;              % Step stride
step_width = PARA.pelvis_width;  % Step width
step_time = 0.5;                 % Step period
L_or_R = 1;                      % First swing foot: 1: Left foot / -1: Right foot

% Disturbance information
Impact_force_x = 0;              % [N]: x-dir impact force
Impact_force_y = 500;            % [N]: y-dir impact force
Impact_duration = 0.05;          % [s]: Impact duration
Impact_timing = 0.3;             % [s]: Timing of impact
Impact_step_number = 3;         

% Flags
flag_HORIZON_CHANGED = 0;       % Set to 1 if the number of MPC horizon is changed
flag_VISUALIZATION = 1;         % Set to 1 for graphic ON
flag_VISUALIZATION_ROBOT = 1;   % Set to 1 to show robot
flag_PLOT = 0;                  % Set to 1 to show plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Get gradients and hessians with CasADI
if flag_HORIZON_CHANGED == 1
    disp(['Now creating gradients and hessians with CasADI for the MPC horizon: ', num2str(PARA.H), '...']);
    getGradHessWithCasadi();
end
%---
%--- Global variables
global flag_EXIT flag_PAUSE;
global dx dy;
dx = 0; dy = 0;
%---

%--- SDB(Step Data Buffer)
p_ref_total = zmpTotal(number_of_step, step_length, step_width, L_or_R);
p_total = zmpTotal(number_of_step, step_length, step_width, L_or_R);

T_step_ref_total = stepTimeTotal(number_of_step, step_time);
T_step_total = stepTimeTotal(number_of_step, step_time);
%---

%--- Initialization
t = 0; t_step = 0;
i = 1; i_max = 1E06;
T_step = T_step_total(:, 1);
T_step_ref = T_step_ref_total(:, 1);
step_phase = 1;
iter_error = 0;
% Flags
flag_STEP_CHANGE = 0;
flag_EXIT = 0;
flag_PAUSE = 0;
flag_ERROR = 0;
% Preview control
[Gi, Gx, Gp] = findPreviewGain(PARA.T_preview, PARA.dt, PARA.zc);
[Gi_MPC, Gx_MPC, Gp_MPC] = findPreviewGain(PARA.T_preview, PARA.dt_MPC, PARA.zc);
p_err_sum_x = 0; p_err_sum_y = 0;
p_err_sum_x_ref = 0; p_err_sum_y_ref = 0;
% ZMP
p_des = p_total(:, 1);
% COM
COM = [0; 0; PARA.zc];
dCOM = [0; 0; 0];
COM_prev_step = COM;
dCOM_prev_step = dCOM;
COM_err = [0; 0; 0];
% COM ref.
COM_ref = [0; 0; PARA.zc];
dCOM_ref = [0; 0; 0];
ddCOM_ref = [0; 0; 0];
% Foot
Foot_state = 2; % DSP
LF = [0;  0.5*PARA.pelvis_width; 0]; LF_prev = LF;
RF = [0; -0.5*PARA.pelvis_width; 0]; RF_prev = RF;
color_LF = [0.9290 0.6940 0.1250];
color_RF = [0.6 0.6 0.6];
% Torso
theta = [0; 0; 0];
w = [0; 0; 0];
%---

%--- Data save
t_stored = zeros(1, i_max);
t_step_stored = zeros(1, i_max);
impact_force_stored = zeros(2, i_max);
T_step_stored = zeros(1, i_max);
COM_stored = zeros(3, i_max);
dCOM_stored = zeros(3, i_max);
COM_ref_stored = zeros(3, i_max);
theta_stored = zeros(3, i_max);

p_stored = zeros(3, i_max); 
mL_stored = zeros(3, i_max); 
fL_stored = zeros(3, i_max); 
mR_stored = zeros(3, i_max); 
fR_stored = zeros(3, i_max);
etaL_stored = zeros(1, i_max);
etaR_stored = zeros(1, i_max);
LF_stored = zeros(3, i_max);
RF_stored = zeros(3, i_max);
ticktock_stored = zeros(1, i_max);

%%

%--- Main Loop
while 1
    tic;

    % Step change
    flag_STEP_CHANGE = checkStepEnd(t_step, T_step);
    if flag_STEP_CHANGE == 1
        % Update step phase
        step_phase = step_phase + 1;
        if step_phase > number_of_step + 3
            flag_EXIT = 1;
        end

        % Update foot state
        Foot_state = (-1)*Foot_state;
        if step_phase == 2
            Foot_state = L_or_R;
        elseif step_phase >= number_of_step + 3
            Foot_state = 2;
        end

        % Update step time
        if step_phase <= number_of_step + 3
            T_step_ref = T_step_ref_total(1, step_phase);
            T_step = T_step_total(:, step_phase);
        end

        % Reset t_step & flag
        t_step = 0;
        flag_STEP_CHANGE = 0; 
    end

    % Exit flag
    if (norm(COM_err) > 0.2)
        disp('Walking fail!!');
        break;
    elseif (flag_EXIT == 1)
        disp('Walking finish!!');
        break;
    end

    % Disturbance
    disturbance_duration = Impact_duration; % [sec]
    disturbance_timing = Impact_timing;
    if (step_phase == Impact_step_number + 1) && ((t_step >= disturbance_timing))
        flag_ERROR = 1;
    end
    if (flag_ERROR == 1) && (iter_error > (Impact_duration/PARA.dt))
        flag_ERROR = 0;
    end
    if flag_ERROR == 1
        disturbance_magnitude = [Impact_force_x; Impact_force_y; 0]; % [N]
        ddCOM_dist = disturbance_magnitude/PARA.m_all;
        dCOM_dist = ddCOM_dist.*PARA.dt;
        COM_dist = dCOM_dist.*PARA.dt + 0.5.*ddCOM_dist.*PARA.dt.*PARA.dt;
        COM = COM + COM_dist;
        dCOM = dCOM + dCOM_dist;
        iter_error = iter_error + 1;
    else     
        disturbance_magnitude = [0; 0; 0]; % [N]
        ddCOM_err = [0; 0; 0];
    end
    
    % Reference    
    [COM_ref_next, dCOM_ref_next, ddCOM_ref_next, p_err_sum_x_ref_next, p_err_sum_y_ref_next] = previewControl(t_step, step_phase, p_ref_total, T_step_ref_total, Gi, Gx, Gp, PARA.A_preview, PARA.B_preview, PARA.C_preview, COM_ref, dCOM_ref, ddCOM_ref, p_err_sum_x_ref, p_err_sum_y_ref);      
    [theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon] = ...
        mpcRefWindow(t_step, step_phase, Foot_state, L_or_R, LF_prev, RF_prev, ...
                     COM_ref, dCOM_ref, ddCOM_ref, ...
                     p_err_sum_x_ref, p_err_sum_y_ref, T_step_ref, p_ref_total, T_step_ref_total, Gi_MPC, Gx_MPC, Gp_MPC);
   
    COM_err = COM - COM_ref;

    % Calc control input
    x0 = [theta; COM; w; dCOM];
    [mL, fL, mR, fR] =  nextState(x0, ...
                                  theta, COM, w, dCOM, ...
                                  theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, ...
                                  rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon);

    contact_wrench_result = [mL; fL; mR; fR];
    rL = LF_prev - COM;
    rR = RF_prev - COM;

    % Plant response
    t_span = [0 PARA.dt];
    y0 = [theta; COM; w; dCOM];
    [t_ode, y_ode] = ode45(@(t_ode, y_ode) odefunc(y_ode, contact_wrench_result, Foot_state, PARA.m_all, PARA.I, PARA.g, rL, rR), t_span, y0);  
    theta_next = [y_ode(end, [1:3])]';
    COM_next   = [y_ode(end, [4:6])]';
    w_next     = [y_ode(end, [7:9])]';
    dCOM_next  = [y_ode(end, [10:12])]';

    % Foot trajectory
    if Foot_state == 2
        LF = LF_prev;
        RF = RF_prev;
    elseif Foot_state ==  1 % LF swing
        LF = footTrajectory(t_step, step_phase, number_of_step, Foot_state, T_step, p_total);
        RF = RF_prev;
    elseif Foot_state == -1 % RF swing
        LF = LF_prev;
        RF = footTrajectory(t_step, step_phase, number_of_step, Foot_state, T_step, p_total);
    end
    if (step_phase == number_of_step + 2)
        LF = LF_prev;
        RF = RF_prev;
    end
    
    % Time ticktock
    ticktock = toc;
    
    % Data save
    t_stored(:, i) = t;
    t_step_stored(:, i) = t_step;
    T_step_stored(:, i) = T_step;
    COM_stored(:, i) = COM;
    dCOM_stored(:, i) = dCOM;
    theta_stored(:,i) = theta;
    COM_ref_stored(:, i) = COM_ref;
    p_stored(:, i) = p_ref_total(:, step_phase);
    mL_stored(:, i) = mL;
    fL_stored(:, i) = fL;
    mR_stored(:, i) = mR;
    fR_stored(:, i) = fR;
    LF_stored(:, i) = LF;
    RF_stored(:, i) = RF;
    etaL_stored(1, i) = etaL_ref_horizon(1);
    etaR_stored(1, i) = etaR_ref_horizon(1);
    ticktock_stored(:, i) = ticktock;

    % Time waits for no one
    t = t + PARA.dt;
    t_step = t_step + PARA.dt;

    COM_ref = COM_ref_next; dCOM_ref = dCOM_ref_next; ddCOM_ref = ddCOM_ref_next;
    p_err_sum_x_ref = p_err_sum_x_ref_next; p_err_sum_y_ref = p_err_sum_y_ref_next;
    
    theta = theta_next; w = w_next;
    COM = COM_next; dCOM = dCOM_next;
    
    LF_prev = LF;
    RF_prev = RF;

    i = i + 1;
end
%---

last_tick = i-1;

%-- Plot
t_stored = t_stored(:, 2:i-1);
t_step_stored = t_step_stored(:, 2:i-1);
impact_force_stored = impact_force_stored(:, 2:i-1);
T_step_stored = T_step_stored(:, 2:i-1);
COM_stored = COM_stored(:, 2:i-1); COM_ref_stored = COM_ref_stored(:, 2:i-1); dCOM_stored = dCOM_stored(:, 2:i-1);
theta_stored = theta_stored (:,2:i-1) 
p_stored = p_stored(:, 2:i-1); 
mL_stored = mL_stored(:, 2:i-1); fL_stored = fL_stored(:, 2:i-1); mR_stored = mR_stored(:, 2:i-1); fR_stored = fR_stored(:, 2:i-1);
etaR_stored = etaR_stored(:, 2:i-1); etaL_stored = etaL_stored(:, 2:i-1);
LF_stored = LF_stored(:, 2:i-1); RF_stored = RF_stored(:, 2:i-1);
ticktock_stored = ticktock_stored(:, 2:i-1);

%% Animation
%--- World generation
if ishandle(10)
    close(10)
end

if flag_VISUALIZATION == 1
    fig_world = figure(10);
    set(fig_world, 'Position', [1000 540 920 455], 'Renderer', 'OpenGL', 'Color',[1,1,1], 'KeyPressFcn', @printfig);
    axe = axes('Parent', fig_world);
    if flag_VISUALIZATION_ROBOT == 1
        ax = 55; ay = 15;
    else
        ax = 0; ay = 89.999;
    end
    view([ax ay]);
    set(axe, 'XLim', [-0.5 0.5], 'YLim', [-0.5 0.5], 'ZLim', [-0.02 1], 'DataAspectRatio', [1 1 1]);
    grid on; grid minor;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
end

i=1;
flag_EXIT= 0;
force_draw_scale = 0.002;

while 1

    COM = COM_stored(:,i);
    COM_ref = COM_ref_stored(:,i);
    theta = theta_stored (:,i);
    LF = LF_stored(:,i);
    RF = RF_stored(:,i);

    if flag_VISUALIZATION == 1
        % Camera control
        ax = ax - dx;
        ay = ay - dy;
        dx = 0; dy = 0;
        view([ax, ay]);
        set(axe,'XLim',[-1.0+COM(1) 1.0+COM(1)],'YLim',[-1.0+COM(2) 1.0+COM(2)],'ZLim',[-0.1 2.0], 'DataAspectRatio', [1 1 1]);
    end

    if i > last_tick
        flag_EXIT = 1
    end

        % Visualization
    if mod(i, 100) == 1  
        if flag_VISUALIZATION == 1
            if flag_VISUALIZATION_ROBOT == 1
                %--- Inverse kinematics - COM
                pCOM = COM; 
                qPEL = mat2quat(rotX_rad(theta(1))*rotY_rad(-theta(2))); rotmPEL = quat2mat(qPEL);
                pLF = LF;   
                qLF = [1, 0, 0, 0];
                pRF = RF;   
                qRF = [1, 0, 0, 0];
            
                x_target = [pCOM; qPEL'; pLF; qLF'; pRF; qRF'];
                [q_target, pPEL] = IK_COM(x_target);
                
                q([1:4],1) = qPEL;
                q([5:7],1) = pPEL;
                q([8:19],1)  = q_target;

                show(tocabi,q,PreservePlot=false,FastUpdate=true);
            end

            % COM
            visual_COM = animatedline('Marker', 'o', 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black');
            addpoints(visual_COM, COM(1), COM(2), COM(3));
    
            visual_COM_ref = animatedline('Marker', 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
            addpoints(visual_COM_ref, COM_ref(1), COM_ref(2), COM_ref(3));
        
            hold on
            % visual_fL = quiver3(LF(1) , LF(2), LF(3), fL(1) * force_draw_scale, fL(2) * force_draw_scale, fL(3) * force_draw_scale, 0, ...
            %                     'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1);
            % visual_fR = quiver3(RF(1), RF(2), RF(3), fR(1) * force_draw_scale, fR(2) * force_draw_scale, fR(3) * force_draw_scale, 0, ...
            %                     'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1);        
            drawnow;
            if flag_PAUSE == 1
                disp('Walking pause!!');
                waitforbuttonpress;
                flag_PAUSE = 0;
            end
            if flag_EXIT == 1
                disp("Walking Finish!!");
                close(10);
                break;
            end

            delete(visual_COM);
            delete(visual_COM_ref);
            % delete(visual_fL); delete(visual_fR);
        end
    end

    i = i + 1;
end


%%
if flag_PLOT == 1
    if ishandle(1)
        close(1)
    end
    fig1 = figure(1);
    set(fig1, 'Position', [100 100 600 800])
    tile = tiledlayout(5, 2); 
    tile.TileSpacing = 'compact';
    tile.Padding = 'compact';
    nexttile([1 2])
    subplot0 = plot(t_stored, impact_force_stored(1, :), 'color', 'k' , 'linewidth', 2.0);
    legend(subplot0, 'Disturbance', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[N]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot1 = plot(t_stored, xi_err_stored(1, :), 'color', 'b' , 'linewidth', 2.0);
    legend(subplot1, '$\xi^{err}_{x}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot2 = plot(t_stored, p_c_stored(1, :), '-', 'color', 'r', 'linewidth', 2.0);
    legend(subplot2, '$p_{\mathrm{c,ZMP},x}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot3 = plot(t_stored, dU_stored(1, :), '-', 'color', [0, 0.5, 0], 'linewidth', 2.0);
    legend(subplot3, '$\delta u_{x}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot4 = plot(t_stored, db_stored(1, :), '-', 'color', 'k', 'linewidth', 2.0);
    legend(subplot4, '$\delta b_{x}$', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot5 = plot(t_stored, dT_stored, '-', 'color', [0.75, 0, 0.75], 'linewidth', 2.0);
    legend(subplot5, '$\delta T$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[s]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot6 = plot(t_stored, ddtheta_result_stored(2, :).*PARA.R2D, 'color', [0, 0.75, 0.75], 'linewidth', 2.0);
    legend(subplot6, '$\ddot{\theta}_{pitch}$', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[deg/s^2]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 2])
    histogram1 = histogram(ticktock_stored.*1000);
    legend(histogram1, 'Solve time', 'interpreter', 'latex', 'location', 'northeast', 'Orientation', 'vertical', 'fontsize', 14);
    grid on; grid minor;
    xlim([0 10])
    ylabel('Frequency')
    xlabel('Solve time [ms]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;    
    
    if ishandle(2)
        close(2)
    end
    fig2 = figure(2);
    set(fig2, 'Position', [700 100 600 800])
    tile = tiledlayout(5, 2); 
    tile.TileSpacing = 'compact';
    tile.Padding = 'compact';
    nexttile([1 2])
    subplot0 = plot(t_stored, impact_force_stored(2, :), 'color', 'k' , 'linewidth', 2.0);
    legend(subplot0, 'Disturbance', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[N]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot1 = plot(t_stored, xi_err_stored(2, :), 'color', 'b' , 'linewidth', 2.0);
    legend(subplot1, '$\xi^{err}_{y}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot2 = plot(t_stored, p_c_stored(2, :), '-', 'color', 'r', 'linewidth', 2.0);
    legend(subplot2, '$p_{\mathrm{c,ZMP},y}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot3 = plot(t_stored, dU_stored(2, :), '-', 'color', [0, 0.5, 0], 'linewidth', 2.0);
    legend(subplot3, '$\delta u_{y}$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot4 = plot(t_stored, db_stored(2, :), '-', 'color', 'k', 'linewidth', 2.0);
    legend(subplot4, '$\delta b_{y}$', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[m]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot5 = plot(t_stored, dT_stored, '-', 'color', [0.75, 0, 0.75], 'linewidth', 2.0);
    legend(subplot5, '$\delta T$', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[s]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 1])
    subplot6 = plot(t_stored, ddtheta_result_stored(1, :).*PARA.R2D, 'color', [0, 0.75, 0.75], 'linewidth', 2.0);
    legend(subplot6, '$\ddot{\theta}_{roll}$', 'interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
    grid on; grid minor;
    ylabel('[deg/s^2]')
    xlabel('Time [s]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    nexttile([1 2])
    histogram1 = histogram(ticktock_stored.*1000);
    legend(histogram1, 'Solve time', 'interpreter', 'latex', 'location', 'northeast', 'Orientation', 'vertical', 'fontsize', 14);
    grid on; grid minor;
    xlim([0 10])
    ylabel('Frequency')
    xlabel('Solve time [ms]')
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
end
%---

%--- Figure function
function printfig(~,evnt)
    global dx dy;
    global flag_EXIT flag_PAUSE;
    
    if double(evnt.Character) == 28
        dx = -5;
    elseif double(evnt.Character) == 29
        dx = 5;
    elseif double(evnt.Character) == 30
        dy = 5;
    elseif double(evnt.Character) == 31
        dy = -5;
    else
        dx = 0; dy = 0;
    end
    
    if evnt.Character == 'e'
        flag_EXIT = 1;
    elseif evnt.Character == 'f'
        flag_PAUSE = 1;        
    end
end

