classdef PARA < handle
    properties (Constant)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H = 20; % : Number of horizons
        dt_MPC = 0.01; % [s] : sampling time of MPC        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % World
        D2R = pi/180;
        R2D = 180/pi;
        g = 9.81;
        dt = 0.001;
        
        % Robot
        pelvis_width = 0.22;
        zc = 0.74;
        Foot_length = 0.30;
        Foot_width  = 0.13;
        Foot_up_max = 0.12;
        J_x = 2;
        J_y = 2;
 
        % MPC
        T_scale  = PARA.dt_MPC/PARA.dt;
        T_window = PARA.dt_MPC * PARA.H;
        N_window = PARA.T_window/PARA.dt;
        state_length = 15; % [p_c; dp_c_x; theta; dtheta; c] in R^{15}
        input_length = 6;  % [F, m] in R^{9}
        
        mu = 0.7;

        dc_x_max =  0.2;
        dc_y_max =  0.2;
        dc_x_min = -0.2;
        dc_y_min = -0.2;
        
        r_x_max = 0.2;
        r_x_min =-0.2;
        r_y_max = 0.2;
        r_y_min =-0.2;

        % Kinematics - from Gazelle`s data
        m_PEL = 15;
        m_TORSO = 5;
        m_LHY = 0.01;
        m_LHR = 3;
        m_LHP = 6;
        m_LKN = 2;
        m_LAP = 0.01;
        m_LAR = 1;
        m_RHY = 0.01;
        m_RHR = 3;
        m_RHP = 6;
        m_RKN = 2;
        m_RAP = 0.01;
        m_RAR = 1;
        m_all = PARA.m_PEL + PARA.m_TORSO + (PARA.m_LHY+PARA.m_LHR+PARA.m_LHP+PARA.m_LKN+PARA.m_LAP+PARA.m_LAR) + (PARA.m_RHY+PARA.m_RHR+PARA.m_RHP+PARA.m_RKN+PARA.m_RAP+PARA.m_RAR);
        
        l0 = 0.15;
        l1 = 0.105;
        l2 = 0.05;
        l3 = 0.06;
        l4 = 0.4;
        l5 = 0.38;
        l6 = 0.07;
        
        c_PEL = [0; 0; 0];
        c_TORSO = [0; 0; 0];
        c_LHY = [0; 0; 0];
        c_LHR = [0; 0.5*PARA.l2; -0.5*PARA.l3];
        c_LHP = [0; 0; -0.5*PARA.l4];
        c_LKN = [0; 0; -0.5*PARA.l5];
        c_LAP = [0; 0; 0];
        c_LAR = [0; 0; -0.5*PARA.l6];
        c_RHY = [0; 0; 0];
        c_RHR = [0; -0.5*PARA.l2; -0.5*PARA.l3];
        c_RHP = [0; 0; -0.5*PARA.l4];
        c_RKN = [0; 0; -0.5*PARA.l5];
        c_RAP = [0; 0; 0];
        c_RAR = [0; 0; -0.5*PARA.l6];
    end
end