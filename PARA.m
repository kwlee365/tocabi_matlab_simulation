classdef PARA < handle
    properties (Constant)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H = 15; % : Number of horizons
        dt_MPC = 0.02; % [s] : sampling time of MPC        
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
        zc = 0.73;
        wn = sqrt(PARA.g/PARA.zc);
        Foot_length = 0.18;
        Foot_width = 0.13;
        Foot_up_max = 0.12;

        I = [18 0 0;
             0 15 0;
             0 0 10];
        
        % Preview control
        T_preview = 1.6;
        NL = PARA.T_preview/PARA.dt;
        A_preview = [1, PARA.dt, (PARA.dt*PARA.dt)/2;
                     0, 1, PARA.dt; 
                     0, 0, 1];
        B_preview = [(PARA.dt*PARA.dt*PARA.dt)/6; 
                     (PARA.dt*PARA.dt)/2; 
                      PARA.dt];
        C_preview = [1, 0, -(1/(PARA.wn^2))];
        
        NL_MPC = PARA.T_preview/PARA.dt_MPC;
        A_preview_MPC = [1, PARA.dt_MPC, (PARA.dt_MPC*PARA.dt_MPC)/2;
                         0, 1, PARA.dt_MPC; 
                         0, 0, 1];
        B_preview_MPC = [(PARA.dt_MPC*PARA.dt_MPC*PARA.dt_MPC)/6; 
                         (PARA.dt_MPC*PARA.dt_MPC)/2; 
                         PARA.dt_MPC];
        C_preview_MPC = [1, 0, -(1/(PARA.wn^2))];
        
        % MPC
        T_scale = PARA.dt_MPC/PARA.dt;
        T_window = PARA.dt_MPC * PARA.H;
        N_window = PARA.T_window/PARA.dt;
        state_length = 12; % [theta, COM, w, dCOM]
        input_length = 12; % [mL, fL, mR, fR];   

        Q_theta = 100000000;
        Q_COM_x = 1000000000;
        Q_COM_y = 1000000000;
        Q_COM_z = 1000000000;
        Q_w     = 10000;
        Q_dCOM  = 10000;
        
        R_mL = 1;
        R_fL = 1;
        R_mR = 1;
        R_fR = 1;

        f_z_max = 500;
        f_z_min = 0;
        mu = 0.7;   % Coulomb's Friction Coefficient
        
        % Kinematics - from TOCABI`s data
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