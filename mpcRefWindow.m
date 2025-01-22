function [theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon] = ...
    mpcRefWindow(t_step, step_phase, Foot_state, L_or_R, LF_prev, RF_prev, COM_ref, dCOM_ref, ddCOM_ref, p_err_sum_x_ref, p_err_sum_y_ref, T_step_ref, p_ref_total, T_step_ref_total, Gi_MPC, Gx_MPC, Gp_MPC)
    
    theta_ref_horizon = zeros(3, PARA.H);    
    COM_ref_horizon = zeros(3, PARA.H);    
    w_ref_horizon = zeros(3, PARA.H);    
    dCOM_ref_horizon = zeros(3, PARA.H);    
    
    rL_ref_horizon = zeros(3, PARA.H);    
    rR_ref_horizon = zeros(3, PARA.H);    
    
    etaL_ref_horizon = ones(PARA.H, 1);    
    etaR_ref_horizon = ones(PARA.H, 1);    

    % etaL_ref_horizon = zeros(PARA.H, 1);    
    % etaR_ref_horizon = zeros(PARA.H, 1);  % why?

    t_step_temp = t_step; T_step_ref_temp = T_step_ref;
    step_phase_temp = step_phase; Foot_state_temp = Foot_state; LF_prev_temp = LF_prev; RF_prev_temp = RF_prev;
    COM_ref_temp = COM_ref; dCOM_ref_temp = dCOM_ref; ddCOM_ref_temp = ddCOM_ref;
    p_err_sum_x_ref_temp = p_err_sum_x_ref; p_err_sum_y_ref_temp = p_err_sum_y_ref;

    for j = 1:PARA.H
        flag_STEP_CHANGE_TEMP = checkStepEnd(t_step_temp, T_step_ref_temp);
       
        if flag_STEP_CHANGE_TEMP == 1
            T_step_ref_temp_prev = T_step_ref_temp;
            
            % Update step phase
            if step_phase_temp < length(p_ref_total)
                step_phase_temp = step_phase_temp + 1;
            end
            
            % Update step time
            T_step_ref_temp = T_step_ref_total(:, step_phase_temp);
            
            % Reset t_step & flag
            t_step_temp = t_step_temp - T_step_ref_temp_prev;
            if t_step_temp < 1E-06
                t_step_temp = 0;
            end
            
            % Update Foot state
            Foot_state_temp = (-1)*Foot_state_temp;
            if step_phase_temp == 2
                Foot_state_temp = L_or_R;
            end

            % Update Support Foot position
            if Foot_state_temp == 2
                LF_prev_temp = LF_prev_temp;
                RF_prev_temp = RF_prev_temp;
            elseif Foot_state_temp ==  1 % LF swing
                LF_prev_temp = LF_prev_temp;
                RF_prev_temp = p_ref_total(:, step_phase_temp);
            elseif Foot_state_temp == -1 % RF swing
                LF_prev_temp = p_ref_total(:, step_phase_temp);
                RF_prev_temp = RF_prev_temp;
            end
        end

        if Foot_state_temp == 2
            etaL_ref_horizon(j) = 1;            
            etaR_ref_horizon(j) = 1;
        elseif Foot_state_temp ==  1 % LF swing
            etaL_ref_horizon(j) = 0;            
            etaR_ref_horizon(j) = 1;
        elseif Foot_state_temp == -1 % RF swing
            etaL_ref_horizon(j) = 1;            
            etaR_ref_horizon(j) = 0;
        end
                
        [COM_ref_temp_next, dCOM_ref_temp_next, ddCOM_ref_temp_next, p_err_sum_x_ref_temp_next, p_err_sum_y_ref_temp_next] = previewControlMPC(t_step_temp, step_phase_temp, p_ref_total, T_step_ref_total, Gi_MPC, Gx_MPC, Gp_MPC, PARA.A_preview_MPC, PARA.B_preview_MPC, PARA.C_preview_MPC, COM_ref_temp, dCOM_ref_temp, ddCOM_ref_temp, p_err_sum_x_ref_temp, p_err_sum_y_ref_temp);  
        
        theta_ref_horizon(:, j) = zeros(3,1);
        COM_ref_horizon(:, j)  = COM_ref_temp;
        
        w_ref_horizon(:, j)     = zeros(3,1);
        dCOM_ref_horizon(:, j) = dCOM_ref_temp;

        if Foot_state_temp == 2
            rL_ref_horizon(:, j) = LF_prev_temp - COM_ref_temp;
            rR_ref_horizon(:, j) = RF_prev_temp - COM_ref_temp;
        elseif Foot_state_temp ==  1 % LF swing
            rL_ref_horizon(:, j) = zeros(3,1);
            rR_ref_horizon(:, j) = RF_prev_temp - COM_ref_temp;
        elseif Foot_state_temp == -1 % RF swing
            rL_ref_horizon(:, j) = LF_prev_temp - COM_ref_temp;
            rL_ref_horizon(:, j) = zeros(3,1);
        end
        
        t_step_temp = t_step_temp + PARA.dt_MPC;
        
        COM_ref_temp = COM_ref_temp_next; dCOM_ref_temp = dCOM_ref_temp_next; ddCOM_ref_temp = ddCOM_ref_temp_next;
        p_err_sum_x_ref_temp = p_err_sum_x_ref_temp_next; p_err_sum_y_ref_temp = p_err_sum_y_ref_temp_next;
    end
end