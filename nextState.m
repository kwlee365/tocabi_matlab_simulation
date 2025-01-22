function [mL, fL, mR, fR] = nextState(x0, ...
                                      theta, COM, w, dCOM, ...
                                      theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, ...
                                      rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon)

[P, c, A, b, G, h] = qpswiftParameters(x0, theta_ref_horizon, COM_ref_horizon, w_ref_horizon, dCOM_ref_horizon, ...
                                           rL_ref_horizon, rR_ref_horizon, etaL_ref_horizon, etaR_ref_horizon);
[v, ~, ~] = qpSWIFT(sparse(P),c,sparse(A),b,sparse(G),h);

if sum(isnan(v)) > 0
    disp(v);
    fprintf("NAN!\n");
    % break;
end

mL = v(PARA.H*PARA.state_length + 1 : PARA.H*PARA.state_length + 3);
fL = v(PARA.H*PARA.state_length + 4 : PARA.H*PARA.state_length + 6);
mR = v(PARA.H*PARA.state_length + 7 : PARA.H*PARA.state_length + 9);
fR = v(PARA.H*PARA.state_length + 10: PARA.H*PARA.state_length + 12);

end