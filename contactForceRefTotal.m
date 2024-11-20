function f_total = contactForceRefTotal(number_of_step, L_step, W_step, L_or_R, m, g)
if number_of_step < 1
    disp('number_of_step should be positive!! (number_of_step < 1)');
    number_of_step = 1;
end

f_total = zeros(6, number_of_step + 3);

for i = 1:number_of_step+3
    if i == 1
        f_total(:, i) = zeros(3, 1);
    elseif i == 2
        f_total(:, i) = [0;0;L_or_R*m*g];
    elseif i == 3
        f_total(:, i) = [0;0;L_or_R*m*g];
    elseif i == number_of_step+2
        f_total(:, i) = [0;0;L_or_R*m*g];
    elseif i == number_of_step+3
        f_total(:, i) = [0;0;L_or_R*m*g];
    else
        f_total(:, i) = [0;0;L_or_R*m*g];
    end
    L_or_R = (-1)*L_or_R;
end
end