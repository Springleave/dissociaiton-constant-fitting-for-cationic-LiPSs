% Define the optimization problem
problem = createOptimProblem('fmincon', 'objective', @residue, 'x0', 0.01489, 'lb', 0.0001152, 'ub', 1, 'options', optimset('Algorithm', 'sqp', 'Disp', 'none'));

% Create a MultiStart object
ms = MultiStart;

% Run the optimization multiple times and obtain the best solution
xm = run(ms, problem, 4);

% Kd_LiTFSI = Li * TFSI / LiTFSI
% Kd_Li3S6 = Li3S6 * TFSI / (Li2S6 * LiTFSI)
% so K2 in article equals (Kd_Li3S6 / Kd_LiTFSI)
K2_value = xm / 6.7622e-4;

% Calculate the residue value using the best solution
fval = double(residue(xm));

% Define the residue function
function rd = residue(input)
    Kd_Li3S6 = input(1);
    v0_Li3S6 = 1.44e4;
    x = [0 0.001 0.002 0.004 0.006 0.008 0.01];
    y = [10.9 12 12.82 16.21 22.4 25.5 28 30.9 32.5 33.1 36.6 39 43.9 46.1];
    y_simu = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    % Calculate simulated values for the given input
    for i = 1:7
        y_simu(i) = cond(x(i), v0_Li3S6, Kd_Li3S6, 0.001);
    end
    
    for i = 8:14
        y_simu(i) = cond(x(i-7), v0_Li3S6, Kd_Li3S6, 0.005);
    end
    
    % Calculate the residue by taking the norm of the difference between
    % the actual and simulated values
    temp = y - y_simu;
    temp = [temp(1:7) 3 * temp(8:14)];
    rd = norm(temp, 1);
end

% Define the calc function
function result = calc(input)
    v0_Li3S6 = input(1);
    Kd_Li3S6 = input(2);
    x = [0 0.001 0.002 0.004 0.006 0.008 0.01];
    calc = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    % Calculate simulated values for the given input
    for i = 1:7
        calc(i) = cond(x(i), v0_Li3S6, Kd_Li3S6, 0.001);
    end
    
    for i = 8:14
        calc(i) = cond(x(i-7), v0_Li3S6, Kd_Li3S6, 0.005);
    end
    
    result = calc;
end

% Define the conductivity function
function p = cond(n_Li2S6, v0_Li3S6, Kd_Li3S6, n_LiTFSI)
    syms LiS6 Li2S6 Li Li3S6 LiTFSI TFSI;
    Kd_Li2S6 = 1.6855e-04;
    Kd_LiTFSI = 6.7622e-04;
    v0_Li2S6 = 1.3812e+04;
    v0_LiTFSI = 2.001e+04;
    
    % Define the equations for the constraints
    eq1 = LiS6 + Li2S6 + Li3S6 == n_Li2S6;
    eq2 = LiTFSI + TFSI == n_LiTFSI;
    eq3 = LiS6 + TFSI == Li3S6 + Li;
    eq4 = Li * LiS6 / Li2S6 == Kd_Li2S6;
    eq5 = Li * TFSI / LiTFSI == Kd_LiTFSI;
    eq6 = Li3S6 * TFSI / (Li2S6 * LiTFSI) == Kd_Li3S6;
    
    % Solve the equations to obtain the values of the variables
    [rLiS6, rLi2S6, rLi, rLi3S6, rLiTFSI, rTFSI] = vpasolve([eq1, eq2, eq3, eq4, eq5, eq6], [LiS6 Li2S6 Li Li3S6 LiTFSI TFSI], [0, 0.3; 0, 0.3; 0, 0.3; 0, 0.3; 0, 0.3; 0, 0.3]);
    
    % Calculate the value of p using the obtained variables
    p = v0_Li2S6 * rLiS6 + v0_Li3S6 * rLi3S6 + v0_LiTFSI * (rTFSI - rLi3S6);
end