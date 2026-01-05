function u_exact = exact_sol(N,points,matrices, fV, continuous)
% This function conputes exact solution of the Time-Order Exponential 
% using matrix expoential in piece wise constant case or 
% ode45 in continuous case
% INPUT:
%     N = number of nodes in the network
%     points = list of time points rescaled to [-1,1]
%     matrices = adjacency matrices in each time point
%     fv = cell with the functions f(t) of interest
%     continuous = True/False parameter that represent if we want classic
%                  piece wise constant or continuous version of TOE

% OUTPUT:
%     u_exact = matrix of the node weights in each time point for each node

% Vector of ones
e = ones(N,1);

% Here we calculate the exact solution depending on the 'continuous' variable
if continuous == 0
    u_exact = zeros(N, length(points));
    u_exact(:,1) = ones(N,1);

    for i = 2:length(points)
        u_exact(:,i) = expm(matrices{i-1} * (points(i) - points(i-1))) * u_exact(:,i-1);
    end
    %u_exact = u_exact(:,2:end);

elseif continuous == 1
    % --- Define the ODE function for ode45
    % The function must have the signature odefun(t, u)
    % It returns the derivative du/dt at that time and state.
    %odefun = @(t, u) (A0 + A1 * fV{2}(t)) * u;
    odefun = @(t, u) matrix_polynomial(matrices, fV, t, u);
    
    % --- Set the time span and call the solver
    sol = ode45(odefun, [points(1), points(end)], e);

    % --- Final solution
    u_ode45 = deval(sol, points);
    u_exact = u_ode45;
end
end

% Helper function that calculates matrix interpolation polynomial
function output = matrix_polynomial(matrices, fV, t, u)
    sum_term = zeros(size(u));
    for j = 1:length(matrices)
        sum_term = sum_term + matrices{j}*(fV{j}(t)*u);
    end

    output = sum_term;
end