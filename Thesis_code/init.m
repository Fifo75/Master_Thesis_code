%% Experiment 1
% The starting adjacency matrix from Experiment 1 section in the thesis
% Contains only two isolated components
matrix = [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %1
          0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %2
          0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %3
          1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %4
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %5
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %6
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %7
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %8
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %9
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %10
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %11
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %12
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %13
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1; %14
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0; %15
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0; %16
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];%17

N = length(matrix); % Number of nodes in the network
n = 10; % number of points to generate
points = linspace(0, n-1, n); % interval [a,b] containing equidistant points
m = 300; % truncation parameter (use max 1500)

% The procedure below creates for each point from 'points' corresponding
% adjacency matrix and stores it to cell called 'matrices'
matrices = cell(1,n);
start_node = 4;
for i = 1:10
    matrix(start_node, start_node+1) = 1;
    matrix(start_node+1, start_node) = 1;
    matrices{i} =  matrix;
    start_node = start_node + 1;
end

% Rescaled interval from [a,b] -> [-1,1]
rescaled_points = (2*points/points(end)) - 1;

% Get the function handles L_i(t) and store them in fV
fV = lagrange_cardinal_basis(rescaled_points); % fV{j} is the handle for L_{j-1}(t)

disp("Lagrangian polynomials calculated");

% In the code bellow we solve the *-Total Communicability centrality for 
% alpha_j = 1 and alpha_j = j
alpha_one_sol = interpolation_solver(N,m,rescaled_points,matrices,fV,0);
alpha_j_sol = interpolation_solver(N,m,rescaled_points,matrices,fV,1);
disp("Interpolation solution calculated");

% Set of nodes we want to study
important_nodes = [4, 8, 9, 10, 14];
%%
% Function that calculates broadcast centralities 
% for Katz and matrix exponential 
% + it creates plots of the temporal evaluation of the nodes
temporal_node_evolution(N, matrices, alpha_one_sol, alpha_j_sol, important_nodes, rescaled_points)
%% Experiment 2
% The matrix of the isolated component with nodes {14,15,16,17} from Example 1 
component_matrix = [0,1,1,1;
                    1,0,0,0;
                    1,0,0,0;
                    1,0,0,0];

N = length(component_matrix); % Number of nodes in the network
n = 10; % number of points to generate
points = linspace(0, n-1, n); % interval [a,b] containing equidistant points
m = 300; % truncation parameter (use max 1500)

% The matrices here represent case 1 and case 2 in Experiment 2
matrices_constant = cell(1,n);
matrices_empty = cell(1,n);

for i = 1:n
    if i == 1
        matrices_empty{i} = component_matrix;
    else
        matrices_empty{i} = zeros(size(component_matrix));
    end
    matrices_constant{i} = component_matrix;
end

% Rescaled interval from [a,b] -> [-1,1]
rescaled_points = (2*points/points(end)) - 1;

% Get the function handles L_i(t) and store them in fV
fV = lagrange_cardinal_basis(rescaled_points); % fV{j} is the handle for L_{j-1}(t)
disp("Lagrangian polynomials calculated");

% In the code bellow we solve the *-Total Communicability for 
% alpha_j = 1 and alpha_j = j in constant case
alpha_one_sol_matrices_constant = interpolation_solver(N,m,rescaled_points,matrices_constant,fV,0);
alpha_j_sol_matrices_constant = interpolation_solver(N,m,rescaled_points,matrices_constant,fV,1);
disp("Interpolation solution calculated for matrices_constant");

% In the code bellow we solve the *-Total Communicability for 
% alpha_j = 1 and alpha_j = j in empty case
alpha_one_sol_matrices_empty = interpolation_solver(N,m,rescaled_points,matrices_empty,fV,0);
alpha_j_sol_matrices_empty = interpolation_solver(N,m,rescaled_points,matrices_empty,fV,1);
disp("Interpolation solution calculated for matrices_empty");

% --- The following part of code is exact from the Experiment 1 ---

% The starting adjacency matrix from Experiment 1 section in the thesis
% Contains only two isolated components
matrix = [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %1
          0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %2
          0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0; %3
          1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %4
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %5
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %6
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %7
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %8
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %9
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %10
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %11
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %12
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; %13
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1; %14
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0; %15
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0; %16
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];%17

% The procedure below creates for each point from 'points' corresponding
% adjacency matrix and stores it to cell called 'matrices'
matrices = cell(1,n);
start_node = 4;
max_spectral_radius = 0;
for i = 1:n
    matrix(start_node, start_node+1) = 1;
    matrix(start_node+1, start_node) = 1;
    matrices{i} =  matrix;
    start_node = start_node + 1;
end

% In the code bellow we solve the *-Total Communicability for 
% alpha_j = 1 and alpha_j = j
alpha_one_sol = interpolation_solver(17,m,rescaled_points,matrices,fV,0);
alpha_j_sol = interpolation_solver(17,m,rescaled_points,matrices,fV,1);
disp("Interpolation solution calculated");

% --- The code from Experiment 1 ends here ---

% Tables with final results

selected_nodes = [14,15,16,17];
alpha_one_sol_selected_rows = alpha_one_sol(selected_nodes,:);
alpha_j_sol_selected_rows = alpha_j_sol(selected_nodes,:);

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = 1) in Example 1\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_one_sol_selected_rows(r,:));
end

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = 1) in constant case\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_one_sol_matrices_constant(r,:));
end

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = j) in Example 1\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_j_sol_selected_rows(r,:));
end

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = j) in constant case\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_j_sol_matrices_constant(r,:));
end

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = j) in empty case\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_one_sol_matrices_empty(r,:));
end

fprintf('\n=================================\n');
fprintf('  Time Evaluation using *-Total Communicability for (alpha_j = j) in empty case\n');
fprintf('---------------------------------\n');
fprintf('  Nodes |  t_0  |  t_1  |  t_2  |  t_3  |  t_4  |  t_5  |  t_6  |  t_7  |  t_8  |  t_9\n');
fprintf('---------------------------------------------------\n');

for r = 1:length(selected_nodes)
    fprintf('  %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f |\n', selected_nodes(r), alpha_j_sol_matrices_empty(r,:));
end

%% Experiment 3
% This function loads the FBsoc dataset
[N, temporal_matrices] = load_FBsoc_dataset("FBsoc.txt");

m = 300; %Max use 1500
n = length(temporal_matrices); %the number of ponits to generate
points = linspace(0,n-1,n); % interval [a,b] containing equidistant points

% Rescaled interval from [a,b] -> [-1,1]
rescaled_points = (2*points/points(end)) - 1;

% Get the function handles L_i(t) and store them in fV
fV = lagrange_cardinal_basis(rescaled_points); % fV{j} is the handle for L_{j-1}(t)

disp("Lagrangian polynomials calculated");

% In the code bellow we solve the *-Total Communicability for 
% alpha_j = 1 and alpha_j = j
alpha_one_sol = interpolation_solver(N,m,rescaled_points,temporal_matrices,fV,0);
alpha_j_sol = interpolation_solver(N,m,rescaled_points,temporal_matrices,fV,1);

disp("Interpolation solution calculated")

% Here we solve the problem using exact solver
exact_solutions_poly = exact_sol(N,rescaled_points, temporal_matrices, fV, 1);
disp("Exact solution calculated")

% Here we calculate relative error of the approximated and exact solution
relative_errors_poly = vecnorm(exact_solutions_poly - alpha_one_sol) ./ vecnorm(exact_solutions_poly);

% Set of nodes we want to study
important_nodes = [9, 12, 67, 431, 557, 1624];

% Function that calculates broadcast centralities 
% for Katz and matrix exponential 
% + it creates plots of the temporal evaluation of the nodes
temporal_node_evolution(N, temporal_matrices, alpha_one_sol, alpha_j_sol, important_nodes, points)

%Visualization of the relative error of exact solution and alpha_one_sol
figure
disp("Relative Errors: ")
disp(relative_errors_poly)
semilogy(relative_errors_poly, LineWidth=2)
title("Relative errors Star approach vs Exponential approach")

% Visualization of the first 100 components of the result vectors
% comparing exact_sol and alpha_one_sol in each time point
figure(2)
for i = 1:n
    % Activate the next subplot
    nexttile; 
    
    % Plot only the i-th component
    %plot(points(2:end), real(interpolation_sol(i,:)), 'b-', 'LineWidth', 2);
    plot(linspace(1, 100, 100), real(alpha_one_sol(1:100,i)), 'b-', 'LineWidth', 2);
    hold on;
    plot(linspace(1, 100, 100), real(exact_solutions_poly(1:100,i)), 'r--', 'LineWidth', 2);
    hold off;
    
    % Add a simple legend and title to THIS subplot
    legend('Approximated sol.', 'Exact sol.');
    title(sprintf('Component %d', i));
    ylabel('Value');
    grid on;
end

% Add a single x-label for the entire figure
xlabel('Time (t)');

