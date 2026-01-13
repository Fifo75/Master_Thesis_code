function temporal_node_evolution(N, matrices, alpha_one_sol, alpha_j_sol, important_nodes, points)
% This function calculates broadcast centralities of 
% Katz and matrix exponential and plots the summary plots 
% representing temporal node evolution for each of the centralities.
% INPUT:
%     N = number of nodes in the network
%     n = number of time points
%     matrices = adjacency matrices in each time point
%     alpha_one_sol = solution of TOE Centrality for alpha_j = 1
%     alpha_j_sol = solution of TOE Centrality for alpha_j = j
%     important_nodes = nodes of the network we want to evaluate
%     points = list of time points rescaled to [-1,1]

% OUTPUT: Plots of the temporal node evolution in each time point




% Matrices where results of the Katz and matrix exponential broadcast centralities
% will be stored
residual_broadcast_results = zeros(N, length(points));
exponential_broadcast_results = zeros(N, length(points));

% Here we calculate the maximum value of alpha that can be used 
% for Katz receiver centrality
max_spectral_radius = 0;
for i = 1:length(points)
    spectral_radius = max(abs(eigs(matrices{i})));

    if spectral_radius > max_spectral_radius
        max_spectral_radius = spectral_radius;
    end
end

% Here we choose at random parameter alpha based on the definition
alpha = 0.00001 + ((1/max_spectral_radius) - 0.00001) * rand(1); % choose 0 < alpha < 1/max_spectral_radius
display("Chosen alpha:")
display(alpha)

% Here we calculate dynamic communicability matrix for both
% Katz and matrix exponential
for i = 1:length(points)
    if i == 1
        Q_residual = inv(eye(N) - alpha*matrices{1});
        Q_exponential = expm(matrices{1});
    else
        Q_residual = Q_residual * (inv(eye(N) - alpha*matrices{i}));
        Q_exponential = Q_exponential * expm(matrices{i});
    end

    residual_broadcast_results(:,i) = Q_residual * ones(N,1);
    exponential_broadcast_results(:,i) = Q_exponential * ones(N,1);
end

% Evolution of the values of the important nodes
% Firstly we create tensor concatenating all the result matrices 
% from the four centralities
nodes_values = zeros(length(important_nodes), length(points), 4);
results = cat(3, alpha_one_sol, alpha_j_sol, residual_broadcast_results, exponential_broadcast_results);

% Here takes only results for the important nodes
[~,~,dimension] = size(nodes_values);
for j = 1:dimension
    nodes_values(:,:,j) = results(important_nodes,:,j);
end


% Plot the results as a grid 2x2 
figure
counter = 1;
for i = 1:2
    for j = 1:2
        subplot(2, 2, counter)
        plot(points, nodes_values(:,:,counter), 'LineWidth', 2)
        legend(string(important_nodes))
        title(['Result #', num2str(counter)])

        ylabel('Values in time');
        xlabel('Time-points');
        ylim([1, max(max(nodes_values(:,:,counter)))]); % Set the limits to be from 1 to N
        counter = counter + 1;
    end
end
end

