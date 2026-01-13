function output = interpolation_solver(N,m,points,matrices,fV,longer_walks)
% This function computes the matrix version of the *-Total Communicability
% based on the section 3.5 using all the theory from chapter 3
% INPUT:
%     N = number of nodes in the network
%     m = truncation parameter
%     points = list of time points rescaled to [-1,1]
%     matrices = adjacency matrices in each time point
%     fv = cell with the functions f(t) of interest
%     longer_walks = True/False parameter that represent if we want classic
%                    alpha_j = 1 or alpha_j = j version of *-Total Communicability

% OUTPUT:
%     output = matrix of the node weights in each time point for each node

tol_linsyst = 1e-12; % Tolerance for the GMRES solver
maxit = 100; %max inner iterstions for GMRES
e = ones(N, 1); % vector of all ones

% Scaling term in case we work with interval that is not [-1,1]
scaling_term = 2/(points(end)-points(1));

% Matrix where we store all the vectors phi
phi_vecs = zeros(length(points),m);

% Here we calculate the matrix representation of the polynomials l_j(t)
fV_H{1} = @(t) 1+0*t; %Special term needed when we don't calculate continuous case
[Fcell, ~] = genCoeffMatrix_SP_interval(fV_H, m, [points(1), points(end)]);
H = Fcell{1};

[Fcell, ~] = genCoeffMatrix_SP_interval(fV, m, [points(1), points(end)]);

% Here we calculate the vector phi(-1)
phim1 = zeros(m,1); % phi(-1)
for i=1:m
    phim1(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);
end

% Right hand side
b = kron(e,phim1);

% Here we solve *-Total Communicability depending if we want alpha_j = 1 or alpha_j = j
% in vectorized way
aprod = @(V) calculate_sum_term(V, matrices, Fcell, N, m);
if longer_walks == 0
    sol = gmres(@(v) aprod(reshape(v,m,N)), b, maxit, tol_linsyst);
    sol = H * reshape(sol,m,N);

elseif longer_walks == 1
    sol1 = gmres(@(v) aprod(reshape(v,m,N)), b, maxit, tol_linsyst);
    sol_vec = gmres(@(v) aprod(reshape(v,m,N)), sol1, maxit, tol_linsyst);

    V_sol = reshape(sol_vec, m, N);
    big_A_times_sol = zeros(m, N);
    
    for k = 1:length(matrices)
        big_A_times_sol = big_A_times_sol + Fcell{k} * V_sol * matrices{k}.';
    end
    
    sol = H * big_A_times_sol;
end

% Here we calculate vectors phi(t_k) for each time point t_k
% In case 'points' is not interval [-1,1] the t_k is rescaled so that it
% is corresponding value in interval [-1,1]
for j = 1:length(points)
    %tau = 2 * (points(j) - points(1)) / (points(end) - points(1)) - 1;
    tau = points(j);
        
    phi_vec = zeros(m,1); % phi(-1)
    for i=0:m-1
        phi_vec(i+1) = legendreP(i, tau) * sqrt((2*i + 1)/2);
    end
    phi_vecs(j,:) = phi_vec.';
end

% Final solution for desired time points
partial_solutions = (phi_vecs * sol * scaling_term).';
%output = partial_solutions(:,2:end);

if longer_walks == 1
    output = partial_solutions + ones(size(partial_solutions));
else
    output = partial_solutions;
end
end


% -- Helper function for calculating aprod--
function output = calculate_sum_term(V, matrices, Fcell, N, m)
    sum_term = reshape(zeros(size(V)), N*m, 1);
    leading_term = reshape(V, N*m, 1);

    for k = 1:length(matrices)
        sum_term = sum_term + reshape(Fcell{k}*V*matrices{k}.', N*m, 1);
    end

    output = leading_term - sum_term;

end
