function L_handles = lagrange_cardinal_basis(points)
% Calculates the function handles for the Lagrange Cardinal Basis L_j(t)
% using Barycentric Interpolation for EQUIDISTANT POINTS.
%
% INPUT:
%      points: Vector of Equidistant nodes (e.g., 0, 1, 2... or rescaled to -1, 1).
% OUTPUT:
%      L_handles: Cell array where L_handles{j} is the @(t) function handle for L_j(t).

    N = length(points);
    L_handles = cell(1, N);
    
    % --- Step A: Calculate Barycentric Weights (W_j) for EQUIDISTANT POINTS ---
    % For equidistant points, weights follow the binomial coefficients:
    % w_j = (-1)^j * (N-1 choose j)
    
    w = zeros(N, 1);
    n = N - 1; % The degree of the polynomial
    
    for j = 0:n
        % Calculate binomial coefficient (n choose j)
        % Note: For very large N (>50), this can overflow, but for N=8 it is perfect.
        binom_coeff = nchoosek(n, j);
        
        % Apply alternating sign
        w(j+1) = (-1)^j * binom_coeff;
    end
    
    % --- Step B: Define the Vectorized Stable Barycentric Evaluation ---
    for j = 1:N
        % This logic remains exactly the same, only the 'w' values changed.
        L_handles{j} = @(t) stable_barycentric_eval(t, points, w, j);
    end
end

function L_t = stable_barycentric_eval(t, points, w, j_index)
% Internal function to evaluate the j_index-th Cardinal Basis polynomial L_{j-1}(t)
% This is the "Second Barycentric Form" - it works for any valid w and points.
    
    % t can be a scalar or a vector
    if isscalar(t)
        t = t(:); % Ensure column vector for matrix operations
    end
    
    N = length(points);
    
    % Create a matrix of differences: Diff(k, i) = t(k) - points(i)
    Diff_matrix = t - points(:)'; 
    
    % Create R_matrix: R(k, i) = w(i) / (t(k) - points(i))
    R_matrix = repmat(w(:)', length(t), 1) ./ Diff_matrix; 
    
    % Sum R_matrix across each row to get the denominator vector: Denom(k) = sum_i R(k, i)
    Denom = sum(R_matrix, 2);
    
    % Numerator vector: Num(k) = R(k, j_index)
    Num = R_matrix(:, j_index);
    
    % Standard case: L_t(k) = Num(k) / Denom(k)
    L_t = Num ./ Denom;
    
    % --- Handle evaluation exactly at nodes (where t - x_j = 0) ---
    % Find indices where t(k) is a node x_i to avoid division by zero (NaN)
    [row_k, col_i] = find(abs(Diff_matrix) < 1e-12);
    
    if ~isempty(row_k)
        for idx = 1:length(row_k)
            k = row_k(idx); % Index in t
            i = col_i(idx); % Index of the node (points(i))
            
            % If t(k) is exactly the j_index-th node, L_j(x_j) = 1
            if i == j_index
                L_t(k) = 1; 
            % If t(k) is another node x_i (i != j_index), L_j(x_i) = 0
            else 
                L_t(k) = 0; 
            end
        end
    end
    
    if isrow(t)
        L_t = L_t';
    end
end