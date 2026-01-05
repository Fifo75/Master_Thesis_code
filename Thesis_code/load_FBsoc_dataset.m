function [N_max, temporal_matrices] = load_FBsoc_dataset(filename)
% This function loads the FBsoc dataset and 
% creates the adjacency matrices according to the time granularity (month, week, day)
% INPUT:
%      filename = name of the file to load (FBsoc.txt)
% OUTPUT:
%       N = number of nodes in network
%       temporal_matrices = cell of the adjacency matrices from the dataset
%filename = 'FBsoc.txt';
fprintf('Loading data from %s...\n', filename);

% Open the text file for reading
fileID = fopen(filename, 'r');
if fileID == -1
    error('Could not open the file. Make sure it is in the correct path.');
end

% Define the format of each line in the file
% %q -> reads a double-quoted string and removes the quotes
% %d -> reads a signed integer
formatSpec = '%q %d %d %d';

% Read the file according to the format
% The output is a cell array where each cell contains a column
data_columns = textscan(fileID, formatSpec);

% Close the file
fclose(fileID);

% The first column is a cell array of date strings. Convert it to a
% MATLAB datetime array, which is much more useful for analysis.
timestamps = datetime(data_columns{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% The other columns are numeric. Concatenate them into a single matrix.
% Column 1: Source, Column 2: Target, Column 3: Weight
network_data = [data_columns{2}, data_columns{3}, data_columns{4}];

% Scan all source and target nodes to find the maximum node ID.
% This ensures all matrices will have the same dimensions.
all_nodes = [network_data(:,1); network_data(:,2)];
N_max = max(all_nodes);
fprintf('Global network size (N_max): %d\n', N_max);

% Find the start date of each month present in the data
month_starts = dateshift(timestamps, 'start', 'month');
unique_months = unique(month_starts);

% Store the final matrices in a cell array for easy access
temporal_matrices = cell(length(unique_months), 1);

fprintf('Generating a sparse matrix for each month...\n');

for i = 1:length(unique_months)
    current_month = unique_months(i);
    
    % Find the start of the next month to define the time interval
    next_month = dateshift(current_month, 'start', 'month', 1);
    
    % --- Filter data for the current month ---
    % Get logical indices for rows that fall within this month's interval
    is_in_month = (timestamps >= current_month) & (timestamps < next_month);
    
    monthly_data = network_data(is_in_month, :);
    
    if isempty(monthly_data)
        % If no data for this month, create an empty sparse matrix
        A_month = sparse(N_max, N_max);
    else
        % --- Efficiently populate the sparse matrix ---
        % accumarray will group all weights for the same (source, target)
        % pair and sum them up by default. This is very fast.
        subs = monthly_data(:, 1:2); % Subscripts (source, target pairs)
        vals = double(monthly_data(:, 3));   % Values (weights)
        
        A_month = accumarray(subs, vals, [N_max, N_max], @sum, 0, true);
    end
    
    % Store the resulting sparse matrix
    %temporal_matrices{i} = A_month/normest(A_month, 1e-15);
    temporal_matrices{i} = A_month/normest(A_month, 1e-15);
    
    % Display progress
    fprintf('  - Created %d x %d sparse matrix for %s with %d non-zero entries.\n', ...
            N_max, N_max, datestr(current_month, 'mmmm yyyy'), nnz(A_month));
end

% You can now access the matrix for the first month with temporal_matrices{1}, etc.
disp('Processing complete.');
end