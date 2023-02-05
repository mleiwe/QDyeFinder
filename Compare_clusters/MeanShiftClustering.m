function [cluster_labels, cluster_centers] = mean_shift(data, bandwidth)
% Inputs:
%   data: NxD matrix of N data points in D dimensions
%   bandwidth: bandwidth parameter for mean-shift
% Outputs:
%   cluster_labels: Nx1 vector indicating the cluster index for each data point
%   cluster_centers: CxD matrix of the C cluster centers

% Initialize variables
[N, D] = size(data);
cluster_labels = zeros(N, 1);
cluster_centers = [];

% For each data point
for i = 1:N
    % Initialize the current mean
    mean = data(i, :);
    
    % Repeat until convergence
    while true
        % Calculate the distances from the current mean
        distances = pdist2(data, mean);
        
        % Calculate the new mean using the weighted average
        weights = exp(-distances.^2 / (2 * bandwidth^2));
        new_mean = sum(bsxfun(@times, data, weights), 1) / sum(weights);
        
        % Check for convergence
        if all(new_mean == mean)
            break;
        end
        
        % Update the mean
        mean = new_mean;
    end
    
    % Check if the mean has already been added to the cluster centers
    match = false;
    for j = 1:size(cluster_centers, 1)
        if all(mean == cluster_centers(j, :))
            match = true;
            break;
        end
    end
    
    % If the mean is not already a cluster center, add it
    if ~match
        cluster_centers = [cluster_centers; mean];
    end
    
    % Assign the data point to the closest cluster center
    distances = pdist2(cluster_centers, mean);
    [~, idx] = min(distances);
    cluster_labels(i) = idx;
end

% Plot the resulting clusters (if 2D data)
if D == 2
    colors = lines(length(cluster_centers));
    figure;
    hold on;
    for i = 1:N
        scatter(data(i,1), data(i,2), [], colors(cluster_labels(i),:), 'filled');
    end
    plot(cluster_centers(:,1), cluster_centers(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
end
end
