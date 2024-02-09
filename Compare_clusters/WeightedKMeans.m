function [cluster_labels, centroids] = WeightedKMeans(data, k, weights)
% Inputs:
%   data: NxD matrix of N data points in D dimensions
%   k: number of clusters
%   weights: Nx1 vector of weights for each data point
% Outputs:
%   cluster_labels: Nx1 vector indicating the cluster index for each data point
%   centroids: kxD matrix of centroids of the k clusters

% Initialize variables
[N, D] = size(data);
cluster_labels = zeros(N, 1);

% Initialize the centroids randomly
centroids = data(randsample(N, k), :);

% Repeat until convergence
while true
    % Assign each data point to the closest centroid
    distances = pdist2(data, centroids);
    [~, cluster_labels] = min(distances, [], 2);
    
    % Check if the centroids have changed
    centroids_prev = centroids;
    
    % Recalculate the centroids
    for i = 1:k
        idx = (cluster_labels == i);
        centroids(i,:) = sum(bsxfun(@times, data(idx,:), weights(idx)), 1) / sum(weights(idx));
    end
    
    % Check for convergence
    if all(all(centroids == centroids_prev))
        break;
    end
end
